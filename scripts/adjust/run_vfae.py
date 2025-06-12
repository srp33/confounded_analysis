from argparse import ArgumentParser
import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
from vfae import VFAE, mmd_rbf_loss 


# Parse command line args --------------------------
parser = ArgumentParser(description="VFAE implementation with MMD for fairness")
parser.add_argument("-i", "--input-file", help="Path to input CSV file", required=True)
parser.add_argument("-o", "--output-file", help="Path to output CSV file for fair reconstructions", required=True)
parser.add_argument("-b", "--batch-col", help="Column name for the sensitive attribute (batch information)", required=True)
parser.add_argument("-l", "--latent-dim", type=int, default=10, help="Dimensionality of the latent space")
parser.add_argument("-hd", "--hidden-dim", type=int, default=128, help="Dimensionality of the hidden layers")
parser.add_argument("-e", "--epochs", type=int, default=100, help="Number of training epochs")
parser.add_argument("-lr", "--learning-rate", type=float, default=1e-3, help="Learning rate for the optimizer")
parser.add_argument("-bs", "--batch-size", type=int, default=64, help="Batch size for training")
parser.add_argument("--w-kl", type=float, default=1.0, help="Weight for the KL divergence loss term")
parser.add_argument("--w-mmd", type=float, default=10.0, help="Weight for the MMD fairness loss term")
parser.add_argument("--mmd-gamma", type=float, default=1.0, help="Gamma parameter for the RBF kernel in MMD")
args = parser.parse_args()

# Set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Load the dataset --------------------------------
if not args.input_file.endswith('.csv'):
    raise ValueError("Input file must be a CSV file.")
try:
    df = pd.read_csv(args.input_file)
except Exception as e:
    raise ValueError(f"Error reading the input file: {args.input_file}\n{e}")

# Prepare data ------------------------------------
if args.batch_col not in df.columns:
    raise ValueError(f"Sensitive attribute column '{args.batch_col}' not found in the dataframe.")

meta_cols = [col for col in df.columns if col.startswith('meta_')]
gene_cols = [col for col in df.columns if col not in meta_cols and col != args.batch_col]
print(f"Found {len(gene_cols)} gene columns.")

genes = df[gene_cols]

# Store min and max for denormalization later
genes_min = genes.min()
genes_max = genes.max()

# Normalize data to [0, 1] range for BCE loss compatibility
print("Normalizing data to [0, 1] range.")
genes_normalized = (genes - genes_min) / (genes_max - genes_min)
genes_normalized = genes_normalized.fillna(0) # Handle potential NaN from columns with no variance

# Prepare sensitive attribute
s_series = df[args.batch_col]
s_codes, s_uniques = pd.factorize(s_series)
num_s_categories = len(s_uniques)
print(f"Found sensitive attribute '{args.batch_col}' with {num_s_categories} unique categories: {s_uniques.tolist()}")

# Convert data to PyTorch tensors
genes_tensor = torch.tensor(genes_normalized.values, dtype=torch.float32)
s_tensor = torch.tensor(s_codes, dtype=torch.long)

# Create DataLoader
dataset = TensorDataset(genes_tensor, s_tensor)
train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True)

# Initialize model and optimizer ------------------
input_dim = len(gene_cols)
model = VFAE(
    input_dim=input_dim,
    num_s_categories=num_s_categories,
    latent_dim=args.latent_dim,
    hidden_dim_enc=args.hidden_dim,
    hidden_dim_dec=args.hidden_dim
).to(device)

optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)

# Training loop -----------------------------------
print("\nStarting VFAE training...")
for epoch in range(args.epochs):
    model.train()
    total_recon_loss = 0
    total_kl_loss = 0
    total_mmd_loss = 0

    for x_batch, s_batch_indices in train_loader:
        x_batch = x_batch.to(device)
        s_batch_indices = s_batch_indices.to(device)

        optimizer.zero_grad()

        x_recon, mu, logvar, z_samples = model(x_batch, s_batch_indices)

        # --- Calculate Loss Components ---

        # 1. Reconstruction Loss (ensures the model can reconstruct the input)
        recon_loss = F.binary_cross_entropy(x_recon, x_batch, reduction='sum')

        # 2. KL Divergence (regularizes the latent space)
        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

        # 3. MMD Loss (enforces fairness)
        mmd_loss = torch.tensor(0.0, device=device)
        unique_s_in_batch = torch.unique(s_batch_indices)
        
        # Only compute MMD if there are multiple sensitive groups in the batch
        if len(unique_s_in_batch) > 1:
            # Compare each group against a sample from the prior (N(0,1))
            # This is a robust way to handle multi-category 's' and batch imbalances.
            prior_samples = torch.randn_like(z_samples)
            mmd_loss = mmd_rbf_loss(z_samples, prior_samples, gamma=args.mmd_gamma)
            
        # --- Combine Losses ---
        # We average losses over the batch size for stable gradients
        batch_size = x_batch.size(0)
        total_loss = (recon_loss / batch_size) + \
                        (args.w_kl * kl_loss / batch_size) + \
                        (args.w_mmd * mmd_loss) # MMD is already a mean

        total_loss.backward()
        optimizer.step()

        total_recon_loss += recon_loss.item()
        total_kl_loss += kl_loss.item()
        total_mmd_loss += mmd_loss.item()

    # Print epoch statistics
    avg_recon_loss = total_recon_loss / len(dataset)
    avg_kl_loss = total_kl_loss / len(dataset)
    avg_mmd_loss = total_mmd_loss / len(train_loader)
    print(f"Epoch {epoch+1}/{args.epochs} | "
            f"Recon Loss: {avg_recon_loss:.4f} | "
            f"KL Loss: {avg_kl_loss:.4f} | "
            f"MMD Loss: {avg_mmd_loss:.4f}")

print("\nTraining finished.")

# Generate predicted representations without batch information ----
print(f"Generating predicted representations without batch info and saving to {args.output_file}")
model.eval()
with torch.no_grad():
    # 1. Get the fair latent representations (mu) for the entire dataset
    full_genes_tensor = genes_tensor.to(device)
    full_s_tensor = s_tensor.to(device)
    s_one_hot_original = F.one_hot(full_s_tensor, num_classes=num_s_categories).float()
    mu, _ = model.encoder(full_genes_tensor, s_one_hot_original)

    # 2. Reconstruct by marginalizing out the sensitive attribute 's'
    all_recons = []
    for s_cat in range(num_s_categories):
        # Create a tensor of this category for all samples
        s_indices_k = torch.full_like(full_s_tensor, fill_value=s_cat)
        s_one_hot_k = F.one_hot(s_indices_k, num_classes=num_s_categories).float().to(device)
        
        # Decode using the latent codes 'mu' but with the modified 's'
        recon_k = model.decoder(mu, s_one_hot_k)
        all_recons.append(recon_k.cpu())

    # Average the reconstructions across all sensitive attribute categories
    fair_recon_normalized = torch.stack(all_recons).mean(dim=0)

    # 3. Denormalize the data to its original scale
    # Convert min/max to tensors for broadcasting
    genes_min_tensor = torch.tensor(genes_min.values, dtype=torch.float32)
    genes_max_tensor = torch.tensor(genes_max.values, dtype=torch.float32)
    
    fair_recon_unnormalized = fair_recon_normalized * (genes_max_tensor - genes_min_tensor) + genes_min_tensor

    # 4. Create and save the output DataFrame
    output_df = pd.DataFrame(fair_recon_unnormalized.numpy(), columns=gene_cols)

    # Add original metadata for context (but not the sensitive attribute itself)
    final_df = pd.concat([df[meta_cols].reset_index(drop=True), output_df.reset_index(drop=True)], axis=1)

    try:
        final_df.to_csv(args.output_file, index=False)
        print(f"Successfully saved the fair reconstructions to {args.output_file}")
    except Exception as e:
        print(f"Error saving the output file: {e}")





# def loss_function(recon_x, x, mean, log_var, class_probs, y_true, loss_fn='cross_entropy'):
#     # Loss is made of several components:
#     # 1. Reconstruction loss, causes the model to learn to reconstruct the input data
#     # 2. KL divergence, causes the model to learn a latent space that is close to a standard normal distribution
#     # 3. Classification loss, which is negated through the encoder portion. This c

#     reconstruction_loss = torch.nn.functional.binary_cross_entropy(recon_x, x, reduction='sum')
#     kl_divergence = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
#     if loss_fn == 'cross_entropy':
#         classification_loss = torch.nn.functional.binary_cross_entropy_with_logits(class_probs, y_true, reduction='sum')
#     elif loss_fn == 'dmi':
#         # -log determinant of the joint probability distribution
#         classification_loss = -torch.logdet(torch.mm(class_probs.t(), class_probs) + 1e-10)


