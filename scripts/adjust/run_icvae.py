from argparse import ArgumentParser
import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
from icvae import VFAE, AuxiliaryClassifier
import sys
import torch.optim as optim
from jaxtyping import Float, Int
from torch import Tensor
from beartype import beartype

parser = ArgumentParser(description="VFAE for learning fair representations using a Mutual Information penalty.")
parser.add_argument("-i", "--input-file", help="Path to input CSV file.", required=True)
parser.add_argument("-o", "--output-file", help="Path to output CSV file for fair reconstructions.", required=True)
parser.add_argument("-b", "--batch-col", help="Column name for the sensitive attribute (batch).", required=True)
parser.add_argument("-l", "--latent-dim", type=int, default=10, help="Dimensionality of the latent space.")
parser.add_argument("-hd", "--hidden-dim", type=int, default=128, help="Dimensionality of hidden layers for VFAE.")
parser.add_argument("-hda", "--hidden-dim-aux", type=int, default=64, help="Dimensionality of hidden layers for Auxiliary Classifier.")
parser.add_argument("-e", "--epochs", type=int, default=100, help="Number of training epochs.")
parser.add_argument("-lr", "--learning-rate", type=float, default=1e-3, help="Learning rate for optimizers.")
parser.add_argument("-bs", "--batch-size", type=int, default=64, help="Batch size for training.")
parser.add_argument("--w-kl", type=float, default=1.0, help="Weight for the KL divergence loss term (beta).")
parser.add_argument("--w-mi-penalty", type=float, default=10.0, help="Weight for the Mutual Information penalty term (gamma).")
args = parser.parse_args()

def print_now(*args, **kwargs):
    print(*args, flush=True, **kwargs)

# --- Setup ---
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print_now(f"Using device: {device}")

# --- Data Loading and Preprocessing ---
print_now(f"Loading data from '{args.input_file}'...")
try:
    df = pd.read_csv(args.input_file)
except Exception as e:
    print_now(f"Error reading the input file: {e}", file=sys.stderr)
    sys.exit(1)

if args.batch_col not in df.columns:
    print_now(f"Error: Batch column '{args.batch_col}' not found.", file=sys.stderr)
    sys.exit(1)

# Separate feature columns from metadata and sensitive attribute
feature_cols = [col for col in df.columns if not col.startswith('meta_') and col != args.batch_col]
print_now(f"Found {len(feature_cols)} feature columns.")

features = df[feature_cols].copy()
features_min = features.min()
features_max = features.max()

# Normalize data to [-1, 1] range for stability with MSE loss
features_normalized = (2 * (features - features_min) / (features_max - features_min)) - 1
features_normalized = features_normalized.fillna(0)  # Handle columns with no variance

batch_series = df[args.batch_col]
batch_codes, unique_batches = pd.factorize(batch_series)
num_batches = len(unique_batches)
print_now(f"Found batch attribute '{args.batch_col}' with {num_batches} categories: {unique_batches.tolist()}")

features_tensor = torch.tensor(features_normalized.values, dtype=torch.float32)
batch_tensor = torch.tensor(batch_codes, dtype=torch.long)

dataset = TensorDataset(features_tensor, batch_tensor)
train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True)

# --- Model Initialization ---
input_dim = len(feature_cols)
vfae_model = VFAE(input_dim, num_batches, args.latent_dim, args.hidden_dim).to(device)
aux_classifier = AuxiliaryClassifier(args.latent_dim, num_batches, args.hidden_dim_aux).to(device)

opt_vfae = optim.Adam(vfae_model.parameters(), lr=args.learning_rate)
opt_aux = optim.Adam(aux_classifier.parameters(), lr=args.learning_rate)

# --- Training Loop ---
print_now("\nStarting training...")
for epoch in range(args.epochs):
    total_recon_loss, total_kl_loss, total_mi_penalty, total_aux_loss = 0, 0, 0, 0

    for x_batch, s_batch in train_loader:
        x_batch = x_batch.to(device)
        s_batch = s_batch.to(device)

        # --- VFAE forward pass ---
        x_recon, mu, logvar, z = vfae_model(x_batch, s_batch)
        
        # --- Update Auxiliary Classifier: q(s|z) ---
        # Train classifier to predict `s` from `z` (gradients do not flow to encoder)
        s_logits_aux = aux_classifier(z.detach())
        loss_aux = F.cross_entropy(s_logits_aux, s_batch)
        opt_aux.zero_grad()
        loss_aux.backward()
        opt_aux.step()

        # --- Update VFAE (Encoder + Decoder) ---
        # 1. Reconstruction Loss (MSE), causes the model to learn to reconstruct the input data
        recon_loss = F.mse_loss(x_recon, x_batch, reduction='mean')
        # 2. KL Divergence, causes the model to learn a latent space that is close to a standard normal distribution
        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

        # 3. Mutual Information Penalty
        # Train encoder to fool the classifier (gradients flow to encoder)
        s_logits_mi: Float[Tensor, "batch num_batches"] = aux_classifier(z)
        log_probs_s: Float[Tensor, "batch num_batches"] = F.log_softmax(s_logits_mi, dim=1)

        s_batch_: Int[Tensor, "batch 1"] = s_batch.unsqueeze(1)
        mi_penalty: Float[Tensor, ""] = log_probs_s.gather(1, s_batch_).mean()


        total_vfae_loss = recon_loss + args.w_kl * kl_loss + args.w_mi_penalty * mi_penalty
        opt_vfae.zero_grad()
        total_vfae_loss.backward()
        opt_vfae.step()
        
        # Accumulate losses for logging
        total_recon_loss += recon_loss.item()
        total_kl_loss += kl_loss.item()
        total_mi_penalty += mi_penalty.item()
        total_aux_loss += loss_aux.item()

    # --- Logging ---
    num_samples = len(dataset)
    avg_recon = total_recon_loss / num_samples
    avg_kl = total_kl_loss / num_samples
    avg_mi = total_mi_penalty / len(train_loader)
    avg_aux = total_aux_loss / len(train_loader)
    
    print_now(f"Epoch {epoch+1}/{args.epochs} | Recon: {avg_recon:.6f} | KL: {avg_kl:.6f} | MI Pen: {avg_mi:.4f} | Aux Loss: {avg_aux:.4f}")

print_now("\nTraining complete.")

# --- Generate and Save Fair Reconstructions ---
print_now(f"Generating fair reconstructions and saving to '{args.output_file}'...")
vfae_model.eval()
with torch.no_grad():
    recon_all = []
    for batch_val in range(num_batches):
        # Batch values here starts from 0, since we used pd.factorize
        batch_indices = torch.full(
            (len(features_tensor),),
            batch_val,
            dtype=torch.long,
            device=device
        )

        # Pass the raw indices to the model. The model will perform the one-hot encoding.
        x_recon, _, _, _ = vfae_model(features_tensor.to(device), batch_indices)

        recon_all.append(x_recon.cpu().numpy())
    
    recon_all_cpu = np.mean(recon_all, axis=0)

# De-normalize the data back to its original scale
recon_denormalized = (recon_all_cpu + 1) / 2 * (features_max.values - features_min.values) + features_min.values

# Create a new DataFrame for the reconstructed data
recon_df = pd.DataFrame(recon_denormalized, columns=feature_cols, index=df.index)

# Combine with metadata if it exists
meta_df = df[[col for col in df.columns if col.startswith('meta_') or col == args.batch_col]]
output_df = pd.concat([meta_df, recon_df], axis=1)

try:
    output_df.to_csv(args.output_file, index=False)
    print_now("Successfully saved reconstructed data.")
except Exception as e:
    print_now(f"Error saving the output file: {e}", file=sys.stderr)
    sys.exit(1)

