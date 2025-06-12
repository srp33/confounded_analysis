import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from jaxtyping import Float
from torch import Tensor
from beartype import beartype 

# --- Architectural Changes based on Moyer et al. (2018) ---
# The methodology requires three key components:
# 1. An Encoder that maps x to a latent space z: q(z|x). It does NOT see the sensitive attribute 's'.
# 2. A Decoder that reconstructs x from z AND s: p(x|z,s).
# 3. An Auxiliary Classifier that estimates the conditional probability of s given z: q(s|z).
#    This is used to compute the mutual information penalty.

class Encoder(nn.Module):
    """
    Encoder module: q(z|x)
    This version is MODIFIED as per the paper's requirements.
    It takes ONLY the input data `x` and does NOT see the sensitive attribute `s`.
    This is a critical change from the original ICVAE implementation.
    """
    def __init__(self, input_dim, latent_dim, hidden_dim=512):
        super(Encoder, self).__init__()
        # Compared to the original ICVAE implementation, the input layer now only depends on `input_dim`, not `input_dim + s_dim`.
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

    def forward(self, x):
        h = F.relu(self.fc1(x))
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar

class Decoder(nn.Module):
    """
    Decoder module: p(x|z,s)
    It takes the latent code `z` and the sensitive attribute `s` to reconstruct `x`.
    This allows the model to explain away variance in `x` due to `s`,
    encouraging `z` to become independent of `s`.
    """
    def __init__(self, latent_dim, s_dim, output_dim, hidden_dim=512):
        super(Decoder, self).__init__()
        # Same as ICVAE, the input layer depends on the concatenated dimensions of z and s.
        self.fc1 = nn.Linear(latent_dim + s_dim, hidden_dim)
        self.fc_out = nn.Linear(hidden_dim, output_dim)

    def forward(self, z, s_one_hot):
        # Concatenate z and s_one_hot to form the input.
        combined_input = torch.cat([z, s_one_hot], dim=1)
        h = F.relu(self.fc1(combined_input))
        x_recon = self.fc_out(h)
        return x_recon


class Treeish(nn.Module):
    """
    Treeish module: p(s|z)
    This is a simple classifier that mimics a single decision tree.
    """
    def __init__(self, latent_dim, s_dim, hidden_dim, ratio_to_keep):
        super(Treeish, self).__init__()
        # Create a binary mask to simulate random feature selection
        self.mask = nn.Parameter(
            (torch.rand(latent_dim) < ratio_to_keep).float(), requires_grad=False
        )
        self.feature_dim = int(self.mask.sum().item())

        # Parameters for linear interactions (like attention, but simpler)
        self.weight_q = nn.Parameter(torch.randn(self.feature_dim) * 0.01)
        self.bias_q = nn.Parameter(torch.zeros(self.feature_dim))

        self.weight_k = nn.Parameter(torch.randn(self.feature_dim) * 0.01)
        self.bias_k = nn.Parameter(torch.zeros(self.feature_dim))

        self.final_layer = nn.Linear(self.feature_dim, s_dim) 

    @beartype
    def forward(self, z: Float[Tensor, "batch latent_dim"]) -> Float[Tensor, "batch s_dim"]:
        # Apply feature mask
        masked_z: Float[Tensor, "batch feature_dim"] = z[:, self.mask.bool()]

        q: Float[Tensor, "batch feature_dim"] = masked_z * self.weight_q + self.bias_q
        k: Float[Tensor, "batch feature_dim"] = masked_z * self.weight_k + self.bias_k

        # Compute interactions (like attention, but simpler)
        interactions: Float[Tensor, "batch feature_dim feature_dim"] = q.unsqueeze(2) * k.unsqueeze(1)
        per_feature_sum: Float[Tensor, "batch feature_dim"] = torch.sum(interactions, dim=1)

        # Final classification layer
        logits: Float[Tensor, "batch s_dim"] = self.final_layer(per_feature_sum)

        return F.log_softmax(logits, dim=1)



class AuxiliaryTreeClassifier(nn.Module):
    """
    Auxiliary Tree Classifier module: q(s|z)
    This network estimates the mutual information I(z;s).
    It takes a latent code `z` and tries to predict the sensitive attribute `s`.
    """
    def __init__(self, latent_dim, s_dim, hidden_dim=128, num_trees=100, ratio_to_keep=0.3):
        super(AuxiliaryClassifier, self).__init__()
        # Mimic a random forest, though trees can help each other.
        self.trees = nn.ModuleList([
            Treeish(latent_dim, s_dim, hidden_dim, ratio_to_keep) for _ in range(num_trees)
        ])
        self.num_trees = num_trees
        self.s_dim = s_dim
        self.hidden_dim = hidden_dim

    def forward(self, z):
        # Pass through each tree and average the outputs
        tree_outputs = [tree(z) for tree in self.trees]
        return sum(tree_outputs) / self.num_trees


class AuxiliaryClassifier(nn.Module):
    """
    Auxiliary Classifier module: q(s|z)
    This network estimates the mutual information I(z;s).
    It takes a latent code `z` and tries to predict the sensitive attribute `s`.
    """
    def __init__(self, latent_dim, s_dim, hidden_dim=128):
        super(AuxiliaryClassifier, self).__init__()
        self.fc_in = nn.Linear(latent_dim, hidden_dim)
        self.fc_mid = nn.ModuleList([
            nn.Linear(hidden_dim, hidden_dim) for _ in range(2)  # Two hidden layers
        ])
        self.fc_out = nn.Linear(hidden_dim, s_dim)

    def forward(self, z):
        h = F.relu(self.fc_in(z))
        for layer in self.fc_mid:
            h = F.relu(layer(h))
        s_logits = self.fc_out(h)
        return s_logits


def loss_function(recon_x, x, mean, log_var, class_probs, y_true, loss_fn='cross_entropy'):
    # Loss is made of several components:
    # 1. Reconstruction loss, causes the model to learn to reconstruct the input data, MSE
    # 2. KL divergence, causes the model to learn a latent space that is close to a standard normal distribution
    # 3. Classification loss, which is negated through the encoder portion. This c

    reconstruction_loss = torch.nn.functional.mse_loss(recon_x, x, reduction='mean')
    kl_divergence = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
    if loss_fn == 'cross_entropy':
        classification_loss = torch.nn.functional.binary_cross_entropy_with_logits(class_probs, y_true, reduction='sum')
    elif loss_fn == 'dmi':
        # -log determinant of the joint probability distribution
        classification_loss = -torch.logdet(torch.mm(class_probs.t(), class_probs) + 1e-10)



class ICVAE(nn.Module):
    """
    Main ICVAE model: Encoder and Decoder.
    """
    def __init__(self, input_dim, num_s_categories, latent_dim, hidden_dim_enc=512, hidden_dim_dec=512):
        super(ICVAE, self).__init__()
        self.latent_dim = latent_dim
        self.num_s_categories = num_s_categories

        # Instantiate the modified Encoder and the standard Decoder
        self.encoder = Encoder(input_dim, latent_dim, hidden_dim_enc)
        self.decoder = Decoder(latent_dim, num_s_categories, input_dim, hidden_dim_dec)

    def sample(self, mu, logvar):
        """Sample from N(mu, var)"""
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x, batch_indices):
        # Convert batch indices to one-hot for the decoder
        s_one_hot = F.one_hot(batch_indices, num_classes=self.num_s_categories).float()

        # 1. Encode x to get parameters of q(z|x)
        mu, logvar = self.encoder(x)

        # 2. Sample z from q(z|x)
        z = self.sample(mu, logvar)

        # 3. Decode (z, s) to reconstruct x
        x_recon = self.decoder(z, s_one_hot)

        return x_recon, mu, logvar, z


from argparse import ArgumentParser
import pandas as pd
from sklearn.model_selection import train_test_split

def print_now(*args, **kwargs):
    print(*args, flush=True, **kwargs)


if __name__ == "__main__":
    parser = ArgumentParser("Auxiliary Classifier test")
    parser.add_argument("-i", "--input-file", help="Path to input CSV file.", required=True)
    parser.add_argument("-b", "--batch-col", help="Column name for the sensitive attribute (batch).", required=True)
    args = parser.parse_args()

    # FIX: Define the device for PyTorch operations.
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print_now(f"Using device: {device}")

    # Load data
    try:
        df = pd.read_csv(args.input_file)
    except FileNotFoundError:
        print_now(f"Error: Input file not found at {args.input_file}")
        sys.exit(1)

    # --- CORRECTED DATA PREPARATION WORKFLOW ---

    # 1. Identify feature columns and sensitive attribute
    feature_cols = [col for col in df.columns if not col.startswith('meta_') and col != args.batch_col]
    if not feature_cols:
        print_now("Error: No feature columns found. Ensure columns don't all start with 'meta_' or match the batch column.")
        sys.exit(1)
    print_now(f"Found {len(feature_cols)} feature columns.")

    # 2. Factorize the sensitive attribute ONCE on the whole dataset for consistent encoding
    batch_series = df[args.batch_col]
    batch_codes, unique_batches = pd.factorize(batch_series)
    df[args.batch_col] = batch_codes # Replace string labels with integer codes
    num_batches = len(unique_batches)
    print_now(f"Found {num_batches} unique batches: {unique_batches.tolist()}")

    # 3. Split data into training and validation sets BEFORE normalization
    # This prevents data leakage from the validation set into the training process.
    train_df, val_df = train_test_split(df, test_size=0.2, random_state=42, stratify=df[args.batch_col])

    # 4. Calculate normalization statistics ONLY from the training data
    train_features = train_df[feature_cols].copy()
    features_min = train_features.min()
    features_max = train_features.max()
    
    # Avoid division by zero for columns with no variance
    range_ = features_max - features_min
    range_[range_ == 0] = 1

    # 5. Apply the training set's normalization to both training and validation sets
    train_features_normalized = 2 * (train_features - features_min) / range_ - 1
    val_features_normalized = 2 * (val_df[feature_cols].copy() - features_min) / range_ - 1

    # Fill any potential NaN values that might result from normalization
    train_features_normalized = train_features_normalized.fillna(0)
    val_features_normalized = val_features_normalized.fillna(0)

    # 6. Prepare PyTorch Tensors
    train_features_tensor = torch.tensor(train_features_normalized.values, dtype=torch.float32)
    train_batch_tensor = torch.tensor(train_df[args.batch_col].values, dtype=torch.long)
    
    val_features_tensor = torch.tensor(val_features_normalized.values, dtype=torch.float32).to(device)
    val_batch_tensor = torch.tensor(val_df[args.batch_col].values, dtype=torch.long).to(device)

    # Create DataLoader for training
    train_dataset = TensorDataset(train_features_tensor, train_batch_tensor)
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)

    # Initialize model and optimizer
    input_dim = train_features_tensor.shape[1]
    # FIX: Removed unused 'num_trees' and 'ratio_to_keep' arguments.
    aux_classifier = AuxiliaryClassifier(latent_dim=input_dim, s_dim=num_batches, hidden_dim=128).to(device)
    optimizer = optim.Adam(aux_classifier.parameters(), lr=0.001)

    # Training loop
    print_now("\nStarting Auxiliary Classifier training...")
    num_epochs = 100
    for epoch in range(num_epochs):
        aux_classifier.train()
        total_loss = 0
        for x_batch, s_batch in train_loader:
            x_batch = x_batch.to(device)
            s_batch = s_batch.to(device)

            optimizer.zero_grad()
            s_logits_aux = aux_classifier(x_batch)
            loss = F.cross_entropy(s_logits_aux, s_batch)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        avg_loss = total_loss / len(train_loader)
        print_now(f"Epoch {epoch+1}/{num_epochs}, Loss: {avg_loss:.4f}")

    # Evaluate accuracy on the validation set
    print_now("\nEvaluating on validation set...")
    aux_classifier.eval()
    with torch.no_grad():
        s_logits_val = aux_classifier(val_features_tensor)
        _, predicted = torch.max(s_logits_val, 1)
        accuracy = (predicted == val_batch_tensor).float().mean().item()
        print_now(f"Validation Accuracy: {accuracy * 100:.2f}%")

