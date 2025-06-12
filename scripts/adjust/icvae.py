import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

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
    This is a critical change from the original VFAE implementation.
    """
    def __init__(self, input_dim, latent_dim, hidden_dim=512):
        super(Encoder, self).__init__()
        # The input layer now only depends on `input_dim`, not `input_dim + s_dim`.
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

    def forward(self, x):
        # The forward pass no longer accepts or uses `s_one_hot`.
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
        # The input layer depends on the concatenated dimensions of z and s.
        self.fc1 = nn.Linear(latent_dim + s_dim, hidden_dim)
        self.fc_out = nn.Linear(hidden_dim, output_dim)

    def forward(self, z, s_one_hot):
        # Concatenate z and s_one_hot to form the input.
        combined_input = torch.cat([z, s_one_hot], dim=1)
        h = F.relu(self.fc1(combined_input))
        # The final sigmoid activation is removed to be compatible with MSE loss.
        x_recon = self.fc_out(h)
        return x_recon


class AuxiliaryClassifier(nn.Module):
    """
    NEW Auxiliary Classifier module: q(s|z)
    This network is introduced to estimate the mutual information I(z;s).
    It takes a latent code `z` and tries to predict the sensitive attribute `s`.
    """
    def __init__(self, latent_dim, s_dim, hidden_dim=128):
        super(AuxiliaryClassifier, self).__init__()
        self.fc1 = nn.Linear(latent_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, s_dim) # Outputs logits for each category of s

    def forward(self, z):
        h = F.relu(self.fc1(z))
        logits_s = self.fc2(h)
        return logits_s



def loss_function(recon_x, x, mean, log_var, class_probs, y_true, loss_fn='cross_entropy'):
    # Loss is made of several components:
    # 1. Reconstruction loss, causes the model to learn to reconstruct the input data
    # 2. KL divergence, causes the model to learn a latent space that is close to a standard normal distribution
    # 3. Classification loss, which is negated through the encoder portion. This c

    reconstruction_loss = torch.nn.functional.binary_cross_entropy(recon_x, x, reduction='sum')
    kl_divergence = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
    if loss_fn == 'cross_entropy':
        classification_loss = torch.nn.functional.binary_cross_entropy_with_logits(class_probs, y_true, reduction='sum')
    elif loss_fn == 'dmi':
        # -log determinant of the joint probability distribution
        classification_loss = -torch.logdet(torch.mm(class_probs.t(), class_probs) + 1e-10)



class VFAE(nn.Module):
    """
    Main VFAE model: Encoder and Decoder.
    """
    def __init__(self, input_dim, num_s_categories, latent_dim, hidden_dim_enc=512, hidden_dim_dec=512):
        super(VFAE, self).__init__()
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


