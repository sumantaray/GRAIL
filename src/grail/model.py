"""Networks and training loop corresponding to Algorithm 1."""
import numpy as np
import torch
from sklearn.neighbors import NearestNeighbors
from torch import nn
from torch.nn import functional as F


def dense_stack(widths, *, linear_last=True):
    layers = []
    for i, (a, b) in enumerate(zip(widths, widths[1:])):
        layers.append(nn.Linear(a, b))
        if i < len(widths) - 2 or not linear_last:
            layers.extend([nn.BatchNorm1d(b), nn.ReLU()])
    return nn.Sequential(*layers)


class Autoencoder(nn.Module):
    def __init__(self, p):
        super().__init__()
        self.encoder = dense_stack([p, 1024, 512, 128, 64], linear_last=False)
        self.decoder = dense_stack([64, 128, 512, 1024, p])

    def forward(self, x):
        return self.decoder(self.encoder(x))


class Generator(nn.Module):
    def __init__(self, p):
        super().__init__()
        self.net = dense_stack([128, 128, 256, 512, p])

    def forward(self, code, noise):
        return self.net(torch.cat([code, noise], dim=1))


def discriminator(p):
    layers = []
    for a, b in zip([p, 512, 256], [512, 256, 64]):
        layers.extend([nn.Linear(a, b), nn.LeakyReLU(0.2)])
    return nn.Sequential(*layers, nn.Linear(64, 1))


def train_models(x_a, x_b_anchors, x_b, cfg, device):
    x_a, x_b_anchors, x_b = [torch.as_tensor(x, device=device) for x in (x_a, x_b_anchors, x_b)]
    p, batch = x_a.shape[1], cfg["batch_size"]
    if batch < 2 or min(cfg["ae_epochs"], cfg["gan_epochs"], cfg["gan_steps_per_epoch"]) < 1:
        raise ValueError("A batch of at least two and positive training durations are required.")
    ae = Autoencoder(p).to(device)
    pooled = torch.cat([x_a, x_b_anchors])
    optimizer = torch.optim.Adam(ae.parameters(), lr=cfg["ae_lr"])
    for _ in range(cfg["ae_epochs"]):
        order = torch.randperm(len(pooled), device=device)
        for group in order.split(batch):
            if len(group) == 1:
                group = torch.cat((group, order[:1]))
            clean = pooled[group]
            optimizer.zero_grad()
            F.mse_loss(ae(clean), clean).backward()
            optimizer.step()
    ae.eval().requires_grad_(False)
    with torch.no_grad():
        z_a = ae.encoder(x_a)
    gen = Generator(p).to(device)
    disc = discriminator(p).to(device)
    adam = dict(lr=cfg["gan_lr"], betas=tuple(cfg["gan_betas"]))
    opt_g, opt_d = torch.optim.Adam(gen.parameters(), **adam), torch.optim.Adam(disc.parameters(), **adam)
    for _ in range(cfg["gan_epochs"] * cfg["gan_steps_per_epoch"]):
        z = z_a[torch.randint(len(z_a), (batch,), device=device)]
        real = x_b[torch.randint(len(x_b), (batch,), device=device)]
        fake = gen(z, torch.randn_like(z))
        real_score = disc(real)
        fake_score = disc(fake.detach())
        loss_d = F.binary_cross_entropy_with_logits(real_score, torch.full_like(real_score, cfg["real_label"]))
        loss_d += F.binary_cross_entropy_with_logits(fake_score, torch.zeros_like(fake_score))
        alpha = torch.rand((batch, 1), device=device)
        between = (alpha * real + (1 - alpha) * fake.detach()).requires_grad_()
        grad = torch.autograd.grad(disc(between).sum(), between, create_graph=True)[0]
        penalty = (grad.norm(2, dim=1) - 1).square().mean()
        opt_d.zero_grad()
        (loss_d + cfg["gradient_penalty"] * penalty).backward()
        opt_d.step()
        disc.requires_grad_(False)
        opt_g.zero_grad()
        # log(1-sigmoid(logit)), the minimax update stated in Algorithm 1.
        (-F.softplus(disc(fake)).mean()).backward()
        opt_g.step()
        disc.requires_grad_(True)
    return ae.cpu().eval(), gen.cpu().eval().requires_grad_(False)


@torch.inference_mode()
def encode(module, X):
    return np.concatenate([module(chunk).numpy() for chunk in torch.from_numpy(X).split(512)])


def pair_codes(za, zb, mode, rng):
    if mode == "centroid":
        return np.repeat(zb.mean(axis=0, keepdims=True), len(za), axis=0), None
    if mode == "nearest":
        indices = NearestNeighbors(n_neighbors=1).fit(zb).kneighbors(za, return_distance=False)[:, 0]
    elif mode == "random":
        indices = rng.integers(len(zb), size=len(za))
    else:
        raise ValueError(f"Unknown pairing: {mode}")
    return zb[indices], indices


@torch.inference_mode()
def generate(gen, codes):
    values = []
    for chunk in torch.from_numpy(np.asarray(codes, dtype=np.float32)).split(512):
        values.append(gen(chunk, torch.zeros_like(chunk)).numpy())
    return np.concatenate(values)
