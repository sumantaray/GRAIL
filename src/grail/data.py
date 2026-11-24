"""Read count matrices and fit the expression feature space to training endpoints."""
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.model_selection import train_test_split


def read_counts(path):
    path = Path(path)
    if path.is_dir():
        manifest = json.loads((path / "manifest.json").read_text())
        panels = []
        for record in manifest["parts"]:
            file = path / record["file"]
            with file.open("rb") as stream:
                if hashlib.file_digest(stream, "sha256").hexdigest() != record["sha256"]:
                    raise ValueError(f"Checksum mismatch: {file}")
            panel = ad.read_h5ad(file)
            if panels and not panels[0].var_names.equals(panel.var_names):
                raise ValueError(f"Gene order differs in {file}")
            panels.append(panel)
        data = ad.concat(panels, merge="same", uns_merge="same")
    else:
        data = ad.read_h5ad(path)
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Cell IDs and gene IDs must be unique.")
    times = pd.to_numeric(data.obs["t"], errors="raise").to_numpy(dtype=float)
    if not np.isfinite(times).all() or (times < 0).any() or (times > 1).any():
        raise ValueError("obs['t'] must contain finite normalized stage coordinates in [0, 1].")
    data.obs["t"] = times
    x = sparse.csr_matrix(data.X, dtype=np.float32)
    x.sum_duplicates()
    x.eliminate_zeros()
    if not np.isfinite(x.data).all() or (x.data < 0).any() or not np.equal(x.data, np.floor(x.data)).all():
        raise ValueError("Supply nonnegative raw integer counts in X.")
    data.X = x
    return data


def feature_space(data, config, seed):
    t = data.obs["t"].to_numpy()
    endpoints = np.flatnonzero((t == 0) | (t == 1))
    if len(np.unique(t[endpoints])) != 2:
        raise ValueError("Both endpoint stages, t=0 and t=1, are required.")
    train, test = train_test_split(endpoints, test_size=config["test_fraction"], stratify=t[endpoints],
                                   random_state=seed)
    passes_qc = np.asarray(data.X.getnnz(axis=1)) >= config["min_genes"]
    train, test = train[passes_qc[train]], test[passes_qc[test]]
    if any(np.sum(t[idx] == s) < 2 for idx in (train, test) for s in (0, 1)):
        raise ValueError("Too few endpoint cells remain in a split after QC.")
    gene_mask = data.X[train].getnnz(axis=0) >= config["min_cells"]
    if gene_mask.sum() < 2:
        raise ValueError("Too few genes pass the training-only expression filter.")
    first = data[:, gene_mask].copy()
    totals = np.asarray(first.X.sum(axis=1)).ravel()
    scale = config["target_sum"] / np.maximum(totals, 1)
    first.X = first.X.multiply(scale[:, None]).tocsr()
    np.log1p(first.X.data, out=first.X.data)
    train_view = first[train].copy()
    if first.n_vars > config["hvg"]:
        # Use the Scanpy Seurat dispersion statistic on normalized/log1p counts.
        sc.pp.highly_variable_genes(train_view, flavor="seurat", n_top_genes=config["hvg"])
        score = np.nan_to_num(train_view.var["dispersions_norm"].to_numpy(), nan=-np.inf)
        chosen = np.sort(np.argsort(-score, kind="stable")[:config["hvg"]])
    else:
        chosen = np.arange(first.n_vars)
    X = first.X[:, chosen].toarray().astype("float32", copy=False)
    specification = dict(normalization_gene_ids=first.var_names.tolist(), target_sum=config["target_sum"],
                         input_genes=first.var_names[chosen].tolist(), min_genes=config["min_genes"])
    return X, train, test, passes_qc, specification


def select_anchors(X, count, rng, method="lsh", tables=20, bits=10):
    count = min(count, len(X))
    if count < 2:
        raise ValueError("At least two anchors per endpoint are required.")
    if method == "random" or count == len(X):
        return rng.choice(len(X), count, replace=False)
    if method != "lsh" or not 1 <= bits <= 30 or tables < 1:
        raise ValueError("Invalid anchor method or hash settings.")
    # One shuffled pass over each table's occupied buckets before recycling them.
    pools = []
    for _ in range(tables):
        proj = rng.standard_normal((X.shape[1], bits), dtype=np.float32)
        codes = ((X @ proj >= 0).astype(np.int64) @ (1 << np.arange(bits)))
        buckets = {}
        for idx in rng.permutation(len(X)):
            buckets.setdefault(codes[idx], []).append(idx)
        order = rng.permutation(list(buckets))
        pools.append([buckets[key] for key in order])
    chosen, used = [], set()
    cursors = [0] * tables
    while len(chosen) < count:
        for j, pool in enumerate(pools):
            bucket = pool[cursors[j] % len(pool)]
            cursors[j] += 1
            while bucket and bucket[-1] in used:
                bucket.pop()
            if bucket:
                idx = bucket.pop()
                used.add(idx)
                chosen.append(idx)
                if len(chosen) == count:
                    break
    return np.asarray(chosen, dtype=int)
