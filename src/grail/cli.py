"""Run and inspect the manuscript-based GRAIL reference implementation."""
import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits

from .audit import evaluate
from .data import feature_space, read_counts, select_anchors
from .model import Generator, encode, generate, pair_codes, train_models


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def profiles(pops, genes, pairs):
    stages = sorted(pops)
    obs = pd.DataFrame([(s, cell) for s in stages for cell in pairs], columns=["t", "source_anchor"])
    obs.index = [f"{s:.12g}:{i}" for s in stages for i in range(len(pairs))]
    result = ad.AnnData(np.vstack([pops[s] for s in stages]), obs=obs, var=pd.DataFrame(index=genes))
    result.uns["expression_scale"] = "model output on log1p normalized-count target scale"
    return result


def load_external(specs, features):
    result = {}
    for spec in specs:
        name, path = spec.split("=", 1)
        if not name or name in result:
            raise ValueError("External method names must be nonempty and unique.")
        ext = ad.read_h5ad(path)
        if not ext.var_names.is_unique or set(features) != set(ext.var_names):
            raise ValueError(f"{name}: gene IDs must match the fitted gene set.")
        ext = ext[:, features]
        times = pd.to_numeric(ext.obs["t"]).to_numpy(dtype=float)
        if not np.isfinite(times).all() or (times < 0).any() or (times > 1).any():
            raise ValueError(f"{name}: invalid interpolation coordinates.")
        X = ext.X.toarray() if hasattr(ext.X, "toarray") else np.asarray(ext.X)
        if not np.isfinite(X).all():
            raise ValueError(f"{name}: non-finite expression value.")
        result[name] = {s: X[times == s] for s in sorted(set(times))}
    return result


def run(args):
    cfg = json.loads(args.config.read_text())
    if args.seed is not None:
        seeds = [args.seed]
    else:
        seeds = cfg["seeds"]
    if not set(cfg["geometry_grid"]).issubset(cfg["grid"]) or len(set(cfg["geometry_grid"])) < 3:
        raise ValueError("geometry_grid must have at least three distinct coordinates in grid.")
    if any(not (0 <= s <= 1) for s in cfg["grid"]) or len(set(seeds)) != len(seeds):
        raise ValueError("Invalid grid or duplicate seeds.")
    if min(cfg["hvg"], cfg["anchors_a"], cfg["anchors_b"], cfg["evaluation_cells"]) < 2:
        raise ValueError("Invalid feature, anchor or evaluation count.")
    if cfg["threads"] < 1:
        raise ValueError("threads must be positive.")
    torch.set_num_threads(cfg["threads"])
    torch.use_deterministic_algorithms(True)
    data = read_counts(args.data)
    t = data.obs["t"].to_numpy()
    stamp_source = args.data if args.data.is_file() else args.data / "manifest.json"
    with stamp_source.open("rb") as f:
        data_hash = hashlib.file_digest(f, "sha256").hexdigest()
    parent = args.out / args.anchors
    parent.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=cfg["threads"]):
        for seed in seeds:
            dest = parent / f"seed-{seed}"
            if dest.exists():
                raise FileExistsError(f"Existing run {dest}. Select a new output directory.")
            with tempfile.TemporaryDirectory(dir=parent) as tmp:
                output = Path(tmp) / "run"
                output.mkdir()
                rng = np.random.default_rng(seed)
                torch.manual_seed(seed)
                X, train, test, qc, features = feature_space(data, cfg, seed)
                ia, ib = [train[t[train] == s] for s in (0, 1)]
                aa = ia[select_anchors(X[ia], cfg["anchors_a"], rng, args.anchors, cfg["lsh_tables"], cfg["lsh_bits"])]
                ab = ib[select_anchors(X[ib], cfg["anchors_b"], rng, args.anchors, cfg["lsh_tables"], cfg["lsh_bits"])]
                ae, gen = train_models(X[aa], X[ab], X[ib], cfg, cfg["device"])
                extra = load_external(args.external, features["input_genes"])
                rows, generated, ids = evaluate(data, X, t, train, test, qc, aa, ab, ae, gen,
                                                cfg, seed, args.anchors, extra)
                for key, records in rows.items():
                    if records:
                        pd.DataFrame(records).assign(seed=seed, anchor_method=args.anchors).to_csv(output/f"{key}.csv", index=False)
                ids.update(train=data.obs_names[train].tolist(), test=data.obs_names[test].tolist(),
                           anchors_a=data.obs_names[aa].tolist(), anchors_b=data.obs_names[ab].tolist())
                dump(output/"cell_ids.json", ids)
                dump(output/"features.json", features)
                dump(output/"config.json", cfg)
                versions = {p: importlib.metadata.version(p) for p in
                            ("grail-scrna", "torch", "numpy", "scipy", "pandas", "scikit-learn", "scanpy", "anndata")}
                dump(output/"run.json", dict(seed=seed, anchor_method=args.anchors, input_sha256=data_hash,
                                             data_name=args.data.name, python=platform.python_version(), versions=versions,
                                             external_inputs=args.external, training_noise="normal(0,I64)",
                                             inference_noise="zeros(64)"))
                torch.save({"ae": ae.state_dict(), "generator": gen.state_dict()}, output/"model.pt")
                za, zb = encode(ae.encoder, X[aa]), encode(ae.encoder, X[ab])
                dest_codes, _ = pair_codes(za, zb, "nearest", rng)
                np.savez_compressed(output/"pairs.npz", za=za, zb=dest_codes)
                profiles(generated, features["input_genes"], ids["anchors_a"]).write_h5ad(output/"generated.h5ad", compression="gzip")
                if args.save_observed:
                    obs = data.obs.iloc[np.flatnonzero(qc)].copy()
                    obs["split"] = "intermediate"
                    obs.loc[data.obs_names[train],"split"] = "train"
                    obs.loc[data.obs_names[test],"split"] = "test"
                    ad.AnnData(X[qc], obs=obs, var=pd.DataFrame(index=features["input_genes"])).write_h5ad(
                        output/"observed.h5ad", compression="gzip")
                os.rename(output, dest)
            print(dest, flush=True)


def regenerate(args):
    ids = json.loads((args.run/"cell_ids.json").read_text())
    genes = json.loads((args.run/"features.json").read_text())["input_genes"]
    model = Generator(len(genes))
    model.load_state_dict(torch.load(args.run/"model.pt", map_location="cpu", weights_only=True)["generator"])
    model.eval()
    if not all(0 <= v <= 1 for v in args.times):
        raise ValueError("All times must lie in [0,1].")
    with np.load(args.run/"pairs.npz") as pairs:
        pops = {s: generate(model, (1 - s)*pairs["za"] + s*pairs["zb"]) for s in args.times}
    if args.out.exists():
        raise FileExistsError(args.out)
    profiles(pops, genes, ids["anchors_a"]).write_h5ad(args.out, compression="gzip")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(required=True)
    fit = sub.add_parser("run", help="train on endpoints, generate profiles, and evaluate")
    fit.add_argument("--data", type=Path, required=True)
    fit.add_argument("--config", type=Path, required=True)
    fit.add_argument("--out", type=Path, required=True)
    fit.add_argument("--seed", type=int)
    fit.add_argument("--anchors", choices=["lsh", "random"], default="lsh")
    fit.add_argument("--external", action="append", default=[], metavar="METHOD=predictions.h5ad")
    fit.add_argument("--save-observed", action="store_true", help="save processed profiles for UMAP/DPT")
    fit.set_defaults(fn=run)
    replay = sub.add_parser("generate", help="use a saved generator and paired endpoint codes")
    replay.add_argument("--run", type=Path, required=True)
    replay.add_argument("--times", type=float, nargs="+", required=True)
    replay.add_argument("--out", type=Path, required=True)
    replay.set_defaults(fn=regenerate)
    args = p.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
