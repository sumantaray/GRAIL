"""Recompute the seed-42 embryoid-body controls with the released generator output.

Run grail generate on results/embryoid_body/seed-42-generator at all five
coordinates first. The same archived cells and pairs enter both distance spaces.
"""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits

from grail.data import feature_space, read_counts
from grail.metrics import wasserstein


ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "datasets" / "embryoid_body"
RESULTS = ROOT / "results" / "embryoid_body"
CHECKPOINT = RESULTS / "seed-42-generator"
TIMES = (0.0, 0.25, 0.5, 0.75, 1.0)
INTERIOR = TIMES[1:-1]
WD_COLUMNS = ["wd_matched", "wd_fixed_a", "wd_other", "wd_real_b"]


def positions(labels, wanted):
    labels = pd.Index(labels)
    if not labels.is_unique or not pd.Index(wanted).is_unique:
        raise ValueError("Expected unique cell or anchor identifiers.")
    found = labels.get_indexer(wanted)
    if (found < 0).any():
        raise ValueError("Missing archived identifier: " + str(wanted[np.flatnonzero(found < 0)[0]]))
    return found


def samples(generated_path, cfg, seed):
    saved = json.loads((CHECKPOINT / "cell_ids.json").read_text())
    data = read_counts(COUNTS)
    x, train, test, qc, features = feature_space(data, cfg, seed)
    if features != json.loads((CHECKPOINT / "features.json").read_text()):
        raise ValueError("Rebuilt feature space differs from the checkpoint.")
    np.testing.assert_array_equal(data.obs_names[train], saved["train"])
    np.testing.assert_array_equal(data.obs_names[test], saved["test"])

    gen = ad.read_h5ad(generated_path)
    np.testing.assert_array_equal(gen.var_names, features["input_genes"])
    np.testing.assert_array_equal(np.unique(gen.obs["t"]), TIMES)
    times = data.obs["t"].to_numpy()
    observed, generated = {}, {}
    for t in TIMES:
        indices = positions(data.obs_names, saved["reference_ids"][str(t)])
        if not (qc[indices] & (times[indices] == t)).all():
            raise ValueError(f"Observed sample disagrees with the archived stage {t}.")
        if t in (0.0, 1.0) and not np.isin(indices, test).all():
            raise ValueError(f"Endpoint sample {t} includes training cells.")
        observed[t] = x[indices]
        block = gen[gen.obs["t"].to_numpy() == t]
        rows = positions(block.obs["source_anchor"], saved["paired_sample_ids"])
        values = block.X[rows]
        generated[t] = values.toarray() if sparse.issparse(values) else np.asarray(values)

    pca = PCA(n_components=10, random_state=seed).fit(x[train])
    return observed, generated, pca


def stage_distances(observed, generated, name, rms):
    records = []
    for t in INTERIOR:
        target = observed[t]
        best_wd, best_t = min(
            (wasserstein(generated[u], target), u) for u in INTERIOR if u != t
        )
        records.append(dict(
            space=name, t=t, n=len(target), paired_rms=rms[t],
            wd_matched=wasserstein(generated[t], target),
            wd_fixed_a=wasserstein(generated[0.0], target),
            wd_other=best_wd, best_other_t=best_t,
            wd_real_b=wasserstein(observed[1.0], target),
        ))
    return pd.DataFrame(records)


def archived(filename, index):
    table = pd.read_csv(RESULTS / filename)
    selection = table.loc[(table["seed"] == 42) & (table["anchor_method"] == "lsh")]
    return selection.set_index(index).sort_index()


def replay(generated_path, cfg, seed):
    real, fake, pca = samples(generated_path, cfg, seed)
    rms = {t: np.sqrt(np.square(fake[t] - fake[0.0]).mean(axis=1)).mean() for t in INTERIOR}
    pc_real = {t: pca.transform(x) for t, x in real.items()}
    pc_fake = {t: pca.transform(x) for t, x in fake.items()}
    pc = stage_distances(pc_real, pc_fake, "pca10_unwhitened", rms)
    measured = pc.set_index("t")
    expected = archived("sensitivity.csv", "t")
    columns = ["paired_rms", *WD_COLUMNS]
    np.testing.assert_array_equal(measured.index, expected.index)
    np.testing.assert_array_equal(measured["n"], expected["n"])
    np.testing.assert_allclose(measured[columns], expected[columns], rtol=1e-4, atol=1e-6)

    endpoints = archived("endpoint_wd.csv", "reference_t")
    np.testing.assert_array_equal(endpoints.index, [0.0, 1.0])
    endpoint_wd = [wasserstein(pc_fake[0.0], pc_real[t]) for t in endpoints.index]
    np.testing.assert_allclose(endpoint_wd, endpoints["wd_generated_t0"], rtol=1e-4, atol=1e-6)

    full = stage_distances(real, fake, "log_expression_2000_genes", rms)
    return pd.concat([pc, full], ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("generated", type=Path, help="H5AD replayed from the seed-42 LSH checkpoint")
    parser.add_argument("--out", type=Path, required=True, help="New CSV file for the six comparisons")
    args = parser.parse_args()
    cfg = json.loads((CHECKPOINT / "config.json").read_text())
    run = json.loads((CHECKPOINT / "run.json").read_text())
    with (COUNTS / "manifest.json").open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != run["input_sha256"]:
            raise ValueError("The count manifest differs from the fitted run.")
    if run["seed"] != 42 or run["anchor_method"] != "lsh":
        raise ValueError("Expected the released LSH fit with seed 42.")
    with threadpool_limits(limits=cfg["threads"]):
        result = replay(args.generated, cfg, run["seed"])
    result.to_csv(args.out, mode="x", index=False, float_format="%.12g")
    print(f"Reproduced the archived 10-PC controls and endpoint distances. Wrote {args.out}")


if __name__ == "__main__":
    main()
