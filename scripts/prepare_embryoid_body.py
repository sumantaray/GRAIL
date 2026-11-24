"""Repack the verified PHATE count data losslessly as annotated sparse matrices."""
import argparse
import hashlib
import json
import shutil
import tempfile
import urllib.request
from pathlib import Path
from zipfile import ZipFile

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread

URL = ("https://data.mendeley.com/public-files/datasets/v6n743h5ng/files/"
       "b1865840-e8df-4381-8866-b04d57309e1d/file_downloaded")
SHA256 = "f0f07b7b5b5bc04e61fbe0e553c1e202e940dd5332cb7ac25f5a27df0170db92"
BATCHES = (("T0_1A", 0), ("T2_3B", 6), ("T4_5C", 12), ("T6_7D", 18), ("T8_9E", 24))


def digest(path):
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def prepare(archive, dest, demo):
    if digest(archive) != SHA256:
        raise ValueError("The archive does not match Mendeley Data v1.")
    dest.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    samples, records = [], []
    with ZipFile(archive) as source:
        for batch, day in BATCHES:
            prefix = f"scRNAseq/{batch}/"
            with source.open(prefix + "matrix.mtx") as fh:
                counts = mmread(fh).T.tocsr().astype("int32")
            with source.open(prefix + "genes.tsv") as fh:
                genes = pd.read_csv(fh, sep="\t", header=None, names=["gene_id", "symbol"],
                                    keep_default_na=False).set_index("gene_id")
            with source.open(prefix + "barcodes.tsv") as fh:
                barcodes = pd.read_csv(fh, header=None, names=["barcode"])
            assert counts.shape == (len(barcodes), len(genes)) and genes.index.is_unique
            obs = barcodes.set_axis(pd.Index(batch + ":" + barcodes.barcode, name="cell_id"))
            obs["batch"] = batch
            obs["collection_window"] = f"Day {day}-{day + 3}"
            obs["t"] = day / 24
            mat = ad.AnnData(counts, obs=obs, var=genes)
            mat.uns["source_doi"] = "10.17632/v6n743h5ng.1"
            # Smaller standalone H5ADs can be opened separately on a laptop.
            for start in range(0, mat.n_obs, 2000):
                stop = min(start + 2000, mat.n_obs)
                part = dest / f"{batch}_{start + 1:04d}-{stop:04d}.h5ad"
                mat[start:stop].copy().write_h5ad(part, compression="gzip")
                records.append(dict(file=part.name, batch=batch, window=obs.collection_window.iloc[0],
                                    t=day / 24, n_cells=stop - start, n_genes=mat.n_vars,
                                    sha256=digest(part)))
            sample = np.sort(rng.choice(mat.n_obs, 96, replace=False))
            samples.append(mat[sample].copy())
    manifest = dict(source_doi="10.17632/v6n743h5ng.1", archive_sha256=SHA256, parts=records)
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    demo.parent.mkdir(parents=True, exist_ok=True)
    example = ad.concat(samples, merge="same", uns_merge="same")
    example.uns["sampling"] = "96 cells per collection window, uniform without replacement, seed 42"
    example.write_h5ad(demo, compression="gzip")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, help="reuse the original scRNAseq.zip")
    p.add_argument("--out", type=Path, default=Path("datasets/embryoid_body"))
    p.add_argument("--demo", type=Path, default=Path("datasets/example/eb_480.h5ad"))
    args = p.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        archive = args.archive or Path(tmp) / "scRNAseq.zip"
        if args.archive is None:
            with urllib.request.urlopen(URL, timeout=120) as response, archive.open("wb") as target:
                shutil.copyfileobj(response, target)
        prepare(archive, args.out, args.demo)
