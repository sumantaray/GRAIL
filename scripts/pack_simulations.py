"""Pack the seven fresh Splatter simulations into cell-by-gene H5ADs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import anndata as ad
import pandas as pd
from scipy.io import mmread


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("--out", type=Path, default=Path("datasets/simulations"))
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for i in range(1, 8):
        directory = args.source / f"P{i}"
        cells = pd.read_csv(directory/"cells.tsv", sep="\t", index_col="cell")
        genes = (directory/"genes.txt").read_text().splitlines()
        counts = mmread(directory/"matrix.mtx").T.tocsr().astype("int32")
        data = ad.AnnData(counts, obs=cells, var=pd.DataFrame(index=genes))
        data.uns["simulation"] = f"P{i} reconstructed Table 1 design"
        data.uns["splatter_parameters"] = (directory/"parameters.R").read_text()
        path = args.out/f"P{i}.h5ad"
        data.write_h5ad(path, compression="gzip")
        with path.open("rb") as f:
            rows.append(dict(file=path.name, n_cells=data.n_obs, n_genes=data.n_vars,
                             sha256=hashlib.file_digest(f,"sha256").hexdigest()))
    (args.out/"source_files.json").write_text(json.dumps(rows,indent=2)+"\n")
    shutil.copyfile(args.source/"R-session.txt", args.out/"R-session.txt")


if __name__ == "__main__":
    main()
