"""Draw four views of one joint embedding from a saved evaluation run."""
import argparse
from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--cap-per-stage", type=int, default=1200)
    args = p.parse_args()
    real = ad.read_h5ad(args.run/"observed.h5ad")
    gen = ad.read_h5ad(args.run/"generated.h5ad")
    if not real.var_names.equals(gen.var_names):
        raise ValueError("Inputs must share ordered genes.")
    rng = np.random.default_rng(args.seed)
    pools = []
    for source, mat in [("Real", real), ("Generated", gen)]:
        keep = []
        for s in sorted(set(mat.obs["t"])):
            ids = np.flatnonzero(mat.obs["t"].to_numpy() == s)
            keep.extend(rng.choice(ids, min(len(ids), args.cap_per_stage), replace=False))
        sub = mat[sorted(keep)].copy()
        sub.obs["source"] = source
        pools.append(sub)
    merged = ad.concat(pools, index_unique=":")
    sc.pp.pca(merged, n_comps=min(30,merged.n_vars-1,merged.n_obs-1), random_state=args.seed)
    sc.pp.neighbors(merged, n_neighbors=min(15, merged.n_obs-1), use_rep="X_pca", random_state=args.seed)
    sc.tl.umap(merged, random_state=args.seed)
    sc.tl.diffmap(merged, random_state=args.seed)
    stage = merged.obs["t"].to_numpy()
    is_real = (merged.obs["source"] == "Real").to_numpy()
    early = np.flatnonzero(is_real & (stage == 0))
    pts = merged.obsm["X_pca"][early]
    merged.uns["iroot"] = int(early[np.argmin(np.linalg.norm(pts - pts.mean(0), axis=1))])
    sc.tl.dpt(merged)
    uv = merged.obsm["X_umap"]
    fig, axs = plt.subplots(2, 2, figsize=(9,8), layout="constrained", sharex=True, sharey=True)
    a,b,c,d = axs.ravel()
    for s in sorted(set(stage[is_real])):
        m = is_real & (stage == s)
        a.scatter(*uv[m].T, s=3, alpha=.65, label=f"t={s:g}")
    a.legend(fontsize=9, markerscale=2)
    for ax, mask, color, label in [(b,~is_real,"#df7d22","Generated"),(c,is_real,"#1767a0","Real")]:
        ax.scatter(*uv[~mask].T, c="#c7c7c7", s=3, alpha=.45)
        ax.scatter(*uv[mask].T, c=color, s=3, alpha=.65, label=f"{label} (n={mask.sum()})")
        ax.legend(fontsize=10, markerscale=2)
    dpt = merged.obs["dpt_pseudotime"].to_numpy()
    finite = np.isfinite(dpt)
    d.scatter(*uv[~finite].T, s=3, color="lightgray")
    points = d.scatter(*uv[finite].T, c=dpt[finite], s=3, cmap="viridis")
    fig.colorbar(points, ax=d, label="DPT")
    for label, ax in zip(["A   Sampling stage","B   Generated highlighted","C   Real highlighted","D   Diffusion pseudotime"],axs.ravel()):
        ax.set(title=label, xlabel="UMAP 1", ylabel="UMAP 2", aspect="equal")
    fig.savefig(args.run/"joint_umap.png",dpi=200)
    coords = merged.obs[["t","source"]].copy()
    coords[["UMAP1","UMAP2"]] = uv
    coords["DPT"] = dpt
    coords.to_csv(args.run/"plot_cells.csv.gz")
    medians = coords.loc[finite].groupby(["source","t"], observed=True).DPT.median().reset_index()
    medians.to_csv(args.run/"dpt_medians.csv", index=False)
    summary = []
    for source, g in medians.groupby("source", observed=True):
        rho = spearmanr(g.t,g.DPT).statistic if len(g)>1 and g.DPT.nunique()>1 else np.nan
        summary.append({"source":source,"stage_order_spearman":rho})
    pd.DataFrame(summary).to_csv(args.run/"dpt_order.csv",index=False)
    print(f"DPT root {merged.obs_names[merged.uns['iroot']]} • {finite.sum()}/{len(finite)} cells reachable")


if __name__ == "__main__":
    main()
