# Inputs, outputs and table mapping

Each seed directory contains its configuration and dependencies (`config.json`, `run.json`), preprocessed feature definitions (`features.json`), training/test/anchor IDs, pairings and evaluation cell IDs (`cell_ids.json`), model weights (`model.pt`), nearest paired codes (`pairs.npz`) and generated profiles (`generated.h5ad`). The row identifier and `source_anchor` column in generated data permit cell-by-cell, time-matched comparisons. Expression values are continuous model outputs. Supply the original input count file to rerun preprocessing.

| Manuscript analysis | Computed output for each seed |
| --- | --- |
| Table 2, centroid path geometry | `geometry.csv` with TL, chord and AC over the common four-point grid |
| Table 3, direct linear DC | `python scripts/directional_check.py` prints both fixed-weight directions, the excluded flat gene and the varying-weight example |
| Figure 3, endpoint classifier | `ccm.csv` with mean predicted Stage B probability by coordinate |
| Table 4, held-out population WD and marker trends | `wd.csv` and `mgtc.csv`, including marker counts |
| Table 5, distinguishability | `population.csv` with five-fold logistic and random-forest AUROC where sample size permits |
| Table 6, population agreement | `population.csv` with HVG correlation, diversity ratio, covariance correlation and optional cell-state JS divergence |
| Table 7, condition sensitivity and source-end agreement | `sensitivity.csv` and `endpoint_wd.csv` with fixed-A, best other intermediate time and real-B controls |
| Table 8, component ablations | Compare nearest, direct expression interpolation and AE-decoder output in `lsh`. Compare `lsh` with `random` using nearest pairing to isolate anchor selection |
| S1, pairing sensitivity | `wd.csv` and `mgtc.csv` for nearest, random and centroid using shared fitted weights |
| S2, branch composition | `branch.csv` for a dataset with annotated terminal branches |
| S3, endpoint reconstruction | `reconstruction.csv` for the training mean, 64-component PCA and autoencoder on held-out A and B |
| S4, encoded intermediate agreement | `latent.csv` with standardized sliced W1 for matched, A, B and closest other intermediate coordinate |
| Figures 4 and 6, joint panels | `scripts/plot_joint.py` writes `joint_umap.png`, plotted cell coordinates and DPT summaries |

Metrics requiring observed intermediate cells are emitted only for data with those cells. The complete EB source lacks cell-type annotations, so cell-state composition requires a separately verified annotation file. Indistinguishability and specific developmental identity are different questions. A low classifier AUROC or a small population distance alone cannot establish the identity of an individual cell trajectory.

The plots jointly fit PCA, a neighbor graph, and a single UMAP on the explicitly selected cells, then reuse their coordinates, limits and aspect ratio for panels B and C. DPT is rooted at the observed Stage A cell nearest the Stage A centroid in the plotting PCA space. Unreachable cells are saved and excluded from the median/rank summaries. Quantitative WD, reconstruction and latent diagnostics are computed independently of the plotting embedding.
