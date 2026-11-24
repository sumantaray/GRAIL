# Datasets

## Human embryoid body

The [embryoid_body](embryoid_body) directory contains the complete **31,161 cells × 33,694 gene features** from Kevin Moon, *Embryoid Body data for PHATE*, Mendeley Data v1 (2018). Each sparse H5AD stores original integer counts in `X`, Ensembl IDs in `var_names`, symbols in `var["symbol"]` and batch-prefixed cell IDs in `obs_names`. Files are partitioned at 2,000 cells for download. The [manifest](embryoid_body/manifest.json) lists every part and checksum.

| Archive batch | Collection window | Manuscript day label | Coordinate | Cells |
| --- | --- | ---: | ---: | ---: |
| T0_1A | Day 0–3 | 0 | 0.00 | 4,649 |
| T2_3B | Day 6–9 | 6 | 0.25 | 7,377 |
| T4_5C | Day 12–15 | 12 | 0.50 | 6,245 |
| T6_7D | Day 18–21 | 18 | 0.75 | 6,563 |
| T8_9E | Day 24–27 | 24 | 1.00 | 6,327 |

Days in the manuscript refer to the start of a collection window. The first and last windows supply training endpoints. The three interior windows supply evaluation data. The [480-cell example](example/eb_480.h5ad) comprises 96 observed cells per window, uniformly sampled without replacement with seed 42 before any model fitting.

**Attribution and changes.** Kevin Moon, Yale University. The data are distributed under CC BY 4.0. We converted the source 10x Matrix Market files to sparse AnnData, preserved row order and counts, added collection-window and coordinate metadata, and supplied a separate subsample. [Original dataset](https://doi.org/10.17632/v6n743h5ng.1) · [Collection-window mapping](https://krishnaswamylab.github.io/tutorial/load_filter_transform/) · [License](https://creativecommons.org/licenses/by/4.0/)

Rebuild in a fresh directory with `python scripts/prepare_embryoid_body.py --out rebuilt/embryoid_body --demo rebuilt/eb_480.h5ad`. Pass `--archive /path/to/scRNAseq.zip` to reuse a saved download. The script verifies the upstream archive SHA-256 `f0f07b7b5b5bc04e61fbe0e553c1e202e940dd5332cb7ac25f5a27df0170db92`. HDF5 byte hashes may depend on the writer version. The checked-in manifest verifies the distributed files.

## Mouse embryonic cortex

The GRAIL manuscript identifies the TrajectoryNet dataset with E12.5 and E17.5 as endpoints, and E14.5 and E16 held out. The corresponding coordinates are 0, 0.40, 0.70 and 1. The cited upstream repository at the inspected revision contains `eb_genes.txt` and `eb_velocity_v5.npz` in its data directory. Both files refer to embryoid body. The TrajectoryNet supplement, section E.2, stated that the cortex data were unavailable when published. The exact GRAIL cortex count matrix, cell identifiers and stage annotations are pending identification. [Upstream data directory](https://github.com/KrishnaswamyLab/TrajectoryNet/tree/810c89b081f95405bc0ad42bfb3cb41038e13700/data) · [TrajectoryNet supplementary methods](https://proceedings.mlr.press/v119/tong20a/tong20a-supp.pdf)

To analyze a verified cortex matrix, use a raw-count H5AD with unique cells and genes and numeric `obs["t"]` as above. Run the same command as for embryoid body with `--data path/to/cortex.h5ad`. Stage labels and expression values on a UMAP figure cannot supply a count matrix.

## Simulated endpoints

[Simulation designs, parameters and commands](simulations/README.md).
