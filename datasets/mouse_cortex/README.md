# Mouse embryonic cortex input

The revised manuscript names four sampling stages from the TrajectoryNet mouse cortex experiment. The source matrix and cell-to-stage assignments used for the manuscript have not been identified in the cited distribution. Figure 4 is available in the [figure gallery](../../docs/results.md). The [dataset guide](../README.md#mouse-embryonic-cortex) links to the inspected upstream repository and its supplementary methods.

| Source stage | Interpolation coordinate | Use |
| --- | ---: | --- |
| E12.5 | 0.00 | Source endpoint |
| E14.5 | 0.40 | Held-out evaluation |
| E16 | 0.70 | Held-out evaluation |
| E17.5 | 1.00 | Target endpoint |

## Upstream source check

The TrajectoryNet supplement reports roughly 20,000 filtered mouse-cortex cells and states that the cortex dataset was unavailable when the paper appeared. [Supplementary Methods, section E.2](https://proceedings.mlr.press/v119/tong20a/tong20a-supp.pdf)

On 24 November 2026, the upstream `data/` directory contained `eb_velocity_v5.npz` and `eb_genes.txt` on all four branch heads (`master`, `anndata-integration`, `l0-reg`, and `v2`) and in tags `0.2.1`, `0.2.3`, and `0.2.4`. The initial upload contained the embryoid-body NPZ alone. The Figure 8 checkpoint was added with the embryoid-body tutorial. An accession or original expression matrix with matching cell-to-stage assignments is still needed to verify the GRAIL cortex rows. [Initial data directory](https://github.com/KrishnaswamyLab/TrajectoryNet/tree/7f6b11559045828b2ad98cae749bdb029d8d578f/data) · [Inspected master data directory](https://github.com/KrishnaswamyLab/TrajectoryNet/tree/810c89b081f95405bc0ad42bfb3cb41038e13700/data) · [Embryoid-body tutorial commit](https://github.com/KrishnaswamyLab/TrajectoryNet/commit/007f3131c3b30ba536a2aa5d4d66dbb78c0b9baf)

## Input format

The reference count workflow expects a cells-by-genes H5AD with integer counts in `X`, unique source cell and gene identifiers, and numeric coordinates in `obs["t"]`. Record the accession or original filename, cell selection and processing history alongside the deposited matrix. If the recovered matrix contains transformed expression only, its preprocessing needs to be matched before running the count workflow.
