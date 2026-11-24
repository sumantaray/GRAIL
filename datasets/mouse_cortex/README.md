# Mouse embryonic cortex input

The revised manuscript names four sampling stages from the TrajectoryNet mouse cortex experiment. The source matrix and cell-to-stage assignments used for the manuscript have not been identified in the cited distribution. Figure 4 is available in the [figure gallery](../../docs/results.md). The [dataset guide](../README.md#mouse-embryonic-cortex) links to the inspected upstream repository and its supplementary methods.

| Source stage | Interpolation coordinate | Use |
| --- | ---: | --- |
| E12.5 | 0.00 | Source endpoint |
| E14.5 | 0.40 | Held-out evaluation |
| E16 | 0.70 | Held-out evaluation |
| E17.5 | 1.00 | Target endpoint |

The reference count workflow expects a cells-by-genes H5AD with integer counts in `X`, unique source cell and gene identifiers, and numeric coordinates in `obs["t"]`. Record the accession or original filename, cell selection and processing history alongside the deposited matrix. If the recovered matrix contains transformed expression only, its preprocessing needs to be matched before running the count workflow.
