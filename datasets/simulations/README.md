# Splatter P1–P7

The seven H5AD files in this directory contain counts for the model-ready endpoints. They were generated with Splatter 1.34.0 and R 4.5.2. Each file records gene identifiers, endpoint labels, simulation parameters and batch. The supplied generator recreates the designs described in Table 1. These are newly generated simulations. The historical simulation matrices and their seeds were unavailable. Starting with simulation seed 42, scenario P1 uses 42, P2 uses 43, and so on. Model training uses an independent seed list.

| Design | Distinguishing setting | Selected A / B cells |
| --- | --- | ---: |
| P1 | Linear ancestry, group probabilities 0.30, 0.30, 0.40 | 300 / 400 |
| P2 | Group probabilities 0.10, 0.10, 0.80 | 100 / 800 |
| P3 | DE probability 0.1, factor location 0.2, factor scale 0.4 | 300 / 400 |
| P4 | DE probability 0.7, factor location 1.0, factor scale 1.2 | 300 / 400 |
| P5 | Dropout midpoint 2.0 | 300 / 400 |
| P6 | Common progenitor and two descendants with probabilities 0.35, 0.35 | 300 / 350 |
| P7 | Three batches with 500, 200, 300 cells and batch factor location 1.5 | 500 / 300 |

Unreported details are explicit choices in this reconstruction. We simulate 3,000 genes, use experiment-wide dropout with midpoint 0 outside P5, and use Splatter's other defaults. For P1–P6 we simulate a pool of 1,400 cells so we can draw the stated endpoint sizes, with Path 1 as A and Path 3 as B for linear scenarios. P6 uses both descendant paths in B. Its realized proportions depend on sampling. P7 uses the first and last batch as the endpoints, giving a deliberately batch-confounded comparison. Only the two designated endpoint populations are exported. The four generated coordinates do not represent measured intermediate populations.

From the repository root, with R, Bioconductor **splatter**, and **Matrix** installed:

```bash
Rscript scripts/simulate.R runs/splatter_input 42
python scripts/pack_simulations.py runs/splatter_input --out runs/rebuilt_simulations
grail run --data runs/rebuilt_simulations/P6.h5ad --config configs/reference.json --out runs/p6
python scripts/summarize.py runs/p6
```

The branch classifier fits annotated terminal training cells. The output `branch.csv` records actual reference and generated counts and totals per branch. Branch proportions follow from the corresponding count divided by its total. Each H5AD retains the Splatter parameters in `uns` and the realized path and batch in `obs`.
