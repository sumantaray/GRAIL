# Measured reference results

These results come from fresh runs of the manuscript-based implementation using seeds 42, 123, and 2026. All fitted models used the [reference configuration](../configs/reference.json). Each run performed 100 autoencoder epochs and 3,000 adversarial updates. The embryoid-body experiment includes both LSH and uniform random anchor selection with an independent fit for each seed and selection method. The P6 experiment uses LSH and a newly generated simulation.

## Embryoid-body expression at held-out stages

Each held-out stage contributes 200 observed cells per seed. The four comparisons use the same observed sample and equally sized comparator populations. The nearest-neighbor anchor pairs and generator weights stay fixed across coordinates. The best other intermediate coordinate is selected within each seed before aggregation. WD is empirical Wasserstein-1 distance with Euclidean costs in ten unwhitened principal components fitted on training endpoints. Smaller distances mean closer distributions.

**Table 7A. Agreement with held-out embryoid-body stages.** Entries give the mean ± sample SD across LSH fits with seeds 42, 123, and 2026. Paired RMS change measures expression change across all 2,000 retained genes.

| Measure | Day 6–9 (t = 0.25) | Day 12–15 (t = 0.50) | Day 18–21 (t = 0.75) |
| --- | ---: | ---: | ---: |
| Paired RMS change from generated t = 0 | 0.052 ± 0.014 | 0.107 ± 0.031 | 0.178 ± 0.068 |
| WD to time-matched generated cells | 17.81 ± 4.21 | 16.69 ± 4.10 | 16.50 ± 3.60 |
| WD to generated cells at fixed t = 0 | 18.17 ± 3.73 | 17.29 ± 3.41 | 17.19 ± 3.08 |
| WD to the best other intermediate coordinate | 17.25 ± 4.13 | 16.59 ± 3.78 | 16.53 ± 3.83 |
| WD to held-out real Stage B cells | **7.93 ± 0.12** | **5.85 ± 0.14** | **5.54 ± 0.10** |

Generated expression changed with the interpolation coordinate. The real Stage B sample had the smallest WD in all nine stage-by-seed comparisons. This differs from the manuscript's Table 7 ordering at Days 6 and 12. In the reference runs, all nine matched-coordinate WDs equal the nearest-neighbor WDs for the same stage and seed in the Table 4 output. [Per-seed Table 7 controls](embryoid_body/sensitivity.csv) · [Table 4 distance output](embryoid_body/wd.csv) · [Three-seed summaries](embryoid_body/summary/sensitivity.csv)

**Table 7B. Agreement of the generated t = 0 population with the real endpoints.** Its mean WD was **20.45 ± 3.54** to held-out Stage A and **17.00 ± 3.13** to held-out Stage B. Stage B was closer in each seed. The LaTeX file contains both embryoid-body panels. [Per-seed endpoint values](embryoid_body/endpoint_wd.csv) · [Table for Overleaf](embryoid_body/table7_reference.tex)

The held-out endpoint audit addresses a separate question. At both endpoints, the autoencoder had lower cell-wise RMSE than the training-mean predictor, and 64-component PCA had the lowest cell-wise RMSE, averaged across seeds. The training-mean predictor had the lowest median gene-wise RMSE at both endpoints. In the encoded space the stage-matched coordinate had the smallest mean sliced Wasserstein distance at the first two held-out stages. At the last held-out stage the closest other intermediate coordinate had a slightly smaller mean distance (0.403) than the matched coordinate (0.406). [Reconstruction](embryoid_body/summary/reconstruction.csv) · [Latent-space distances](embryoid_body/summary/latent.csv) · [Marker trend correlations](embryoid_body/summary/mgtc.csv)

### Does the distance space change the seed-42 result?

We replayed the frozen seed-42 generator and recovered the saved real cell IDs, anchor pairs, training split and ordered genes. We evaluated all four controls first in the ten released principal components and then directly across the 2,000 retained genes on the model's log-normalized expression scale. Each comparison used the same 200 cells or profiles per population. Both calculations used Euclidean costs and uniform weights. Distances are reported in the units of each space.

| Observed window | 10 PCs, matched | 10 PCs, real Stage B | 2,000 genes, matched | 2,000 genes, real Stage B |
| --- | ---: | ---: | ---: | ---: |
| Day 6–9 | 21.71 | 8.06 | 34.36 | 20.41 |
| Day 12–15 | 20.35 | 6.00 | 34.54 | 20.49 |
| Day 18–21 | 19.51 | 5.56 | 34.08 | 20.64 |

Real Stage B was closer than the time-matched generated population at all three stages in both spaces. The full-gene comparison covers seed 42, the fit with a published generator checkpoint. The replay checks the four archived 10-PC distances and the paired RMS at each stage, plus both t = 0 endpoint distances, with relative tolerance `1e-4` and absolute tolerance `1e-6`. [All six rows and controls](embryoid_body/seed-42-distance-spaces.csv) · [Verification run](https://github.com/sumantaray/GRAIL/actions/runs/45414419675)

To reproduce the comparison, use the command in [Replay the embryoid-body generator](#replay-the-embryoid-body-generator) to create `eb_seed42_predictions.h5ad`. From the repository root, write the comparison to a new CSV file.

```bash
python scripts/replay_table7.py eb_seed42_predictions.h5ad --out eb_seed42_distances.csv
```

## Simulated P6 branch composition

The newly simulated Stage B sample contained 187 of 350 cells on branch 1 (53.4%) and 163 of 350 on branch 2 (46.6%). A classifier fitted to annotated Stage B training cells assigned branch membership to 200 generated terminal profiles per seed.

| Training seed | Generated branch 1 | Generated branch 2 |
| ---: | ---: | ---: |
| 42 | 187/200 (93.5%) | 13/200 (6.5%) |
| 123 | 200/200 (100.0%) | 0/200 (0.0%) |
| 2026 | 200/200 (100.0%) | 0/200 (0.0%) |

These are predicted branch assignments in generated expression space. The reference counts use the observed simulator labels. The [per-seed counts](P6/branch.csv) retain both denominators, and the [P6 dataset description](../datasets/simulations/README.md) specifies how the endpoints were sampled.

## Replay the embryoid-body generator

The [seed-42 generator](embryoid_body/seed-42-generator) includes the fitted inference weights, paired endpoint codes, ordered features, cell and anchor IDs, and configuration from the full LSH run. Run the following command from the repository root after installing the package.

```bash
grail generate --run results/embryoid_body/seed-42-generator --times 0 0.25 0.50 0.75 1 --out eb_seed42_predictions.h5ad
```

On the tested CPU environment, this command reproduced all 1,000 generated profiles over 2,000 genes, with identical expression values, cell IDs and gene ordering. The training command builds both autoencoder and generator checkpoints for any seed. [Training commands](../docs/quickstart.md) · [Tested dependencies](../requirements-tested.txt)

## File guide

Each dataset directory contains one CSV per metric and a `summary/` directory with per-variant means, sample standard deviations and contributing seed counts. Each compressed `run_records.json.gz` contains the fitted settings, input hash, software versions, split and reference IDs, marker and gene IDs, and per-run output file hashes. A reconstruction row already summarizes errors across genes or cells within one seed. The aggregated tables then summarize those seed-level values. [Embryoid-body files](embryoid_body) · [P6 files](P6) · [Metrics and manuscript table mapping](../docs/api.md)
