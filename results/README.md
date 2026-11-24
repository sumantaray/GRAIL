# Measured reference results

These results come from fresh runs of the manuscript-based implementation using seeds 42, 123, and 2026. All fitted models used the [reference configuration](../configs/reference.json). Each run performed 100 autoencoder epochs and 3,000 adversarial updates. The embryoid-body experiment includes both LSH and uniform random anchor selection with an independent fit for each seed and selection method. The P6 experiment uses LSH and a newly generated simulation.

## Embryoid-body expression at held-out stages

For each observed intermediate population, the table reports its distance to 200 time-matched generated profiles and to 200 held-out real Stage B cells. The same observed intermediate cells are used for both comparisons within each seed. Distances are empirical Wasserstein-1 in ten unwhitened principal components fitted on training endpoints. Smaller values indicate closer distributions. Entries give the three-seed mean and sample SD for LSH anchor selection and nearest-neighbor pairing.

| Observed collection window | Coordinate | Distance to time-matched GRAIL | Distance to real Stage B | RMS expression change from generated t = 0 |
| --- | ---: | ---: | ---: | ---: |
| Day 6–9 | 0.25 | 17.81 ± 4.21 | 7.93 ± 0.12 | 0.052 ± 0.014 |
| Day 12–15 | 0.50 | 16.69 ± 4.10 | 5.85 ± 0.14 | 0.107 ± 0.031 |
| Day 18–21 | 0.75 | 16.50 ± 3.60 | 5.54 ± 0.10 | 0.178 ± 0.068 |

The generated profiles changed as the input coordinate moved. Held-out real Stage B cells were closer to the observations in all nine stage-by-seed comparisons. In particular, this run does not reproduce the Table 7 ordering at Days 6 and 12. The original analysis matrices, preprocessing choices, distance implementation and model checkpoints would be needed to reconcile the results. These distances are specific to the feature space and units defined above. [Per-seed expression comparisons](embryoid_body/sensitivity.csv) · [Aggregates, including the fixed-A and other-time controls](embryoid_body/summary/sensitivity.csv)

The held-out endpoint audit addresses a separate question. At both endpoints, the autoencoder had lower cell-wise RMSE than the training-mean predictor, and 64-component PCA had the lowest cell-wise RMSE, averaged across seeds. The training-mean predictor had the lowest median gene-wise RMSE at both endpoints. In the encoded space the stage-matched coordinate had the smallest mean sliced Wasserstein distance at the first two held-out stages. At the last held-out stage the closest other intermediate coordinate had a slightly smaller mean distance (0.403) than the matched coordinate (0.406). [Reconstruction](embryoid_body/summary/reconstruction.csv) · [Latent-space distances](embryoid_body/summary/latent.csv) · [Marker trend correlations](embryoid_body/summary/mgtc.csv)

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
