# Running GRAIL

Run commands from a clone of the repository root, using Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python scripts/directional_check.py
grail run --data datasets/example/eb_480.h5ad --config configs/quick.json --out runs/example --seed 42 --save-observed
grail generate --run runs/example/lsh/seed-42 --times 0 0.25 0.50 0.75 1 --out runs/example/replayed.h5ad
python scripts/plot_joint.py runs/example/lsh/seed-42
```

The short run trains for two autoencoder epochs and three GAN updates. It checks the software path. For an independent analysis of the full source counts, use the reference configuration below.

## Full embryoid-body analysis and ablations

```bash
grail run --data datasets/embryoid_body --config configs/reference.json --out runs/eb
grail run --data datasets/embryoid_body --config configs/reference.json --anchors random --out runs/eb
python scripts/summarize.py runs/eb
```

Each command trains and evaluates seeds **42, 123, and 2026**. Use `--seed 42` to run a single seed, and `python scripts/summarize.py runs/eb --seeds 42` for that case. Sample SD is undefined for a single seed. Existing seed folders are protected against accidental replacement. Interrupted fits leave completed seed folders intact. Resume by running the missing seeds individually. Set `--save-observed` during fitting to retain processed real cells for the optional joint plot.

Nearest-neighbor, random and Stage B centroid pairing use the same fitted weights and source anchors within a seed. Direct expression interpolation and AE-decoder readout also reuse the matched pairs. The second command retrains with random anchor selection at the same counts for the no-LSH comparison. Hold all other configuration fields fixed when comparing variants. See the [output guide](api.md).

## Comparing external methods

External predictions from TrajectoryNet or Waddington-OT can be evaluated with the same population metrics. Obtain the external methods' own outputs using endpoint-only fitting. Align them to the seed's ordered `features.json` on the same log-normalized scale and provide H5ADs with `obs["t"]` covering the evaluation grid, including 0 and 1. The command uses a fresh directory so existing results are preserved.

```bash
grail run --data datasets/embryoid_body --config configs/reference.json --seed 42 \
  --out runs/eb_external --external TrajectoryNet=path/to/tn_seed42.h5ad \
  --external WaddingtonOT=path/to/wot_seed42.h5ad
```

The authors' external-model checkpoints, exact historical solver configurations and original baseline predictions have not been deposited in this reference package. The import path evaluates externally supplied predictions. [TrajectoryNet](https://github.com/KrishnaswamyLab/TrajectoryNet) · [Waddington-OT](https://github.com/broadinstitute/wot)

## Other datasets

Use a single H5AD with nonnegative integer counts in `X` and unique cell and gene IDs. Set the normalized numeric stage in `obs["t"]`, with 0 and 1 for the endpoints. For a verified mouse cortex input, the intermediate values are 0.40 and 0.70. Additional observed intermediate stages are evaluated automatically. For synthetic endpoint-only data, the pipeline reports geometric and endpoint diagnostics. Add real cell-type annotations in `obs["cell_type"]` to enable composition analysis, or terminal branch labels in `obs["branch"]` for the P6 branch comparison.
