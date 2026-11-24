# GRAIL

**Generative Reconstruction of Unobserved Dynamics of Single-Cell Transcriptomic Trajectories**
Sumanta Ray, Sukhen Das Mandal, Snehalika Lall, and Saumyadipta Pyne

GRAIL takes single-cell expression counts from two endpoint populations and generates candidate expression profiles between them. The model selects endpoint anchors, encodes them with a shared autoencoder, and applies a separately trained generator to interpolated pairs of codes. Observed intermediate cells are reserved for evaluation.

![GRAIL workflow, Figure 1 from the revised manuscript](docs/figures/figure-1-workflow.png)

## Run the example

The example contains 480 **observed human embryoid-body cells** with their original counts, sampled across the five collection windows. The short configuration exercises training and evaluation on a CPU. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python scripts/directional_check.py
grail run --data datasets/example/eb_480.h5ad --config configs/quick.json --out runs/example --seed 42 --save-observed
python scripts/plot_joint.py runs/example/lsh/seed-42
```

Open `runs/example/lsh/seed-42/generated.h5ad` for the generated profiles and `joint_umap.png` for the four-panel figure. An [example notebook](notebooks/quickstart.ipynb) reads the cells and runs the same workflow.

## Data, figures and experiments

| Resource | Location |
| --- | --- |
| Full embryoid-body counts, collection windows and file checksums | [datasets/embryoid_body](datasets/embryoid_body) |
| Dataset provenance, source DOI and mouse cortex availability | [Dataset guide](datasets/README.md) |
| Mouse cortex input status and required stage metadata | [Mouse cortex](datasets/mouse_cortex/README.md) |
| Seven Splatter simulation designs and regeneration commands | [P1–P7](datasets/simulations/README.md) |
| All six figures from the revised manuscript | [Figure gallery](docs/results.md) |
| Measured three-seed results, full embryoid-body dataset and branching simulation | [Results and per-seed tables](results/README.md) |
| Trained embryoid-body generator, seed 42, with paired codes and feature IDs | [Replay the checkpoint](results/README.md#replay-the-embryoid-body-generator) |
| Full training, three seeds, ablations and held-out-stage evaluation | [Experiment commands](docs/quickstart.md) |
| Exact model inputs and implementation choices | [Algorithm and settings](docs/theory.md) |
| Checkpoints, evaluation outputs and table mapping | [Output guide](docs/api.md) |

## Version and provenance

This is a reference implementation reconstructed from the revised Methods and Algorithm 1. It records the configuration, seed, selected genes, cell IDs, anchor pairs and software versions for each new run. The original training checkpoints, historical per-seed outputs, and exact mouse cortex input were unavailable for comparison. The manuscript figures are supplied as figure exports. The comparison tables produced by the commands above contain measurements from the corresponding new run. See the [validation guide](docs/results.md) for interpretation.
