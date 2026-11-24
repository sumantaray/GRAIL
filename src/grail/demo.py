"""Convenience entry point for the included real-data example from a checkout."""
import argparse
from pathlib import Path
from .cli import run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("runs/example"))
    args = p.parse_args()
    args.config = Path("configs/quick.json")
    args.data = Path("datasets/example/eb_480.h5ad")
    args.seed, args.anchors, args.external, args.save_observed = 42, "lsh", [], False
    run(args)


if __name__ == "__main__":
    main()
