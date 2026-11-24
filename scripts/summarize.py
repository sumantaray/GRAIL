"""Aggregate per-seed results for one dataset/configuration; SD uses ddof=1."""
import argparse
import json
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("runs", type=Path)
    p.add_argument("--seeds", type=int, nargs="+", default=[42, 123, 2026])
    args = p.parse_args()
    if len(set(args.seeds)) != len(args.seeds):
        p.error("Supply distinct seed IDs.")
    groups = sorted(q for q in args.runs.iterdir() if q.is_dir() and list(q.glob("seed-*/run.json")))
    if not groups:
        p.error("No completed runs found.")
    out = args.runs/"summary"
    data_ids, settings = set(), set()
    todo = []
    for group in groups:
        dirs = [group/f"seed-{seed}" for seed in args.seeds]
        for d,seed in zip(dirs,args.seeds):
            rec = json.loads((d/"run.json").read_text())
            if rec["seed"] != seed or rec["anchor_method"] != group.name:
                raise ValueError(f"Mismatched run identity in {d}")
            data_ids.add(rec["input_sha256"])
            settings.add((d/"config.json").read_text())
        metric_names = {"reconstruction", "wd", "population", "sensitivity", "latent", "endpoint_wd",
                        "mgtc", "geometry", "ccm", "branch"}
        names = [{f.name for f in d.glob("*.csv") if f.stem in metric_names} for d in dirs]
        if any(files != names[0] for files in names):
            raise ValueError("Inconsistent metric outputs across seeds.")
        if len(settings) != 1 or len(data_ids) != 1:
            raise ValueError("Runs must use the same input and configuration.")
        for name in sorted(names[0]):
            df = pd.concat([pd.read_csv(d/name) for d in dirs], ignore_index=True)
            keys = [k for k in ("anchor_method","stage","t","reference_t","method","branch") if k in df]
            vals = [c for c in df.select_dtypes("number") if c not in keys+["seed"]]
            table = df.groupby(keys, dropna=False)[vals].agg(["mean","std","count"])
            table.columns = ["_".join(c) for c in table.columns]
            todo.append((group.name, name, table))
    for group, name, table in todo:
        (out/group).mkdir(parents=True, exist_ok=True)
        table.to_csv(out/group/name)
    print(out)


if __name__ == "__main__":
    main()
