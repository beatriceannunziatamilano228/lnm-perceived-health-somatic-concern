#!/usr/bin/env python3
"""Query Harvard-Oxford cortical and subcortical atlases at cluster peaks."""
from __future__ import annotations
import argparse
from pathlib import Path
import csv
import shutil
import subprocess

ATLASES = [
    "Harvard-Oxford Cortical Structural Atlas",
    "Harvard-Oxford Subcortical Structural Atlas",
]

def query(atlasq: str, atlas: str, x: str, y: str, z: str) -> str:
    cmd = [atlasq, "query", atlas, "-c", x, y, z]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return " ".join(completed.stdout.strip().split())

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cluster-dir", type=Path, required=True)
    args = p.parse_args()
    atlasq = shutil.which("atlasq")
    if atlasq is None:
        raise SystemExit("FSL atlasq was not found on PATH.")

    for table in sorted(args.cluster_dir.glob("*_clusters_frozen.tsv")):
        rows = list(csv.DictReader(table.open(), delimiter="\t"))
        output = table.with_name(table.stem.replace("_frozen", "_atlas_labels") + ".tsv")
        fields = list(rows[0].keys()) + ["harvard_oxford_cortical", "harvard_oxford_subcortical"]
        with output.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in rows:
                labels = [query(atlasq, atlas, row["max_x_mm"], row["max_y_mm"], row["max_z_mm"]) for atlas in ATLASES]
                row["harvard_oxford_cortical"], row["harvard_oxford_subcortical"] = labels
                writer.writerow(row)
        print(output)

if __name__ == "__main__":
    main()
