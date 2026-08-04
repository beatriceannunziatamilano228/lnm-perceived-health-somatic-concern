#!/usr/bin/env python3
"""Run the validated analysis workflow from the repository root."""

from pathlib import Path
import argparse
import subprocess
import sys


def run(command):
    print("\n$ " + " ".join(map(str, command)))
    subprocess.run(command, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--anatomy-permutations", type=int, default=100_000)
    parser.add_argument("--maxstat-permutations", type=int, default=10_000)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    scripts = root / "scripts"
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    python = sys.executable
    data = args.data_dir.resolve()

    run([python, scripts / "01_primary_maps.py", "--data-dir", data, "--output-dir", out])
    run([python, scripts / "02_supplementary_figure_S1.py", "--input-dir", out, "--output-dir", out])
    run([python, scripts / "05_cross_dataset_prediction.py", "--data-dir", data, "--output-dir", out])
    run([python, scripts / "06_nrs_item_level_comparison.py", "--data-dir", data, "--output-dir", out])
    run([
        python, scripts / "03_anatomy_conditioned_permutation.py",
        "--data-dir", data,
        "--output-dir", out / "anatomy_conditioned_permutation",
        "--n-permutations", str(args.anatomy_permutations),
    ])
    run([
        python, scripts / "04_combined_map_maxstat.py",
        "--data-dir", data,
        "--output-dir", out,
        "--n-permutations", str(args.maxstat_permutations),
    ])


if __name__ == "__main__":
    main()
