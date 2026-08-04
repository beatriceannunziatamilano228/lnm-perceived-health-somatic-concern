#!/usr/bin/env python3
"""
Supplementary Figure S1: agreement between maps estimated with and without covariates.

This script reproduces the *numerical comparison* underlying Supplementary Figure S1
from the supplied connectivity matrices and behavioral tables.

The exact historical brain-slice rendering cannot be reconstructed from these four
files alone because the voxel-to-NIfTI mask/template used by the MATLAB pipeline is
not included. If that mask is recovered, the generated .npy vectors can be exported
to NIfTI without changing the statistical analysis.
"""

from pathlib import Path
import argparse
import json
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from analysis_utils import spatial_correlation


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    c_with = np.load(args.input_dir / "corbetta_primary_partial_r.npy")
    c_without = np.load(args.input_dir / "corbetta_no_covariates_r.npy")
    g_with = np.load(args.input_dir / "grafman_primary_partial_r.npy")
    g_without = np.load(args.input_dir / "grafman_no_covariates_r.npy")

    c_r, c_nvox = spatial_correlation(c_with, c_without)
    g_r, g_nvox = spatial_correlation(g_with, g_without)

    summary = {
        "corbetta_spatial_r": c_r,
        "grafman_spatial_r": g_r,
        "corbetta_voxels": c_nvox,
        "grafman_voxels": g_nvox,
    }
    with open(args.output_dir / "supplementary_figure_S1_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Two independent publication-quality agreement plots.
    for label, with_cov, no_cov, r in [
        ("Dataset 1 (Corbetta)", c_with, c_without, c_r),
        ("Dataset 2 (Grafman)", g_with, g_without, g_r),
    ]:
        valid = np.isfinite(with_cov) & np.isfinite(no_cov)

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.hexbin(no_cov[valid], with_cov[valid], gridsize=80, mincnt=1)
        ax.set_xlabel("Voxelwise r, no covariates")
        ax.set_ylabel("Voxelwise partial r, primary model")
        ax.set_title(f"{label}: spatial r = {r:.3f}")
        fig.tight_layout()

        stem = "S1_corbetta" if "Corbetta" in label else "S1_grafman"
        fig.savefig(args.output_dir / f"{stem}.png", dpi=300)
        fig.savefig(args.output_dir / f"{stem}.pdf")
        plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
