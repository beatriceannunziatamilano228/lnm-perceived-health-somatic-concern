#!/usr/bin/env python3
"""Reconstruct primary and no-covariate lesion network maps."""

from pathlib import Path
import argparse
import json
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from analysis_utils import (
    load_inputs,
    complete_case,
    partial_correlation_map,
    correlation_map,
    spatial_correlation,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--chunk-size", type=int, default=10_000)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    c_maps, c_beh, g_maps, g_beh = load_inputs(args.data_dir)

    c_valid, c_y, c_cov = complete_case(c_beh)
    g_valid, g_y, g_cov = complete_case(g_beh)

    c_with = partial_correlation_map(
        c_maps[:, c_valid], c_y, c_cov, args.chunk_size
    )
    g_with = partial_correlation_map(
        g_maps[:, g_valid], g_y, g_cov, args.chunk_size
    )

    # No-covariate sensitivity maps use all participants with a finite primary outcome.
    c_outcome_valid = np.isfinite(c_beh[:, 0])
    g_outcome_valid = np.isfinite(g_beh[:, 0])

    c_without = correlation_map(
        c_maps[:, c_outcome_valid], c_beh[c_outcome_valid, 0], args.chunk_size
    )
    g_without = correlation_map(
        g_maps[:, g_outcome_valid], g_beh[g_outcome_valid, 0], args.chunk_size
    )

    c_r, c_vox = spatial_correlation(c_with, c_without)
    g_r, g_vox = spatial_correlation(g_with, g_without)

    c_r_z, _ = spatial_correlation(c_with, c_without, fisher_z=True)
    g_r_z, _ = spatial_correlation(g_with, g_without, fisher_z=True)

    np.save(args.output_dir / "corbetta_primary_partial_r.npy", c_with)
    np.save(args.output_dir / "corbetta_no_covariates_r.npy", c_without)
    np.save(args.output_dir / "grafman_primary_partial_r.npy", g_with)
    np.save(args.output_dir / "grafman_no_covariates_r.npy", g_without)

    summary = {
        "corbetta_n_primary_complete_case": int(c_valid.sum()),
        "corbetta_n_no_covariates": int(c_outcome_valid.sum()),
        "grafman_n_primary_complete_case": int(g_valid.sum()),
        "grafman_n_no_covariates": int(g_outcome_valid.sum()),
        "corbetta_with_vs_without_covariates_spatial_r": c_r,
        "grafman_with_vs_without_covariates_spatial_r": g_r,
        "corbetta_with_vs_without_covariates_fisher_z_spatial_r": c_r_z,
        "grafman_with_vs_without_covariates_fisher_z_spatial_r": g_r_z,
        "voxels_compared_corbetta": c_vox,
        "voxels_compared_grafman": g_vox,
    }

    with open(args.output_dir / "primary_map_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
