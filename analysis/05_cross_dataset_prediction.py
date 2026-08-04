#!/usr/bin/env python3
"""Cross-dataset prediction analyses for Figures 3 and the historical reverse check."""

from pathlib import Path
import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, rankdata, t as tdist

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from analysis_utils import load_inputs, complete_case, partial_correlation_map


def corr_each_column(vector: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Pearson spatial correlation between one voxel vector and each matrix column."""
    vector = np.asarray(vector, dtype=np.float64)
    matrix = np.asarray(matrix, dtype=np.float64)
    vector_centered = vector - np.mean(vector)
    matrix_centered = matrix - np.mean(matrix, axis=0, keepdims=True)
    return (vector_centered @ matrix_centered) / (
        np.linalg.norm(vector_centered) * np.linalg.norm(matrix_centered, axis=0)
    )


def partial_correlation(x, y, covariates, rank=False):
    """Partial Pearson or rank-based partial correlation with two-sided p-value."""
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(covariates).all(axis=1)
    x = np.asarray(x[valid], dtype=float)
    y = np.asarray(y[valid], dtype=float)
    covariates = np.asarray(covariates[valid], dtype=float)

    if rank:
        x = rankdata(x)
        y = rankdata(y)
        covariates = np.column_stack(
            [rankdata(covariates[:, j]) for j in range(covariates.shape[1])]
        )

    design = np.column_stack([np.ones(len(y)), covariates])
    residual_maker = np.eye(len(y)) - design @ np.linalg.pinv(design)
    xr = residual_maker @ x
    yr = residual_maker @ y
    r = float((xr @ yr) / (np.linalg.norm(xr) * np.linalg.norm(yr)))
    df = len(y) - covariates.shape[1] - 2
    t_value = r * np.sqrt(df / (1 - r * r))
    p_value = float(2 * tdist.sf(abs(t_value), df))
    return r, p_value, len(y), df


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

    c_primary = partial_correlation_map(
        c_maps[:, c_valid], c_y, c_cov, args.chunk_size
    )
    g_primary = partial_correlation_map(
        g_maps[:, g_valid], g_y, g_cov, args.chunk_size
    )

    # Figure 3: Dataset 1 map evaluated in all 181 Dataset 2 participants.
    predictor_g = corr_each_column(c_primary, g_maps)
    predictor_g_z = np.arctanh(np.clip(predictor_g, -0.999999, 0.999999))
    nrs02 = g_beh[:, 0]
    forward_s = spearmanr(predictor_g_z, nrs02)
    forward_p = pearsonr(predictor_g_z, nrs02)

    forward = {
        "cohort_n": int(len(nrs02)),
        "analytic_n": int(np.isfinite(nrs02).sum()),
        "experimental_unit": "participant-level lesion-connectivity map",
        "spearman_rho": float(forward_s.statistic),
        "spearman_p_two_sided": float(forward_s.pvalue),
        "pearson_r_sensitivity": float(forward_p.statistic),
        "pearson_p_two_sided": float(forward_p.pvalue),
    }

    # Historical reverse analysis; not confirmatory and not reported in current Results.
    predictor_c = corr_each_column(g_primary, c_maps)
    predictor_c_z = np.arctanh(np.clip(predictor_c, -0.999999, 0.999999))
    reverse_pearson = partial_correlation(
        predictor_c_z, c_beh[:, 0], c_beh[:, 1:], rank=False
    )
    reverse_spearman = partial_correlation(
        predictor_c_z, c_beh[:, 0], c_beh[:, 1:], rank=True
    )

    reverse = {
        "cohort_n": int(c_beh.shape[0]),
        "analytic_n": int(reverse_pearson[2]),
        "partial_pearson_r": reverse_pearson[0],
        "partial_pearson_p_two_sided": reverse_pearson[1],
        "partial_spearman_rho": reverse_spearman[0],
        "partial_spearman_p_two_sided": reverse_spearman[1],
        "recommendation": (
            "Historical exploratory analysis; do not present as confirmatory validation."
        ),
    }

    pd.DataFrame([forward]).to_csv(
        args.output_dir / "figure3_prediction.csv", index=False
    )
    pd.DataFrame([reverse]).to_csv(
        args.output_dir / "reverse_prediction.csv", index=False
    )

    with open(args.output_dir / "cross_dataset_prediction_summary.json", "w") as f:
        json.dump({"forward": forward, "reverse": reverse}, f, indent=2)

    print(json.dumps({"forward": forward, "reverse": reverse}, indent=2))


if __name__ == "__main__":
    main()
