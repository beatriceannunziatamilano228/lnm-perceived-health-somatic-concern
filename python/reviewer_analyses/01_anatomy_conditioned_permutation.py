#!/usr/bin/env python3
"""
reviewer-requested anatomy-conditioned cross-dataset permutation test.

Inputs expected in the same folder as this script:
- Corbettaconnectivitymaps101.mat
- CorbettaSF02_controls.mat
- GrafmanmapsnoNaN.mat
- NBRwithoutNaN_NBR02column1.mat

Data structure:
- Corbetta behavioral matrix: column 1 = SF-02; columns 2:end = covariates,
  including lesion size.
- Grafman behavioral matrix: column 1 = NBR-02; columns 2:end = covariates,
  including lesion size.

The analysis:
1. Reconstructs each voxel-wise partial-correlation map by residualizing both
   connectivity and the outcome for the supplied covariates.
2. Computes the observed spatial Pearson correlation between the two maps.
3. Keeps the Corbetta map fixed.
4. Keeps Grafman lesion connectivity maps, lesion locations, and covariates
   (including lesion size) fixed.
5. Permutes the residualized Grafman outcome and recomputes the cross-dataset
   spatial correlation efficiently.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt


FILE_NAMES = {
    "corbetta_maps": "Corbettaconnectivitymaps101.mat",
    "corbetta_behavior": "CorbettaSF02_controls.mat",
    "grafman_maps": "GrafmanmapsnoNaN.mat",
    "grafman_behavior": "NBRwithoutNaN_NBR02column1.mat",
}


def load_inputs(data_dir: Path):
    """Load the four MATLAB files and return numeric arrays."""
    c_maps_mat = sio.loadmat(data_dir / FILE_NAMES["corbetta_maps"])
    c_beh_mat = sio.loadmat(data_dir / FILE_NAMES["corbetta_behavior"])
    g_maps_mat = sio.loadmat(
        data_dir / FILE_NAMES["grafman_maps"],
        squeeze_me=True,
        struct_as_record=False,
    )
    g_beh_mat = sio.loadmat(data_dir / FILE_NAMES["grafman_behavior"])

    corbetta_maps = np.asarray(c_maps_mat["CorbettaFinal"])
    corbetta_behavior = np.asarray(c_beh_mat["question2_controls"], dtype=float)

    grafman_struct = g_maps_mat["Grafman2"]
    grafman_maps = np.asarray(grafman_struct.mapsCopy)

    grafman_behavior = np.asarray(
        g_beh_mat["NBRwithoutNaN_NBR02column1"],
        dtype=float,
    )

    return (
        corbetta_maps,
        corbetta_behavior,
        grafman_maps,
        grafman_behavior,
    )


def residual_projection(behavior: np.ndarray):
    """
    Return:
    - unit-norm outcome residuals;
    - raw outcome residuals;
    - residual-maker matrix for connectivity;
    - valid participant mask;
    - rank of the design matrix.
    """
    valid = np.isfinite(behavior).all(axis=1)
    b = behavior[valid]

    y = b[:, 0]
    covariates = b[:, 1:]
    design = np.column_stack([np.ones(len(y)), covariates])

    beta_y, *_ = np.linalg.lstsq(design, y, rcond=None)
    y_residual = y - design @ beta_y

    norm = np.linalg.norm(y_residual)
    if not np.isfinite(norm) or norm <= np.finfo(float).eps:
        raise ValueError("Residualized outcome has zero or invalid variance.")

    y_unit = y_residual / norm

    projection = np.eye(len(y)) - design @ np.linalg.pinv(design)
    rank = np.linalg.matrix_rank(design)

    return y_unit, y_residual, projection, valid, rank


def observed_partial_map(
    maps: np.ndarray,
    projection: np.ndarray,
    y_unit: np.ndarray,
    chunk_size: int = 10_000,
) -> np.ndarray:
    """Compute a voxel-wise partial-correlation map in chunks."""
    n_voxels = maps.shape[0]
    result = np.full(n_voxels, np.nan, dtype=np.float64)

    for start in range(0, n_voxels, chunk_size):
        stop = min(start + chunk_size, n_voxels)
        block = np.asarray(maps[start:stop, :], dtype=np.float64)

        residualized = block @ projection
        norms = np.linalg.norm(residualized, axis=1)
        good = np.isfinite(norms) & (norms > np.finfo(float).eps)

        result[start:stop][good] = (
            residualized[good] @ y_unit
        ) / norms[good]

    return result


def grafman_sufficient_statistics(
    grafman_maps: np.ndarray,
    projection: np.ndarray,
    fixed_corbetta_map: np.ndarray,
    chunk_size: int = 10_000,
):
    """
    Accumulate sufficient statistics allowing fast permutation testing without
    rebuilding and retaining a 285k-voxel map for every permutation.
    """
    fixed_centered = fixed_corbetta_map - np.mean(fixed_corbetta_map)
    fixed_ss = float(fixed_centered @ fixed_centered)

    n_subjects = grafman_maps.shape[1]
    sum_u = np.zeros(n_subjects, dtype=np.float64)
    utu = np.zeros((n_subjects, n_subjects), dtype=np.float64)
    utc = np.zeros(n_subjects, dtype=np.float64)
    count = 0

    for start in range(0, grafman_maps.shape[0], chunk_size):
        stop = min(start + chunk_size, grafman_maps.shape[0])
        block = np.asarray(grafman_maps[start:stop, :], dtype=np.float64)

        residualized = block @ projection
        norms = np.linalg.norm(residualized, axis=1)
        good = np.isfinite(norms) & (norms > np.finfo(float).eps)

        u = residualized[good] / norms[good, None]
        c = fixed_centered[start:stop][good]

        sum_u += u.sum(axis=0)
        utu += u.T @ u
        utc += u.T @ c
        count += int(good.sum())

    gram_centered = utu - np.outer(sum_u, sum_u) / count

    return utc, gram_centered, fixed_ss, count


def spatial_r_from_stats(
    outcome_unit: np.ndarray,
    cross_vector: np.ndarray,
    gram_matrix: np.ndarray,
    fixed_ss: float,
) -> float:
    numerator = float(cross_vector @ outcome_unit)
    moving_ss = float(outcome_unit @ gram_matrix @ outcome_unit)
    denominator = np.sqrt(fixed_ss * moving_ss)

    if not np.isfinite(denominator) or denominator <= 0:
        return np.nan

    return numerator / denominator


def run_permutations(
    outcome_residual: np.ndarray,
    cross_vector: np.ndarray,
    gram_matrix: np.ndarray,
    fixed_ss: float,
    n_permutations: int,
    seed: int,
    batch_size: int = 1_000,
) -> np.ndarray:
    """Run residual permutations in vectorized batches."""
    rng = np.random.default_rng(seed)
    norm = np.linalg.norm(outcome_residual)
    n_subjects = len(outcome_residual)
    null = np.empty(n_permutations, dtype=np.float64)

    for start in range(0, n_permutations, batch_size):
        current_batch = min(batch_size, n_permutations - start)

        permuted = np.column_stack(
            [
                outcome_residual[rng.permutation(n_subjects)]
                for _ in range(current_batch)
            ]
        ) / norm

        numerators = cross_vector @ permuted
        moving_ss = np.sum(permuted * (gram_matrix @ permuted), axis=0)
        denominators = np.sqrt(fixed_ss * moving_ss)

        null[start:start + current_batch] = numerators / denominators

    return null


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Folder containing the four .mat files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output folder. Default: <data-dir>/anatomy_conditioned_results.",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=100_000,
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260724,
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000,
    )
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else data_dir / "anatomy_conditioned_results"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    for filename in FILE_NAMES.values():
        path = data_dir / filename
        if not path.is_file():
            raise FileNotFoundError(f"Missing required file: {path}")

    (
        corbetta_maps,
        corbetta_behavior,
        grafman_maps,
        grafman_behavior,
    ) = load_inputs(data_dir)

    c_y_unit, c_y_raw, c_projection, c_valid, c_rank = residual_projection(
        corbetta_behavior
    )
    g_y_unit, g_y_raw, g_projection, g_valid, g_rank = residual_projection(
        grafman_behavior
    )

    corbetta_maps = corbetta_maps[:, c_valid]
    grafman_maps = grafman_maps[:, g_valid]

    if corbetta_maps.shape[0] != grafman_maps.shape[0]:
        raise ValueError("The two datasets do not have the same voxel count.")

    print("Reconstructing observed Corbetta map...")
    corbetta_observed = observed_partial_map(
        corbetta_maps,
        c_projection,
        c_y_unit,
        args.chunk_size,
    )

    print("Reconstructing observed Grafman map...")
    grafman_observed = observed_partial_map(
        grafman_maps,
        g_projection,
        g_y_unit,
        args.chunk_size,
    )

    finite = np.isfinite(corbetta_observed) & np.isfinite(grafman_observed)
    corbetta_observed = corbetta_observed[finite]
    grafman_observed = grafman_observed[finite]
    grafman_maps = grafman_maps[finite]

    observed_r = float(
        np.corrcoef(corbetta_observed, grafman_observed)[0, 1]
    )

    # This Fisher-z value is reported only as a reproducibility check against
    # the historical MATLAB implementation. The permutation statistic below
    # uses the raw partial-correlation maps consistently for observed and null.
    fisher_r = float(
        np.corrcoef(
            np.arctanh(np.clip(corbetta_observed, -0.999999, 0.999999)),
            np.arctanh(np.clip(grafman_observed, -0.999999, 0.999999)),
        )[0, 1]
    )

    print("Preparing fast permutation statistics...")
    cross, gram, fixed_ss, n_voxels = grafman_sufficient_statistics(
        grafman_maps,
        g_projection,
        corbetta_observed,
        args.chunk_size,
    )

    observed_check = spatial_r_from_stats(
        g_y_unit,
        cross,
        gram,
        fixed_ss,
    )
    if not np.isclose(observed_r, observed_check, atol=1e-10):
        raise RuntimeError(
            "Internal consistency check failed for observed spatial r."
        )

    print(f"Running {args.n_permutations:,} permutations...")
    null = run_permutations(
        g_y_raw,
        cross,
        gram,
        fixed_ss,
        args.n_permutations,
        args.seed,
    )

    valid_null = np.isfinite(null)
    exceedances_one_sided = int(np.sum(null[valid_null] >= observed_r))
    exceedances_two_sided = int(
        np.sum(np.abs(null[valid_null]) >= abs(observed_r))
    )
    p_one_sided = (
        1 + exceedances_one_sided
    ) / (1 + int(valid_null.sum()))
    p_two_sided = (
        1 + exceedances_two_sided
    ) / (1 + int(valid_null.sum()))

    summary = {
        "corbetta_n": int(c_valid.sum()),
        "grafman_n_complete_case": int(g_valid.sum()),
        "n_voxels": int(n_voxels),
        "corbetta_design_rank": int(c_rank),
        "grafman_design_rank": int(g_rank),
        "observed_spatial_r_raw_maps": observed_r,
        "observed_spatial_r_fisher_z_check": fisher_r,
        "n_permutations": int(args.n_permutations),
        "seed": int(args.seed),
        "one_sided_permutation_p": float(p_one_sided),
        "two_sided_permutation_p": float(p_two_sided),
        "null_mean": float(np.mean(null[valid_null])),
        "null_sd": float(np.std(null[valid_null], ddof=1)),
        "null_2_5_percentile": float(np.quantile(null[valid_null], 0.025)),
        "null_97_5_percentile": float(np.quantile(null[valid_null], 0.975)),
        "null_95_percentile": float(np.quantile(null[valid_null], 0.95)),
        "one_sided_exceedances": exceedances_one_sided,
        "two_sided_exceedances": exceedances_two_sided,
    }

    pd.DataFrame([summary]).to_csv(
        output_dir / "anatomy_conditioned_summary.csv",
        index=False,
    )
    np.save(output_dir / "anatomy_conditioned_null_distribution.npy", null)

    with open(output_dir / "anatomy_conditioned_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    plt.figure(figsize=(8, 5))
    plt.hist(null[valid_null], bins=60)
    plt.axvline(observed_r, linewidth=2)
    plt.xlabel("Spatial correlation with the fixed Corbetta map")
    plt.ylabel("Permutations")
    plt.title(
        f"Anatomy-conditioned permutation test\n"
        f"Observed r = {observed_r:.3f}, two-sided p = {p_two_sided:.4f}"
    )
    plt.tight_layout()
    plt.savefig(
        output_dir / "anatomy_conditioned_null_distribution.png",
        dpi=300,
    )
    plt.close()

    print("\nFINAL RESULTS")
    print(f"Observed raw-map spatial r: {observed_r:.6f}")
    print(f"Fisher-z reproducibility check: {fisher_r:.6f}")
    print(f"One-sided permutation p: {p_one_sided:.6f}")
    print(f"Two-sided permutation p: {p_two_sided:.6f}")
    print(f"Complete-case Grafman n: {int(g_valid.sum())}")
    print(f"Results saved in: {output_dir}")


if __name__ == "__main__":
    main()
