#!/usr/bin/env python3
"""Combined lesion-network map and whole-brain maximum-statistic correction.

The script reconstructs the two covariate-adjusted voxelwise partial-correlation
maps, combines them using sample-size weights, and applies a whole-brain
family-wise error correction using a Freedman-Lane residual-permutation scheme.

Default inferential choices
---------------------------
- Complete-case sample size is used for weighting (Dataset 1 n=101; Dataset 2
  n=173 with the currently supplied covariate table).
- The combined statistic is the sample-size-weighted mean of the two raw
  partial-correlation maps, matching the manuscript's stated weighted-mean
  approach.
- A two-sided maximum statistic, max(abs(combined map)), controls FWE across
  the analysed brain mask.
- A plus-one empirical p-value correction is used.

The output vectors remain in the supplied 285,903-voxel MATLAB mask order.
Use `07_export_vectors_to_nifti.py` with the supplied historical group-level
reference NIfTI to export reconstructed vectors into the 2-mm MNI grid.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json
import shutil
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio


FILES = {
    "corbetta_maps": "Corbettaconnectivitymaps101.mat",
    "corbetta_behavior": "CorbettaSF02_controls.mat",
    "grafman_maps": "GrafmanmapsnoNaN.mat",
    "grafman_behavior": "NBRwithoutNaN_NBR02column1.mat",
}


def load_behavior(data_dir: Path):
    c = np.asarray(
        sio.loadmat(data_dir / FILES["corbetta_behavior"])["question2_controls"],
        dtype=np.float64,
    )
    g = np.asarray(
        sio.loadmat(data_dir / FILES["grafman_behavior"])[
            "NBRwithoutNaN_NBR02column1"
        ],
        dtype=np.float64,
    )
    return c, g


def prepare_behavior(behavior: np.ndarray):
    valid = np.isfinite(behavior).all(axis=1)
    b = behavior[valid]
    outcome = b[:, 0]
    covariates = b[:, 1:]
    design = np.column_stack([np.ones(len(outcome)), covariates])
    residual_maker = np.eye(len(outcome)) - design @ np.linalg.pinv(design)
    outcome_residual = residual_maker @ outcome
    residual_norm = np.linalg.norm(outcome_residual)
    if not np.isfinite(residual_norm) or residual_norm <= np.finfo(float).eps:
        raise ValueError("Residualized outcome has zero or invalid variance.")
    outcome_unit = (outcome_residual / residual_norm).astype(np.float32)
    return {
        "valid": valid,
        "outcome": outcome,
        "design_rank": int(np.linalg.matrix_rank(design)),
        "residual_maker": residual_maker,
        "outcome_residual": outcome_residual,
        "outcome_unit": outcome_unit,
    }


def load_maps(path: Path, dataset: str):
    if dataset == "corbetta":
        return sio.loadmat(path)["CorbettaFinal"]
    mat = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    return mat["Grafman2"].mapsCopy


def write_unit_residualized_maps(
    maps_path: Path,
    dataset: str,
    valid: np.ndarray,
    residual_maker: np.ndarray,
    output_path: Path,
    chunk_size: int,
):
    maps = load_maps(maps_path, dataset)
    n_voxels = maps.shape[0]
    n_participants = int(valid.sum())
    unit_maps = np.memmap(
        output_path,
        dtype=np.float32,
        mode="w+",
        shape=(n_voxels, n_participants),
    )
    valid_voxels = np.ones(n_voxels, dtype=bool)

    for start in range(0, n_voxels, chunk_size):
        stop = min(start + chunk_size, n_voxels)
        block = np.asarray(maps[start:stop, :][:, valid], dtype=np.float64)
        residualized = block @ residual_maker
        norms = np.linalg.norm(residualized, axis=1)
        good = np.isfinite(norms) & (norms > np.finfo(float).eps)
        unit_maps[start:stop] = 0
        unit_maps[start:stop][good] = (
            residualized[good] / norms[good, None]
        ).astype(np.float32)
        valid_voxels[start:stop] = good

    unit_maps.flush()
    del maps
    return unit_maps, valid_voxels


def permutation_outcomes(
    rng: np.random.Generator,
    outcome_residual: np.ndarray,
    residual_maker: np.ndarray,
    batch_size: int,
):
    """True Freedman-Lane residual permutations for a reduced covariate model."""
    n = len(outcome_residual)
    permuted = np.column_stack(
        [outcome_residual[rng.permutation(n)] for _ in range(batch_size)]
    )
    # Residualizing fitted + permuted residuals against the reduced model is
    # algebraically equivalent to applying M to the permuted residuals.
    projected = residual_maker @ permuted
    norms = np.linalg.norm(projected, axis=0, keepdims=True)
    return (projected / norms).astype(np.float32)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/maxstat"))
    parser.add_argument("--n-permutations", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20_260_731)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--chunk-size", type=int, default=5_000)
    parser.add_argument(
        "--weighting",
        choices=["effective", "nominal"],
        default="effective",
        help="Use complete-case n or nominal cohort n (101 and 181).",
    )
    parser.add_argument("--keep-cache", action="store_true")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for name in FILES.values():
        if not (data_dir / name).is_file():
            raise FileNotFoundError(data_dir / name)

    c_behavior, g_behavior = load_behavior(data_dir)
    c = prepare_behavior(c_behavior)
    g = prepare_behavior(g_behavior)

    cache_dir = Path(tempfile.mkdtemp(prefix="lnm_maxstat_", dir=output_dir))
    try:
        print("Residualizing Dataset 1 connectivity maps...")
        c_unit, c_good = write_unit_residualized_maps(
            data_dir / FILES["corbetta_maps"],
            "corbetta",
            c["valid"],
            c["residual_maker"],
            cache_dir / "corbetta_unit_maps.dat",
            args.chunk_size,
        )
        print("Residualizing Dataset 2 connectivity maps...")
        g_unit, g_good = write_unit_residualized_maps(
            data_dir / FILES["grafman_maps"],
            "grafman",
            g["valid"],
            g["residual_maker"],
            cache_dir / "grafman_unit_maps.dat",
            args.chunk_size,
        )

        common = c_good & g_good
        if not common.all():
            c_unit = np.asarray(c_unit[common])
            g_unit = np.asarray(g_unit[common])

        c_map = np.asarray(c_unit @ c["outcome_unit"], dtype=np.float64)
        g_map = np.asarray(g_unit @ g["outcome_unit"], dtype=np.float64)

        n1 = int(c["valid"].sum())
        n2_effective = int(g["valid"].sum())
        n2_weight = n2_effective if args.weighting == "effective" else 181
        combined = (n1 * c_map + n2_weight * g_map) / (n1 + n2_weight)

        # Fisher-z weighted sensitivity map; not used for the primary max-stat test.
        c_z = np.arctanh(np.clip(c_map, -0.999999, 0.999999))
        g_z = np.arctanh(np.clip(g_map, -0.999999, 0.999999))
        fisher_combined = np.tanh(
            ((n1 - 3) * c_z + (n2_effective - 3) * g_z)
            / ((n1 - 3) + (n2_effective - 3))
        )

        observed_max = float(np.max(np.abs(combined)))
        observed_positive_max = float(np.max(combined))
        null_two_sided = np.empty(args.n_permutations, dtype=np.float32)
        null_one_sided = np.empty(args.n_permutations, dtype=np.float32)
        rng = np.random.default_rng(args.seed)

        print(f"Running {args.n_permutations:,} Freedman-Lane permutations...")
        index = 0
        while index < args.n_permutations:
            current = min(args.batch_size, args.n_permutations - index)
            c_perm = permutation_outcomes(
                rng,
                c["outcome_residual"],
                c["residual_maker"],
                current,
            )
            g_perm = permutation_outcomes(
                rng,
                g["outcome_residual"],
                g["residual_maker"],
                current,
            )
            c_permuted_maps = c_unit @ c_perm
            g_permuted_maps = g_unit @ g_perm
            combined_permuted = (
                n1 * c_permuted_maps + n2_weight * g_permuted_maps
            ) / (n1 + n2_weight)

            null_two_sided[index:index + current] = np.max(
                np.abs(combined_permuted), axis=0
            )
            null_one_sided[index:index + current] = np.max(
                combined_permuted, axis=0
            )
            index += current
            if index % 1_000 == 0 or index == args.n_permutations:
                print(f"  {index:,}/{args.n_permutations:,}")

        threshold_two_sided = float(np.quantile(null_two_sided, 0.95))
        threshold_one_sided = float(np.quantile(null_one_sided, 0.95))
        p_two_sided = float(
            (1 + np.sum(null_two_sided >= observed_max))
            / (args.n_permutations + 1)
        )
        p_one_sided = float(
            (1 + np.sum(null_one_sided >= observed_positive_max))
            / (args.n_permutations + 1)
        )

        significant = np.abs(combined) >= threshold_two_sided
        thresholded = np.where(significant, combined, 0)

        np.save(output_dir / "dataset1_partial_r.npy", c_map)
        np.save(output_dir / "dataset2_partial_r.npy", g_map)
        np.save(output_dir / "combined_weighted_raw_r.npy", combined)
        np.save(output_dir / "combined_fisher_z_sensitivity.npy", fisher_combined)
        np.save(output_dir / "maxstat_null_two_sided.npy", null_two_sided)
        np.save(output_dir / "maxstat_null_one_sided.npy", null_one_sided)
        np.save(output_dir / "maxstat_significant_mask.npy", significant)
        np.save(output_dir / "maxstat_thresholded_combined_map.npy", thresholded)

        summary = {
            "dataset1_complete_case_n": n1,
            "dataset2_complete_case_n": n2_effective,
            "dataset2_weight_used": n2_weight,
            "weighting": args.weighting,
            "n_voxels": int(len(combined)),
            "cross_dataset_fisher_z_spatial_r": float(
                np.corrcoef(c_z, g_z)[0, 1]
            ),
            "observed_combined_positive_max": observed_positive_max,
            "observed_combined_absolute_max": observed_max,
            "n_permutations": args.n_permutations,
            "seed": args.seed,
            "two_sided_FWE_threshold": threshold_two_sided,
            "two_sided_global_pFWE": p_two_sided,
            "one_sided_FWE_threshold": threshold_one_sided,
            "one_sided_global_pFWE": p_one_sided,
            "significant_voxels_two_sided": int(significant.sum()),
            "significant_positive_voxels": int(
                np.sum(combined >= threshold_two_sided)
            ),
            "significant_negative_voxels": int(
                np.sum(combined <= -threshold_two_sided)
            ),
            "dataset1_design_rank": c["design_rank"],
            "dataset2_design_rank": g["design_rank"],
        }

        with open(output_dir / "maxstat_summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        pd.DataFrame([summary]).to_csv(
            output_dir / "maxstat_summary.csv", index=False
        )

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(null_two_sided, bins=60)
        ax.axvline(observed_max, linewidth=2)
        ax.set_xlabel("Maximum absolute weighted partial correlation")
        ax.set_ylabel("Permutations")
        ax.set_title(
            "Whole-brain maximum-statistic correction\n"
            f"observed max = {observed_max:.3f}, pFWE = {p_two_sided:.4f}"
        )
        fig.tight_layout()
        fig.savefig(output_dir / "maxstat_null_distribution.png", dpi=300)
        fig.savefig(output_dir / "maxstat_null_distribution.pdf")
        plt.close(fig)

        print(json.dumps(summary, indent=2))

    finally:
        if args.keep_cache:
            print(f"Cache retained at: {cache_dir}")
        else:
            shutil.rmtree(cache_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
