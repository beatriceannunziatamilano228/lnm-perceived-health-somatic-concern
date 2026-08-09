#!/usr/bin/env python3
"""Compare the historical combined map with the Kletenik cross-modal map.

Observed-only mode requires public NIfTI files. Optional permutation mode also
requires the four restricted participant-level analysis matrices.
"""
from __future__ import annotations
import argparse
import json
import shutil
import tempfile
from pathlib import Path

import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
import pandas as pd
import scipy.io as sio

FILES = {
    "c_maps": "Corbettaconnectivitymaps101.mat",
    "c_behavior": "CorbettaSF02_controls.mat",
    "g_maps": "GrafmanmapsnoNaN.mat",
    "g_behavior": "NBRwithoutNaN_NBR02column1.mat",
}

def vector_from_reference(image: nib.spatialimages.SpatialImage, reference: nib.spatialimages.SpatialImage):
    ref_data = np.asanyarray(reference.dataobj)
    if image.shape != reference.shape or not np.allclose(image.affine, reference.affine, atol=1e-4):
        image = resample_from_to(image, reference, order=1)
    data = np.asanyarray(image.dataobj, dtype=np.float64)
    ref_flat = ref_data.ravel(order="F")
    mask = np.isfinite(ref_flat) & (ref_flat != 0)
    return data.ravel(order="F")[mask], mask

def pearson(a: np.ndarray, b: np.ndarray) -> float:
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3:
        return float("nan")
    return float(np.corrcoef(a[valid], b[valid])[0, 1])

def load_behavior(data_dir: Path):
    c = np.asarray(sio.loadmat(data_dir / FILES["c_behavior"])["question2_controls"], dtype=float)
    g = np.asarray(sio.loadmat(data_dir / FILES["g_behavior"])["NBRwithoutNaN_NBR02column1"], dtype=float)
    return c, g

def prepare_behavior(behavior: np.ndarray):
    valid = np.isfinite(behavior).all(axis=1)
    b = behavior[valid]
    y, cov = b[:, 0], b[:, 1:]
    design = np.column_stack([np.ones(len(y)), cov])
    residual_maker = np.eye(len(y)) - design @ np.linalg.pinv(design)
    y_res = residual_maker @ y
    return valid, residual_maker, y_res / np.linalg.norm(y_res), y_res

def load_maps(path: Path, dataset: str):
    if dataset == "c":
        return sio.loadmat(path)["CorbettaFinal"]
    return sio.loadmat(path, squeeze_me=True, struct_as_record=False)["Grafman2"].mapsCopy

def unit_residualized_maps(path: Path, dataset: str, valid: np.ndarray, M: np.ndarray, output: Path, chunk: int):
    maps = load_maps(path, dataset)
    out = np.memmap(output, dtype=np.float32, mode="w+", shape=(maps.shape[0], int(valid.sum())))
    good = np.ones(maps.shape[0], dtype=bool)
    for start in range(0, maps.shape[0], chunk):
        stop = min(start + chunk, maps.shape[0])
        block = np.asarray(maps[start:stop, :][:, valid], dtype=np.float64)
        residualized = block @ M
        norms = np.linalg.norm(residualized, axis=1)
        ok = np.isfinite(norms) & (norms > np.finfo(float).eps)
        out[start:stop] = 0
        out[start:stop][ok] = (residualized[ok] / norms[ok, None]).astype(np.float32)
        good[start:stop] = ok
    out.flush()
    return out, good

def permuted_unit_outcomes(rng, residual, M, batch):
    n = len(residual)
    values = np.column_stack([residual[rng.permutation(n)] for _ in range(batch)])
    values = M @ values
    return values / np.linalg.norm(values, axis=0, keepdims=True)

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--combined-map", type=Path, required=True)
    p.add_argument("--anosognosia-map", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=Path("outputs/anosognosia_comparison"))
    p.add_argument("--data-dir", type=Path, default=None)
    p.add_argument("--n-permutations", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=20260806)
    p.add_argument("--batch-size", type=int, default=25)
    p.add_argument("--chunk-size", type=int, default=5_000)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    combined_img = nib.load(str(args.combined_map))
    combined_vec, mask = vector_from_reference(combined_img, combined_img)
    external_img = nib.load(str(args.anosognosia_map))
    external_vec, _ = vector_from_reference(external_img, combined_img)
    observed_r = pearson(combined_vec, external_vec)

    summary = {
        "external_map": str(args.anosognosia_map),
        "neurovault_image_id": 795501,
        "n_analysis_voxels": int(mask.sum()),
        "observed_spatial_pearson_r": observed_r,
        "permutation_method": None,
        "n_permutations": 0,
        "two_sided_empirical_p": None,
    }

    if args.data_dir is not None:
        data_dir = args.data_dir.resolve()
        for filename in FILES.values():
            if not (data_dir / filename).is_file():
                raise FileNotFoundError(data_dir / filename)
        c_behavior, g_behavior = load_behavior(data_dir)
        c_valid, c_M, c_y, c_res = prepare_behavior(c_behavior)
        g_valid, g_M, g_y, g_res = prepare_behavior(g_behavior)
        cache = Path(tempfile.mkdtemp(prefix="anosognosia_", dir=args.output_dir))
        try:
            c_unit, c_good = unit_residualized_maps(data_dir / FILES["c_maps"], "c", c_valid, c_M, cache / "c.dat", args.chunk_size)
            g_unit, g_good = unit_residualized_maps(data_dir / FILES["g_maps"], "g", g_valid, g_M, cache / "g.dat", args.chunk_size)
            common = c_good & g_good
            external = external_vec[common]
            external_centered = external - np.nanmean(external)
            ext_ss = float(external_centered @ external_centered)
            c_unit = np.asarray(c_unit[common])
            g_unit = np.asarray(g_unit[common])

            # Verify the public historical map definition.
            c_map = c_unit @ c_y
            g_map = g_unit @ g_y
            observed_recomputed = (101 * np.arctanh(np.clip(c_map, -.999999, .999999)) + 181 * np.arctanh(np.clip(g_map, -.999999, .999999))) / 282
            recomputed_r = pearson(observed_recomputed, external)

            rng = np.random.default_rng(args.seed)
            null = np.empty(args.n_permutations, dtype=np.float32)
            index = 0
            while index < args.n_permutations:
                current = min(args.batch_size, args.n_permutations - index)
                c_perm_y = permuted_unit_outcomes(rng, c_res, c_M, current)
                g_perm_y = permuted_unit_outcomes(rng, g_res, g_M, current)
                c_maps = c_unit @ c_perm_y
                g_maps = g_unit @ g_perm_y
                combined = (101 * np.arctanh(np.clip(c_maps, -.999999, .999999)) + 181 * np.arctanh(np.clip(g_maps, -.999999, .999999))) / 282
                combined -= combined.mean(axis=0, keepdims=True)
                numerator = external_centered @ combined
                denominator = np.sqrt(ext_ss * np.sum(combined * combined, axis=0))
                null[index:index + current] = numerator / denominator
                index += current

            p_two = float((1 + np.sum(np.abs(null) >= abs(recomputed_r))) / (args.n_permutations + 1))
            np.save(args.output_dir / "anosognosia_null_two_sided.npy", null)
            summary.update({
                "observed_spatial_pearson_r_recomputed_from_restricted_inputs": float(recomputed_r),
                "permutation_method": "Freedman-Lane residual permutation; historical nominal 101/181 Fisher-z weighted combined map",
                "n_permutations": args.n_permutations,
                "seed": args.seed,
                "two_sided_empirical_p": p_two,
            })
        finally:
            shutil.rmtree(cache, ignore_errors=True)

    (args.output_dir / "anosognosia_comparison_summary.json").write_text(json.dumps(summary, indent=2))
    pd.DataFrame([summary]).to_csv(args.output_dir / "anosognosia_comparison_summary.csv", index=False)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
