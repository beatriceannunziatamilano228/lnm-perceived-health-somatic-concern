"""Shared utilities for lesion network mapping analyses."""

from __future__ import annotations
from pathlib import Path
import numpy as np
import scipy.io as sio


EXPECTED_FILES = {
    "corbetta_maps": "Corbettaconnectivitymaps101.mat",
    "corbetta_behavior": "CorbettaSF02_controls.mat",
    "grafman_maps": "GrafmanmapsnoNaN.mat",
    "grafman_behavior": "NBRwithoutNaN_NBR02column1.mat",
}


def load_inputs(data_dir: Path):
    """Load the four analysis inputs."""
    data_dir = Path(data_dir)

    c_maps_mat = sio.loadmat(data_dir / EXPECTED_FILES["corbetta_maps"])
    c_beh_mat = sio.loadmat(data_dir / EXPECTED_FILES["corbetta_behavior"])
    g_maps_mat = sio.loadmat(
        data_dir / EXPECTED_FILES["grafman_maps"],
        squeeze_me=True,
        struct_as_record=False,
    )
    g_beh_mat = sio.loadmat(data_dir / EXPECTED_FILES["grafman_behavior"])

    corbetta_maps = np.asarray(c_maps_mat["CorbettaFinal"])
    corbetta_behavior = np.asarray(c_beh_mat["question2_controls"], dtype=float)

    grafman_struct = g_maps_mat["Grafman2"]
    grafman_maps = np.asarray(grafman_struct.mapsCopy)

    grafman_behavior = np.asarray(
        g_beh_mat["NBRwithoutNaN_NBR02column1"], dtype=float
    )

    return corbetta_maps, corbetta_behavior, grafman_maps, grafman_behavior


def complete_case(behavior: np.ndarray):
    """Return complete-case mask, outcome and covariates."""
    valid = np.isfinite(behavior).all(axis=1)
    b = np.asarray(behavior[valid], dtype=np.float64)
    y = b[:, 0]
    covariates = b[:, 1:]
    return valid, y, covariates


def residual_maker(covariates: np.ndarray):
    """Residual-maker matrix for an intercept + covariates model."""
    x = np.column_stack([np.ones(covariates.shape[0]), covariates])
    return np.eye(x.shape[0]) - x @ np.linalg.pinv(x)


def partial_correlation_map(
    maps: np.ndarray,
    outcome: np.ndarray,
    covariates: np.ndarray,
    chunk_size: int = 10_000,
):
    """Voxelwise partial correlation between connectivity and outcome."""
    maps = np.asarray(maps)
    outcome = np.asarray(outcome, dtype=np.float64)
    m = residual_maker(np.asarray(covariates, dtype=np.float64))

    y_res = m @ outcome
    y_norm = np.linalg.norm(y_res)

    result = np.full(maps.shape[0], np.nan, dtype=np.float64)

    for start in range(0, maps.shape[0], chunk_size):
        stop = min(start + chunk_size, maps.shape[0])
        block = np.asarray(maps[start:stop], dtype=np.float64)
        map_res = block @ m
        denom = np.linalg.norm(map_res, axis=1) * y_norm
        good = np.isfinite(denom) & (denom > np.finfo(float).eps)
        result[start:stop][good] = (map_res[good] @ y_res) / denom[good]

    return result


def correlation_map(
    maps: np.ndarray,
    outcome: np.ndarray,
    chunk_size: int = 10_000,
):
    """Voxelwise Pearson correlation without covariates."""
    maps = np.asarray(maps)
    outcome = np.asarray(outcome, dtype=np.float64)
    y = outcome - np.mean(outcome)
    y_norm = np.linalg.norm(y)

    result = np.full(maps.shape[0], np.nan, dtype=np.float64)

    for start in range(0, maps.shape[0], chunk_size):
        stop = min(start + chunk_size, maps.shape[0])
        block = np.asarray(maps[start:stop], dtype=np.float64)
        block_centered = block - np.mean(block, axis=1, keepdims=True)
        denom = np.linalg.norm(block_centered, axis=1) * y_norm
        good = np.isfinite(denom) & (denom > np.finfo(float).eps)
        result[start:stop][good] = (block_centered[good] @ y) / denom[good]

    return result


def spatial_correlation(a: np.ndarray, b: np.ndarray, fisher_z: bool = False):
    """Pearson spatial correlation between two voxelwise maps."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    valid = np.isfinite(a) & np.isfinite(b)
    a = a[valid]
    b = b[valid]

    if fisher_z:
        a = np.arctanh(np.clip(a, -0.999999, 0.999999))
        b = np.arctanh(np.clip(b, -0.999999, 0.999999))

    return float(np.corrcoef(a, b)[0, 1]), int(valid.sum())
