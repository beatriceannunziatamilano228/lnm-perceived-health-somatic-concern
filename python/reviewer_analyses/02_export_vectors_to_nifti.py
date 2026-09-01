#!/usr/bin/env python3
"""Export a 285,903-element group vector to the historical 2-mm MNI grid.

The non-zero voxels of the supplied historical combined-map NIfTI define the
analysis mask. Vector assignment follows MATLAB column-major linear indexing,
which preserves the ordering used by the restricted connectivity matrices.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np


def export_vector(vector_file: Path, reference_file: Path, output_file: Path) -> None:
    vector = np.asarray(np.load(vector_file), dtype=np.float32).reshape(-1)
    reference_img = nib.load(str(reference_file))
    reference = np.asanyarray(reference_img.dataobj)

    reference_flat = reference.ravel(order="F")
    mask_flat = np.isfinite(reference_flat) & (reference_flat != 0)
    if int(mask_flat.sum()) != int(vector.size):
        raise ValueError(
            f"Reference mask has {int(mask_flat.sum())} voxels, "
            f"but vector has {int(vector.size)} values."
        )

    volume_flat = np.zeros(reference_flat.size, dtype=np.float32)
    volume_flat[mask_flat] = vector
    volume = volume_flat.reshape(reference.shape, order="F")

    header = reference_img.header.copy()
    header.set_data_dtype(np.float32)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(volume, reference_img.affine, header), str(output_file))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vector", type=Path, required=True, help="Input .npy vector.")
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Historical combined-map NIfTI defining the 285,903-voxel mask.",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output .nii or .nii.gz.")
    args = parser.parse_args()
    export_vector(args.vector.resolve(), args.reference.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
