#!/usr/bin/env python3
"""Export a MATLAB-mask-ordered vector to a NIfTI volume.

The supplied historical combined map has 285,903 nonzero voxels and provides the
2-mm MNI affine and mask. MATLAB linear indexing is column-major, so vectors are
inserted into the reference mask using Fortran-order flattening.
"""

from pathlib import Path
import argparse
import numpy as np
import nibabel as nib


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vector", type=Path, required=True, help="Input .npy vector")
    parser.add_argument(
        "--reference-nifti",
        type=Path,
        default=Path("reference_maps/combinedmap_historical.nii"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    vector = np.asarray(np.load(args.vector), dtype=np.float32).reshape(-1)
    reference = nib.load(str(args.reference_nifti))
    reference_data = np.asarray(reference.get_fdata())

    mask_flat = reference_data.ravel(order="F") != 0
    if int(mask_flat.sum()) != vector.size:
        raise ValueError(
            f"Reference mask has {int(mask_flat.sum())} nonzero voxels, "
            f"but vector has {vector.size} values."
        )

    output_flat = np.zeros(reference_data.size, dtype=np.float32)
    output_flat[mask_flat] = vector
    output_data = output_flat.reshape(reference_data.shape, order="F")

    image = nib.Nifti1Image(output_data, reference.affine, reference.header.copy())
    image.set_data_dtype(np.float32)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    nib.save(image, str(args.output))
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
