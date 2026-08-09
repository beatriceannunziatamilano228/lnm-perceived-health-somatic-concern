# Analyses introduced during peer review

## `01_anatomy_conditioned_permutation.py`

Reviewer-requested test of whether the observed Dataset 1–Dataset 2 spatial correspondence could arise from the fixed Dataset 2 lesion-connectivity/anatomical sampling structure. Dataset 1 is held fixed; Dataset 2 connectivity maps and covariates remain fixed; residualized Dataset 2 outcome values are permuted. The manuscript uses the two-sided absolute empirical result.

Frozen result: observed raw-map spatial `r = 0.5902`, 100,000 permutations, two-sided `p = 0.04050`.

## `02_corrected_maximum_statistic.py`

Two-sided whole-brain maximum-statistic correction using Freedman-Lane residual permutations and complete-case sample-size weighting. This is the authoritative revised multiple-comparisons analysis, distinct from the recovered historical Westfall-Young working script.

Frozen result: 10,000 permutations; threshold `|r| = 0.28064`; global two-sided `pFWE = 0.00580`; 35 positive significant voxels and no negative significant voxels.

Both scripts require the four restricted matrices described in `../../data/README.md`. Frozen non-identifiable null distributions and outputs are included under `../../results/revision/`.

## `03_export_vectors_to_nifti.py`

Exports a 285,903-element group vector from NumPy format to the historical 2-mm MNI grid using MATLAB-compatible column-major mask ordering. The supplied historical combined map can be used as the reference mask.
