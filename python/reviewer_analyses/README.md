# Analyses introduced during peer review

## `01_anatomy_conditioned_permutation.py`

Reviewer-requested test of whether the observed Dataset 1–Dataset 2 spatial correspondence could arise from the fixed Dataset 2 lesion-connectivity/anatomical sampling structure. Dataset 1 is held fixed; Dataset 2 connectivity maps and covariates remain fixed; residualized Dataset 2 outcome values are permuted. The manuscript uses the two-sided absolute empirical result.

Frozen result: observed raw-map spatial `r = 0.5902`, 100,000 permutations, two-sided `p = 0.04050`.


## `02_export_vectors_to_nifti.py`

Exports a 285,903-element group vector from NumPy format to the historical 2-mm MNI grid using MATLAB-compatible column-major mask ordering. The supplied historical combined map can be used as the reference mask.
