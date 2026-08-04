# Lesion network mapping of perceived health change and somatic concern

##

This repository is the public code release accompanying the revised manuscript. It contains:

- the cleaned Python code used to validate the reported analyses;
- frozen, non-identifiable outputs generated from the validated study inputs;
- group-level NIfTI maps and permutation null distributions;
- a clear crosswalk between manuscript claims, analysis scripts, and output files.

**The authors do not need to rerun the analyses locally for resubmission.** The validated
outputs used for the revision are already included under `results/`. The code is provided so
editors and readers can inspect the implementation and, where authorized data access is
available, reproduce the analysis.

## Headline validated results

- Dataset 1 primary analytic sample: **n = 101**.
- Dataset 2 cohort: **n = 181**; complete-case primary model: **n = 173**.
- Cross-dataset spatial correspondence: **r = 0.590** on raw partial-correlation maps
  and **r = 0.589** after Fisher transformation.
- Dataset 1 network predicting Dataset 2 NRS-02: **Spearman rho = 0.233,
  p = 0.00162, n = 181**.
- NRS-02 showed the largest adjusted item-level map correspondence among the 27 NRS items.
- Anatomy-conditioned permutation analysis: **100,000 permutations,
  two-sided empirical p = 0.0405**.
- Whole-brain maximum-statistic correction: **10,000 permutations,
  pFWE = 0.0058**, threshold **|r| = 0.28064**, with **35 significant voxels**
  in the medial orbitofrontal cluster; historical peak **MNI [28, 40, -22]**.

## Repository contents

```text
analysis/       Cleaned Python analysis code
src/            Shared numerical functions
results/        Frozen validated tables, maps, figures and null distributions
docs/           Methods, code/data statements and manuscript-output crosswalk
data/           Input schema only; participant-level inputs are not distributed
provenance/     Sanitized historical MATLAB working code
```

## Scope

The public release reproduces the analyses that can be validated from the available final
connectivity matrices and behavioral tables:

1. primary voxelwise partial-correlation maps;
2. cross-dataset spatial correspondence;
3. Dataset 1-to-Dataset 2 prediction;
4. NRS item-level map comparison;
5. no-covariate sensitivity analysis underlying Supplementary Figure S1;
6. anatomy-conditioned permutation analysis;
7. sample-size-weighted combined map;
8. two-sided whole-brain maximum-statistic family-wise error correction;
9. NIfTI export of group-level vectors.

TFCE, split-half validation and five-fold cross-validation were exploratory historical analyses
and are not part of the revised manuscript.

## Data access

Participant-level and participant-derived input files are not included. Access is subject to
the policies, legal requirements and ethical approvals of the institutions responsible for the
two original cohorts. Frozen non-identifiable outputs are supplied for transparent review.

## Citation and versioning

Use the exact tagged release archived for the journal revision. Update `CITATION.cff` with the
final repository DOI after Zenodo archiving.

## License

MIT License applies to code only. Data, historical group-level maps and third-party resources
remain subject to their respective permissions.
