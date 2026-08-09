# Manuscript alignment notes

The revised manuscript should use the following verified values and terminology:

- Dataset 1 primary analytic N: **101**.
- Dataset 2 cohort N: **181**; complete-case primary map N: **173**.
- Cross-dataset spatial correspondence: **r = 0.59**; describe as moderate and as partial correspondence, not direct replication.
- Reviewer-requested anatomy-conditioned test: **two-sided p = 0.0405**, 100,000 permutations.
- Forward participant-level prediction: **Spearman rho = 0.233, p = 0.00162, N = 181**.
- Figure 4: NRS-02 is rank 1 among 27 adjusted item maps; describe as a descriptive item-level comparison rather than a formal test of specificity.
- Supplementary Figure S1: with-versus-without-covariate spatial r values are **0.826** (Dataset 1) and **0.833** (Dataset 2) from the current final matrices.
- Corrected two-sided maximum-statistic result: threshold **|r| = 0.28064**, global **pFWE = 0.00580**, 35 positive significant voxels.

Reconstructed Dataset 1 additional-control similarities, comparing base and augmented maps within the same complete-case subset:

- NIHSS: N=89, r=0.99988
- FIM: N=101, r=0.99484
- FAM: N=101, r=0.99635
- RNLI: N=100, r=0.99489
- GDS: N=100, r=0.99876
- SIP domains jointly: N=101, r=0.95377
- Sex: N=101, r=0.99992

These verified values supersede earlier approximate working-document values if the corresponding controls are retained in the final manuscript.


## Additional alignment points introduced in this release

- Historical FSL cluster thresholds were Dataset 1 `r = 0.310`, Dataset 2 `r = 0.235`, and combined map `r = 0.190`, with 26-connectivity and a 40-voxel (320 mm³) extent criterion.
- In the lesion-size-only NRS item analysis, NRS-02 is not the highest-ranked item; the claim of strongest item-level correspondence applies to the fully adjusted maps. The manuscript should not describe the lesion-size-only ranking as confirming specificity.
- The Kletenik cross-modal comparator is NeuroVault image 795501. Its source and the corresponding anosognosia spatial-comparison workflow are documented in the repository.
