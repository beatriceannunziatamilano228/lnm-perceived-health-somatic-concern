# Statistical implementation summary

Voxelwise lesion-connectivity associations were estimated using partial Pearson correlation
after adjustment for the prespecified covariates. Cross-dataset prediction was assessed with
Spearman correlation. Spatial correspondence between the independently estimated maps was
quantified with Pearson correlation, with Fisher-transformed map correlation retained as a
reproducibility check.

The reviewer-requested anatomy-conditioned analysis held the Dataset 1 map fixed and
preserved Dataset 2 lesion-connectivity maps, lesion distribution, lesion size and covariates
while permuting the residualized Dataset 2 primary outcome. The reported two-sided empirical
p-value used 100,000 permutations and the plus-one correction.

Whole-brain family-wise error correction used a two-sided maximum-absolute-statistic
Freedman-Lane procedure with 10,000 permutations. The maximum absolute sample-size-
weighted partial-correlation coefficient was retained at every iteration.
