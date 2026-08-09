# Frozen non-identifiable results

## `group_maps/`

- `dataset1_primary_partial_r.nii.gz`: covariate-adjusted Dataset 1 map.
- `dataset2_primary_partial_r.nii.gz`: covariate-adjusted Dataset 2 map.
- `dataset1_no_covariates_r.nii.gz` and `dataset2_no_covariates_r.nii.gz`: Supplementary Figure S1 sensitivity maps.
- `combined_lesion_associated_map_historical.nii.gz`: historical combined group map used as input to target identification.
- `target_connectivity_map.nii.gz`: group-level output of the Center for Brain Circuit Therapeutics normative-connectome targeting pipeline.

## `cluster_localization/`

Frozen cluster geometry tables reconstructed with the exact archived thresholds, 26-connectivity and 40-voxel extent criterion. These tables reproduce the cluster sizes and peak coordinates visible in the historical FSL terminal output.

## `verified/`

Summary tables and arrays independently reconstructed from the currently available final analysis matrices. The item-level table includes both fully adjusted and lesion-size-only comparisons. No participant identifiers or participant-level clinical records are included.

## `revision/`

Reviewer-requested and corrected inferential outputs, including frozen permutation null distributions and corrected whole-brain maps.
