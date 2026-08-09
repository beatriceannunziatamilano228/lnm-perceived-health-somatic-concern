# Manuscript protocol coverage

| Manuscript component | Public implementation/output | Coverage status |
|---|---|---|
| Dataset 1 SF-02 map | `matlab/manuscript_analysis/run_primary_maps.m`; public group map | Covered; restricted participant inputs required for rerun |
| Dataset 2 NRS-02 map | same | Covered; complete-case N documented |
| Cross-dataset spatial correspondence | `run_cross_dataset_correspondence.m` | Observed statistic covered; historical null seed/output not preserved |
| Dataset 1 → Dataset 2 prediction | `run_cross_dataset_prediction.m` | Covered and frozen result supplied |
| Figure 4 adjusted item maps | `run_nrs_item_level_analysis.m` | Covered; all other NRS items and lesion size included |
| Figure 4 sensitivity without remaining NRS items | same | Covered as lesion-size-only item maps; verified summary table supplied |
| Historical combined map | `run_historical_combined_map.m` | Covered; nominal 101/181 Fisher-z weighting |
| Cluster localization | `cluster_localization/` | Covered using exact archived thresholds, 26-connectivity and 40-voxel extent |
| Harvard-Oxford labels | `query_harvard_oxford_labels.py` | Covered; requires local FSL atlases |
| Historical Westfall–Young code | `matlab/historical/` | Preserved for provenance; revised inference is separate |
| Anatomy-conditioned reviewer test | `python/reviewer_analyses/01_anatomy_conditioned_permutation.py` | Covered; frozen 100,000-permutation null supplied |
| Corrected maximum statistic | `python/reviewer_analyses/02_corrected_maximum_statistic.py` | Covered; frozen 10,000-permutation null and thresholded map supplied |
| Additional Dataset 1 controls | `run_dataset1_additional_controls.m` | Covered; de-identified restricted input required |
| No-covariate sensitivity | `run_supplementary_no_covariates.m` | Covered; group maps supplied |
| Anosognosia comparison | `external_maps/`; `python/external_map_comparisons/` | Public comparator and code covered; historical p=0.93 not independently verified |
| Target-connectivity map | public input/output maps and provenance statement | Output covered; upstream CBCT connectome infrastructure not redistributed |
| SimNIBS Figure 6 | `simnibs/` | Placement matrices and workflow covered; original software logs are not published because they contained personal paths |
