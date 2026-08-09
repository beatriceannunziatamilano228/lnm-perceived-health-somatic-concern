# MATLAB manuscript-analysis implementation

These modular MATLAB R2023a functions restate the final analysis logic from the recovered working scripts and the restricted input schemas. They are designed for documented reruns by investigators who hold authorized access to the participant-level inputs.

## Entry point

```matlab
dataDir = "/secure/path/to/restricted_inputs";
outDir = fullfile(pwd,"outputs");
run_all_core_analyses(dataDir,outDir,10000,20260724);
```

## Analyses

1. `run_primary_maps.m`
2. `run_cross_dataset_correspondence.m`
3. `run_cross_dataset_prediction.m`
4. `run_nrs_item_level_analysis.m`
   - maps each NRS item while controlling for the other 26 items and lesion size;
   - also maps each item while controlling for lesion size only, matching the manuscript sensitivity analysis performed without the remaining NRS items.
5. `run_historical_combined_map.m`
6. `run_dataset1_additional_controls.m`
7. `run_supplementary_no_covariates.m`

Cluster localization, external-map comparison, reviewer-requested Python analyses, and SimNIBS modelling are documented in separate top-level directories.

## Important boundary

The scripts begin from precomputed participant-level lesion-connectivity maps. The upstream normative-connectome infrastructure and participant-derived inputs are not redistributed.
