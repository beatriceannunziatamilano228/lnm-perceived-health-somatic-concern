# Restricted input data schema

No participant-level data are distributed in this repository.

The public scripts expect the following files in a user-specified restricted-data directory.

## 1. `Corbettaconnectivitymaps101.mat`

- Variable: `CorbettaFinal`
- Shape: `285903 × 101`
- Rows: voxels in the common analysis mask
- Columns: Dataset 1 participants, row-aligned with `question2_controls`

## 2. `CorbettaSF02_controls.mat`

- Variable: `question2_controls`
- Shape: `101 × 28`
- Column 1: SF-02 outcome
- Columns 2–11: SF-03 to SF-12
- Columns 12–27: SF-17 to SF-32
- Column 28: lesion size

The final primary model therefore includes 27 covariates: 26 SF-36 items and lesion size.

## 3. `GrafmanmapsnoNaN.mat`

- Variable: `Grafman2.mapsCopy`
- Shape: `285903 × 181`
- Rows: voxels in the same common analysis mask
- Columns: Dataset 2 participants, row-aligned with `NBRwithoutNaN_NBR02column1`

## 4. `NBRwithoutNaN_NBR02column1.mat`

- Variable: `NBRwithoutNaN_NBR02column1`
- Shape: `181 × 28`
- Column 1: NRS-02 somatic concern
- Column 2: NRS-01
- Columns 3–27: NRS-03 to NRS-27
- Column 28: lesion size

For the primary NRS-02 map, columns 2–28 are covariates. For the item-level analysis, the first 27 columns are restored to NRS-01…NRS-27 order, and each item is analysed while adjusting for the other 26 items and lesion size.

## 5. Optional `Corbetta_additional_controls.mat`

This restricted, identifier-free file is expected only by `run_dataset1_additional_controls.m`. Variables:

- `nih_total`
- `fim_total`
- `fam_tot`
- `rnltotal`
- `gdss_score`
- `gender_code`
- `sip_body`, `sip_social`, `sip_mob`, `sip_com`, `sip_emo`, `sip_house`, `sip_alert`, `sip_amb`, `sip_psychosoc`, `sip_physical`

Each vector must contain 101 rows in the same participant order as `CorbettaFinal` and `question2_controls`. Participant identifiers must not be included in the public repository.

## Group-level NIfTI reference

`results/group_maps/combined_lesion_associated_map_historical.nii.gz` provides the 91 × 109 × 91, 2-mm MNI grid and the 285,903-voxel analysis mask used by the modular MATLAB NIfTI-writing helper.
