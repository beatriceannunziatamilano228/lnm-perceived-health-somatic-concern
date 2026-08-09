# FSL cluster localization workflow

This directory formalizes the terminal workflow used to identify and report positive clusters from the three group maps.

## Historical parameters

The maps use a 2 mm isotropic grid, so each voxel has a volume of 8 mm³. The manuscript minimum extent of greater than 320 mm³ corresponds to 40 voxels. The archived terminal output establishes the map-specific positive r thresholds and a 26-voxel connectivity definition:

| Map | r threshold | Minimum extent |
|---|---:|---:|
| Dataset 1 | 0.310 | 40 voxels / 320 mm³ |
| Dataset 2 | 0.235 | 40 voxels / 320 mm³ |
| Historical combined map | 0.190 | 40 voxels / 320 mm³ |

The minimum values shown on figure colour bars can be slightly higher than the command threshold because the bars reflect values remaining after thresholding and extent filtering.

## Run

FSL 6.0.6 or later is recommended. From the repository root:

```bash
bash cluster_localization/run_fsl_cluster_localization.sh
python cluster_localization/query_harvard_oxford_labels.py   --cluster-dir results/cluster_localization
```

The shell script calls `fsl-cluster` directly when available, uses 26-connectivity, reports millimetre coordinates, and writes cluster-index, thresholded-map and local-maximum outputs. The labelling script calls FSL `atlasq` for both the Harvard-Oxford cortical and subcortical structural atlases.

## Frozen tables

`results/cluster_localization/` contains coordinate tables reconstructed from the public group maps using the exact historical thresholds, 26-connectivity and 40-voxel extent criterion. These tables reproduce the cluster sizes and peaks visible in the archived terminal output. Atlas labels should be interpreted as probabilistic anatomical aids and were reviewed against the images used for the manuscript.
