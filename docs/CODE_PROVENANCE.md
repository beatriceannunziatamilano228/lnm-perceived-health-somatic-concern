# Code provenance

## Original analysis language

The original analyses were conducted in MATLAB R2023a. The repository includes the recovered master working script (`SF3_32_Final.m`), the recovered NRS item-level script (`predictinggrafman.m`), and the recovered Westfall-Young working script. Public copies are path-sanitized and retained under `matlab/historical/`.

## MATLAB manuscript-analysis implementation

`matlab/manuscript_analysis/` restates the final manuscript models without relying on an interactive workspace. The Dataset 1 primary model uses SF-02 as outcome and 27 covariates: SF-03–SF-12, SF-17–SF-32, and lesion size. The Dataset 2 primary model uses NRS-02 as outcome and the remaining 26 NRS items plus lesion size.

## Peer-review analyses

The anatomy-conditioned permutation and corrected maximum-statistic analyses were added during peer review and implemented in Python. They are separated from the historical MATLAB code and accompanied by frozen outputs.

## Stimulation-target map

The authors supplied the final combined group-level map as input to a normative-connectome targeting pipeline maintained by the Center for Brain Circuit Therapeutics. The repository provides the input combined map, the resulting group-level target-connectivity map, target coordinates, and SimNIBS modelling parameters. The upstream normative connectome and associated processing infrastructure are not redistributed.

## Participant identifiers and additional controls

The 101 Dataset 1 analysis rows were internally aligned to the clinical database using exact 27-variable SF-36 response profiles. That internal matching table and all participant identifiers are excluded. The public repository contains only code, input schemas and aggregate map-similarity results for FIM, FAM, NIHSS, RNLI, GDS, SIP domains and sex.

## Verification boundary

The modular MATLAB refactor was checked against the recovered working scripts, the final four-matrix input schema, and independently reconstructed frozen outputs. The modular MATLAB files were not re-executed in the release-building environment because MATLAB and the restricted participant-level inputs were not available there. This repository does not claim that the refactored files are byte-identical to the historical interactive workspace; instead, it preserves the historical scripts separately and documents the final analysis logic in modular form.


## FSL cluster localization

The original peak localization was performed interactively in the terminal with FSL `cluster`. Archived terminal output established the exact positive thresholds (0.310, 0.235 and 0.190), 26-connectivity and 40-voxel minimum extent. This release converts that workflow into a versioned shell script and supplies frozen geometry tables matching the terminal output.

## Kletenik external map

The cross-modal awareness/anosognosia comparator is sourced from NeuroVault collection 13792, image 795501, linked to DOI 10.1002/ana.26709. The authoritative image is downloaded from NeuroVault rather than duplicated. The historical p-value was not accompanied by a preserved null distribution; a revised, explicit permutation workflow is therefore supplied separately.
