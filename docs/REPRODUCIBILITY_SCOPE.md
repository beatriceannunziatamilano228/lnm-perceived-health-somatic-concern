# Reproducibility scope

## Authoritative analysis

The Python files in `analysis/` are the authoritative revision code. Frozen outputs generated
from those analyses are under `results/`.

## Historical code

The MATLAB file in `provenance/historical_matlab/` is retained only to document the
historical workflow. It is not used as evidence for the final corrected inference.

## Analyses not included

The following historical exploratory analyses are not reported in the revised manuscript and
are excluded from the validated pipeline: TFCE, split-half validation, and five-fold
cross-validation.

## Controls requiring restricted inputs

Sex and disability-control analyses and the external anosognosia comparison require
additional restricted or externally sourced inputs. Their availability and reporting should be
documented separately in the final manuscript and repository release notes.
