#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/results/cluster_localization/fsl_outputs"
mkdir -p "$OUT"

if [[ -n "${FSLDIR:-}" && -x "$FSLDIR/bin/fsl-cluster" ]]; then
  CLUSTER="$FSLDIR/bin/fsl-cluster"
elif command -v fsl-cluster >/dev/null 2>&1; then
  CLUSTER="$(command -v fsl-cluster)"
elif command -v cluster >/dev/null 2>&1; then
  CLUSTER="$(command -v cluster)"
else
  echo "FSL cluster/fsl-cluster was not found." >&2
  exit 1
fi

run_cluster () {
  local id="$1" input="$2" threshold="$3"
  "$CLUSTER" \
    --in="$input" \
    --thresh="$threshold" \
    --connectivity=26 \
    --minextent=40 \
    --mm \
    --oindex="$OUT/${id}_cluster_index.nii.gz" \
    --othresh="$OUT/${id}_extent_thresholded.nii.gz" \
    --olmax="$OUT/${id}_local_maxima.tsv" \
    > "$OUT/${id}_clusters.tsv"
}

run_cluster dataset1 \
  "$ROOT/results/group_maps/dataset1_primary_partial_r.nii.gz" 0.310
run_cluster dataset2 \
  "$ROOT/results/group_maps/dataset2_primary_partial_r.nii.gz" 0.235
run_cluster combined \
  "$ROOT/results/group_maps/combined_lesion_associated_map_historical.nii.gz" 0.190

echo "FSL cluster outputs written to: $OUT"
