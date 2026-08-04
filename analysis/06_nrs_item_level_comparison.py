#!/usr/bin/env python3
"""Reconstruct the covariate-adjusted and unadjusted NRS item-level comparisons."""

from pathlib import Path
import argparse
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from analysis_utils import (
    load_inputs,
    complete_case,
    partial_correlation_map,
    correlation_map,
    spatial_correlation,
)


ITEM_NAMES = {
    1: "comprehension deficit",
    2: "somatic concern",
    3: "guilt feelings",
    4: "inaccurate insight and self-appraisal",
    5: "expressive deficit",
    6: "tension",
    7: "conceptual disorganization",
    8: "inattention/reduced alertness",
    9: "blunted affect",
    10: "disorientation",
    11: "speech articulation defect",
    12: "lability of mood",
    13: "agitation",
    14: "hostility/uncooperativeness",
    15: "depressive mood",
    16: "emotional withdrawal",
    17: "unusual thought content",
    18: "hallucinatory behaviour",
    19: "motor retardation",
    20: "fatigability",
    21: "memory deficit",
    22: "excitement",
    23: "poor planning",
    24: "anxiety",
    25: "suspiciousness",
    26: "disinhibition",
    27: "decreased initiative/motivation",
}

# Supplied Dataset 2 table is ordered as NRS-02, NRS-01, NRS-03 ... NRS-27,
# followed by lesion size.
COLUMN_TO_ITEM = [2, 1] + list(range(3, 28))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--chunk-size", type=int, default=10_000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    c_maps, c_beh, g_maps, g_beh = load_inputs(args.data_dir)
    c_valid, c_y, c_cov = complete_case(c_beh)
    c_primary = partial_correlation_map(
        c_maps[:, c_valid], c_y, c_cov, args.chunk_size
    )

    rows = []
    for column_index, item_number in enumerate(COLUMN_TO_ITEM):
        outcome = g_beh[:, column_index]
        other_nrs_columns = [j for j in range(27) if j != column_index]
        covariates = g_beh[:, other_nrs_columns + [27]]

        adjusted_valid = np.isfinite(outcome) & np.isfinite(covariates).all(axis=1)
        adjusted_map = partial_correlation_map(
            g_maps[:, adjusted_valid],
            outcome[adjusted_valid],
            covariates[adjusted_valid],
            args.chunk_size,
        )
        adjusted_n = int(adjusted_valid.sum())
        adjusted_r, adjusted_voxels = spatial_correlation(
            c_primary, adjusted_map, fisher_z=True
        )

        unadjusted_valid = np.isfinite(outcome)
        unadjusted_map = correlation_map(
            g_maps[:, unadjusted_valid], outcome[unadjusted_valid], args.chunk_size
        )
        unadjusted_n = int(unadjusted_valid.sum())
        unadjusted_r, unadjusted_voxels = spatial_correlation(
            c_primary, unadjusted_map, fisher_z=True
        )

        rows.append(
            {
                "NRS_item_number": item_number,
                "NRS_item_name": ITEM_NAMES[item_number],
                "adjusted_complete_case_n": adjusted_n,
                "adjusted_spatial_r": adjusted_r,
                "adjusted_voxels": adjusted_voxels,
                "unadjusted_available_n": unadjusted_n,
                "unadjusted_spatial_r": unadjusted_r,
                "unadjusted_voxels": unadjusted_voxels,
                "is_primary_outcome": item_number == 2,
            }
        )

    table = pd.DataFrame(rows)
    table["adjusted_rank_descending"] = (
        table["adjusted_spatial_r"].rank(method="min", ascending=False).astype(int)
    )
    table = table.sort_values("adjusted_spatial_r", ascending=False)
    table.to_csv(args.output_dir / "nrs_item_level_correlations.csv", index=False)

    primary = table.loc[table["NRS_item_number"] == 2].iloc[0]
    summary = {
        "common_adjusted_complete_case_n": int(primary["adjusted_complete_case_n"]),
        "primary_item": "NRS-02 somatic concern",
        "primary_adjusted_spatial_r": float(primary["adjusted_spatial_r"]),
        "primary_adjusted_rank": int(primary["adjusted_rank_descending"]),
        "interpretation": (
            "Descriptive item-level ranking. NRS-02 has the largest adjusted spatial "
            "correlation among the 27 item maps; no inferential test of differences "
            "between item correlations is performed by this script."
        ),
    }
    with open(args.output_dir / "nrs_item_level_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
