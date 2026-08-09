#!/usr/bin/env python3
"""Fail if the public tree contains obvious private inputs or personal paths."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = {".xlsx", ".xls", ".mat", ".log"}
FORBIDDEN_TEXT = ["/" + "Users/", "\\" + "Users\\", "MATLAB" + "-Drive"]
FORBIDDEN_NAME_TOKENS = {
    "internal",
    "depressionscores",
    "fcs_demographics",
    "participant_id",
    "subject_id",
}
errors: list[str] = []

for path in ROOT.rglob("*"):
    if not path.is_file() or ".git" in path.parts or path.resolve() == Path(__file__).resolve():
        continue
    rel = path.relative_to(ROOT)
    rel_lower = rel.as_posix().lower()

    if path.suffix.lower() in FORBIDDEN_SUFFIXES:
        errors.append(f"Forbidden file type: {rel}")
    for token in FORBIDDEN_NAME_TOKENS:
        if token in rel_lower:
            errors.append(f"Forbidden private-name token: {rel}")

    # NIfTI and NumPy files are allowed only as explicitly public group/frozen results.
    if rel_lower.endswith((".nii", ".nii.gz", ".npy")) and rel.parts[0] not in {"results", "external_maps"}:
        errors.append(f"Imaging/array file outside results/: {rel}")

    if path.suffix.lower() in {".md", ".txt", ".csv", ".py", ".m", ".yml", ".yaml", ".cff", ".json"}:
        text = path.read_text(errors="ignore")
        for token in FORBIDDEN_TEXT:
            if token in text:
                errors.append(f"Forbidden personal-path token: {rel}")

if errors:
    print("\n".join(sorted(set(errors))))
    sys.exit(1)
print("Repository privacy/path audit passed.")
