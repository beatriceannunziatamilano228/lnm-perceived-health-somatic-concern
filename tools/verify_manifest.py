#!/usr/bin/env python3
"""Verify all SHA-256 entries in MANIFEST.sha256."""
from pathlib import Path
import hashlib
import sys

root = Path(__file__).resolve().parents[1]
manifest = root / "MANIFEST.sha256"
if not manifest.is_file():
    raise SystemExit("MANIFEST.sha256 is missing.")

errors = []
listed = set()
for line in manifest.read_text().splitlines():
    if not line.strip():
        continue
    digest, rel = line.split("  ", 1)
    path = root / rel
    listed.add(rel)
    if not path.is_file():
        errors.append(f"Missing file: {rel}")
        continue
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != digest:
        errors.append(f"Checksum mismatch: {rel}")

actual = {
    p.relative_to(root).as_posix()
    for p in root.rglob("*")
    if p.is_file() and ".git" not in p.parts and p.name != "MANIFEST.sha256"
}
for rel in sorted(actual - listed):
    errors.append(f"Unlisted file: {rel}")
for rel in sorted(listed - actual):
    errors.append(f"Manifest-only file: {rel}")

if errors:
    print("\n".join(errors))
    sys.exit(1)
print(f"Manifest verified for {len(listed)} files.")
