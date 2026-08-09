#!/usr/bin/env python3
"""Download an authoritative NeuroVault image and save its API metadata."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse
import requests

ALLOWED = {795500, 795501, 795502}

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--image-id", type=int, default=795501, choices=sorted(ALLOWED))
    p.add_argument("--output-dir", type=Path, default=Path("external_maps/downloaded"))
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    api_url = f"https://neurovault.org/api/images/{args.image_id}/"
    response = requests.get(api_url, headers={"Accept": "application/json"}, timeout=60)
    response.raise_for_status()
    metadata = response.json()
    file_url = str(metadata["file"]).replace("http://", "https://", 1)
    filename = Path(urlparse(file_url).path).name or f"neurovault_{args.image_id}.nii.gz"

    image_response = requests.get(file_url, stream=True, timeout=120)
    image_response.raise_for_status()
    output = args.output_dir / filename
    sha = hashlib.sha256()
    with output.open("wb") as f:
        for chunk in image_response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)
                sha.update(chunk)

    metadata["downloaded_file"] = output.name
    metadata["sha256"] = sha.hexdigest()
    metadata["persistent_identifier"] = f"https://identifiers.org/neurovault.image:{args.image_id}"
    (args.output_dir / f"neurovault_{args.image_id}_metadata.json").write_text(
        json.dumps(metadata, indent=2)
    )
    print(output)
    print(f"SHA-256: {sha.hexdigest()}")

if __name__ == "__main__":
    main()
