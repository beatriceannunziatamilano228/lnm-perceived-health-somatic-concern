#!/usr/bin/env python3
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
lines=[]
for p in sorted(root.rglob('*')):
    if p.is_file() and '.git' not in p.parts and p.name!='MANIFEST.sha256':
        h=hashlib.sha256(p.read_bytes()).hexdigest()
        lines.append(f"{h}  {p.relative_to(root).as_posix()}")
(root/'MANIFEST.sha256').write_text('\n'.join(lines)+'\n')
