#!/usr/bin/env python3
"""Create a SHA-256 manifest for release files."""

import argparse
import hashlib
import os
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


parser = argparse.ArgumentParser(
    description="Write SHA-256 fingerprints for files and directory contents."
)
parser.add_argument("paths", nargs="+", type=Path)
parser.add_argument("--output", type=Path, default=Path("SHA256SUMS.txt"))
args = parser.parse_args()

files = []
for path in args.paths:
    files.extend(item for item in path.rglob("*") if item.is_file()) if path.is_dir() else files.append(path)

with args.output.open("w", encoding="utf-8", newline="\n") as handle:
    for path in sorted(set(item.resolve() for item in files)):
        handle.write(
            f"{sha256(path)}  {Path(os.path.relpath(path, Path.cwd())).as_posix()}\n"
        )

print(f"Wrote {len(set(files))} checksums to {args.output}")
