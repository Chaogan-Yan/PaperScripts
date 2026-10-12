"""Create a streaming SHA-256 manifest for reproducibility artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--base", type=Path, default=Path.cwd())
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base = args.base.resolve()
    output = args.output.resolve()
    files: list[Path] = []
    for supplied in args.paths:
        path = supplied.resolve()
        if path.is_dir():
            files.extend(sorted(item for item in path.rglob("*") if item.is_file()))
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(path)

    # Never hash the manifest into itself when an existing output file is
    # located inside one of the supplied directories.
    files = [path for path in files if path != output]
    output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, Path]] = []
    for path in sorted(set(files)):
        try:
            display = path.relative_to(base).as_posix()
        except ValueError:
            display = path.name
        rows.append((display, path))
    display_names = [display for display, _ in rows]
    if len(display_names) != len(set(display_names)):
        raise ValueError(
            "Manifest paths are not unique relative to --base; choose a common base directory"
        )

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        for display, path in rows:
            writer.writerow({"path": display, "bytes": path.stat().st_size, "sha256": sha256(path)})


if __name__ == "__main__":
    main()
