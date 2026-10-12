from __future__ import annotations

import csv
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / "tests" / "_synthetic_runtime"


class HashManifestTests(unittest.TestCase):
    def tearDown(self) -> None:
        if RUNTIME.exists():
            for path in sorted(RUNTIME.rglob("*"), reverse=True):
                if path.is_file():
                    path.unlink()
                elif path.is_dir():
                    path.rmdir()
            RUNTIME.rmdir()

    def test_output_file_is_not_hashed_into_itself(self) -> None:
        RUNTIME.mkdir(parents=True, exist_ok=True)
        artifact = RUNTIME / "artifact.txt"
        artifact.write_text("synthetic artifact\n", encoding="utf-8")
        manifest = RUNTIME / "manifest.csv"
        manifest.write_text("old manifest\n", encoding="utf-8")

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "code" / "hash_manifest.py"),
                str(RUNTIME),
                "--base",
                str(RUNTIME),
                "--output",
                str(manifest),
            ],
            cwd=ROOT,
            check=True,
        )

        with manifest.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([row["path"] for row in rows], ["artifact.txt"])
        self.assertEqual(len(rows[0]["sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
