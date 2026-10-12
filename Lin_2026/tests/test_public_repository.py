from __future__ import annotations

import csv
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_DIRECTORIES = {".git", "__pycache__", "_synthetic_runtime"}
PRIVATE_EXTENSIONS = {
    ".docx",
    ".xls",
    ".xlsx",
    ".nii",
    ".gz",
    ".mat",
    ".npy",
    ".npz",
    ".pkl",
    ".sav",
    ".pdf",
    ".png",
    ".jpg",
    ".jpeg",
    ".bmp",
    ".tif",
    ".tiff",
}
TEXT_EXTENSIONS = {"", ".csv", ".gitignore", ".m", ".md", ".py", ".txt", ".yml", ".yaml"}
USER_ROOT_NAMES = ("Users", "home")
USER_ROOT_PATTERN = "|".join(USER_ROOT_NAMES)
PERSONAL_PATH = re.compile(
    rf"(?:[A-Za-z]:[\\/](?:{USER_ROOT_PATTERN})[\\/]|/(?:{USER_ROOT_PATTERN})/)"
)
LOCAL_MOUNT_PATH = re.compile(r"/(?:mnt|media)/[^/\s]+/", re.IGNORECASE)
WINDOWS_ABSOLUTE_PATH = re.compile(r"(?<![A-Za-z0-9_])[A-Za-z]:[\\/]")


def public_files() -> list[Path]:
    return sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file() and not any(part in EXCLUDED_DIRECTORIES for part in path.parts)
    )


class PublicRepositoryTests(unittest.TestCase):
    def test_fc_seed_manifest_is_fixed_to_the_reported_four_seeds(self) -> None:
        manifest_path = ROOT / "config" / "fc_seed_manifest.csv"
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(
            [row["SeedName"] for row in rows],
            [
                "Right_Basal_nucleus",
                "Left_Accessory_Basal_nucleus",
                "Right_Accessory_Basal_nucleus",
                "Right_Cortical_nucleus",
            ],
        )
        self.assertEqual({row["LabelID"] for row in rows}, {"7008"})
        self.assertTrue(all(not Path(row["StatFile"]).is_absolute() for row in rows))

    def test_no_private_data_or_generated_results(self) -> None:
        forbidden_directories = {"data", "derived", "private_data", "results", "analysis"}
        violations = []
        for path in public_files():
            relative = path.relative_to(ROOT)
            if any(part.lower() in forbidden_directories for part in relative.parts[:-1]):
                violations.append(str(relative))
            if path.suffix.lower() in PRIVATE_EXTENSIONS:
                violations.append(str(relative))
            allowed_public_csv = relative.as_posix() == "config/fc_seed_manifest.csv" or (
                relative.parts[0] == "config" and path.name.endswith(".example.csv")
            )
            if path.suffix.lower() in {".csv", ".tsv"} and not allowed_public_csv:
                violations.append(str(relative))
            if path.stat().st_size > 1_000_000:
                violations.append(f"{relative} (>1 MB)")
        self.assertEqual(violations, [])

    def test_no_personal_absolute_paths_in_text_files(self) -> None:
        violations = []
        for path in public_files():
            if path.suffix.lower() not in TEXT_EXTENSIONS:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            if (
                PERSONAL_PATH.search(text)
                or LOCAL_MOUNT_PATH.search(text)
                or WINDOWS_ABSOLUTE_PATH.search(text)
            ):
                violations.append(str(path.relative_to(ROOT)))
        self.assertEqual(violations, [])


if __name__ == "__main__":
    unittest.main()
