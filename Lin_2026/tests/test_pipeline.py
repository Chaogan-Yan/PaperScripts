from __future__ import annotations

import subprocess
import shutil
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from roi_analysis import (  # noqa: E402
    AMYGDALA_ROIS,
    FALFF_NAME_MAP,
    HIPPOCAMPAL_ROIS,
)


class SyntheticPipelineTests(unittest.TestCase):
    def test_end_to_end_cli_and_matching_glm_schema(self) -> None:
        rng = np.random.default_rng(20260830)
        n_hc = 40
        n_mdd = 40
        n = n_hc + n_mdd
        subtype = np.r_[np.ones(n_hc, dtype=int), np.full(n_mdd, 2, dtype=int)]
        subject_ids = [f"synthetic-{i:03d}" for i in range(n)]
        age = rng.uniform(18, 65, n)
        sex = rng.choice([1, 2], n)
        etiv = rng.normal(1_450_000, 120_000, n)
        site = rng.choice(["site-A", "site-B"], n)
        shared = pd.DataFrame(
            {
                "SubID": subject_ids,
                "Sex": sex,
                "Age": age,
                "subtype": subtype,
                "EstimatedTotalIntraCranialVol": etiv,
                "Site": site,
            }
        )

        hippocampus = shared.copy()
        for index, roi in enumerate(HIPPOCAMPAL_ROIS):
            hippocampus[roi] = 100 + index + 0.5 * age + rng.normal(0, 2, n)

        amygdala = shared.copy()
        for index, roi in enumerate(AMYGDALA_ROIS):
            amygdala[roi] = 50 + index + 0.2 * age + rng.normal(0, 1, n)

        mdd_ids = subject_ids[n_hc:]
        clinical = pd.DataFrame(
            {
                "SubID": mdd_ids,
                "age_at_onset_years": np.tile([18.0, 30.0], n_mdd // 2),
                "illness_duration_months": rng.uniform(6, 120, n_mdd),
                "first_episode": np.tile([0, 1], n_mdd // 2),
                "currently_medicated": np.tile([1, 0], n_mdd // 2),
            }
        )

        motion = pd.DataFrame(
            {
                "SubID": subject_ids,
                "Age": age,
                "Sex": sex,
                "Motion": rng.uniform(0.01, 0.20, n),
            }
        )
        falff = pd.DataFrame({"SubID": subject_ids})
        for index, column in enumerate(FALFF_NAME_MAP.values()):
            falff[column] = 0.5 + index / 100 + rng.normal(0, 0.02, n)

        temporary = ROOT / "tests" / "_synthetic_runtime"
        if temporary.exists():
            shutil.rmtree(temporary)
        temporary.mkdir()
        try:
            hippocampal_file = temporary / "hippocampal.csv"
            amygdala_file = temporary / "amygdala.csv"
            clinical_file = temporary / "clinical.csv"
            motion_file = temporary / "motion.csv"
            falff_file = temporary / "falff.csv"
            output = temporary / "output"
            hippocampus.to_csv(hippocampal_file, index=False)
            amygdala.to_csv(amygdala_file, index=False)
            clinical.to_csv(clinical_file, index=False)
            motion.to_csv(motion_file, index=False)
            falff.to_csv(falff_file, index=False)

            completed = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "code" / "roi_analysis.py"),
                    "--hippocampal-input",
                    str(hippocampal_file),
                    "--amygdala-input",
                    str(amygdala_file),
                    "--clinical-input",
                    str(clinical_file),
                    "--falff-input",
                    str(falff_file),
                    "--falff-covariates",
                    str(motion_file),
                    "--output",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("Analysis outputs written to", completed.stdout)

            primary = pd.read_csv(output / "volume_primary.csv")
            falff_result = pd.read_csv(output / "falff_mdd_vs_hc.csv")
            direct = pd.read_csv(output / "volume_within_mdd.csv")
            subgroup = pd.read_csv(output / "volume_subgroup_vs_hc.csv")
            self.assertEqual(len(primary), 42)
            self.assertEqual(len(falff_result), 13)
            self.assertEqual(len(direct), 39)
            self.assertEqual(len(subgroup), 78)
            self.assertEqual(list(falff_result.columns), list(primary.columns))
            self.assertEqual(set(primary["n"]), {n})
            self.assertEqual(set(falff_result["n"]), {n})
        finally:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
