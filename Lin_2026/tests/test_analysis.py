from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from roi_analysis import fdr_by_family, fit_ols  # noqa: E402


class AnalysisUnitTests(unittest.TestCase):
    def synthetic_frame(self) -> pd.DataFrame:
        group = np.tile([0.0, 1.0], 20)
        age = np.linspace(20.0, 60.0, group.size)
        sex = np.tile([1.0, 2.0, 2.0, 1.0], 10)
        covariate = np.linspace(-1.0, 1.0, group.size) ** 2
        outcome = 3.0 * group + 0.2 * age - 0.5 * sex + 0.3 * covariate
        outcome += np.sin(np.arange(group.size)) * 0.1
        return pd.DataFrame(
            {
                "outcome": outcome,
                "group": group,
                "age": age,
                "sex": sex,
                "covariate": covariate,
            }
        )

    def test_complete_case_fit_and_output_schema(self) -> None:
        frame = self.synthetic_frame()
        frame.loc[0, "outcome"] = np.nan
        result = fit_ols(frame, "outcome", "group", ["age", "sex", "covariate"])

        self.assertEqual(result["n"], 39)
        self.assertEqual(result["n_group_1"], 20)
        self.assertEqual(result["n_group_0"], 19)
        self.assertEqual(
            set(result),
            {
                "n",
                "n_group_1",
                "n_group_0",
                "beta",
                "residual_sd",
                "adjusted_cohens_d",
                "ci95_low",
                "ci95_high",
                "t_value",
                "df_residual",
                "p_value",
            },
        )

    def test_rank_deficient_design_fails(self) -> None:
        frame = self.synthetic_frame()
        frame["duplicate_age"] = frame["age"]
        with self.assertRaisesRegex(ValueError, "rank deficient"):
            fit_ols(frame, "outcome", "group", ["age", "duplicate_age"])

    def test_zero_residual_variance_fails(self) -> None:
        frame = self.synthetic_frame()
        frame["outcome"] = 3.0 * frame["group"] + 0.2 * frame["age"]
        with self.assertRaisesRegex(ValueError, "residual standard deviation"):
            fit_ols(frame, "outcome", "group", ["age"])

    def test_fdr_is_applied_separately_by_family(self) -> None:
        frame = pd.DataFrame(
            {
                "family": ["hippocampus", "hippocampus", "amygdala", "amygdala"],
                "p_value": [0.01, 0.04, 0.02, 0.90],
            }
        )
        result = fdr_by_family(frame)
        expected = [0.02, 0.04, 0.04, 0.90]
        np.testing.assert_allclose(result["q_fdr_bh"], expected)


if __name__ == "__main__":
    unittest.main()
