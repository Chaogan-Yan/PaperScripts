"""Reproduce the manuscript's ROI-level volume and optional fALFF GLMs.

Volume analyses run from the harmonized workbooks. fALFF analysis is enabled
only when a final post-ComBat ROI table is supplied explicitly.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


GROUP_COLS = ["Age", "Sex", "EstimatedTotalIntraCranialVol"]
CLINICAL_ONSET = "age_at_onset_years"
CLINICAL_DURATION = "illness_duration_months"
CLINICAL_EPISODE = "first_episode"
CLINICAL_MEDICATION = "currently_medicated"

LEGACY_CLINICAL_COLUMNS = {
    "ID": "SubID",
    "*发病年龄（岁）（缺失=留空)": CLINICAL_ONSET,
    "*病程（月）（缺失=留空)": CLINICAL_DURATION,
    "*是否首发(是=1，否=0，缺失=留空)": CLINICAL_EPISODE,
    "*是否正在用药(是=1，否=0，缺失=留空)": CLINICAL_MEDICATION,
}

HIPPOCAMPAL_ROIS = [
    "Left-parasubiculum",
    "Left-Hippocampal_tail",
    "Left-hippocampal_fissure",
    "Left-subiculum",
    "Left-CA1",
    "Left-presubiculum",
    "Left-molecular_layer",
    "Left-GC_ML_DG",
    "Left-CA3",
    "Left-CA4",
    "Left-fimbria",
    "Left-HATA",
    "Right-Hippocampal_tail",
    "Right-subiculum",
    "Right-hippocampal_fissure",
    "Right-CA1",
    "Right-presubiculum",
    "Right-parasubiculum",
    "Right-molecular_layer",
    "Right-GC_ML_DG",
    "Right-CA3",
    "Right-CA4",
    "Right-fimbria",
    "Right-HATA",
]

AMYGDALA_ROIS = [
    "Left-Lateral_nucleus",
    "Left-Basal_nucleus",
    "Left-Accessory_Basal_nucleus",
    "Left-Anterior_amygdaloid_area_AAA",
    "Left-Central_nucleus",
    "Left-Medial_nucleus",
    "Left-Cortical_nucleus",
    "Left-Corticoamygdaloid_transitio",
    "Left-Paralaminar_nucleus",
    "Right-Lateral_nucleus",
    "Right-Basal_nucleus",
    "Right-Accessory_Basal_nucleus",
    "Right-Anterior_amygdaloid_area_AAA",
    "Right-Central_nucleus",
    "Right-Medial_nucleus",
    "Right-Cortical_nucleus",
    "Right-Corticoamygdaloid_transitio",
    "Right-Paralaminar_nucleus",
]

STRUCTURAL_HIPPOCAMPAL_ROIS = [
    "Right-presubiculum",
    "Right-parasubiculum",
    "Right-HATA",
]

STRUCTURAL_AMYGDALA_ROIS = [
    "Left-Accessory_Basal_nucleus",
    "Left-Anterior_amygdaloid_area_AAA",
    "Left-Medial_nucleus",
    "Left-Cortical_nucleus",
    "Left-Corticoamygdaloid_transitio",
    "Right-Basal_nucleus",
    "Right-Accessory_Basal_nucleus",
    "Right-Anterior_amygdaloid_area_AAA",
    "Right-Cortical_nucleus",
    "Right-Corticoamygdaloid_transitio",
]

FALFF_NAME_MAP = {
    "Right-presubiculum": "Right_presubiculum",
    "Right-parasubiculum": "Right_parasubiculum",
    "Right-HATA": "Right_HATA",
    **{name: name.replace("-", "_") for name in STRUCTURAL_AMYGDALA_ROIS},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--hippocampal-input", type=Path, required=True)
    parser.add_argument("--amygdala-input", type=Path, required=True)
    parser.add_argument("--clinical-input", type=Path, required=True)
    parser.add_argument("--clinical-sheet", default=None)
    parser.add_argument("--hippocampal-column-map", type=Path, default=None)
    parser.add_argument("--amygdala-column-map", type=Path, default=None)
    parser.add_argument("--clinical-column-map", type=Path, default=None)
    parser.add_argument(
        "--site-regex",
        default=None,
        help="Optional regular expression with one capture group used to derive Site from SubID.",
    )
    parser.add_argument(
        "--duration-derived-onset-site",
        action="append",
        default=[],
        help="Site code for which onset is derived as age minus duration; repeat as needed.",
    )
    parser.add_argument(
        "--falff-input",
        type=Path,
        default=None,
        help="Final post-ComBat ROI CSV. Omit to run volume analyses only.",
    )
    parser.add_argument(
        "--falff-covariates",
        type=Path,
        default=None,
        help="fALFF covariate CSV or Excel table.",
    )
    parser.add_argument("--falff-column-map", type=Path, default=None)
    parser.add_argument("--falff-covariate-column-map", type=Path, default=None)
    return parser.parse_args()


def require_columns(frame: pd.DataFrame, columns: Iterable[str], source: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def apply_column_map(frame: pd.DataFrame, map_path: Path | None) -> pd.DataFrame:
    if map_path is None:
        return frame
    mapping = pd.read_csv(map_path)
    require_columns(mapping, ["source_name", "canonical_name"], str(map_path))
    if mapping["source_name"].duplicated().any():
        raise ValueError(f"{map_path} contains duplicate source_name values")
    if mapping["canonical_name"].duplicated().any():
        raise ValueError(f"{map_path} contains duplicate canonical_name values")
    rename = dict(zip(mapping["source_name"], mapping["canonical_name"], strict=True))
    return frame.rename(columns=rename)


def read_derived_table(path: Path, sheet_name: str | int | None = None) -> pd.DataFrame:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".xlsx":
        return pd.read_excel(path, sheet_name=0 if sheet_name is None else sheet_name)
    raise ValueError(f"Unsupported derived-table format: {path.suffix}")


def fdr_by_family(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["q_fdr_bh"] = np.nan
    for _, index in out.groupby("family", sort=False).groups.items():
        out.loc[index, "q_fdr_bh"] = multipletests(
            out.loc[index, "p_value"].to_numpy(), method="fdr_bh"
        )[1]
    out["significant_fdr_0_05"] = out["q_fdr_bh"] < 0.05
    return out


def fit_ols(
    frame: pd.DataFrame,
    outcome: str,
    predictor: str,
    covariates: list[str],
) -> dict:
    columns = [outcome, predictor, *covariates]
    data = frame.loc[:, columns].replace([np.inf, -np.inf], np.nan).dropna()
    if data[predictor].nunique() != 2:
        raise ValueError(f"{outcome}: predictor {predictor} does not have two levels")

    design = sm.add_constant(data[[predictor, *covariates]].astype(float), has_constant="add")
    if len(data) <= design.shape[1]:
        raise ValueError(f"{outcome}: insufficient observations for the requested model")
    if np.linalg.matrix_rank(design.to_numpy()) != design.shape[1]:
        raise ValueError(f"{outcome}: design matrix is rank deficient")
    model = sm.OLS(data[outcome].astype(float), design).fit()
    residual_sd = float(np.sqrt(model.mse_resid))
    outcome_scale = max(1.0, float(np.max(np.abs(data[outcome].astype(float)))))
    numerical_zero = np.finfo(float).eps * outcome_scale * len(data)
    if not np.isfinite(residual_sd) or residual_sd <= numerical_zero:
        raise ValueError(f"{outcome}: residual standard deviation is not positive")
    beta = float(model.params[predictor])
    ci_low, ci_high = model.conf_int(alpha=0.05).loc[predictor].astype(float)

    return {
        "n": int(model.nobs),
        "n_group_1": int((data[predictor] == 1).sum()),
        "n_group_0": int((data[predictor] == 0).sum()),
        "beta": beta,
        "residual_sd": residual_sd,
        "adjusted_cohens_d": beta / residual_sd,
        "ci95_low": float(ci_low / residual_sd),
        "ci95_high": float(ci_high / residual_sd),
        "t_value": float(model.tvalues[predictor]),
        "df_residual": float(model.df_resid),
        "p_value": float(model.pvalues[predictor]),
    }


def load_volumes(
    root: Path,
    hippocampal_input: Path,
    amygdala_input: Path,
    hippocampal_column_map: Path | None = None,
    amygdala_column_map: Path | None = None,
    site_regex: str | None = None,
) -> pd.DataFrame:
    del root
    hippocampal_path = hippocampal_input.resolve()
    amygdala_path = amygdala_input.resolve()
    hippocampus = apply_column_map(read_derived_table(hippocampal_path), hippocampal_column_map)
    amygdala = apply_column_map(read_derived_table(amygdala_path), amygdala_column_map)
    shared = ["SubID", "Sex", "Age", "subtype", "EstimatedTotalIntraCranialVol"]
    require_columns(hippocampus, [*shared, *HIPPOCAMPAL_ROIS], str(hippocampal_path))
    require_columns(amygdala, [*shared, *AMYGDALA_ROIS], str(amygdala_path))
    if hippocampus["SubID"].duplicated().any():
        raise ValueError(f"{hippocampal_path} contains duplicate subject identifiers")
    if amygdala["SubID"].duplicated().any():
        raise ValueError(f"{amygdala_path} contains duplicate subject identifiers")

    check = hippocampus[shared].merge(amygdala[shared], on="SubID", suffixes=("_h", "_a"))
    for column in shared[1:]:
        if not check[f"{column}_h"].equals(check[f"{column}_a"]):
            raise ValueError(f"Volume workbooks disagree in {column}")

    if "Site" in hippocampus.columns and "Site" in amygdala.columns:
        site_check = hippocampus[["SubID", "Site"]].merge(
            amygdala[["SubID", "Site"]], on="SubID", suffixes=("_h", "_a")
        )
        if not site_check["Site_h"].equals(site_check["Site_a"]):
            raise ValueError("Volume workbooks disagree in Site")
    merged = hippocampus.merge(
        amygdala.drop(columns=[*shared[1:], "Site"], errors="ignore"),
        on="SubID",
        validate="one_to_one",
    )
    if "Site" not in merged.columns:
        if site_regex is None:
            raise ValueError("Provide a canonical Site column or --site-regex")
        merged["Site"] = merged["SubID"].str.extract(site_regex, expand=False)
    if merged["Site"].isna().any():
        raise ValueError("Site is missing for at least one participant")
    return merged


def primary_volume_glm(volumes: pd.DataFrame) -> pd.DataFrame:
    data = volumes.loc[volumes["subtype"].isin([1, 2])].copy()
    data["mdd"] = (data["subtype"] == 2).astype(int)
    rows = []
    for family, rois in (("hippocampus", HIPPOCAMPAL_ROIS), ("amygdala", AMYGDALA_ROIS)):
        for roi in rois:
            rows.append(
                {
                    "analysis": "MDD_minus_HC",
                    "family": family,
                    "roi": roi,
                    **fit_ols(data, roi, "mdd", GROUP_COLS),
                }
            )
    return fdr_by_family(pd.DataFrame(rows))


def prepare_clinical(
    volumes: pd.DataFrame,
    root: Path,
    clinical_input: Path,
    clinical_sheet: str | None = None,
    clinical_column_map: Path | None = None,
    duration_derived_onset_sites: Iterable[str] = (),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    del root
    clinical_path = clinical_input.resolve()
    clinical = read_derived_table(clinical_path, clinical_sheet)
    clinical = apply_column_map(clinical, clinical_column_map).rename(
        columns=LEGACY_CLINICAL_COLUMNS
    )
    require_columns(
        clinical,
        ["SubID", CLINICAL_ONSET, CLINICAL_DURATION, CLINICAL_EPISODE, CLINICAL_MEDICATION],
        str(clinical_path),
    )
    mdd = volumes.loc[volumes["subtype"] == 2].merge(
        clinical, on="SubID", how="left", validate="one_to_one"
    )

    onset_reported = pd.to_numeric(mdd[CLINICAL_ONSET], errors="coerce")
    duration_months = pd.to_numeric(mdd[CLINICAL_DURATION], errors="coerce")
    onset_derived = mdd["Age"] - duration_months / 12.0

    derivation_sites = set(duration_derived_onset_sites)
    derive_from_duration = (
        onset_reported.isna()
        & mdd["Site"].isin(derivation_sites)
        & duration_months.notna()
    )
    mdd["age_at_onset"] = onset_reported.mask(derive_from_duration, onset_derived)
    mdd["onset_source"] = np.select(
        [onset_reported.notna(), derive_from_duration],
        ["reported", "derived_age_minus_duration"],
        default="missing",
    )
    mdd["early_onset"] = np.where(
        mdd["age_at_onset"].notna(), (mdd["age_at_onset"] <= 21).astype(float), np.nan
    )

    episode = pd.to_numeric(mdd[CLINICAL_EPISODE], errors="coerce")
    mdd["recurrent"] = np.where(episode.isin([0, 1]), (episode == 0).astype(float), np.nan)
    medication = pd.to_numeric(mdd[CLINICAL_MEDICATION], errors="coerce")
    mdd["medicated"] = np.where(
        medication.isin([0, 1]), (medication == 1).astype(float), np.nan
    )

    audit = pd.DataFrame(
        [
            {
                "analysis": "age_at_onset",
                "included": int(mdd["early_onset"].notna().sum()),
                "group_1": int((mdd["early_onset"] == 1).sum()),
                "group_0": int((mdd["early_onset"] == 0).sum()),
                "derived_from_duration": int(derive_from_duration.sum()),
            },
            {
                "analysis": "episode_status",
                "included": int(mdd["recurrent"].notna().sum()),
                "group_1": int((mdd["recurrent"] == 1).sum()),
                "group_0": int((mdd["recurrent"] == 0).sum()),
                "derived_from_duration": 0,
            },
            {
                "analysis": "current_medication",
                "included": int(mdd["medicated"].notna().sum()),
                "group_1": int((mdd["medicated"] == 1).sum()),
                "group_0": int((mdd["medicated"] == 0).sum()),
                "derived_from_duration": 0,
            },
        ]
    )
    return mdd, audit


def make_sample_reports(
    volumes: pd.DataFrame,
    clinical_audit: pd.DataFrame,
    motion_covariate_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create de-identified analytic-sample summaries for Table 1 and flow QC."""
    primary = volumes.loc[volumes["subtype"].isin([1, 2])].copy()
    primary["Motion"] = np.nan
    if motion_covariate_path is not None and motion_covariate_path.is_file():
        motion = read_derived_table(motion_covariate_path)
        require_columns(motion, ["SubID", "Motion"], str(motion_covariate_path))
        if motion["SubID"].duplicated().any():
            raise ValueError(f"{motion_covariate_path} contains duplicate subject identifiers")
        primary = primary.drop(columns="Motion").merge(
            motion[["SubID", "Motion"]], on="SubID", how="left", validate="one_to_one"
        )
    summary_rows = []
    for code, label in ((1, "HC"), (2, "MDD")):
        group = primary.loc[primary["subtype"].eq(code)]
        sex = pd.to_numeric(group["Sex"], errors="coerce")
        age = pd.to_numeric(group["Age"], errors="coerce")
        etiv = pd.to_numeric(group["EstimatedTotalIntraCranialVol"], errors="coerce")
        motion = pd.to_numeric(group["Motion"], errors="coerce")
        summary_rows.append(
            {
                "group": label,
                "n": int(len(group)),
                "age_n": int(age.notna().sum()),
                "age_mean": float(age.mean()),
                "age_sd": float(age.std(ddof=1)),
                "sex_1_n": int((sex == 1).sum()),
                "sex_2_n": int((sex == 2).sum()),
                "sex_other_or_missing_n": int((~sex.isin([1, 2])).sum()),
                "etiv_n": int(etiv.notna().sum()),
                "etiv_mean": float(etiv.mean()),
                "etiv_sd": float(etiv.std(ddof=1)),
                "mean_fd_n": int(motion.notna().sum()),
                "mean_fd_mean": float(motion.mean()),
                "mean_fd_sd": float(motion.std(ddof=1)),
            }
        )
    characteristics = pd.DataFrame(summary_rows)

    site_counts = (
        primary.assign(group=np.where(primary["subtype"].eq(2), "MDD", "HC"))
        .groupby(["Site", "group"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    for group in ("MDD", "HC"):
        if group not in site_counts:
            site_counts[group] = 0
    site_counts = site_counts[["Site", "MDD", "HC"]]
    site_counts["total"] = site_counts["MDD"] + site_counts["HC"]

    flow_rows = [
        {
            "stage": "matched_harmonized_volume_workbooks",
            "n": int(len(volumes)),
            "note": "Rows present after merging hippocampal and amygdala workbooks",
        },
        {
            "stage": "primary_MDD_HC_volume_analysis",
            "n": int(len(primary)),
            "note": "subtype codes 1 (HC) and 2 (MDD)",
        },
        {
            "stage": "excluded_other_subtype_codes",
            "n": int((~volumes["subtype"].isin([1, 2])).sum()),
            "note": "subtype codes outside the manuscript MDD-HC comparison",
        },
    ]
    for row in clinical_audit.to_dict(orient="records"):
        flow_rows.append(
            {
                "stage": f"within_MDD_{row['analysis']}",
                "n": int(row["included"]),
                "note": "Complete cases for the corresponding clinical variable",
            }
        )
    return characteristics, site_counts, pd.DataFrame(flow_rows)


def clinical_volume_glms(
    volumes: pd.DataFrame, mdd: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    roi_families = (
        ("hippocampus", STRUCTURAL_HIPPOCAMPAL_ROIS),
        ("amygdala", STRUCTURAL_AMYGDALA_ROIS),
    )
    comparisons = [
        ("early_minus_late", "early_onset"),
        ("recurrent_minus_first", "recurrent"),
        ("medicated_minus_unmedicated", "medicated"),
    ]

    direct_rows = []
    for analysis, predictor in comparisons:
        for family, rois in roi_families:
            for roi in rois:
                direct_rows.append(
                    {
                        "analysis": analysis,
                        "family": family,
                        "roi": roi,
                        **fit_ols(mdd, roi, predictor, GROUP_COLS),
                    }
                )
    direct = pd.concat(
        [fdr_by_family(group) for _, group in pd.DataFrame(direct_rows).groupby("analysis", sort=False)],
        ignore_index=True,
    )

    hc = volumes.loc[volumes["subtype"] == 1].copy()
    subgroup_rows = []
    subgroup_specs = [
        ("early_MDD_minus_HC", "early_onset", 1),
        ("late_MDD_minus_HC", "early_onset", 0),
        ("recurrent_MDD_minus_HC", "recurrent", 1),
        ("first_episode_MDD_minus_HC", "recurrent", 0),
        ("medicated_MDD_minus_HC", "medicated", 1),
        ("unmedicated_MDD_minus_HC", "medicated", 0),
    ]
    for analysis, predictor, level in subgroup_specs:
        patient_group = mdd.loc[mdd[predictor] == level, volumes.columns].copy()
        patient_group["patient_subgroup"] = 1
        controls = hc.copy()
        controls["patient_subgroup"] = 0
        comparison = pd.concat([controls, patient_group], ignore_index=True)
        for family, rois in roi_families:
            for roi in rois:
                subgroup_rows.append(
                    {
                        "analysis": analysis,
                        "family": family,
                        "roi": roi,
                        **fit_ols(comparison, roi, "patient_subgroup", GROUP_COLS),
                    }
                )
    subgroup = pd.concat(
        [fdr_by_family(group) for _, group in pd.DataFrame(subgroup_rows).groupby("analysis", sort=False)],
        ignore_index=True,
    )
    return subgroup, direct


def falff_glm(
    falff_path: Path,
    covariate_path: Path,
    volumes: pd.DataFrame,
    falff_column_map: Path | None = None,
    covariate_column_map: Path | None = None,
) -> pd.DataFrame:
    falff = apply_column_map(read_derived_table(falff_path), falff_column_map).rename(
        columns={"SubjectID": "SubID"}
    )
    covariates = apply_column_map(
        read_derived_table(covariate_path), covariate_column_map
    )
    require_columns(covariates, ["SubID", "Motion", "Age", "Sex"], str(covariate_path))
    require_columns(falff, ["SubID", *FALFF_NAME_MAP.values()], str(falff_path))
    if falff["SubID"].duplicated().any():
        raise ValueError(f"{falff_path} contains duplicate subject identifiers")
    if covariates["SubID"].duplicated().any():
        raise ValueError(f"{covariate_path} contains duplicate subject identifiers")

    diagnosis = volumes[["SubID", "subtype"]]
    data = falff.merge(covariates, on="SubID", validate="one_to_one").merge(
        diagnosis, on="SubID", validate="one_to_one"
    )
    data = data.loc[data["subtype"].isin([1, 2])].copy()
    data["mdd"] = (data["subtype"] == 2).astype(int)

    rows = []
    families = (
        ("hippocampus", STRUCTURAL_HIPPOCAMPAL_ROIS),
        ("amygdala", STRUCTURAL_AMYGDALA_ROIS),
    )
    for family, structural_names in families:
        for structural_name in structural_names:
            falff_name = FALFF_NAME_MAP[structural_name]
            rows.append(
                {
                    "analysis": "MDD_minus_HC",
                    "family": family,
                    "roi": structural_name,
                    **fit_ols(data, falff_name, "mdd", ["Age", "Sex", "Motion"]),
                }
            )
    return fdr_by_family(pd.DataFrame(rows))


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    output = (args.output or (root / "results" / "reproduced")).resolve()
    output.mkdir(parents=True, exist_ok=True)

    volumes = load_volumes(
        root,
        args.hippocampal_input,
        args.amygdala_input,
        args.hippocampal_column_map,
        args.amygdala_column_map,
        args.site_regex,
    )
    primary = primary_volume_glm(volumes)
    mdd, clinical_audit = prepare_clinical(
        volumes,
        root,
        args.clinical_input,
        args.clinical_sheet,
        args.clinical_column_map,
        args.duration_derived_onset_site,
    )
    motion_path = args.falff_covariates.resolve() if args.falff_covariates else None
    characteristics, site_counts, analytic_flow = make_sample_reports(
        volumes,
        clinical_audit,
        motion_path if motion_path is not None and motion_path.is_file() else None,
    )
    subgroup, direct = clinical_volume_glms(volumes, mdd)
    falff = None
    if args.falff_input is not None:
        falff_path = args.falff_input.resolve()
        if args.falff_covariates is None:
            raise ValueError("--falff-covariates is required with --falff-input")
        covariate_path = args.falff_covariates.resolve()
        if not falff_path.is_file():
            raise FileNotFoundError(f"fALFF input not found: {falff_path}")
        if not covariate_path.is_file():
            raise FileNotFoundError(f"fALFF covariates not found: {covariate_path}")
        falff = falff_glm(
            falff_path,
            covariate_path,
            volumes,
            args.falff_column_map,
            args.falff_covariate_column_map,
        )
    primary.to_csv(output / "volume_primary.csv", index=False)
    subgroup.to_csv(output / "volume_subgroup_vs_hc.csv", index=False)
    direct.to_csv(output / "volume_within_mdd.csv", index=False)
    if falff is not None:
        falff.to_csv(output / "falff_mdd_vs_hc.csv", index=False)
    clinical_audit.to_csv(output / "clinical_sample_audit.csv", index=False)
    characteristics.to_csv(output / "sample_characteristics_primary.csv", index=False)
    site_counts.to_csv(output / "site_counts_primary.csv", index=False)
    analytic_flow.to_csv(output / "analytic_sample_flow.csv", index=False)
    print(f"Analysis outputs written to: {output}")


if __name__ == "__main__":
    main()
