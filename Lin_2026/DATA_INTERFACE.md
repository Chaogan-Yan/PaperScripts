# Data interface

DIRECT distributes the underlying neuroimaging data. The participant-level
spreadsheets and CSV files used during local analysis are study-derived
artifacts produced after preprocessing, ROI extraction, quality control, and
ComBat harmonization; they are not DIRECT-distributed files and are not included
in this repository.

Source imaging access: <https://rfmri.org/DIRECTPhase2> and
<https://doi.org/10.57760/sciencedb.psych.00689>. Access is governed by the
DIRECT Phase II Data Use Agreement described on the official sharing page.

The analysis code consumes canonical derived tables, independent of local file
names. Supply file paths with the command-line options in `code/roi_analysis.py`.
If an upstream tool uses different headers, provide a two-column mapping CSV
with `source_name` and `canonical_name`; see
`config/column_map.example.csv`. Separate maps may be supplied for hippocampal,
amygdala, clinical, fALFF, and fALFF-covariate tables.

## Canonical identifiers and covariates

- `SubID`: unique participant identifier shared across derived tables
- `Site`: scanning-site identifier; alternatively provide `--site-regex`
- `subtype`: diagnosis code (`1 = HC`, `2 = MDD` for the manuscript contrasts)
- `Age`, `Sex`, `EstimatedTotalIntraCranialVol`; `Sex` uses the study's binary
  codes `1` and `2`
- fALFF covariates: `SubID`, `Age`, `Sex`, `Motion`
- clinical fields: `SubID`, `age_at_onset_years`,
  `illness_duration_months`, `first_episode`, `currently_medicated`

Clinical binary coding is `first_episode: 1 = yes, 0 = no` and
`currently_medicated: 1 = yes, 0 = no`. Other values are treated as missing.

ROI columns use the canonical names defined in `code/roi_analysis.py`. Column
maps allow differently named upstream outputs to be used without editing the
analysis code.

## Output tables

All ROI GLM tables share the columns `analysis`, `family`, `roi`, `n`,
`n_group_1`, `n_group_0`, `beta`, `residual_sd`, `adjusted_cohens_d`,
`ci95_low`, `ci95_high`, `t_value`, `df_residual`, `p_value`, `q_fdr_bh`, and
`significant_fdr_0_05`. Generated outputs are written under `results/`, which is
excluded from Git.

## Provenance boundary

The public repository should contain code, canonical schemas, field-mapping
templates, software versions, and DPABI settings. A hash manifest may be public
only when its displayed paths and artifact labels contain no participant-level
identifiers and release is authorized; detailed participant-level manifests stay
private. The repository should not contain local participant-level derived
tables unless their release is separately authorized. The Methods and Data
Availability sections should distinguish DIRECT source imaging data from these
derived analysis inputs and describe how each derived artifact was produced.
Personal absolute paths, execution dates, and analyst identities are not needed
by the public code. Runtime roots are supplied locally, while public manifests
use paths relative to those roots.
