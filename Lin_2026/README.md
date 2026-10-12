# MDD hippocampal/amygdala multimodal analysis code

This repository contains the custom ROI-level and FC cluster-reporting code
described in the manuscript. It is designed to run without embedded personal
paths or participant-level data.

## Inputs

DIRECT provides the underlying neuroimaging data, not the participant-level
tables used by these GLMs. The GLM inputs are study-derived artifacts created
after preprocessing, ROI extraction, quality control, and ComBat harmonization;
they are not distributed in this repository.

DIRECT Phase II access and its Data Use Agreement are described on the
[official data-sharing page](https://rfmri.org/DIRECTPhase2). The dataset record
is available at <https://doi.org/10.57760/sciencedb.psych.00689>.

File names and source column names are configurable. See `DATA_INTERFACE.md`
and `config/column_map.example.csv` for the canonical schema and mapping
interface. The diagnosis code is `1 = HC`, `2 = MDD`; codes 3 and 4 are excluded
from the manuscript MDD-HC analyses.

## ROI statistics

Use Python 3.11-3.13. Create an isolated environment, install
`requirements.txt`, and run:

```text
python code/roi_analysis.py \
  --hippocampal-input derived/hippocampal_volumes.csv \
  --amygdala-input derived/amygdala_volumes.csv \
  --clinical-input derived/clinical_variables.csv
python -m unittest discover -s tests -v
```

For differently named derived inputs, pass `--hippocampal-input`,
`--amygdala-input`, `--clinical-input`, and the corresponding optional
`--*-column-map` files. Run `python code/roi_analysis.py --help` for the complete
interface. Where specified by the study protocol, use
`--duration-derived-onset-site <site-code>` to derive onset as scan age minus
illness duration for the relevant site; repeat the option for multiple sites.

Outputs are written to `results/reproduced/`. The script covers:

- all 42 MDD-HC volume GLMs;
- early/late, first/recurrent, and medicated/unmedicated subgroup-vs-HC GLMs;
- the three direct within-MDD volume contrasts;
- adjusted effect sizes, 95% CIs, t statistics, raw P values, and family-wise
  Benjamini-Hochberg q values.

It also creates de-identified primary sample characteristics, site counts, and
an analytic-sample flow table to support Table 1 and the analysis flow.

The clinical analysis implements the manuscript-defined age-at-onset, episode,
and medication groupings and records aggregate analysis counts at run time.

## fALFF extraction and ComBat

`code/extract_falff_individual_rois.m` extracts mean values from individual ROI
masks after nearest-neighbour reslicing to each fALFF map. Start from
`config/falff_manifest.example.csv` and `config/roi_definitions.example.csv`.

From MATLAB, with `DPABI_ROOT` set locally, the entry points are:

```matlab
combat_falff_maps('config/falff_manifest.csv', ...
    fullfile(getenv('DPABI_ROOT'), 'Templates', 'BrainMask_05_97x115x97.nii'), ...
    'results/falff_combat', getenv('DPABI_ROOT'));

extract_falff_individual_rois('config/falff_postcombat_manifest.csv', ...
    'config/roi_definitions.csv', 'results/falff_roi_values.csv', ...
    getenv('DPABI_ROOT'));
```

The manifest paths are resolved relative to their CSV files; no absolute
workstation path is required in a public configuration. The extraction
manifest must reference the final post-ComBat maps produced by the reported
DPABI GUI workflow (or by the reference implementation above).

The reported analysis used DPABI Harmonization through the graphical interface.
`code/combat_falff_maps.m` is a portable reference implementation of the same
parametric model: site is the batch variable, and age, sex, mean FD, and
diagnosis are preserved. It is not a replacement for the filled GUI provenance
record. Both MATLAB scripts require DPABI on the MATLAB path.

Run the final fALFF GLMs explicitly after the final post-ComBat table is ready:

```text
python code/roi_analysis.py --falff-input derived/falff_roi_postcombat.csv \
  --falff-covariates derived/motion_covariates.csv \
  --hippocampal-input derived/hippocampal_volumes.csv \
  --amygdala-input derived/amygdala_volumes.csv \
  --clinical-input derived/clinical_variables.csv
```

The script applies the manuscript GLM and FDR definitions to the selected
hippocampal and amygdala ROI families.

## DPABI GUI operations

GUI-based DPABI stages remain valid parts of the workflow, but their settings
must be versioned alongside the code. `workflow/DPABI_GUI_OPERATIONS.md` records
the required provenance for structural/fALFF/FC ComBat, FC calculation, the FC
group GLM, and GRF thresholding. Copy
`config/dpabi_gui_parameters.yml` when adapting the workflow. The
published configuration uses logical artifact names and relative paths; it does
not require an execution date, analyst identity, or personal workstation path.

Use `code/hash_manifest.py` to create SHA-256 manifests for the final analysis
inputs and outputs so a GUI rerun can be tied to exact files. Detailed manifests
whose paths contain participant identifiers remain private; only publish a
sanitized artifact-level manifest when release is authorized.

## FC cluster and peak reporting

`ana_DPABIViewerFixed.m` reads the four statistic maps listed in a seed
manifest, reproduces the two-tailed GRF cluster threshold, reports peak
coordinates, and annotates clusters with DPABI atlases. Configure it with
`config/fc_config_example.m` or these environment variables:

- `MDD_FC_RESULTS_ROOT`
- `DPABI_TEMPLATE_DIR`
- `MDD_FC_SEED_MANIFEST`

With the first two variables set, run the supplied configuration from the
repository root:

```matlab
run('config/fc_config_example.m')
```

The four manuscript-defined seeds are the right basal nucleus, bilateral
accessory basal nuclei, and right cortical nucleus. They use DIRECT Phase II
individualized MNI-space labels selected by DPABI wildcard matching with label
value 7008. The analysis mask is DPABI Templates
`BrainMask_05_97x115x97.nii`, and the group contrast is MDD minus HC.

Positive and negative tails are labelled independently to prevent adjacent
opposite-sign voxels from being merged. The script requires the statistic maps,
their DF/smoothness header metadata, analysis mask, and atlas resources.

## Tests

The public test suite contains no study data. It verifies complete-case model
fitting, output structure, rank checks, and family-wise Benjamini-Hochberg
correction, and it runs the complete Python CLI on synthetic volume, clinical,
motion, and fALFF inputs. A release guard also rejects participant-level file
types, generated result directories, oversized files, and personal absolute
paths.

## Software requirements

- Python 3.11-3.13 and the pinned packages in `requirements.txt`.
- MATLAB with Image Processing Toolbox and Statistics and Machine Learning
  Toolbox.
- DPABI/DPABISurf, including the atlas resources and helper functions listed in
  `workflow/DPABI_GUI_OPERATIONS.md`.

The reported analysis environment used MATLAB R2020a on Ubuntu 18.04. The
postprocessing configuration records DPABI V9.0_250415 and DPABISurf 3.3, with
DPABI releases checked through April 2026. Versions used to create the shared
DIRECT preprocessing derivatives are recorded separately in the manuscript
Methods and GUI parameter record.

The public configuration also records the mask, atlas, and template resources
used by the workflow.
Third-party dependencies and their licenses are summarized in
`THIRD_PARTY_NOTICES.md`; none of those packages or resources is redistributed
here.
