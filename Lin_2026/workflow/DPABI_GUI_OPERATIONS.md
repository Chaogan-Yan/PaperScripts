# DPABI GUI operations and provenance

GUI operations are part of the reproducible workflow. They are not replaced by
pseudocode; instead, the repository records the analysis-relevant settings,
software versions, and operation order. The public parameter record is
`config/dpabi_gui_parameters.yml`.

For GUI navigation, the free DPABI/DPABISurf/DPARSF course is available at
<https://edu.deepbrain.com/RFMRICourse>. Its modules cover DPABI statistical
analysis, DPABI Viewer, individualized seed-based connectivity, and DPABI
Harmonization. The course is a supplementary operating tutorial; the filled
configuration in this repository remains the record of the study-specific
settings.

DIRECT source imaging files and locally derived analysis tables are distinct
provenance levels. Each GUI record must identify its DIRECT-derived imaging
inputs and the derived artifact it creates; local participant-level tables are
not assumed to be DIRECT-distributed files.

## Public reproducibility record

Record all of the following:

1. DPABI/DPABISurf version or release date, MATLAB release, and operating system.
2. DPABI module and method selected in the GUI.
3. Logical input artifact, mask, and subject-order rule.
4. Batch variable, preserved covariates, categorical coding, and missing-data rule.
5. Every numeric threshold and checkbox state that changes the analysis.
6. Logical output artifact and output naming rule.

Execution dates, analyst identities, personal workstation paths, completion
logs, and detailed participant-level manifests are not part of the public Code
Availability repository. Runtime roots are supplied locally; the FC manifest
uses paths relative to the locally supplied results root.

## Operation order

### 1. Structural volume ComBat

- Module: DPABI Harmonization.
- Method: ComBat, parametric empirical Bayes.
- Batch: scanning site.
- Preserved variables for the primary MDD-HC analysis: age, sex, eTIV, diagnosis.
- Input: 24 hippocampal and 18 amygdala volumes plus identifiers/covariates.
- Output: harmonized volume matrix used by `code/roi_analysis.py`.
- QC: unique subject identifiers, identical subject order across matrices, no
  missing volume/covariate values in the primary comparison, and no constant ROI.

### 2. fALFF ComBat

- Module: DPABI Harmonization.
- Method: ComBat, parametric empirical Bayes.
- Batch: scanning site.
- Preserved variables: age, sex, mean FD, diagnosis.
- State whether ComBat was applied to voxelwise fALFF maps or to an organized
  ROI matrix. Deposit the exact post-ComBat artifact used downstream.
- Record the analysis mask, interpolation, subject order, and output suffix.

### 3. Individual-specific fALFF extraction

- Code: `code/extract_falff_individual_rois.m`.
- Inputs: post-ComBat fALFF map and participant-specific ROI mask.
- Reslicing: nearest-neighbour for label masks.
- Summary: mean over finite in-mask voxels.
- Retain the companion voxel-count QC table.

### 4. Seed-based FC calculation

- Module: DPABI seed-based functional connectivity.
- Seeds: the four structurally altered ROIs with FDR-significant fALFF results.
- Seed identities: right basal nucleus, bilateral accessory basal nuclei, and
  right cortical nucleus.
- Seed labels: DIRECT Phase II individualized labels transformed to MNI space,
  located with DPABI FC wildcard matching; label value 7008.
- GUI input branch: `FunVoluGlobal`, based on `FunSurf`; the `global` suffix is
  DPABI's global-signal-regression branch.
- One functional session; zero GUI parallel workers; 2 mm isotropic volume
  normalization.
- Functional Connectivity, Extract ROI time courses, Normalize Derivatives,
  and Smooth Derivatives selected. Direct input smoothing was not selected;
  derivative volume smoothing used 6 x 6 x 6 mm FWHM.
- Nuisance Covariates Regression, ICA-AROMA Noise Remove, Filter, Scrubbing,
  and head-motion scrubbing regressors were not selected in this FC GUI run.
  The input was already preprocessed as described in the manuscript Methods.
- No additional initial time points were removed in this downstream GUI run.
- Extract the mean seed time series, calculate voxelwise Pearson correlations,
  and apply Fisher r-to-z transformation.

### 5. FC ComBat

- Module: DPABI Harmonization.
- Method: ComBat, parametric empirical Bayes.
- Batch: scanning site.
- Preserved variables: age, sex, mean FD, diagnosis.
- Apply the same subject ordering and analysis mask used for group statistics.
- The manuscript Methods must explicitly state this step if it was used.

### 6. FC MDD-HC group GLM

- Module: DPABI Statistical Analysis.
- Outcome: Fisher-z FC map.
- Variable of interest: diagnosis, with the sign of the contrast recorded.
- Contrast direction: MDD minus HC.
- Covariates: age, sex, mean FD.
- Record group coding, design-matrix column order, contrast vector, analysis
  mask, and contrast direction.
- Save the unthresholded T map with degrees of freedom and smoothness metadata.

### 7. GRF correction and peak reporting

- Correction: voxel P < 0.001, cluster P < 0.05, two-tailed, 26-connectivity.
- Correction is performed separately for each seed.
- Code: `ana_DPABIViewerFixed.m`.
- Seed mapping: `config/fc_seed_manifest.csv`, containing exactly four rows and
  statistic-map paths relative to the locally supplied results root.
- Analysis mask: `BrainMask_05_97x115x97.nii` from the DPABI Templates folder.
- Inputs: unthresholded `T2.nii`, brain mask, DPABI atlas files, and NIfTI header
  containing DF and dLh/FWHM, or explicitly recorded manual values.
- The script labels positive and negative connected components independently.

## Recommended release artifacts

- Filled GUI parameter YAML.
- Software-version text file.
- Optional sanitized artifact-level SHA-256 manifests, if their release is
  separately authorized; participant-level manifests stay private.
- De-identified GUI settings screenshots or sanitized human-readable parameter
  exports. Do not publish binary settings files that have not been audited.
- A sanitized DPABI console transcript saved as text, with participant labels
  and workstation paths removed, if one is retained.
- Full ROI result tables, not only significant rows.
