# Yan_2025

This repository contains scripts for statistical analysis, meta-analysis, and visualization of cortical measures in Major Depressive Disorder (MDD).

---

## 1. System Requirements

- **Operating System**
  - Tested on: **Ubuntu 18.04 LTS**

- **Software Dependencies**
  - **MATLAB R2022a** (required)
    - SPM12
    - DPABISurf
  - **Python 3.6.13**
    - Standard scientific libraries: `numpy`, `scipy`, `pandas`, `matplotlib`
  - **R 3.6.1**
    - Standard statistical libraries (e.g., `meta`, `metafor`)

- **Hardware Requirements**
  - Actual hardware requirements depend on the size of the dataset
  - Demo and typical analyses run in a short time (< several minutes)

---

## 2. Installation Guide

1. **MATLAB**
   - Install MATLAB R2022a or later.
   - Add SPM12 and DPABISurf to MATLAB path.

2. **Python**
   - Install Python 3.6.13.
   - Create a virtual environment and install dependencies:
     ```bash
     conda create -n yan2025 python=3.6.13
     conda activate yan2025
     pip install numpy scipy pandas matplotlib
     ```

3. **R**
   - Install R 3.6.1.
   - Install required packages from CRAN:
     ```R
     install.packages(c("meta", "metafor"))
     ```

**Typical installation time**: a few minutes on a standard desktop computer.

---

## 3. Demo

A small demo dataset is provided in `4_DemoData/`, including two sites of MDD vs. healthy controls, and a `ModelSummaryInfo.txt`.

### Steps

1. **Site-level statistics**
   - Run scripts in `1_StatsAtSites/` on demo data.
   - Output: statistical results for each site.

2. **Meta-analysis**
   - Use site-level results as input for scripts in `2_MetaAnalysis/`.
   - Output: aggregated meta-analysis results.

3. **Visualization**
   - Run scripts in `3_plot/` to generate final figures.
   - Output: brain maps, forest plots, and summary statistics.

**Expected runtime**: short (a few minutes) on a standard desktop.

---

## 4. Instructions for Site-level Statistics

### Running the Site-level Statistics on Your Data

1. **Prepare the environment**
   - Install SPM12 and DPABISurf.
   - Alternatively, install FreeSurfer and GNU Parallel.

2. **Set up data paths in MATLAB**
   ```matlab
   DataDir='/data/SurfData/Results'; %Data Dir for stats. Should have AnatSurfLH AnatSurfRH under this dir.
   CovariatesFile='/1_StatsAtSites/Covariates_Data.xlsx'; %Covariates File (.xlsx).  SHOULD have a column of eTIV!!!
   FreesurferSubjectsDir='/opt/freesurfer/subjects'; %Freesurfer dir for smoothing. Should have 'fsaverage' subject under this dir.
   ```

3. **Optional parameters**
   ```matlab
   MinimumSubjectNumber = 10; % Default: 10
   ParallelWorkersNumber = 1; % Default: 1
   IsDPABISurfDataStyle = 0;  % Default: 0 (set 1 if data is from DPABISurf)
   ```

4. **Run the main script**
   ```matlab
   ENIGMA_RMD_Vertex_Stats(DataDir, CovariatesFile, FreesurferSubjectsDir, ...
                           MinimumSubjectNumber, ParallelWorkersNumber, IsDPABISurfDataStyle);
   ```

---

## 5. Folder Contents

### `1_StatsAtSites/`
| File | Description |
|------|-------------|
| `Dockerfile` | Defines a containerized environment. |
| `Commands.txt` | Example shell/Matlab commands for running site-level statistics. Serves as quick-start instruction. |
| `Covariates_Data.xlsx` | Example covariates file. **Must include an "eTIV" column.** |
| `*.m` (statistical code) | MATLAB scripts for computing site-specific statistics (patient vs. control comparisons). Output: site-level statistical summaries. |

### `2_MetaAnalysis/`
| File | Description |
|------|-------------|
| `DIRECT_Meta_Analysis_onD_run.m` | Performs the main meta-analysis across sites. |
| `DIRECT_Meta_Regression_Model1.m` | Conducts regression analysis (e.g., effect of covariates across sites). |

### `3_plot/`
| File | Description |
|------|-------------|
| `cluster_report_deal.py` | Extracts peak coordinates from cluster reports. |
| `deal_sort.py` | Filters clusters by vertex size threshold and sorts them. |
| `peak_extract.m` | MATLAB script: uses peak indices to extract site-specific values. |
| `R_Cal_Meta_peak.R` | R script: computes meta-analysis forest plots. |
| `mergeenigmastatcsv.py` | Merges cluster report contents with statistical result tables. |

---
