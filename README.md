# FuSD
Functional Signed Directionality (FuSD) is a distribution-free, nonparametric framework for quantile-based inference in function spaces. This repository provides code to compute FuSD, set-valued functional quantiles, quantile-based spread and skewness summaries, and tight probability bands, with reproducible environmental case studies. For reproducibility, we include processed functional datasets for the Monterrey metropolitan area (Nuevo León, Mexico), including Landsat 8 Land Surface Temperature (LST) and Normalized Difference Vegetation Index (NDVI) curves from 2018 to 2024, and scripts to reproduce the empirical illustrations reported in the paper and Supplementary Material.


# FuSD — LST/NDVI analysis (2018–2024)

This repository provides the computational support for the paper:

**Nonparametric functional quantiles, skewness, and probability bands via Functional Signed Directionality**

It contains a minimal, script-based implementation to compute FuSD quantities and generate the main analysis outputs for monthly **LST** and **NDVI** functional data.

## Contents

- **`LSTNDVI18to24.RData`**  
  Preprocessed functional datasets:
  - **LST**: monthly *maxima* (2018–2024)
  - **NDVI**: monthly *means* (2018–2024)

- **`FuSDPackage060126.R`**  
  Core FuSD functions (“the engine”): FuSD computation (with variants depending on the chosen order function), FuSD quantiles, tight probability bands, and related utilities.

- **`FuSDAnalysis.R`**  
  Main analysis script: loads data, runs the FuSD pipeline, computes evaluation summaries (e.g., \(F_{\beta}^{\gamma}\) and its components), FuSD coverage, spread/skewness, and generates figures.

- **`GraphsAux.R`**  
  Plot helpers (themes, legend presets, small utilities used by the analysis figures).

## Requirements

- R (recommended: RStudio)
- Packages used heavily are in the tidyverse ecosystem:
  - `dplyr`, `tidyr`, `ggplot2` (or simply `tidyverse`)

> Note: The scripts may use additional packages depending on your setup. If something is missing, R will error with the package name—install it and re-run.

## Quick start

1. Put the four files in the same folder and set that folder as your working directory.
2. Run the analysis script.

This workflow reproduces the empirical results and figures exactly as reported in the reference paper and its Supplementary Material (for the included preprocessed datasets and the default settings in `FuSDAnalysis.R

Example:

```r
# From the repo root (folder containing the 4 files)
# setwd("path/to/repo")  # optional if you open the folder directly

source("FuSDPackage060126.R")
source("GraphsAux.R")

load("LSTNDVI18to24.RData")

# Run the full pipeline (produces outputs + figures)
source("FuSDAnalysis.R")

````


## Switching data (LST vs NDVI) and study period

`FuSDAnalysis.R` is designed to be the single “control panel”. The dataset and time window are selected by subsetting the monthly index.

In the current version, the study period is controlled by a block like:

```r
J = 6  # 2024
per = (J*12 + 1):((J + 1)*12)
lst = Dlst[per, ]
Fdata = t(lst)
````

Minimal guidance:

- `J` selects the year **relative to the first year in the dataset** (e.g., if the data start in 2018, then `J = 0` corresponds to 2018, `J = 1` to 2019, …, `J = 6` to 2024).
- `per` is the vector of **monthly indices** extracted from the full series.
- `Dlst` is used for **LST**. To run **NDVI**, replace it with the NDVI object (e.g., `Dndvi`) in the same block and keep the remaining lines unchanged.
- `Fdata = t(lst)` converts the selected year block into the expected orientation for downstream functions.
- To analyze a different period, adjust `J` (or directly edit `per` to select any custom set of months).

## Outputs

Running `FuSDAnalysis.R` typically produces:

- Functional FuSD quantiles plot
- Probability bands
- Functional skewness plot
- Spread comparison
- \(F_{\beta}^{\gamma}\) comparison and components
- FuSD coverage comparison

> If you want outputs to go into a dedicated folder (e.g., `fig/`), create it and update the `ggsave()` paths inside `FuSDAnalysis.R`.

## Notes

- This is research code intended to reproduce FuSD analyses and figures.
- The preprocessed `.RData` object names are loaded into the workspace by `load(...)`. If you are unsure what was loaded, run:
  ```r
  load("LSTNDVI18to24.RData")
  ls()

## Authors

- Emmanuel Ambriz — Statistics Program, King Abdullah University of Science and Technology (KAUST), Thuwal 23955-6900, Saudi Arabia — emmanuel.ambrizlobato@kaust.edu.sa  
- Graciela González Farías — Centro de Investigación en Matemáticas (CIMAT), Guanajuato 36023, Mexico — farias@cimat.mx  
- Emilien Joly — Centro de Investigación en Matemáticas (CIMAT), Guanajuato 36023, Mexico — emilien.joly@cimat.mx  
- Maurizio Filippone — Statistics Program, King Abdullah University of Science and Technology (KAUST), Thuwal 23955-6900, Saudi Arabia — maurizio.filippone@kaust.edu.sa  
- Hernando Ombao — Statistics Program, King Abdullah University of Science and Technology (KAUST), Thuwal 23955-6900, Saudi Arabia — hernando.ombao@kaust.edu.sa  




