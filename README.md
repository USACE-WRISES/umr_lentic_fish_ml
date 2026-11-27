# Fish Occurrence Models within Backwater Habitats in the Upper Mississippi River

This repository contains the R code used to model the occurrence of 33 fish species 
in backwater (contiguous floodplain lake) habitats of Pools 4, 8, and 13 of the 
Upper Mississippi River, to compare spatial and non-spatial cross-validation 
strategies, and to identify species groups with similar habitat affinities.

The code represents the analysis conducted for a manuscript entitled "Predicting 
Upper Mississippi River backwater fish presence from habitat conditions managers 
can forecast for restoration" a manuscript by Edward Stowe, Kristen Bouska, 
Nicole Ward, Sebastian Teran-Hidalgo, and Kyle McKay currently in peer review.

---

## Overview

**Goal.** Predict species occurrence in backwater habitats from environmental variables 
that can reasonably be estimated for a site *before* a restoration project is built, 
and quantify how well those predictions transfer to new locations.

**Approach.**

- Random forest classification models (`ranger` via `tidymodels`) fit separately for each of 33 species.
- Two nested predictor sets: the **Bluegill** variables (depth, dissolved oxygen, 
temperature, current velocity, plus sampling period) and the **restoration** variables 
(Bluegill variables plus vegetation density, substrate type, and snag presence).
- Three cross-validation schemes of increasing spatial independence: random, 
lake-cluster, and pool-level spatial. Comparing them quantifies how much spatial 
autocorrelation inflates apparent model performance.
- Class-probability thresholds chosen per species by minimizing the Euclidean 
distance to perfect classification on the ROC curve.
- Variable importance and partial dependence for the primary model (pool-level 
spatial CV, restoration variables), followed by clustering of species by habitat-affinity profile.

**Six model configurations** are fit per species (3 CV schemes x 2 predictor sets).

---

## Code structure

```
.
├── config.R                       # Shared settings used by all scripts
├── 01_data_prep.R                 # Build analysis data and the three CV fold sets
├── 02_run_models.R                # Tune random forests for all six model configurations
├── 03_thresholds_diagnostics.R    # Thresholds, performance metrics, calibration diagnostics
├── 04_importance_pdp.R            # Final model fits, permutation importance, partial dependence
├── 05_synthesis.R                 # CV comparisons, figures
├── 06_save_shiny_objects.R        # Assemble objects needed by the Shiny app
├── 07_site_map.R                  # Create map of project area
├── 08_clustering.R                # Clustering species based on model results
└── output/
    ├── models/                    # Fold objects and tuning results (one file per model config)
    ├── results/                   # Metrics, diagnostics, importance, PDPs, synthesis objects
    ├── figures/                   # Manuscript and supplemental figures
    └── shiny/                     # Final objects loaded by a separate Shiny app
```

## Data

### Source

Fish data are from the **Long Term Resource Monitoring (LTRM) element of the Upper Mississippi River 
Restoration Program**, collected by the U.S. Geological Survey Upper Midwest Environmental Sciences Center 
and partner state agencies. Backwater habitat was identified using the LTRM **Aquatic Areas** spatial database.

- LTRM data portal: https://umesc.usgs.gov/data_library/fisheries/fish_page.html
- Aquatic Areas database: https://doi.org/10.5066/F7Z89BMP

### Archive

The specific fish data used in this analysis is archived on Zenodo (https://doi.org/10.5281/zenodo.24604233) 
