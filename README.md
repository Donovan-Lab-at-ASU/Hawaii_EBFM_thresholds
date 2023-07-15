[![DOI](https://zenodo.org/badge/666594763.svg)](https://zenodo.org/badge/latestdoi/666594763)

This repository includes data and analysis scripts to accompany:

# Utility of indicator thresholds across spatial gradients for applications to ecosystem-based fisheries management
*Under consideration for publication*

### Author of analysis and code: Shannon M Hennessey
-----

### Description:
This work analyzes data from the Hawaii Monitoring and Reporting Collaborative (HIMARC) to assess the utility of a multi-model approach to identifying thresholds between ecosystem indicators and environmental drivers.

### Contents:
#### Scripts:
* **1_gradient_forests.R:** R script for fitting a gradient forest model to the indicators and drivers.
* **2_gamm_fitting.R:** R script that fits Generalized Additive Mixed Models (GAMMs) for indicator-driver pairs.
* **3_gamm_prediction.R:** R script that reads in output of 2_gamm_fitting.R and generates predicted GAMM relationships for indicator-driver pairs.
* **4_gamm_bootstrapping.R:** R script that reads in output of 3_gamm_prediction.R and generates bootstrapped GAMM predicted fits for indicator-driver pairs.
* **5_gamm_CI_derivatives.R:** R script that reads in outputs of 3_gamm_prediction.R and 4_gamm_bootstrapping.R and calculates confidence intervals and 1st and 2nd derivatives of fitted GAMM relationships.
* **6_gamm_threshold_locations.R:** R script that reads in outputs of 5_gamm_CI_derivatives.R and identifies threshold ranges and point estimates.
* **7_gamm_threshold_plots.R:** R script that reads in outputs of 6_gamm_threshold_locations.R and plots fitted GAMMs with bootstrapped confidence intervals, and threshold range(s) and point estimate(s).
* **8_spp_correlations.R:** R script with analyses of indicator correlations.
* **functions.R:** R script with functions to extract GAMM diagnostics and plot fitted GAMMs with threshold locations.

#### Data:
A data file with fish indicator biomass and environmental driver values can be found [here](https://github.com/Donovan-Lab-at-ASU/Hawaii_herbivore_management/blob/main/data/donovan_etal_HIMARCdata_for_herbivore_publication.csv).
