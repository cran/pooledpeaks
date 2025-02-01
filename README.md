<!-- badges: start -->
  [![R-CMD-check](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kmkuesters/pooledpeaks/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Package Overview

pooledpeaks is designed for analyzing genetic data obtained from Fragment 
Analysis output files (.fsa) of pooled biological samples. It provides functions
for a comprehensive analysis pipeline from processing .fsa files, to cleaning 
the peak data, and conducting population genetic analyses. Some features are 
listed below and a usage example of the entire pipeline is included as a vignette.

### Features
  * **Peak Scoring:** Process .fsa files and score peaks contained therein.
  * **Data Manipulation:** Clean and prepare peak data for downstream analyses.
  * **Population Genetics Analysis:** 
    * Calculate Gene Identity Matrix and Genetic Distance Matrix
    * Calculate diversity indices
    * Calculate differentiation indices
    * Perform cluster analysis
  * **Visualization:** Visualize the peak scoring and genetic analysis results.

The pooledpeaks package was developed by the Blanton Lab as part of Kathleen 
Kuesters' dissertation.

# Installation Instructions

You can install the package directly from GitHub using the following instructions:

Open R and copy the following code into your console

**Install devtools and pooledpeaks from GitHub**

  * install.packages("devtools")

  * devtools::install_github("kmkuesters/pooledpeaks")

**Install pooledpeaks directly from CRAN**

  * install.packages("pooledpeaks")

# References:

* Covarrubias-Pazaran et al. (2016) <doi:10.1186/s12863-016-0365-6> 
* Long et al. (2022) <doi:10.1038/s41598-022-04776-0> 
* Jost (2008) <doi:10.1111/j.1365-294x.2008.03887.x> 
* Nei (1973) <doi:10.1073/pnas.70.12.3321>
* Foulley et al. (2006) <doi:10.1016/j.livprodsci.2005.10.021>
* Chao et al. (2008) <doi:10.1111/j.1541-0420.2008.01010.x>

