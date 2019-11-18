# Purpose

This repository contains all the data, functions, scripts to run simulations and analysis, and scripts to generate plots for the paper
"Covariance-based sample selection for heterogenous data: Applications to gene expression and autism risk gene detection".

# Installation

This package can be installed through `devtools` in R.

```{r}
library("devtools")
devtools::install_github("linnylin92/covarianceSelection", subdir = "covarianceSelection")
```
The package and downstream analysis and simulations depends on several packages, some of which originate
from BioConductor. 
These include `DBI`, `dequer`, `devtools`, `doMC`, `foreach`, `glmnet`, `hash`, `huge`, `igraph`, `MASS`, `Matrix`, and `org.Hs.eg.db`. 

Warning: The `doMC` package does not really work on Windows, as it does not seem to actually parallelize the code.
See: http://stackoverflow.com/questions/16453625/package-domc-not-available-for-r-version-3-0-0-warning-in-install-packages

# Data 

The two major datasets used in this article are also included in the repository.
The first dataset is the BrainSpan microarray measurements collected by Kang et al. (2011). While the original dataset 
is publicly available on GEO (\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25219}),
we provide a locally preprocessed dataset, which was created to be amendable for our analysis in R.
This dataset is a 105 MB `.RData` file, and is synced onto GitHub using the Git Large File Storage system (https://git-lfs.github.com/). Please
install this system before proceeding.


The second dataset is the p-value risk scores (called TADA scores in the paper) for the genes obtained applying
 the TADA framework (He et al., 2013) to the data available in De Rubeis et al. (2014). The full citations are given in the paper.
 
Supplementary datasets were also used to manage the data and assess the results. These are all documented appropriately under the file `covarianceSelection/R/data.R`.
All data used in this entire project are either publicly available or processed by our lab on publicly available data.

# Reproducing the results

## Note

All the code below were run on a server with 10 cores. If you do not have 10 cores, be sure to change the value of `ncores` appropriately in the following files:
`main/main.R` (for the data analysis) or `simulation/gaussian.R`, 
`simulation/gaussian_beta.R`, `simulation/nonparanormal.R`, `simulation/nonparanormal_beta.R`, 
`simulation/nonparanormal_goodness.R`, `simulation/nonparanormal_accelerated.R` and
`simulation/nonparanormal_accelerated_beta.R` (for the simulations).

All results produced are automatically placed in the `results/` folder as `.RData` files.

The simulations rely on a custom-made simulation engine developed in https://github.com/linnylin92/simulation. 
To download this (before running any simulation), use the code chunk below in R.

```{r}
library("devtools")
devtools::install_github("linnylin92/simulation", subdir = "simulation")
```

## Running the simulations

To reproduce the simulations (Section 5 of our paper), navigate to the `simulation` folder. All
simulations in this folder (documented below) take about 8 hours to finish.
From this location, run the following line in the command window to reproduce the 
results for Figures 6, 7, and 8A.

```
R CMD BATCH nonparanormal.R # this runs the simulation
R CMD BATCH nonparanormal_postprocess.R # this generates Figure 6
R CMD BATCH nonparanormal_postprocess_2.R # this generates Figure 7
R CMD BATCH nonparanormal_postprocess_3.R # this generates Figure 8A
```

Run the following line in the command window to reproduce the results for Figure 8B. 

```
R CMD BATCH nonparanormal_beta.R # this runs the simulation
R CMD BATCH nonparanormal_beta_postprocess.R # this generates Figure 8B
```

Run the following line in the command window to reproduce the results for Figures S.2 and S.3.

```
R CMD BATCH nonparanormal_goodness.R # this runs the simulation
R CMD BATCH nonparanormal_goodness_postprocess.R # this generates Figure 8B
```

Run the following line in the command window to reproduce the results for Figures S.4 and S.5A.

```
R CMD BATCH gaussian.R # this runs the simulation
R CMD BATCH gaussian_postprocess.R # this generates Figure S.4
R CMD BATCH gaussian_postprocess_2.R # this generates Figure S.5A
```

Run the following line in the command window to reproduce the results for Figure S.5B.

```
R CMD BATCH gaussian_beta.R # this runs the simulation
R CMD BATCH gaussian_beta_postprocess.R # this generates Figure S.5B
```

Run the following line in the command window to reproduce the results for Figures S.6 and S.7A.

```
R CMD BATCH nonparanormal_accelerated.R # this runs the simulation
R CMD BATCH nonparanormal_accelerated_postprocess.R # this generates Figure S.6
R CMD BATCH nonparanormal_accelerated_postprocess_2.R # this generates Figure S.7A
```

Run the following line in the command window to reproduce the results for Figure S.7B.

```
R CMD BATCH nonparanormal_accelerated_beta.R # this runs the simulation
R CMD BATCH nonparanormal_accelerated_beta_postprocess.R # this generates Figure S.7B
```


## Running the analysis

To reproduce the analysis (Section 6 of our paper), navigate to the `main` folder. From this location, run the following lines in the command window.

```
R CMD BATCH main.R
```

These took 15 hours on our machine respecitvely. 
This runs a sequence of steps in the analysis pipeline, which we briefly describe here.
* `step0_loading.R` loads the BrainSpan dataset and TADA dataset (i.e., matching the genes in both
datasets, resolving gene synonyms, removing genes not expressed in the brain).
* `step1_screening.R` screens the genes according to Liu et al. (2015). This is reported in Section 6.1.
* `step2_nodawn_analysis.R` detects the risk genes only based on the TADA dataset.
* `step3_pfc35_analysis.R` implicates risk genes in the DAWN framework using the Window 1B partitions. This is reported in Section 6.4.
* `step4_subjectselection.R` uses COBS and is the most computational-intensive part of our procedure. 
This selects the 24 partitions in BrainSpan that we report in our paper. This is reported in Section 6.2.
* `step5_our_analysis.R` implicates risk genes in the DAWN framework using our 24 partitions selected by COBS.
This is reported in Section 6.4.
* `step6_our_analysis_robustness.R` performs the robustness analysis. This is reported in Section 6.4.
* `step7_goodness.R` performs the goodness of fit diagnostic. This is reported in Section 3.2 and 6.2.
* `step8_results.R` collects all the key results from all the above sections. This is reported in Section 6.4.
* `step9_figures.R` generates all the figures related to the analysis sections in our paper.

All the resulting figures will be placed in the `figures` folder.
