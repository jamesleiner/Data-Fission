# Data-Fission

This repository contains code to replicate all figures and experiments found in our paper: [Data Fission: splitting a single data point](https://arxiv.org/abs/2112.11079).

## Pre-requisites
Prior to running the code, please install and download the [trendfiltering](https://capolitsch.github.io/trendfiltering/) package by  [Collin Politsch](https://collinpolitsch.com/) as follows.

```
install.packages("remotes")
remotes::install_github("capolitsch/trendfiltering")
```

The remainder of the pre-requistie packages can be download through the CRAN repository. 

## Instructions
* The experiments for the interactive multiple testing application (section 3) need to be run on a cluster. To do this, please
  * Create a results folder in your working directory to store output
  * Run the script `interactive_experiments_batch.sh` on a cluster (currently written assuming a slurm scheduler)
  * Run the script `Interactive Testing -- Figures.R` to reproduce the plots in the document. 
* The experiments for sections 4,5, and 6 may be run on a personal computer but can be sped up considerably if using a a machine with multiple cores. To reproduce the results, please:
  * Create a results folder in your working directory to store output
  * Run the script `Regression Experiments.R`
  * Run the script `Regression Experiments -- Figures.R` to reproduce the plots in the document (excluding section 6.3). 
  * Run the script `Spectroscopy Examples.R` to reproduce the plots for section 6.3 in the document (the spectroscopy live data example). 
* An additional file, `Reproduce Miscellaneous Figures.R` includes additional figures that are not core results but are used for illustrative purposes in the paper. 



## Acknolwedgements
* The dataset for the spectroscopy example (section 6.3) was graciously provided to us by Collin Politsch, compiled from the Baryon Oscillation Spectroscopic Survey. 
* All of the base underlying functions used in the interactive testing experiments are taken from the [public repository](https://github.com/lihualei71/STAR) for the paper [STAR: A general interactive framework for FDR control under structural constraints](https://arxiv.org/pdf/1710.02776.pdf). Thank you to Lihua Lei and coauthors.
