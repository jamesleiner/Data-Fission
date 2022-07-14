###########################################################################
# Load required packages and scripts for regression simulations

# Author(s): Boyan Duan, James Leiner
#############################################################################
suppressPackageStartupMessages({
  library(glmnet)
  library(clubSandwich)
  library(sandwich)
  library(MASS)
  library(rootSolve)
  library(genlasso)
  library(splines)
  library(cplm)
  library(quantreg)
  library(parallel)
  library(ggplot2)
  library(dplyr)
  library(magrittr)
  library(tidyr)
  library(rlang)
  library(purrr)
  library(parallel)
  library(zoo)
})

source("regression_code/regression.R")
source("regression_code/sim_dat_regression.R")
