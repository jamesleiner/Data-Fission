###################################################################################################
# Load in required packages and files 

###################################################################################################


source("regression_code/regression.R")
source("regression_code/sim_dat_regression.R")
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
  library(matrixStats)
  library(latex2exp)
  library(ggforce)
  library(VGAM)
  library(sn)
  
  #library(glmgen)
  #library(trendfiltering)
})