source("regression.R")
source("sim_dat_regression.R")
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
})