###################################################################################################
# Default parameters to load into regression simulations
#
# Parameters:
#   - n (integer): Number of observations.
#   - prob (numeric): Parameter tau
#   - beta (numeric vector): Coefficients for the covariates.
#   - type (character): Type of data generation for the design matrix ("independent" or "correlated").
#   - slope (numeric): Slope parameter for time series trend construction.
#   - p (numeric): Probability of changing the trend for time series trend construction.
#   - sigma (numeric): Standard deviation of the noise for time series trend construction.
#   - degree (numeric): Degree parameter to use for Toeplitz matrix for covariance matrix construction.
#   - scale (numeric): Scale parameter for regression simulations.
#   - alpha (numeric): Significance level for constructing confidence intervals.
#   - methods (character vector): Methods to use in simulations ("masking" - data fission, "full" - full dataset reuse, "split" - data splitting).
#   - CI_type (character): Type of confidence intervals for trend filtering uncertainty quantification ("uniform" or "pointwise").
#   - error_type (character): Type of error distribution ("gaussian").
#   - R (integer): Number of iterations in simulations.
#
###################################################################################################

n = 1000
prob = 0.5
beta = c(1, 0, rep(1,5), rep(0,10), rep(2,3))
type = "independent"

slope = 0.5
p = 0.99
sigma = 1
degree = 0


scale = 1
alpha = 0.2
methods = c("masking", "full", "split")
CI_type = "uniform"
error_type = "gaussian"
R = 100