###################################################################################################
# Default parameters to load into regression simulations
#
# Author(s): Boyan Duan, James Leiner
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
R = 100