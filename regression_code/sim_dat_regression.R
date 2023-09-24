###################################################################################################
# Code for generating synthetic datasets to be used in large-scale testing of experiments 
# Includes code for Gaussian, Poisson, Trend Filtering, and and Binomial Regression synthetic datasets. 
#
###################################################################################################



###################################################################################################
# Generate Poisson-distributed data
# Generates synthetic data with Poisson-distributed response variable.
# Inputs:
#   - n (integer): Number of observations.
#   - p (integer): Number of covariates.
#   - beta (numeric vector): Coefficients for the covariates.
#   - type: Type of data generation for the design matrix ("independent" or "correlated"). If the data
#           is correlated, then the covariates are drawn from a multivariate normal with covariance matrix
#           designed as Toeplitz matrix with correlation parameter rho.
#   - rho (numeric): Correlation parameter for correlated data.
#   - add_influential (numeric vector): Additional influential covariate to add into the design matrix. By default
#           the user's vector will be interpreted as a multiple to apply to the largest existing entry in the design matrix.
# Outputs:
#   - X (matrix): Design matrix.
#   - Y (numeric vector): Poisson-distributed response.
#   - cluster (factor): Cluster information for use in some methods that assumed heteroscedastic responses.
#     This functionality is not required for these experiments, so only a vector of 1s is returned.
###################################################################################################

generate_poisson = function(n, p, beta, type, rho, add_influential = c()){
  if (type == "independent") {
    X = cbind(
      # c(rep(0, n/2), rep(1, n/2)), #a1
      # c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
      matrix(rbinom(n = n*2, size = 1, prob = 0.5), ncol = 2),
      matrix(rnorm(n = n*(p-2), sd = 1), ncol = p-2)) #a3 
  } else {
    p_size = p/5
    #Sigma = toeplitz((p_size:1)^degree/p_size^degree)#could be other decreasing rate
    Sigma = toeplitz(rho^(0:(p_size-1)))
    Sigma = bdiag(lapply(1:5, function(x) Sigma))
    X = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
  }
  if(length(add_influential)>0) {
    baseline = apply(X,2,max)
    for(i in 1:length(add_influential)){
      mult = add_influential[i]
      X=rbind(X,mult*baseline)
    }
  }
  cluster = factor(1:nrow(X))
  Y = sapply(exp(X%*%beta), function(lambda){rpois(1, lambda)})
  return(list(X = X, Y = Y, cluster = cluster))
}




###################################################################################################
# generate_linear: Generates synthetic data with Gaussian-distributed response variable.
# Inputs:
#   - n (integer): Number of observations.
#   - p (integer): Number of covariates.
#   - beta (numeric vector): Coefficients for the covariates.
#   - type: Type of data generation for the design matrix ("independent" or "correlated"). If the data
#           is correlated, then the covariates are drawn from a multivariate normal with covariance matrix
#           designed as Toeplitz matrix with correlation parameter rho.
#   - rho (numeric): Correlation parameter for correlated data.
#   - add_influential (numeric vector): Additional influential covariates (optional).
#   - error_type (character): Type of error distribution to use ("gaussian" - Gaussian distrubted errors, "t"
#      - t-distributed errors, "laplace" -- Laplace distributed errors, "sn" - skewed normal distributed errors).
# Outputs:
#   - X (matrix): Design matrix.
#   - Y (numeric vector): Gaussian-distributed response variable.
#   - cluster (factor): Cluster information for use in some methods that assumed heteroscedastic responses.
#     This functionality is not required for these experiments, so only a vector of 1s is returned.
#   - Sigma (matrix): Covariance matrix for covariance (if correlated data).
#   - sd (numeric): Standard deviation of the response variable.
###################################################################################################

generate_linear = function(n, p, beta, type, rho, add_influential = c(),error_type="gaussian"){
  if (type == "independent") {
    # X = cbind(
    #   # c(rep(0, n/2), rep(1, n/2)), #a1
    #   # c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
    #   matrix(rbinom(n = n*2, size = 1, prob = 0.5), ncol = 2),
    #   matrix(rnorm(n = n*(p-2), sd = 1), ncol = p-2)) #a3 
    X = matrix(rnorm(n = n*p, sd = 1), ncol = p)
    Sigma = diag(rep(1, p))
  } else {
    p_size = p/5
    #Sigma = toeplitz((p_size:1)^degree/p_size^degree)#could be other decreasing rate
    Sigma = toeplitz(rho^(0:(p_size-1)))
    Sigma = bdiag(lapply(1:5, function(x) Sigma))
    X = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
  }
  
  if(length(add_influential)>0) {
    baseline = apply(X,2,max)
    for(i in 1:length(add_influential)){
      mult = add_influential[i]
      X=rbind(X,mult*baseline)
    }
  }
  if(error_type == "gaussian"){
    Y = rnorm(n+length(add_influential)) + X%*%beta
    sd = 1 
  }
  if(error_type == "t") {
    print("T")
    df = 5
    error = rt(n+length(add_influential),df) 
    sd = sqrt(df/(df-2))
    Y = error/sd + X%*%beta
  }
  if(error_type == "laplace") {
    print("LAPLACE")
    Y = rlaplace(n+length(add_influential),1/sqrt(2)) + X%*%beta
    sd=1
  }
  if(error_type == "sn") {
    print("Skewed Normal")
    omega = 1
    alpha = 5
    error = rsn(n+length(add_influential),0,omega,alpha)
    mu = omega*(alpha/(sqrt(1+alpha**2)))*sqrt(2/pi)
    sd = sqrt(omega**2 *(1-2*(alpha**2)/(1+alpha**2)/pi))
    Y = (error-mu)/sd + X%*%beta
  }
  cluster = factor(1:(n+length(add_influential)))
  
  return(list(X = X, Y = Y, cluster = cluster, Sigma = Sigma,sd=sd))
}


###################################################################################################
# generate_logistic: Generates synthetic data with Binomial-distributed response variable.
# Inputs:
#   - n (integer): Number of observations.
#   - p (integer): Number of covariates.
#   - beta (numeric vector): Coefficients for the covariates.
#   - type (character): Type of data generation for the design matrix ("independent" or "correlated").
#                       If the data is correlated, then the covariates are drawn from a multivariate normal
#                       with a covariance matrix designed as a Toeplitz matrix with correlation parameter rho.
#   - add_influential (numeric vector): Additional influential covariate to add into the design matrix.
#                                      By default, the user's vector will be interpreted as a multiple to apply
#                                      to the largest existing entry in the design matrix.
# Outputs:
#   - X (matrix): Design matrix.
#   - Y (numeric vector): Binomial-distributed response.
#   - exp_y (numeric vector): Expected probabilities.
#   - cluster (factor): Cluster information for use in some methods that assume heteroscedastic responses.
#     This functionality is not required for these experiments, so only a vector of 1s is returned.
###################################################################################################


generate_logistic = function(n, p, beta, type = "independent",add_influential = c()){
  if (type == "independent") {
    if (p <= 2) {
      X = matrix(rnorm(n = n*p, sd = 1), ncol = p)
      cluster = factor(1:n)
    } else {
      X = cbind(
        # c(rep(0, n/2), rep(1, n/2)), #a1
        # c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
        matrix(rbinom(n = n*2, size = 1, prob = 0.5), ncol = 2),
        matrix(rnorm(n = n*(p-2), sd = 1), ncol = p-2)) #a3 
      cluster = factor(1:n)
    }
  } else {
    p_size = p/5
    cluster = factor(1:n)
    #Sigma = toeplitz((p_size:1)^degree/p_size^degree)#could be other decreasing rate
    Sigma = toeplitz(rho^(0:(p_size-1)))
    Sigma = bdiag(lapply(1:5, function(x) Sigma))
    X = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
  }
  if(length(add_influential)>0) {
    baseline = apply(X,2,max)
    for(i in 1:length(add_influential)){
      mult = add_influential[i]
      X=rbind(X,mult*baseline)
    }
  }
  exp_y = 1/(1 + exp(-X%*%beta))
  Y = sapply(exp_y, function(prob){rbinom(1, 1, prob)})
  cluster = factor(1:nrow(X))
  return(list(X = X, Y = Y, exp_y = exp_y, cluster = cluster))
}


###################################################################################################
# time_seq: Construct time series trend to evaluate trend filtering on, as a function of number of knots, slope, and noise. 
# Inputs:
#   - n (integer): Number of observations.
#   - slope (numeric): Slope of the trend at change points.
#   - p (numeric): Probability of changing the trend.
#   - sigma (numeric): Standard deviation of the noise.
# Outputs:
#   - X (numeric vector): True trend.
#   - Y (numeric vector): Observed time series with noise.
#   - knots (integer vector): Indices of trend change points.
###################################################################################################
time_seq = function(n, slope, p, sigma){
  x = rep(0, n); v = rep(0, n)
  x[1] = 0; v[1] = 0; knots = c()
    #runif(1, min = -0.5, max = 0.5)
  for(i in 2:n){
    ind = rbinom(1, 1, 1 - p)
    v[i] = v[i-1] + ind*(runif(1, min = -slope, max = slope) - v[i-1])
    x[i] = x[i-1] + v[i-1]
    if (ind) {knots = c(knots, i)}
  }
  y = x + rnorm(n, sd = sigma)
  return(list(X = x, Y = y, knots = knots))
}

###################################################################################################
# true_para: Helper function to identify where the KL minimizer is by finding the root of a given score equation. 
# Inputs:
#   - X (matrix): Design matrix.
#   - beta (numeric vector): True coefficients used for data generation.
#   - selected (integer vector): Indices of selected covariates.
#   - type (character): Type of data generation ("poisson", "binomial_mask" - binomial data when fission is used, 
#       "binomial_split" - binomial data when data spliitting is used, "binomial").
#   - g_Y (numeric): For binomial data with fission, the fissioned data generated.
#   - prob (numeric): For binomial data with fission, the parameter tau used with data fissioned.
#   - split_ind (numeric vector): For binomial data with split, the indicator indicating which points are used to generate f(Y).
# Outputs:
#   - ss (numeric vector): KL minimized parameters for the selected model
###################################################################################################

true_para = function(X, beta, selected, type, g_Y = NA, prob = NA, split_ind = NA){
  # exp_omit = rowSums(X[,selected] %*% beta[selected]) +
  #   colSums(beta[-selected]*(mu[-selected] +
  #                          Sigma[-selected,selected]%*%solve(Sigma[selected,selected])%*%(t(X[,selected]) - mu[selected])))
  # sigma_omit = 0.5*t(beta[-selected]) %*%
  #   (Sigma[-selected,-selected] - Sigma[-selected,selected]%*%solve(Sigma[selected,selected])%*%Sigma[selected,-selected])%*%
  #   beta[-selected]
  # exp_y = exp(exp_omit + as.numeric(sigma_omit))
  if (type == "poisson") {
    exp_y = exp(rowSums(X %*% beta))
    model = function(para){
      return(t(X[,selected])%*%exp_y - t(X[,selected])%*%exp(X[,selected]%*%para))
    }
    ss <- multiroot(f = model, start = beta[selected])$root
  } else if (type == "binomial_mask") {
    exp_y = 1/(1 + exp(-X[,-1]%*%beta))
    cond_exp_y = exp_y/(exp_y + (1 - exp_y)*(prob/(1 - prob))^(2*g_Y - 1))
    sl = c(1, selected + 1)
    model = function(para){
      return(t(X[,sl])%*%cond_exp_y - t(X[,sl])%*%(1/(1 + exp(-X[,sl]%*%para))))
    }
    # if (length(selected) > 1) {
    # } else {
    #   model = function(para){
    #     return(t(X[,selected])%*%cond_exp_y - t(X[,selected])%*%(1/(1 + exp(-X[,selected]*para))))
    #   }
    # }
    ss <- multiroot(f = model, start = c(0, beta[selected]))$root
} else if (type == "binomial_split") {
  exp_y = 1/(1 + exp(-X[split_ind == 0, -1]%*%beta))
  sl = c(1, selected + 1)
  model = function(para){
    return(t(X[split_ind == 0, sl])%*%exp_y -
             t(X[split_ind == 0,sl])%*%(1/(1 + exp(-X[split_ind == 0,sl]%*%para))))
  }
  ss <- multiroot(f = model, start = c(0, beta[selected]))$root
} else if (type == "binomial") {
  exp_y = 1/(1 + exp(-X[,-1]%*%beta))
  sl = c(1, selected + 1)
  model = function(para){
    return(t(X[,sl])%*%exp_y - t(X[,sl])%*%(1/(1 + exp(-X[,sl]%*%para))))
  }
  # if (length(selected) > 1) {
  #   
  # } else {
  #   model = function(para){
  #     return(t(X[,selected])%*%exp_y - t(X[,selected])%*%(1/(1 + exp(-X[,selected]*para))))
  #   }
  # }
  ss <- multiroot(f = model, start = c(0, beta[selected]))$root
}
  return(ss)
}
