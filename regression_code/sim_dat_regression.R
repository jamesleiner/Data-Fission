###########################################################################
# Simulate synthetic datasets for regression experiments

# Author(s): Boyan Duan, James Leiner
#############################################################################

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




generate_linear = function(n, p, beta, type, rho, add_influential = c()){
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
  Y = rnorm(n+length(add_influential)) + X%*%beta
  cluster = factor(1:(n+length(add_influential)))
  return(list(X = X, Y = Y, cluster = cluster, Sigma = Sigma))
}


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
