###################################################################################################
# Helper functions to execute interactive multiple testing simulations
#
###################################################################################################


###################################################################################################
# Generate Gaussian noise with specific correlation structure (rho) using Toeplitz matrix.
# Inputs:
#   - n (numeric): Number of observations.
#   - mu (numeric): Mean vector.
#   - rho (numeric): Correlation parameter.
#   - type (numeric): Type of correlation structure (1, 2, or 3).
# Outputs:
#   - Sigma (matrix): Covariance matrix with the specified correlation structure.
###################################################################################################

sig.gen <- function(n, mu, rho, type){
  if (type == 1){
    new.rho <- rho / n
    Sigma <- diag(rep(1 - new.rho, n)) + new.rho
    z <- mvrnorm(1, mu, Sigma)
  } else if (type == 2){
    m <- ceiling(n / 2)
    ind1 <- sample(n, m)
    ind2 <- (1:n)[-ind1]
    Sigma <- diag(rep(0, n))
    Sigma[ind1, ind1] <- abs(rho)
    Sigma[ind2, ind2] <- abs(rho)
    Sigma[ind1, ind2] <- -abs(rho)
    Sigma[ind2, ind1] <- -abs(rho)
    diag(Sigma) <- 1            
    z <- mvrnorm(1, mu, Sigma)
  } else if (type == 3){
    Sigma <- diag(rep(0, n))
    for (i in 1:n){
      Sigma[i, ] <- rho * (1 - 2 * abs(rho))^abs(1:n - i)
    }
    diag(Sigma) <- 1
  }
  return(Sigma)
}


###################################################################################################
# Takes in p-values and returns rejection sequences for pre-determined set of alphas using BH.
# Inputs:
#   - pvals (numeric vector): Vector of p-values.
#   - H0 (logical vector): Logical vector indicating true null hypotheses.
#   - values (numeric vector): Vector of values associated with location at each point on the grid.
#   - mu (numeric vector): Mean vector associated with the entry at each point on the grid.
#   - variance (numeric): Variance of the error term.
#   - alpha.list (numeric vector): Vector of significance levels to evaluate BH on.
# Outputs:
#   - df (data frame): Data frame containing number of rejections, FDP, power, and mean of rejected region
###################################################################################################

summary.BH <- function(pvals, H0,values, mu,variance,
                       alpha.list = seq(0.01, 0.3, 0.01)){
  n <- length(pvals)
  nfrej <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    sum(pvals[!H0] < alpha, na.rm = TRUE)
  })
  ntrej <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    sum(pvals[H0] < alpha, na.rm = TRUE)
  })
  true_mean <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    mean(mu[pvals < alpha], na.rm = TRUE)
  })
  measured_mean <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    mean(values[pvals < alpha], na.rm = TRUE)
  })
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0), 1)
  sd <- sqrt(variance/nrej)
  df <- data.frame(nrej = nrej, FDP = FDP, power = power, mean=measured_mean, true_mean=true_mean,sd=sd)
  return(df)
}


###################################################################################################
# Takes in p-values and returns rejection sequences for pre-determined set of alphas using AdaPT.
# Inputs:
#   - adapt (AdaPT object): AdaPT object.
#   - pvals (numeric vector): Vector of p-values.
#   - H0 (logical vector): Logical vector indicating true null hypotheses.
#   - values (numeric vector): Vector of values associated with location at each point on the grid.
#   - mu (numeric vector): Mean vector associated with the entry at each point on the grid.
#   - variance (numeric): Variance of the error term.
# Outputs:
#   - df (data frame): Data frame containing number of rejections, FDP, power, and mean of rejected region
###################################################################################################

summary.AdaPT <- function(adapt, H0, pvals,vals,mu,variance){
  nfrej <- apply(adapt$s, 2, function(s){
    tmp <- (pvals <= s)
    sum(tmp[!H0], na.rm = TRUE)
  })
  ntrej <- apply(adapt$s, 2, function(s){
    tmp <- (pvals <= s)        
    sum(tmp[H0], na.rm = TRUE)
  })
  true_mean <- apply(adapt$s, 2, function(s){
    tmp <- (pvals <= s)        
    mean(mu[tmp], na.rm = TRUE)
  })
  measured_mean <- apply(adapt$s, 2, function(s){
    tmp <- (pvals <= s)        
    mean(vals[tmp], na.rm = TRUE)
  })
  
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0),1)
  sd <- sqrt(variance/nrej)
  df <- data.frame(nrej = nrej, FDP = FDP, power = power, mean=measured_mean, true_mean =true_mean,sd=sd)
  return(df)
}


###################################################################################################
# Takes in p-values and returns rejection sequences for pre-determined set of alphas using STAR.
# Inputs:
#   - STAR.obj (STAR object): STAR object.
#   - H0 (logical vector): Logical vector indicating true null hypotheses.
#   - values (numeric vector): Vector of values associated with location at each point on the grid.
#   - mu (numeric vector): Mean vector associated with the entry at each point on the grid.
#   - variance (numeric): Variance of the error term.
# Outputs:
#   - df (data frame): Data frame containing number of rejections, FDP, power, and mean of rejected region
###################################################################################################

summary.STAR <- function(STAR.obj, H0,vals,mu,variance){
  nfrej <- apply(STAR.obj$mask, 2, function(x){
    sum(x[!H0], na.rm = TRUE)
  })
  ntrej <- apply(STAR.obj$mask, 2, function(x){
    sum(x[H0], na.rm = TRUE)
  })
  true_mean <- apply(STAR.obj$mask, 2, function(x){
    mean(mu[x])
  })
  measured_mean <- apply(STAR.obj$mask, 2, function(x){
    mean(vals[x])
  })
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0),1)
  sd <- sqrt(variance/nrej)
  df <- data.frame(nrej = nrej, FDP = FDP, power = power, mean=measured_mean, true_mean =true_mean,sd=sd)
  return(df)
}

###################################################################################################
# Run interactive testing experiments for a single trial run. Generating data according to the specified distribution type 
# and running multiple statistical tests (BH, AdaPT, STAR) on the data.
# Inputs:
#   - x (numeric): Input parameter representing locations on grid.
#   - mu (numeric vector): Mean vector representing value at each parameter.
#   - tau (numeric): Tau chosen for splitting the data.
#   - alpha.list (numeric vector): Vector of significance levels to evaluate BH on.
#   - num.steps.update.score (numeric): Number of steps to update the score.
#   - scope.params (numeric): Scope parameters.
#   - alt (numeric): Value to assign to grid for locations following an alternative hypothesis.
#   - null (numeric): Value to assign to grid for locations following a null hypothesis (typically 0).
#   - type (character): Type of distribution to run simulation on ("normal" or "poisson").
#   - filename (character): Optional filename for saving results.
#   - seed (numeric): Random seed for reproducibility.
# Outputs:
#   - List of experimental results, one for each trial run.
###################################################################################################

experiment_masked <- function(x,mu,tau,alpha.list,num.steps.update.score,scope.params,alt,null,type="normal",filename=NULL,seed=1){
  set.seed(seed)
  if(type == "normal"){
    
    #fission the data
    var <- sig.gen(n,mu,0,1)
    Y <- mvrnorm(1,mu,var)
    Z <- mvrnorm(1,rep(0,n),var)
    f_Y <- Y+tau*Z
    g_Y <- Y-(1/tau)*Z
    
    #return p-values for both fissioned and non-fissioned data
    pvals <- 1 - pnorm(Y)
    pvals_mask <- 1- pnorm(f_Y,sd=sqrt((1+tau**2)))
  }
  if(type == "poisson"){
    Y <- rpois(n,mu)
    Z <- rbinom(n,Y,tau)
    g_Y = Y - Z
    f_Y = Z
    
    #masking scheme to make p-values continuous
    Y_minus = Y- 0.1
    f_Y_minus = f_Y - 0.1
    U1 <- runif(length(mu),0,1)
    U2 <- runif(length(mu),0,1)
    rand <- U1*ppois(Y, null) + (1-U1)* ppois(Y_minus, null)
    rand_F <- U2*ppois(f_Y, null*tau) + (1-U2)* ppois(f_Y_minus, null*tau)
    pvals <- 1-rand
    pvals_mask <- 1 - rand_F
  }
  
  
  #Reject with BH, AdaPT, STAR for full data
  STAR.obj1 <- STAR.convex(pvals, x, 
                           alpha.list = alpha.list,
                           type = "model-assist",
                           update.mask.params = list(dir = "min"),
                           num.steps.update.score = num.steps.update.score,
                           score.params = score.params)
  AdaPT.obj1 <- AdaPT(x, pvals, cov.formula = cov.formula,
                      q.list = alpha.list, plot.quiet = TRUE)
  
  
  BH.result_full <- summary.BH(pvals, H0,Y,mu,1, alpha.list = alpha.list)
  STAR.result_full <- summary.STAR(STAR.obj1,H0,Y,mu,1)
  adapt.result_full <- summary.AdaPT(AdaPT.obj1, H0, pvals,Y,mu,1)

  
  #Reject with BH, AdaPT, STAR for fissioned data
  STAR.obj2 <- STAR.convex(pvals_mask, x, 
                           alpha.list = alpha.list,
                           type = "model-assist",
                           update.mask.params = list(dir = "min"),
                           num.steps.update.score = num.steps.update.score,
                           score.params = score.params)
  AdaPT.obj2 <- AdaPT(x, pvals_mask, cov.formula = cov.formula,
                      q.list = alpha.list, plot.quiet = TRUE)
  
  var_mask <- 1+ (1/(tau**2))
  BH.result_mask <- summary.BH(pvals_mask, H0,g_Y,mu,var_mask,alpha.list = alpha.list)
  STAR.result_mask <- summary.STAR(STAR.obj2,H0,g_Y,mu,var_mask)
  adapt.result_mask <- summary.AdaPT(AdaPT.obj2, H0, pvals,g_Y,mu,var_mask)
  return(list(x=x,mu=mu,tau=tau,H0=H0,mu=mu,alt=alt,null=null,type=type,var=var,var_mask=var_mask,
              Y=Y, f_Y=f_Y,g_Y=g_Y,pvals=pvals,pvals_mask=pvals_mask,
              BH.result.mask=BH.result_mask, BH.result.full=BH.result_full,
              adapt.result.mask=adapt.result_mask,adapt.result.full=adapt.result_full,
              STAR.result.mask=STAR.result_mask, STAR.result.full = STAR.result_full,
              STAR.full.object=STAR.obj1,AdaPT.full.object=AdaPT.obj1,
              STAR.mask.object=STAR.obj2,AdaPT.mask.object=AdaPT.obj2))
  if(is.null(filename)){
    save(file = paste(filename,tau,type,format(Sys.time(), "%Y-%m-%d_%H%M"),sep="_"), results)
  }
}
