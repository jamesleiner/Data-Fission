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
  #library(glmgen)
  library(trendfiltering)
  library(dplyr)
  library(magrittr)
  library(tidyr)
  library(rlang)
  library(purrr)
  library(parallel)
  library(matrixStats)
  library(zoo)
})
source("regression.R")
source("sim_dat_regression.R")
p = 20
influ=1
para_vary = list(list(name = "scale", value = scale),
                 list(name = "R", value = 100),
                 list(name = "n", value = 25),
                 list(name = "p", value = p),
                 list(name = "prob",value = 0.2),
                 list(name = "beta", 
                      value = c(1, rep(0,15),1,-1,1,0)),
                 list(name="add_influential",value=c(influ)))

h = experiment_poisson(para_vary)

source("input_regression.R", local = TRUE)
for (single_para_vary in para_vary) {
  assign(single_para_vary$name, single_para_vary$value)
}

CIs = list()
selected_list = list()
projected = list()

dat = generate_poisson(n = n, p = p, beta = beta*scale, type = type, rho = rho)

if ("masking" %in% methods) {
  g_Y = sapply(dat$Y, function(x){rbinom(1, size = x, prob = prob)})
  h_Y = dat$Y - g_Y
  
  select_model = cv.glmnet(dat$X, g_Y, family = "poisson")
  selected = which(coef(select_model, s = 'lambda.1se')[-1] != 0) 
  
  CIs[["masking"]]$standard = matrix(NA, nrow = p, ncol = 2)
  CIs[["masking"]]$sandwich = matrix(NA, nrow = p, ncol = 2)
  CIs[["masking"]]$fahrmeir = matrix(NA, nrow = p, ncol = 2)
  if (length(selected) > 0) {
    
    selected_list[["masking"]] = selected
    off = rep(offset(log(1/(1-prob))),length(dat$Y))
    infer_model = glm(h_Y ~ dat$X[,selected], family = "poisson",x=TRUE,y=TRUE, offset = off)
    temp = tryCatch({confint(infer_model, level = 1 - alpha)[-1,] },
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 2) ) })
    CIs[["masking"]]$standard[selected,] = temp
    temp = tryCatch({conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]},
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 6) ) })
    CIs[["masking"]]$sandwich[selected,] = cbind(temp[,5],temp[,6])
    temp = tryCatch({ fahrmeir_CIs(infer_model,alpha)[-1,]},
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 2) ) })
    CIs[["masking"]]$fahrmeir[selected,] = temp
    projected[["masking"]] = glm(exp(dat$X%*%beta) ~ dat$X[,selected], family = "poisson")$coefficients[-1]
  }
  
}
if ("full" %in% methods) {
  select_model = cv.glmnet(dat$X, dat$Y, family = "poisson")
  selected = which(coef(select_model, s = 'lambda.1se')[-1] != 0) 
  CIs[["full"]]$standard = matrix(NA, nrow = p, ncol = 2)
  CIs[["full"]]$sandwich = matrix(NA, nrow = p, ncol = 2)
  CIs[["full"]]$fahrmeir = matrix(NA, nrow = p, ncol = 2)
  
  if (length(selected) > 0) {
    selected_list[["full"]] = selected
    infer_model = glm(dat$Y ~ dat$X[,selected], family = "poisson",x=TRUE,y=TRUE)
    temp = tryCatch({confint(infer_model, level = 1 - alpha)[-1,] },
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 2) ) })
    CIs[["full"]]$standard[selected,] = temp
    temp = tryCatch({conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]},
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 6) ) })
    CIs[["full"]]$sandwich[selected,] = cbind(temp[,5],temp[,6])
    temp = tryCatch({ fahrmeir_CIs(infer_model,alpha)[-1,]},
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 2) ) })
    CIs[["full"]]$fahrmeir[selected,] = temp
    projected[["full"]] = glm(exp(dat$X%*%beta) ~ dat$X[,selected], family = "poisson")$coefficients[-1]
  }
}
if ('split' %in% methods) {
  split_ind = sample(1:n, size = n/2)
  select_model = cv.glmnet(dat$X[split_ind,], dat$Y[split_ind], family = "poisson")
  selected = which(coef(select_model, s = 'lambda.1se')[-1] != 0) 
  CIs[["split"]]$standard = matrix(NA, nrow = p, ncol = 2)
  CIs[["split"]]$sandwich = matrix(NA, nrow = p, ncol = 2)
  CIs[["split"]]$fahrmeir = matrix(NA, nrow = p, ncol = 2)
  if (length(selected) > 0) {
    selected_list["split"] = selected
    infer_model = glm(dat$Y[-split_ind] ~ dat$X[-split_ind,selected],x=TRUE,y=TRUE,family = "poisson")
    temp = tryCatch({confint(infer_model, level = 1 - alpha)[-1,] },
             error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 2) ) })
    CIs[["split"]]$standard[selected,] = temp
    temp = tryCatch({conf_int(infer_model, vcov = "CR2", cluster = factor(1:length(dat$Y[-split_ind])), level = 1 - alpha)[-1,]},
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 6) ) })
    CIs[["split"]]$sandwich[selected,] = cbind(temp[,5],temp[,6])
    temp = tryCatch({ fahrmeir_CIs(infer_model,alpha)[-1,]},
                    error=function(cond){ return(matrix(NA, nrow = length(selected), ncol = 2) ) })
    CIs[["split"]]$fahrmeir[selected,] = temp
    projected[["split"]] = glm(exp(dat$X%*%beta) ~ dat$X[,selected], family = "poisson")$coefficients[-1]
  }
}
return(list(CIs = CIs, projected=projected,selected=selected_list,beta=beta))






scale = 2
p = 20
para_vary = list(list(name = "scale", value = scale),
                 list(name = "prob", value = 0.2),
                 list(name = "R", value = 2),
                 list(name = "n", value = 30),
                 list(name = "p", value = p),
                 list(name = "beta", 
                      value = c(1, rep(0,15),1,-1,1,0)),
                 list(name="add_influential",value=c(influ)))

scale=0.5
p=100
n=1000
para_vary = list(list(name = "scale", value = scale),
                 list(name = "prob", value = 0.2),
                 list(name = "R", value = 100),
                 list(name = "p", value = p),
                 list(name = "beta", 
                      value = c(1, 0, rep(1,20), rep(0, p - 31), rep(2,9))
                      # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                      # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                 )) 

source("input_regression.R", local = TRUE)
for (single_para_vary in para_vary) {
  assign(single_para_vary$name, single_para_vary$value)
}



CIs = list(); selected_list = list();  projected = list();

dat = generate_logistic(n = n, p = p, beta = beta*scale)

if ("masking" %in% methods) {
  Z = sapply(dat$Y, function(x){rbinom(1, size = 1, prob = prob)})
  g_Y = (1 - dat$Y)*Z + dat$Y*(1 - Z)
  
  select_model = cv.glmnet(dat$X, g_Y, family = "binomial")
  selected = which(coef(select_model, s = 'lambda.1se')[-1] != 0)
  
  CIs[["masking"]]$standard = matrix(NA, nrow = p, ncol = 2)
  CIs[["masking"]]$sandwich = matrix(NA, nrow = p, ncol = 2)
  CIs[["masking"]]$fahrmeir = matrix(NA, nrow = p, ncol = 2)

  if (length(selected) > 0) {
    selected_list[["masking"]] = selected
    
    infer_model = glm(dat$Y ~ dat$X[,selected], family = "binomial",x=TRUE,y=TRUE)
    CIs[["masking"]]$standard[selected,] = confint(infer_model, level = 1 - alpha)[-1,]
    temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]
    CIs[["masking"]]$sandwich[selected,] = cbind(temp[,5],temp[,6])
    CIs[["masking"]]$fahrmeir[selected,] = fahrmeir_CIs(infer_model,alpha)[-1,]
    
    cond_exp_y = dat$exp_y/(dat$exp_y + (1 - dat$exp_y)*(prob/(1 - prob))^(2*g_Y - 1))
    X = cbind(rep(1, nrow(dat$X)), dat$X)
    projected[["masking"]] = projected_beta(
      exp_y = cond_exp_y, X = X[,selected], beta = scale*beta[selected])
  }
  
}

X = cbind(rep(1, nrow(dat$X)), dat$X)
projected[["masking"]] = projected_beta(
  exp_y = cond_exp_y, X = dat$X[,selected], beta = scale*beta[selected])










experiment_logistic = function(para_vary){
  source("input_regression.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    CIs = list(); selected = list(); nCIs = list(); projected = list(); ll_ratio = list()
    
    dat = generate_logistic(n = n, p = p, beta = beta*scale)
    sim_dat = generate_logistic(n = n, p = p, beta = beta*scale)
    X = cbind(rep(1, n), dat$X); sim_X = cbind(rep(1,n), sim_dat$X); beta = c(0, beta)
    if ("masking" %in% methods) {
      Z = sapply(dat$Y, function(x){rbinom(1, size = 1, prob = prob)})
      g_Y = (1 - dat$Y)*Z + dat$Y*(1 - Z)
      
      select_model = cv.glmnet(dat$X, g_Y, family = "binomial")
      selected[["masking"]] = which(coef(select_model, s = 'lambda.1se') != 0)
      if (length(selected[["masking"]]) > 0) {
        infer_model = glm(dat$Y ~ X[,selected[["masking"]]] -1, family = "binomial")
        temp = conf_int(infer_model, vcov = "CR0", cluster = dat$cluster, level = 1 - alpha)
        CIs[["masking"]] = cbind(temp[,5],temp[,6])
        # cov_vec[["masking"]] = temp[,2]
        
        sim_infer_model = glm(sim_dat$Y~ sim_X[,selected[["masking"]]] -1, family = "binomial")
        temp = conf_int(sim_infer_model, vcov = "CR0", cluster = dat$cluster, level = 1 - alpha)
        nCIs[["masking"]] = cbind(temp[,5],temp[,6])
        
        # B = diag(bread(infer_model)); M = diag(meatHC(infer_model))
        # sim_B = diag(bread(sim_infer_model)); sim_M = diag(meatHC(sim_infer_model))
        ratio_B = diag(bread(infer_model))/diag(bread(sim_infer_model))
        ratio_M = diag(meatHC(infer_model))/diag(meatHC(sim_infer_model))
        
        cond_exp_y = dat$exp_y/(dat$exp_y + (1 - dat$exp_y)*(prob/(1 - prob))^(2*g_Y - 1))
        projected[["masking"]] = projected_beta(
          exp_y = cond_exp_y, X = X[,selected[["masking"]]], beta = scale*beta[selected[["masking"]]])
        if (length(selected[["masking"]]) == 1) {
          projected_exp = 1/(1 + exp(-X[,selected[["masking"]]]*projected[["masking"]]))
          fi_selected = X[,selected[["masking"]]] *
            as.vector(projected_exp/(1 + projected_exp)^2) * X[,selected[["masking"]]]
          fi_sim = sim_X[,selected[["masking"]]] * 
            as.vector(sim_dat$exp_y/(1 + sim_dat$exp_y)^2) * sim_X[,selected[["masking"]]]
        } else {
          projected_exp = 1/(1 + exp(-X[,selected[["masking"]]]%*%projected[["masking"]]))
          fi_selected = diag(t(X[,selected[["masking"]]]) %*%
                               diag(as.vector(projected_exp/(1 + projected_exp)^2)) %*% X[,selected[["masking"]]])
          fi_sim = diag(t(sim_X[,selected[["masking"]]]) %*%
                          diag(as.vector(sim_dat$exp_y/(1 + sim_dat$exp_y)^2)) %*% sim_X[,selected[["masking"]]])
        }
        ratio_fi = fi_selected/fi_sim
        ll_ratio[["masking"]] = sum(dat$Y*projected_exp/dat$exp_y +
                                      (1 - dat$Y)*(1 - projected_exp)/(1 - dat$exp_y))
      } else {CIs[["masking"]] = NA; nCIs[["masking"]] = NA; projected[["masking"]] = NA;
      ratio_B = NA; ratio_M = NA; projected[["masking"]] = NA; ratio_fi = NA}
    }
    if ("full" %in% methods) {
      select_model = cv.glmnet(dat$X, dat$Y, family = "binomial")
      selected[["full"]] = which(coef(select_model, s = 'lambda.1se') != 0)
      if (length(selected[["full"]]) > 0) {
        infer_model = glm(dat$Y ~ X[,selected[["full"]]] -1, family = "binomial")
        temp = conf_int(infer_model, vcov = "CR0", cluster = dat$cluster, level = 1 - alpha)
        CIs[["full"]] = cbind(temp[,5],temp[,6])
        # cov_vec[["full"]] = temp[,2]
        
        sim_infer_model = glm(sim_dat$Y~ sim_X[,selected[["full"]]] -1, family = "binomial")
        temp = conf_int(sim_infer_model, vcov = "CR0", cluster = dat$cluster, level = 1 - alpha)
        nCIs[["full"]] = cbind(temp[,5],temp[,6])
        
        projected[["full"]] = projected_beta(
          exp_y = dat$exp_y, X = X[,selected[["full"]]], beta = scale*beta[selected[["full"]]])
        projected[["sim"]] = projected_beta(
          exp_y = sim_dat$exp_y, X = sim_X[,selected[["full"]]], beta = scale*beta[selected[["full"]]])
        if (length(selected[["full"]]) == 1) {
          projected_exp = 1/(1 + exp(-X[,selected[["full"]]]*projected[["full"]]))
        } else {
          projected_exp = 1/(1 + exp(-X[,selected[["full"]]]%*%projected[["full"]]))
        }
        ll_ratio[["full"]] = sum(dat$Y*projected_exp/dat$exp_y +
                                   (1 - dat$Y)*(1 - projected_exp)/(1 - dat$exp_y))
      } else {CIs[["full"]] = NA; nCIs[["full"]] = NA; projected[["full"]] = NA; projected[["sim"]] = NA}
    }
    if ('split' %in% methods) {
      # split_ind = sample(1:n, size = n/2)
      split_ind = rbinom(n, 1, 1/2)
      select_model = cv.glmnet(dat$X[split_ind == 1,], dat$Y[split_ind == 1], family = "binomial")
      selected[["split"]] = which(coef(select_model, s = 'lambda.1se') != 0)
      if (length(selected[["split"]]) > 0) {
        infer_model = glm(dat$Y[split_ind == 0] ~ X[split_ind == 0,selected[["split"]]] -1,
                          family = "binomial")
        temp = conf_int(infer_model, vcov = "CR0", cluster = factor(1:(n - sum(split_ind))),
                        level = 1 - alpha)
        CIs[["split"]] = cbind(temp[,5],temp[,6])
        # cov_vec[["split"]] = temp[,2]
        
        sim_infer_model = glm(sim_dat$Y~ sim_X[,selected[["split"]]] -1, family = "binomial")
        temp = conf_int(sim_infer_model, vcov = "CR0", cluster = dat$cluster, level = 1 - alpha)
        nCIs[["split"]] = cbind(temp[,5],temp[,6])
        
        projected[["split"]] = projected_beta( 
          exp_y = dat$exp_y[split_ind == 0], X = X[split_ind == 0,selected[["split"]]],
          beta = scale*beta[selected[["split"]]])
        if (length(selected[["split"]]) == 1) {
          projected_exp = 1/(1 + exp(-X[,selected[["split"]]]*projected[["split"]]))
        } else {
          projected_exp = 1/(1 + exp(-X[,selected[["split"]]]%*%projected[["split"]]))
        }
        ll_ratio[["split"]] = sum(dat$Y*projected_exp/dat$exp_y +
                                    (1 - dat$Y)*(1 - projected_exp)/(1 - dat$exp_y))
      } else {CIs[["split"]] = NA; nCIs[["split"]] = NA; projected[["split"]] = NA}
    }
    if ('test' %in% methods) {
      infer_model = glm(dat$Y~ dat$X, family = "binomial")
      # temp = conf_int(infer_model, vcov = "CR2", cluster = factor(1:n), level = 1 - alpha)[-1,]
      # CIs[["test"]] = cbind(temp[,5],temp[,6])
      CIs[["test"]] = confint(infer_model, level = 1 - alpha)[-1,]
      sd = summary(infer_model)$coefficients[2,2]
      est = summary(infer_model)$coefficients[2,1]
      g_Y = NA
    }
    return(list(CIs = CIs, nCIs = nCIs, projected = projected, ll_ratio = ll_ratio,
                ratio_B = ratio_B, ratio_M = ratio_M, ratio_fi = ratio_fi, selected = selected
                #sd = sd, est = est, X = cbind(rep(1, n), dat$X), split_ind = split_ind, g_Y = g_Y, Y = dat$Y
    ))
  }
  # result_binomial = lapply(1:R, wrapper_func)
  result_binomial = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  return(result_binomial)
}



