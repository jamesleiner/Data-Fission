source("input_regression.R", local = TRUE)
for (single_para_vary in para_vary) {
  assign(single_para_vary$name, single_para_vary$value)
}

para_vary = list(list(name = "scale", value = scale),
                 list(name = "R", value = 10),
                 list(name = "p", value = p),
                 list(name="rho", value = 0),
                 list(name = "beta", 
                      value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))))
CIs = list()
selected = list()
projected = list()

dat = generate_linear(n = n, p = p, beta = beta*scale, type = type, rho = rho)

if ("masking" %in% methods) {
  noise = rnorm(n, sd = 1)
  g_Y = dat$Y + noise
  h_Y = dat$Y - noise
  
  select_model = cv.glmnet(dat$X, g_Y, family = "gaussian")
  selected[["masking"]] = which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
  if (length(selected[["masking"]]) > 0) {
    infer_model = lm(h_Y ~ dat$X[,selected[["masking"]]])
    temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]
    CIs[["masking"]] = cbind(temp[,4],temp[,5])
    
    projected[["masking"]] = beta[selected[["masking"]]]*scale +
      scale*beta[-selected[["masking"]]] %*% 
      dat$Sigma[-selected[["masking"]], selected[["masking"]]] %*%
      solve(dat$Sigma[selected[["masking"]], selected[["masking"]]]) 
  } else {CIs[["masking"]] = NA;  projected[["masking"]] = NA; selected[["masking"]] = NA}
}





experiment_linear = function(para_vary){
  source("input_regression.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    CIs = list()
    selected = list()
    projected = list()
    
    dat = generate_linear(n = n, p = p, beta = beta*scale, type = type, rho = rho)
    if ("masking" %in% methods) {
      noise = rnorm(n, sd = 1)
      g_Y = dat$Y + noise
      h_Y = dat$Y - noise
      
      select_model = cv.glmnet(dat$X, g_Y, family = "gaussian")
      selected[["masking"]] = which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
      if (length(selected[["masking"]]) > 0) {
        infer_model = lm(h_Y ~ dat$X[,selected[["masking"]]])
        temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]
        CIs[["masking"]] = cbind(temp[,4],temp[,5])
        
        projected[["masking"]] = beta[selected[["masking"]]]*scale +
          scale*beta[-selected[["masking"]]] %*% 
          dat$Sigma[-selected[["masking"]], selected[["masking"]]] %*%
          solve(dat$Sigma[selected[["masking"]], selected[["masking"]]]) 
      } else {CIs[["masking"]] = NA;  projected[["masking"]] = NA; selected[["masking"]] = NA}
    }
    if ("full" %in% methods) {
      select_model = cv.glmnet(dat$X, dat$Y, family = "gaussian")
      selected[["full"]] = which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
      if (length(selected[["full"]]) > 0) {
        infer_model = lm(dat$Y ~ dat$X[,selected[["full"]]])
        temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]
        CIs[["full"]] = cbind(temp[,4],temp[,5])
        
        projected[["full"]] = beta[selected[["full"]]]*scale +
          scale*beta[-selected[["full"]]] %*% 
          dat$Sigma[-selected[["full"]], selected[["full"]]] %*%
          solve(dat$Sigma[selected[["full"]], selected[["full"]]])
      } else {CIs[["full"]] = NA;  projected[["full"]] = NA; selected[["full"]] = NA}
    }
    if ('split' %in% methods) {
      split_ind = rbinom(n, 1, 1/2)
      select_model = cv.glmnet(dat$X[split_ind == 1,], dat$Y[split_ind == 1], family = "gaussian")
      selected[["split"]] = which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
      if (length(selected[["split"]]) > 0) {
        infer_model = lm(dat$Y[split_ind == 0] ~ dat$X[split_ind == 0,selected[["split"]]])
        temp = conf_int(infer_model, vcov = "CR2", cluster = factor(1:(n - sum(split_ind))),
                        level = 1 - alpha)[-1,]
        CIs[["split"]] = cbind(temp[,4],temp[,5])
        
        projected[["split"]] = beta[selected[["split"]]]*scale +
          scale*beta[-selected[["split"]]] %*% 
          dat$Sigma[-selected[["split"]], selected[["split"]]] %*%
          solve(dat$Sigma[selected[["split"]], selected[["split"]]])
      } else {CIs[["split"]] = NA; selected[["split"]] = NA; projected[["split"]] = NA}
    }
    if ('test' %in% methods) {
      infer_model = lm(Y~ X)
      CIs[["test"]] = confint(infer_model, level = 1 - alpha)[-1,]
    }
    return(list(CIs = CIs, selected = selected, projected = projected))
  }
  # result_linear = lapply(1:R, wrapper_func)
  result_linear = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  return(result_linear)
}


t = experiment_linear(para_vary)
