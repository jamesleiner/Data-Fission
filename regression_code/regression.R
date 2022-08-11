###################################################################################################
# Helper functions to execute experiments (Gaussian, Poisson, Logistic Regression and Trend Filtering)
#
# Author(s): Boyan Duan, James Leiner
###################################################################################################

fahrmeir_CIs <- function(infer_model,alpha){
  resid = infer_model$fitted.values - infer_model$y
  Hn_inv <- vcov(infer_model)
  
  
  if(infer_model$family$family == "poisson"){
    V = infer_model$fitted.values
  }
  if(infer_model$family$family == "binomial"){
    V = infer_model$fitted.values *(1-infer_model$fitted.values)
  }
  
  Vn <- t(infer_model$x) %*% diag(resid**2) %*% infer_model$x
  vcov <- Hn_inv  %*% Vn %*% Hn_inv 
  se <- sqrt(diag(vcov))
  CI <- cbind(infer_model$coefficients + qnorm(alpha/2)*se,infer_model$coefficients- qnorm(alpha/2)*se)
  return(CI) 
}


experiment_poisson = function(para_vary){
  source("regression_code/input_regression.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  
  wrapper_func = function(i){
    print(i)
    CIs = list()
    selected_list = list()
    projected = list()
    
    beta= beta*scale
    dat = generate_poisson(n = n, p = p, beta = beta, type = type, rho = rho)
    
    if ("masking" %in% methods) {
      g_Y = sapply(dat$Y, function(x){rbinom(1, size = x, prob = prob)})
      h_Y = dat$Y - g_Y
      
      select_model = cv.glmnet(dat$X, g_Y, family = "poisson")
      selected = which(coef(select_model, s = 'lambda.1se')[-1] != 0) 
      #selected = which(beta !=0)
      
      CIs[["masking"]]$standard = matrix(NA, nrow = p, ncol = 2)
      CIs[["masking"]]$sandwich = matrix(NA, nrow = p, ncol = 2)
      CIs[["masking"]]$fahrmeir = matrix(NA, nrow = p, ncol = 2)
      #CIs[["masking"]]$fahrmeir2 = matrix(NA, nrow = p, ncol = 2)
      if (length(selected) > 0) {
        
        selected_list[["masking"]] = selected
        off = rep(offset(log(1-prob)),length(dat$Y))
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
        #sum1 <- gee(h_Y ~ dat$X[,selected], id=dat$cluster, family=poisson, scale.fix=TRUE)
        #se <- sqrt(diag(sum1$robust.variance))
        #CIs[["masking"]]$fahrmeir2[selected,] <- cbind(infer_model$coefficients + qnorm(alpha/2)*se,infer_model$coefficients- qnorm(alpha/2)*se)[-1,]
        projected[["masking"]] = glm(exp(dat$X%*%beta) ~ dat$X[,selected], family = "poisson")$coefficients[-1]
      }
      
    }
    if ("full" %in% methods) {
      select_model = cv.glmnet(dat$X, dat$Y, family = "poisson")
      selected = which(coef(select_model, s = 'lambda.1se')[-1] != 0) 
      CIs[["full"]]$standard = matrix(NA, nrow = p, ncol = 2)
      CIs[["full"]]$sandwich = matrix(NA, nrow = p, ncol = 2)
      CIs[["full"]]$fahrmeir = matrix(NA, nrow = p, ncol = 2)
      #CIs[["full"]]$fahrmeir2 = matrix(NA, nrow = p, ncol = 2)
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
        #sum1 <- gee(dat$Y ~ dat$X[,selected], id=dat$cluster, family=poisson, scale.fix=TRUE)
        #se <- sqrt(diag(sum1$robust.variance))
        #CIs[["full"]]$fahrmeir2[selected,] <- cbind(infer_model$coefficients + qnorm(alpha/2)*se,infer_model$coefficients- qnorm(alpha/2)*se)[-1,]
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
      #CIs[["split"]]$fahrmeir2 = matrix(NA, nrow = p, ncol = 2)
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
        #sum1 <- gee(dat$Y[-split_ind] ~ dat$X[-split_ind,selected],id=factor(1:length(dat$Y[-split_ind])),  family=poisson,scale.fix=TRUE)
        #se <- sqrt(diag(sum1$robust.variance))
        #CIs[["split"]]$fahrmeir2[selected,] <- cbind(infer_model$coefficients + qnorm(alpha/2)*se,infer_model$coefficients- qnorm(alpha/2)*se)[-1,]
        
        
        projected[["split"]] = glm(exp(dat$X%*%beta) ~ dat$X[,selected], family = "poisson")$coefficients[-1]
      }
    }
    return(list(CIs = CIs, projected=projected,selected=selected_list,beta=beta))
  }
  
  #result_poisson = lapply(1:R, wrapper_func)
  result_poisson = mclapply(1:R, wrapper_func, mc.cores = detectCores(),mc.preschedule = FALSE)
  return(result_poisson)
}



experiment_linear = function(para_vary){
  source("regression_code/input_regression.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    CIs = list()
    selected = list()
    projected = list()
    if(exists("add_influential")){
      dat = generate_linear(n = n, p = p, beta = beta*scale, type = type, rho = rho,add_influential=add_influential)
      n= n + length(add_influential)
    }
    else{
      dat = generate_linear(n = n, p = p, beta = beta*scale, type = type, rho = rho)
    }
    print(dat$cluster)
    if ("masking" %in% methods) {
      noise = rnorm(n, sd = 1)
      g_Y = dat$Y + noise
      h_Y = dat$Y - noise
      
      select_model = cv.glmnet(dat$X, g_Y, family = "gaussian")
      selected[["masking"]] = which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
      if (length(selected[["masking"]]) > 0) {
        infer_model = lm(h_Y ~ dat$X[,selected[["masking"]]])

        if(CI_type == "normal"){
          temp <- confint(infer_model)[-1,]
          CIs[["masking"]] = temp
        }
        else{
          temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]
          CIs[["masking"]] = cbind(temp[,5],temp[,6])
        }
        
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

        
        if(CI_type == "normal"){
          temp <- confint(infer_model)[-1,]
          CIs[["full"]] = temp
        }
        else{
          temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)[-1,]
          CIs[["full"]] = cbind(temp[,5],temp[,6])
        }
        
        
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
        
        if(CI_type == "normal"){
          temp <- confint(infer_model,level=1-alpha)[-1,]
          CIs[["full"]] = temp
        }
        else{
          temp = conf_int(infer_model, vcov = "CR2", cluster = factor(1:(n - sum(split_ind))),
                          level = 1 - alpha)[-1,]
          CIs[["split"]] = cbind(temp[,5],temp[,6])
        }
    
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


projected_beta = function(exp_y, X, beta) {
  if (is.matrix(X)) {
    model = function(para){
      return(t(X)%*%exp_y - t(X)%*%(1/(1 + exp(-X%*%para))))
    }
  } else {
    model = function(para){
      return(sum(X*exp_y) - sum(X*(1/(1 + exp(-X*para)))))
    }
  }
  ss <- tryCatch(multiroot(f = model, start = beta)$root,
                 warning = function(w) {rep(NA, length(beta))})
  return(ss)
}

experiment_logistic = function(para_vary){
  source("regression_code/input_regression.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    CIs = list(); selected_list = list(); projected = list();
    if(exists("add_influential")){
      dat = generate_logistic(n = n, p = p, beta = beta*scale,add_influential=add_influential)
      n=n+length(add_influential)
    }
    else {
      dat = generate_logistic(n = n, p = p, beta = beta*scale) 
    }
    X = cbind(rep(1, n), dat$X)
    beta = c(0, beta)
    if ("masking" %in% methods) {
      CIs[["masking"]]$standard = matrix(NA, nrow = ncol(X), ncol = 2)
      CIs[["masking"]]$sandwich = matrix(NA, nrow = ncol(X), ncol = 2)
      CIs[["masking"]]$fahrmeir = matrix(NA, nrow = ncol(X), ncol = 2)
      
      Z = sapply(dat$Y, function(x){rbinom(1, size = 1, prob = prob)})
      g_Y = (1 - dat$Y)*Z + dat$Y*(1 - Z)
      
      select_model = cv.glmnet(dat$X, g_Y, family = "binomial")
      selected = which(coef(select_model, s = 'lambda.1se') != 0)
      selected_list[["masking"]] = selected
      if (length(selected) > 0) {
        infer_model = glm(dat$Y ~ X[,selected_list[["masking"]]] -1, family = "binomial")
        CIs[["masking"]]$standard[selected,] = confint(infer_model, level = 1 - alpha)
        temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)
        CIs[["masking"]]$sandwich[selected,] = cbind(temp[,5],temp[,6])
        temp = conf_int(infer_model, vcov = "CR0", cluster = dat$cluster, level = 1 - alpha)
        CIs[["masking"]]$fahrmeir[selected,] = cbind(temp[,5],temp[,6])
        
        
        cond_exp_y = dat$exp_y/(dat$exp_y + (1 - dat$exp_y)*(prob/(1 - prob))^(2*g_Y - 1))
        projected[["masking"]] = projected_beta(exp_y = cond_exp_y, X = X[,selected_list[["masking"]]], beta = scale*beta[selected_list[["masking"]]])
      } else {projected[["masking"]] = NA}
    }
    if ("full" %in% methods) {
      CIs[["full"]]$standard = matrix(NA, nrow = ncol(X), ncol = 2)
      CIs[["full"]]$sandwich = matrix(NA, nrow = ncol(X), ncol = 2)
      CIs[["full"]]$fahrmeir = matrix(NA, nrow = ncol(X), ncol = 2)
      
      select_model = cv.glmnet(dat$X, dat$Y, family = "binomial")
      selected = which(coef(select_model, s = 'lambda.1se') != 0)
      selected_list[["full"]] = selected
      if (length(selected) > 0) {
        infer_model = glm(dat$Y ~ X[,selected_list[["full"]]] -1, family = "binomial")
        CIs[["full"]]$standard[selected,] = confint(infer_model, level = 1 - alpha)
        temp = conf_int(infer_model, vcov = "CR2", cluster = dat$cluster, level = 1 - alpha)
        CIs[["full"]]$sandwich[selected,] = cbind(temp[,5],temp[,6])
        temp = conf_int(infer_model, vcov = "CR0", cluster = dat$cluster, level = 1 - alpha)
        CIs[["full"]]$fahrmeir[selected,] = cbind(temp[,5],temp[,6])
        
        projected[["full"]] = projected_beta(
          exp_y = dat$exp_y, X = X[,selected_list[["full"]]], beta = scale*beta[selected_list[["full"]]])
      } else {projected[["full"]] = NA}
    }
    if ('split' %in% methods) {
      CIs[["split"]]$standard = matrix(NA, nrow = ncol(X), ncol = 2)
      CIs[["split"]]$sandwich = matrix(NA, nrow = ncol(X), ncol = 2)
      CIs[["split"]]$fahrmeir = matrix(NA, nrow = ncol(X), ncol = 2)
      # split_ind = sample(1:n, size = n/2)
      split_ind = rbinom(n, 1, 1/2)
      select_model = cv.glmnet(dat$X[split_ind == 1,], dat$Y[split_ind == 1], family = "binomial")
      selected = which(coef(select_model, s = 'lambda.1se') != 0)
      selected_list[["split"]] = selected
      if (length(selected) > 0) {
        infer_model = glm(dat$Y[split_ind == 0] ~ X[split_ind == 0,selected_list[["split"]]] -1,
                          family = "binomial")
        CIs[["split"]]$standard[selected,] = confint(infer_model, level = 1 - alpha)
        temp = conf_int(infer_model, vcov = "CR2", cluster = factor(1:(n - sum(split_ind))),
                        level = 1 - alpha)
  
        CIs[["split"]]$sandwich[selected,] = cbind(temp[,5],temp[,6])
        temp = conf_int(infer_model, vcov = "CR0", cluster = factor(1:(n - sum(split_ind))),
                        level = 1 - alpha)
        CIs[["split"]]$fahrmeir[selected,] = cbind(temp[,5],temp[,6])
  
        
        projected[["split"]] = projected_beta( 
          exp_y = dat$exp_y[split_ind == 0], X = X[split_ind == 0,selected_list[["split"]]],
          beta = scale*beta[selected_list[["split"]]])
       
      } else {projected[["split"]] = NA}
    }
    return(list(CIs = CIs, projected = projected, selected = selected_list,beta = scale*beta))
  }
  #result_binomial = lapply(1:R, wrapper_func)
  result_binomial = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  return(result_binomial)
}




unif_CI = function(n, knots, Y, X, alpha,degree=1){
  k = length(knots) + degree
  if (length(knots) == 0) {
    bs_X = cbind(rep(1, n), 1:n)
  } else if (k + 1 < n) {
    basis = tp(1:n, knots = knots, degree=degree, k = k)
    bs_X = cbind(rep(1, n), basis$X, basis$Z)
  } else {
    stop("too many knots")
  }
  E <- eigen(solve(t(bs_X) %*% bs_X))
  B <- E$vectors %*% diag(sqrt(pmax(0,E$values))) %*% t(E$vectors)
  BX1 <- B %*% t(bs_X[-1, ])
  BX1 <- BX1/sqrt(apply(BX1^2, 2, sum))
  BX0 <- B %*% t(bs_X[-n, ])
  BX0 <- BX0/sqrt(apply(BX0^2, 2, sum))
  kappa <- sum(sqrt(apply((BX1 - BX0)^2, 2, sum)))
  v = n - k - 1
  cvu <- uniroot(function(x) {kappa*(1 + x^2/v)^(-v/2)/pi + 2*(1 - pt(x, df = v)) - alpha},
                 c(0, max(2*kappa/alpha/pi, 10)))$root
  
  true_trend = lm(X ~ bs_X)$fitted.values
  temp = predict(lm(Y ~ bs_X), se.fit = TRUE, level = 1 - alpha)
  width = cvu*temp$se.fit
  return(list(predicted = cbind(temp$fit, temp$fit - width, temp$fit + width, true_trend),
              c = cvu, mean_se = mean(temp$se.fit)))
}


find_knots <- function(fit_y,x,k=1,tol=0.01){
  first_deriv = diff(fit_y)/diff(x)
  second_deriv = diff(first_deriv)/diff(x)[2:(length(x)-1)]
  third_deriv = diff(second_deriv)/diff(x)[3:(length(x)-1)]
  if(k==1) {
    knots = which(abs(second_deriv)>tol) +2
  }
  if(k==2){
    knots = which(abs(third_deriv)>tol) +3
  }  
  return(knots)
}



experiment_trendfilter = function(para_vary){
  source("regression_code/input_regression.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    CI_bend = list()
    project_trend = list()
    c = list()
    mean_se = list()
    nk_selected = list()
    
    dat = time_seq(n = n, slope = slope, p = p, sigma = sigma)
    if ("masking" %in% methods) {
      if(estimate_var){
        sigma_hat = sqrt(mean((dat$Y[-1] - dat$Y[-n])^2)/2)
      } else {
        sigma_hat = sigma 
      }
      noise = rnorm(n, sd = sigma_hat)
      g_Y = dat$Y + noise
      h_Y = dat$Y - noise
      
      if(type=="SURE"){
        tf = sure_trendfilter(1:n,g_Y, weights =rep(1,n),k=deg)
        lambda_min = tf$lambda_min
        lambda_1se= tf$lambda_1se
      } else{
        tf = cv_trendfilter(1:n,g_Y, weights =rep(1,n),k=deg)
        lambda_min = tf$lambda_min[3]
        lambda_1se= tf$lambda_1se[3]
      }
      if(onesd_rule) {
        fit_y = tf$fitted_values[,which(tf$lambda == lambda_1se)]
      }else{
        fit_y = tf$fitted_values[,which(tf$lambda == lambda_min)]
      }
      knots = find_knots(fit_y,1:n,k=deg)
      nk_selected[["masking"]] = length(knots)
      if (CI_type == "uniform") {
        res_masking = tryCatch(unif_CI(n = n, knots = knots, Y = h_Y, X = dat$X, alpha = alpha,degree=deg),
                               error = function(e){
                                 print(e)
                                 return(list(predicted = matrix(NA, nrow = n, ncol = 4),
                                             c = NA, mean_se = NA))})
        CI_bend[["masking"]] = res_masking$predicted[,1:3]
        project_trend[["masking"]] = res_masking$predicted[,4]
        c[["masking"]] = res_masking$c
        mean_se[["masking"]] = res_masking$mean_se
      } else if (CI_type == "pointwise") {
        k= length(knots)+deg
        basis = tp(1:n, knots = knots, degree=deg,k=k)
        bs_X = cbind(rep(1, n), basis$X, basis$Z)
        n_x = ncol(bs_X) 
        temp = predict(lm(h_Y ~ bs_X), interval="confidence", se.fit = TRUE, level = 1 - alpha)
        CI_bend[["masking"]] = temp$fit
        project_trend[["masking"]] = lm(dat$X ~ bs_X)$fitted.values
        c[["masking"]] = qt(1 - alpha/2, n - n_x)
        mean_se[["masking"]] = mean(temp$se.fit)
      }
    }
    
    if ("full" %in% methods) {
      if(estimate_var){
        sigma_hat = sqrt(mean((dat$Y[-1] - dat$Y[-n])^2)/2)
      } else {
        sigma_hat = sigma 
      }
      
      if(type=="SURE"){
        tf = sure_trendfilter(1:n,dat$Y, weights =rep(1,n),k=deg)
        lambda_min = tf$lambda_min
        lambda_1se= tf$lambda_1se
      } else{
        tf = cv_trendfilter(1:n,dat$Y, weights =rep(1,n),k=deg)
        lambda_min = tf$lambda_min[3]
        lambda_1se= tf$lambda_1se[3]
      }
      if(onesd_rule) {
        fit_y = tf$fitted_values[,which(tf$lambda == lambda_1se)]
      }else{
        fit_y = tf$fitted_values[,which(tf$lambda == lambda_min)]
      }
      knots = find_knots(fit_y,1:n,k=deg)
      nk_selected[["full"]] = length(knots)
      if (CI_type == "uniform") {
        res_masking = tryCatch(unif_CI(n = n, knots = knots, Y = dat$Y, X = dat$X, alpha = alpha,degree=deg),
                               error = function(e){
                                 print(e)
                                 return(list(predicted = matrix(NA, nrow = n, ncol = 4),
                                             c = NA, mean_se = NA))})
        CI_bend[["full"]] = res_masking$predicted[,1:3]
        project_trend[["full"]] = res_masking$predicted[,4]
        c[["full"]] = res_masking$c
        mean_se[["full"]] = res_masking$mean_se
      } else if (CI_type == "pointwise") {
        k= length(knots)+deg
        basis = tp(1:n, knots = knots, degree=deg,k=k)
        bs_X = cbind(rep(1, n), basis$X, basis$Z)
        n_x = ncol(bs_X) 
        temp = predict(lm(dat$Y ~ bs_X), interval="confidence", se.fit = TRUE, level = 1 - alpha)
        CI_bend[["full"]] = temp$fit
        project_trend[["full"]] = lm(dat$X ~ bs_X)$fitted.values
        c[["full"]] = qt(1 - alpha/2, n - n_x)
        mean_se[["full"]] = mean(temp$se.fit)
      }
    }
    if ("split" %in% methods) {
      if(estimate_var){
        sigma_hat = sqrt(mean((dat$Y[-1] - dat$Y[-n])^2)/2)
      } else {
        sigma_hat = sigma 
      }
      
      ind_g = seq(from=2,to=length(dat$Y),by=2)
      ind_h = seq(from=1,to=length(dat$Y),by=2)
      
      g_Y = dat$Y[ind_g]
      
      h_Y = dat$Y
      h_Y[ind_g] = NA
      h_Y = na.approx(h_Y)
      
      
      if(type=="SURE"){
        tf = sure_trendfilter(ind_g,g_Y, weights =rep(1,n/2),k=deg)
        lambda_min = tf$lambda_min
        lambda_1se= tf$lambda_1se
      } else{
        tf = cv_trendfilter(ind_g,g_Y, weights =rep(1,n/2),k=deg)
        lambda_min = tf$lambda_min[3]
        lambda_1se= tf$lambda_1se[3]
      }
      if(onesd_rule) {
        fit_y = tf$fitted_values[,which(tf$lambda == lambda_1se)]
      }else{
        fit_y = tf$fitted_values[,which(tf$lambda == lambda_min)]
      }
      knots = ind_g[find_knots(fit_y,ind_g,k=deg)]
      nk_selected[["split"]] = length(knots)
      
      if (CI_type == "uniform") {
        res_masking = tryCatch(unif_CI(n = length(h_Y), knots = knots, Y = h_Y, X = dat$X[1:length(h_Y)], alpha = alpha,degree=deg),
                               error = function(e){
                                 print(e)
                                 return(list(predicted = matrix(NA, nrow = n, ncol = 4),
                                             c = NA, mean_se = NA))})
        CI_bend[["split"]] = res_masking$predicted[ind_h,1:3]
        project_trend[["split"]] = res_masking$predicted[,4][ind_h]
        c[["split"]] = res_masking$c
        mean_se[["split"]] = res_masking$mean_se
      } else if (CI_type == "pointwise") {
        k= length(knots)+deg
        basis = tp(1:length(h_Y), knots = knots, degree=deg,k=k)
        bs_X = cbind(rep(1,length(h_Y) ), basis$X, basis$Z)
        n_x = ncol(bs_X) 
        temp = predict(lm(h_Y ~ bs_X), interval="confidence", se.fit = TRUE, level = 1 - alpha)
        CI_bend[["split"]] = temp$fit[ind_h,]
        project_trend[["split"]] = lm(dat$X[1:length(h_Y)] ~ bs_X)$fitted.values[ind_h]
        c[["split"]] = qt(1 - alpha/2, n - n_x)
        mean_se[["split"]] = mean(temp$se.fit)
      }
    }
    return(list(CI_bend = CI_bend, project_trend = project_trend,real_trend = dat$Y,
                c = c, mean_se = mean_se, sigma_hat = sigma_hat,
                nk_selected = nk_selected, nk_true = length(dat$knots)))
  }
  #result_trendfilter = lapply(1:R, wrapper_func)
  result_trendfilter = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  return(result_trendfilter)
}

