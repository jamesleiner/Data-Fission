
source("trendfiltering-stats-scripts/trendfilter.R")
source("trendfiltering-stats-scripts/utility-functions.R")
source("trendfiltering-stats-scripts/hyperparameter-tuning.R")




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


p=0.99
para_vary = list(list(name = "n", value = 200),
                 list(name = "p", value = prob),
                 list(name = "slope", value = 0.5),
                 list(name = "sigma", value = 0.5),
                 list(name = "CI_type", value = "uniform"),
                 list(name = "R", value = 2),
                 list(name = "estimate_var", value = 0),
                 list(name="onesd_rule",value = 0),
                 list(name="type",value="CV"),
                 list(name="deg",value=1)) 

test = experiment_trendfilter(para_vary)


CI_bend = list()
project_trend = list()
c = list()
mean_se = list()
nk_selected = list()
source("input_regression.R", local = TRUE)
for (single_para_vary in para_vary) {
  assign(single_para_vary$name, single_para_vary$value)
}



experiment_trendfilter = function(para_vary){
  source("input_regression.R", local = TRUE)
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



