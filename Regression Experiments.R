
source("regression_code/set_up_regression.R")

#######################################################################################################################
# Linear Regression Experiments
#######################################################################################################################


scale_seq = seq(0, 0.2, length.out = 5)
p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent.Rdata",sep = ""))
  
dat <- generate_linear(500,100, c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9)),"dependent",0.3)
select_model = cv.glmnet(dat$X, dat$Y, family = "gaussian")

rho_seq = seq(-0.5, 0.5, length.out = 5)
p = 100
result = list()
for (rho in rho_seq) { 
  print(rho)
  para_vary = list(list(name = "rho", value = rho),
                   list(name = "scale", value = 0.1),
                   list(name = "R", value = 500),
                   list(name = "p", value = 100),
                   list(name = "type", value = "dependent"),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = rep(c(1, 0, 0, 1, 0, 0, 1, 0, 0, 0), 10)
                   ))
  result[[as.character(rho)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_dependent.Rdata",sep = ""))



influ_seq = seq(2, 6, length.out = 5)
result=list()
for (influ in influ_seq) { 
  para_vary = list(list(name = "rho", value = 0),
                   list(name = "scale", value = 1),
                   list(name = "R", value = 500),
                   list(name = "n", value = 15),
                   list(name = "p", value = 20),
                   list(name = "type", value = "independent"),
                   list(name = "beta", 
                        value = c(1, rep(0,15),1,-1,1,0)),
                   list(name="add_influential",value=c(influ))
  )
  experiment_linear(para_vary)
  
  result[[as.character(influ)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_influential.Rdata",sep = ""))


#######################################################################################################################
# Linear Regression with Misspecified Models
#######################################################################################################################

scale_seq = seq(0, 0.2, length.out = 5)
p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "est_var", value = TRUE),
                   list(name = "error_type", value = "t"),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent_t_estvar.Rdata",sep = ""))

p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "est_var", value = FALSE),
                   list(name = "error_type", value = "t"),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent_t.Rdata",sep = ""))



p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "est_var", value = TRUE),
                   list(name = "error_type", value = "sn"),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent_sn_estvar.Rdata",sep = ""))

p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "est_var", value = FALSE),
                   list(name = "error_type", value = "sn"),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent_sn.Rdata",sep = ""))

p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "est_var", value = TRUE),
                   list(name = "error_type", value = "laplace"),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent_laplace_estvar.Rdata",sep = ""))

p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "est_var", value = FALSE),
                   list(name = "error_type", value = "laplace"),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent_laplace.Rdata",sep = ""))



scale_seq = seq(0, 0.2, length.out = 5)
p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "est_var", value = TRUE),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_linear(para_vary)
}
save(result, file=paste("results/regression_linear_independent_estvar.Rdata",sep = ""))


#######################################################################################################################
# Poisson Regression Experiments
#######################################################################################################################

influ_seq = seq(1, 6, length.out = 6)
result=list()
for (influ in influ_seq) { 
  para_vary = list(list(name = "scale", value = 1),
                   list(name = "R", value = 500),
                   list(name = "n", value = 15),
                   list(name = "p", value = 20),
                   list(name = "prob",value = 0.2),
                   list(name = "beta", 
                        value = c(1, rep(0,15),1,-1,1,0)),
                   list(name="add_influential",value=c(influ)))
  
  result[[as.character(influ)]] = experiment_poisson(para_vary)
}

save(result, file=paste("results/regression_poisson_influential.Rdata",sep = ""))


scale_seq = seq(0, 0.5, length.out = 5)
p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "n", value = 1000),
                   list(name = "R", value = 500),
                   list(name = 'prob',value=0.5),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 30), rep(2,8))
                   )) 
  result[[as.character(scale)]] = experiment_poisson(para_vary)
}
save(result, file=paste("results/regression_poisson_varyscale.Rdata",sep = ""))



rho_seq = seq(-0.5, 0.5, length.out = 5)
result = list()
for (rho in rho_seq) { 
  p=100
  para_vary = list(list(name = "rho", value = rho),
                   list(name = "scale", value = 0.25),
                   list(name = "n", value = 1000),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = 'prob',value=0.2),
                   list(name = "type", value = "dependent"),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 30), rep(2,8))
                   ))
  result[[as.character(rho)]] = experiment_poisson(para_vary)
}
save(result, file=paste("results/regression_poisson_varyRho.Rdata",sep = ""))


if(0){
  n_seq = c(10,100,1000,10000,100000)
  p = 100
  result = list()
  for (n in n_seq) { 
    print(n)
    para_vary = list(list(name = "scale", value = 0.5),
                     list(name = "n", value = n),
                     list(name = "R", value = 500),
                     list(name = 'prob',value=0.5),
                     list(name = "p", value = p),
                     list(name = "beta", 
                          value = c(1, 0, rep(1,20), rep(0, p - 30), rep(2,8))
                     )) 
    result[[as.character(log(n,base=10))]] = experiment_poisson(para_vary)
  }
  save(result, file=paste("results/regression_poisson_varynum.Rdata",sep = ""))
}

#######################################################################################################################
# Logistic Regression Experiments
#######################################################################################################################




scale_seq = seq(0, 0.5, length.out = 6)
p = 100
result = list()
for (scale in scale_seq) { 
  print(scale)
  para_vary = list(list(name = "scale", value = scale),
                   list(name = "prob", value = 0.2),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(2,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(scale)]] = experiment_logistic(para_vary)
}
save(result, file=paste("results/regression_logistic_varyscale.Rdata",sep = ""))




prob_seq = seq(0.05, 0.45, 0.05)
p = 100
result = list()
for (prob in prob_seq) { 
  print(prob)
  para_vary = list(list(name = "scale", value = 0.5),
                   list(name = "prob", value = prob),
                   list(name = "R", value = 500),
                   list(name = "p", value = p),
                   list(name = "beta", 
                        value = c(1, 0, rep(1,20), rep(0, p - 31), rep(2,9))
                        # value = c(1, 0, rep(1,5), rep(0, p - 10), rep(2,3))
                        # value = rep(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 10))
                   )) 
  result[[as.character(prob)]] = experiment_logistic(para_vary)
}
save(result, file=paste("results/regression_logistic_varyprob.Rdata",sep = ""))


influ_seq = seq(2, 6, length.out = 5)
result=list()
for (influ in influ_seq) { 
  para_vary = list(list(name = "scale", value = 1),
                   list(name = "R", value = 5),
                   list(name = "n", value = 15),
                   list(name = "p", value = 20),
                   list(name = "prob",value = 0.2),
                   list(name = "beta", 
                        value = c(1, rep(0,15),1,-1,1,0)),
                   list(name="add_influential",value=c(influ)))
  
  result[[as.character(influ)]] = experiment_logistic(para_vary)
}

save(result, file=paste("results/regression_logistic_influential.Rdata",sep = ""))


#######################################################################################################################
# Trend Filter Experiments
#######################################################################################################################


sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
result = list()
for (sigma in sigma_seq) {
  result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    print(prob)
    para_vary = list(list(name = "n", value = 200),
                     list(name = "p", value = prob),
                     list(name = "slope", value = 0.5),
                     list(name = "sigma", value = sigma),
                     list(name = "CI_type", value = "uniform"),
                     list(name = "R", value = 100),
                     list(name = "estimate_var", value = 0),
                     list(name="onesd_rule",value = 0),
                     list(name="type",value="SURE"),
                     list(name="deg",value=1)) 
    result[[as.character(sigma)]][[as.character(prob)]] = experiment_trendfilter(para_vary)
  }
}
save(result, file=paste("results/regression_trendfilter_varySigmaProb_uniform_SURE.Rdata",sep = ""))

sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
result = list()
for (sigma in sigma_seq) {
  result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    print(prob)
    para_vary = list(list(name = "n", value = 200),
                     list(name = "p", value = prob),
                     list(name = "slope", value = 0.5),
                     list(name = "sigma", value = sigma),
                     list(name = "CI_type", value = "pointwise"),
                     list(name = "R", value = 100),
                     list(name = "estimate_var", value = 0),
                     list(name="onesd_rule",value = 0),
                     list(name="type",value="SURE"),
                     list(name="deg",value=1)) 
    result[[as.character(sigma)]][[as.character(prob)]] = experiment_trendfilter(para_vary)
  }
}
save(result, file=paste("results/regression_trendfilter_varySigmaProb_pointwise_SURE.Rdata",sep = ""))


sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
result = list()
for (sigma in sigma_seq) {
  result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    print(prob)
    para_vary = list(list(name = "n", value = 200),
                     list(name = "p", value = prob),
                     list(name = "slope", value = 0.5),
                     list(name = "sigma", value = sigma),
                     list(name = "CI_type", value = "uniform"),
                     list(name = "R", value = 100),
                     list(name = "estimate_var", value = 0),
                     list(name="onesd_rule",value = 0),
                     list(name="type",value="CV"),
                     list(name="deg",value=1)) 
    result[[as.character(sigma)]][[as.character(prob)]] = experiment_trendfilter(para_vary)
  }
}
save(result, file=paste("results/regression_trendfilter_varySigmaProb_uniform_CV.Rdata",sep = ""))

sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
result = list()
for (sigma in sigma_seq) {
  result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    print(prob)
    para_vary = list(list(name = "n", value = 200),
                     list(name = "p", value = prob),
                     list(name = "slope", value = 0.5),
                     list(name = "sigma", value = sigma),
                     list(name = "CI_type", value = "pointwise"),
                     list(name = "R", value = 100),
                     list(name = "estimate_var", value = 0),
                     list(name="onesd_rule",value = 0),
                     list(name="type",value="CV"),
                     list(name="deg",value=1)) 
    result[[as.character(sigma)]][[as.character(prob)]] = experiment_trendfilter(para_vary)
  }
}
save(result, file=paste("results/regression_trendfilter_varySigmaProb_pointwise_CV.Rdata",sep = ""))



sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
result = list()
for (sigma in sigma_seq) {
  result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    print(prob)
    para_vary = list(list(name = "n", value = 200),
                     list(name = "p", value = prob),
                     list(name = "slope", value = 0.5),
                     list(name = "sigma", value = sigma),
                     list(name = "CI_type", value = "uniform"),
                     list(name = "R", value = 100),
                     list(name = "estimate_var", value = 1),
                     list(name="onesd_rule",value = 0),
                     list(name="type",value="SURE"),
                     list(name="deg",value=1)) 
    result[[as.character(sigma)]][[as.character(prob)]] = experiment_trendfilter(para_vary)
  }
}
save(result, file=paste("results/regression_trendfilter_varySigmaProb_uniform_SURE_estimate.Rdata",sep = ""))

sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
result = list()
for (sigma in sigma_seq) {
  result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    print(prob)
    para_vary = list(list(name = "n", value = 200),
                     list(name = "p", value = prob),
                     list(name = "slope", value = 0.5),
                     list(name = "sigma", value = sigma),
                     list(name = "CI_type", value = "pointwise"),
                     list(name = "R", value = 100),
                     list(name = "estimate_var", value = 1),
                     list(name="onesd_rule",value = 0),
                     list(name="type",value="SURE"),
                     list(name="deg",value=1)) 
    result[[as.character(sigma)]][[as.character(prob)]] = experiment_trendfilter(para_vary)
  }
}
save(result, file=paste("results/regression_trendfilter_varySigmaProb_pointwise_SURE_estimate.Rdata",sep = ""))

