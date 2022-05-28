library("MASS")
library("mgcv")
library("RCurl")

## AdaPT
source("STAR_code/AdaPT.R")
source("STAR_code/AdaPT_gam.R")
source('STAR_code/summarize_methods.R')
source("STAR_code/expr_func.R")
source('STAR_code/STAR_convex.R')

wd = paste(getwd(),"/results/",sep="")

iter_run <- function(x,mu,tau,alpha.list,num.steps.update.score,scope.params,filename,repeat.times=500,base=1,type="normal"){
  results <- list()
  for(i in 1:repeat.times){
    print(i)
    results[[i]] <- tryCatch(experiment_masked(x,mu,tau,alpha.list,num.steps.update.score,scope.params,base=base,type=type) ,
                             error = function(e){
                               print(e)
                               return(list())})
    save(file = filename, results)
  }
}


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


experiment_masked <- function(x,mu,tau,alpha.list,num.steps.update.score,scope.params,base=1,type="normal"){
  if(type == "normal"){
    var <- sig.gen(n,mu,0,1)
    Y <- mvrnorm(1,mu,var)
    Z <- mvrnorm(1,rep(0,n),var)
    f_Y <- Y+tau*Z
    g_Y <- Y-(1/tau)*Z
    pvals <- 1 - pnorm(Y)
    pvals_mask <- 1- pnorm(f_Y,sd=sqrt((1+tau**2)))
  }
  if(type == "poisson"){
    Y <- rpois(length(mu),base+mu)
    Z <- rbinom(length(Y),Y,tau)
    g_Y = Y - Z
    f_Y = Z
    Y_minus = sapply(Y,function(x) max(x-1,0))
    f_Y_minus = sapply(f_Y,function(x) max(x-1,0))
    U1 <- runif(length(mu),0,1)
    U2 <- runif(length(mu),0,1)
    rand <- U1*ppois(Y, base) + (1-U1)* ppois(Y_minus, base)
    rand_F <- U2*ppois(f_Y, base*tau) + (1-U2)* ppois(f_Y_minus, base*tau)
    pvals <- 1-rand
    pvals_mask <- 1 - rand_F
  }
  
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
  
  browser()
  
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
  return(list(x=x,mu=mu,tau=tau,H0=H0,mu=mu,var=var,var_mask=var_mask,Y=Y, f_Y=f_Y,g_Y=g_Y,pvals=pvals,pvals_mask=pvals_mask,
              BH.result.mask=BH.result_mask, BH.result.full=BH.result_full,
              adapt.result.mask=adapt.result_mask,adapt.result.full=adapt.result_full,
              STAR.result.mask=STAR.result_mask, STAR.result.full = STAR.result_full,
              STAR.full.object=STAR.obj1,AdaPT.full.object=AdaPT.obj1,
              STAR.mask.object=STAR.obj2,AdaPT.mask.object=AdaPT.obj2))
}


rho <- 0
type <- 0
repeat.times <- 500
seed <- 1
set.seed(seed)
n <- 25*25
x1 <- seq(-100, 100, length.out =25)
x2 <- seq(-100, 100, length.out =25)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
cov.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 0.3, 0.01)
num.steps.update.score <- 10
num.steps.gam <- 5
score.params <- list(cov.formula = cov.formula,
                     num.steps = num.steps.gam)


## Case 1: a circle in the center
H0 <- apply(x, 1, function(coord){sum(coord^2) < 1600})


mu <- ifelse(H0, 10, 0)
experiment_masked(x,mu,0.1,alpha.list,num.steps.update.score,scope.params,base=5,type="poisson")

iter_run(x,mu,0.1,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.1',sep=""),base=5,type="poisson")
iter_run(x,mu,0.2,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.2',sep=""),base=5,type="poisson")
iter_run(x,mu,0.3,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.3',sep=""),base=5,type="poisson")
iter_run(x,mu,0.4,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.4',sep=""),base=5,type="poisson")
iter_run(x,mu,0.5,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.5',sep=""),base=5,type="poisson")
iter_run(x,mu,0.6,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.6',sep=""),base=5,type="poisson")
iter_run(x,mu,0.7,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.7',sep=""),base=5,type="poisson")
iter_run(x,mu,0.8,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.8',sep=""),base=5,type="poisson")
iter_run(x,mu,0.9,alpha.list,num.steps.update.score,scope.params,paste(wd,'hyptest_poisson_10_0.9',sep=""),base=5,type="poisson")


mu <- ifelse(H0, 1, 0)
iter_run(x,mu,0.1,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.1')
iter_run(x,mu,0.2,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.2')
iter_run(x,mu,0.3,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.3')
iter_run(x,mu,0.4,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.4')
iter_run(x,mu,0.5,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.5')
iter_run(x,mu,0.6,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.6')
iter_run(x,mu,0.7,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.7')
iter_run(x,mu,0.8,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.8')
iter_run(x,mu,0.9,alpha.list,num.steps.update.score,scope.params,'experimentsref_1_0.9')

mu <- ifelse(H0, 2, 0)
iter_run(x,mu,0.1,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.1')
iter_run(x,mu,0.2,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.2')
iter_run(x,mu,0.3,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.3')
iter_run(x,mu,0.4,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.4')
iter_run(x,mu,0.5,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.5')
iter_run(x,mu,0.6,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.6')
iter_run(x,mu,0.7,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.7')
iter_run(x,mu,0.8,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.8')
iter_run(x,mu,0.9,alpha.list,num.steps.update.score,scope.params,'experimentsref_2_0.9')

mu <- ifelse(H0, 3, 0)
iter_run(x,mu,0.2,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.1')
iter_run(x,mu,0.4,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.2')
iter_run(x,mu,0.6,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.3')
iter_run(x,mu,0.8,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.4')
iter_run(x,mu,0.2,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.5')
iter_run(x,mu,0.4,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.6')
iter_run(x,mu,0.6,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.7')
iter_run(x,mu,0.8,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.8')
iter_run(x,mu,0.9,alpha.list,num.steps.update.score,scope.params,'experimentsref_3_0.9')



rho <- 0
type <- 0
repeat.times <- 500
seed <- 1
set.seed(seed)
n <- 50*50
x1 <- seq(-100, 100, length.out =50)
x2 <- seq(-100, 100, length.out =50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
cov.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 0.3, 0.01)
num.steps.update.score <- 10
num.steps.gam <- 5
score.params <- list(cov.formula = cov.formula,
                     num.steps = num.steps.gam)


## Case 1: a circle in the center
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)
tau <- 0.1

t = experiment_masked(x,mu,tau,alpha.list,num.steps.update.score,scope.params) 



pdf(file="convex_example_run2.pdf")
par(mfrow=c(2,3))
khat <- max(c(0,which(sort(t$pvals)<=0.1*(1:n)/n)))
alpha <- 0.1 * khat / n
plot.rej.convex(x, ifelse(t$pvals < alpha,1,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "BH (Full)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$pvals < t$AdaPT.full.object$s,1,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "AdaPT (Full)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$STAR.full.object$mask,2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "STAR (Full)", xaxt = "n", yaxt = "n")
khat <- max(c(0,which(sort(t$pvals_mask)<=0.1*(1:n)/n)))
alpha <- 0.1 * khat / n
plot.rej.convex(x, ifelse(t$pvals_mask < alpha,2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "BH (Fission)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$pvals_mask < t$AdaPT.mask.object$s,2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "AdaPT (Fission)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$STAR.mask.object$mask,100,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "STAR (Fission)", xaxt = "n", yaxt = "n")
dev.off()

mask <- ifelse(t$pvals < t$AdaPT.full.object$s,1,0)
color <- ifelse(mask, col.fg, col.bg)

plot(x[, 1], x[, 2], type = "n", xlab = "", ylab = "")
points(x[, 1], x[, 2], col = mask, cex =0.5)

khat <- max(c(0,which(sort(t$pvals)<=0.1*(1:n)/n)))
alpha <- 0.1 * khat / n
t$pvals < alpha


