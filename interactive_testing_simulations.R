library("MASS")
library("mgcv")
library("RCurl")
library("parallel")

## AdaPT
source("STAR_code/AdaPT.R")
source("STAR_code/AdaPT_gam.R")
source('STAR_code/summarize_methods.R')
source("STAR_code/expr_func.R")
source('STAR_code/STAR_convex.R')
source('Interactive Testing Scripts.R')

wd = paste(getwd(),"/results/",sep="")

args = commandArgs(trailingOnly=TRUE)
type = args[1]
tau  = as.numeric(args[2])
alt = as.numeric(args[3])
null = as.numeric(args[4])
grid_size = as.numeric(args[5])
repeat.times = as.numeric(args[6])
seed = as.numeric(args[7])*1000

rho <- 0
n <- grid_size*grid_size
x1 <- seq(-100, 100, length.out =grid_size)
x2 <- seq(-100, 100, length.out =grid_size)
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
mu <- ifelse(H0, alt, null)

wrapper_func <- function(x,mu,tau,alpha.list,num.steps.update.score,scope.params,alt,null,type=type,seed=seed) {
  res <- tryCatch(experiment_masked(x,mu,tau,alpha.list,num.steps.update.score,scope.params,alt,null,type=type,seed=seed) ,
           error = function(e){
             print(e)
             return(list())
             })
}

results <- mclapply(1:repeat.times, function(r) {return(wrapper_func(x,mu,tau,alpha.list,num.steps.update.score,scope.params,alt,null,type=type,seed=(seed+r)))}, mc.cores = detectCores())
filename = paste(wd,paste("interactive",type,tau,alt,null,grid_size,repeat.times,sep="_"),sep="")
save(file = filename, results)




