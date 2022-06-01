library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(gridExtra)
library(patchwork)
library(latex2exp)

agg_stats <- function(df,alpha.FCR) {
  q <- nrow(df$x)
  pvals <- df$pvals
  pvals_mask <- df$pvals_mask
  tau <- df$tau
  type <- df$type
  
  BH_full_CI <- sapply(alpha.list,function(alpha) {
    n <- length(pvals)
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    x <- which(pvals < alpha)
    alpha.2 <- alpha.FCR*length(x)/q
    nrej <- length(x)
    if(type == "poisson"){
      m <- sum(df$Y[x])
      true_mean <- sum(df$mu[x])
      low = qchisq(alpha.2/2, 2*m)/2
      up = qchisq(1-alpha.2/2, 2*(m+1))/2
      
      cov <- true_mean <= up & (true_mean >= low)
      if(nrej > 0 ){
        CI_length = (up - low)/nrej
      }
      else{
        CI_length = NA
      }
      ret <- c(true_mean,m,sum(cov),length(cov),CI_length)
    }
    else{
      sd = 1
      up <- df$Y[x] + qnorm(1-alpha.2/2,mean=0,sd=sd)
      low <- df$Y[x] - qnorm(1-alpha.2/2,mean=0,sd=sd)
      cov <- df$mu[x]  <= up & (df$mu[x] >= low)
      ret<- c(mean(df$mu[x]),mean(df$Y[x]),sum(cov),length(cov),mean(up-low))
    }
    return(ret)
  })
  
  ada_full_CI <- apply(df$AdaPT.full.object$s,2,
                       function(s) {
                         tmp <- (df$pvals <= s)  
                         x <- (df$pvals_mask <= s)  
                         nrej <- sum(x)
                         alpha.2 <- alpha.FCR*sum(tmp)/q
                         if(type == "poisson"){
                           m <- sum(df$Y[x])
                           true_mean <- sum(df$mu[x])
                           low = qchisq(alpha.2/2, 2*m)/2
                           up = qchisq(1-alpha.2/2, 2*(m+1))/2
                           
                           cov <- true_mean <= up & (true_mean >= low)
                           if(nrej > 0 ){
                             CI_length = (up - low)/nrej
                           }
                           else{
                             CI_length = NA
                           }
                           ret <- c(true_mean,m,sum(cov),length(cov),CI_length)
                         }
                         else{
                           sd = 1
                           up <- df$Y[tmp] + qnorm(1-alpha.2/2,mean=0,sd=sd)
                           low <- df$Y[tmp] - qnorm(1-alpha.2/2,mean=0,sd=sd)
                           cov <- df$mu[tmp] <= up & (df$mu[tmp] >= low)
                           
                           ret<- c(mean(df$mu[tmp]),mean(df$Y[tmp]),sum(cov),length(cov),mean(up-low))
                         }
                         return(ret)
                       })
  
  STAR_full_CI <- apply(df$STAR.full.object$mask,2,
                        function(x) {
                          nrej <- sum(x)
                          alpha.2 <- alpha.FCR* sum(x)/q
                          if(type == "poisson"){
                            m <- sum(df$Y[x])
                            true_mean <- sum(df$mu[x])
                            low = qchisq(alpha.2/2, 2*m)/2
                            up = qchisq(1-alpha.2/2, 2*(m+1))/2
                            
                            cov <- true_mean <= up & (true_mean >= low)
                            if(nrej > 0 ){
                              CI_length = (up - low)/nrej
                            }
                            else{
                              CI_length = NA
                            }
                            ret <- c(true_mean,m,sum(cov),length(cov),CI_length)
                          }
                          else{
                            sd = 1
                            up <- df$Y[x] + qnorm(1-alpha.2/2,mean=0,sd=sd)
                            low <- df$Y[x] - qnorm(1-alpha.2/2,mean=0,sd=sd)
                            cov <- df$mu[x] <= up & (df$mu[x] >= low)
                            ret<- c(mean(df$mu[x]),mean(df$Y[x]),sum(cov),length(cov),mean(up-low))
                          }

                          return(ret)
                        })
  
  BH_mask_CI <- sapply(alpha.list,function(alpha) {
    n <- length(pvals_mask)
    khat <- max(c(0,which(sort(pvals_mask)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    x <- which(pvals_mask < alpha)
    nrej <- length(x)
    if(type == "poisson"){
      m <- sum(df$g_Y[x])
      true_mean <- sum(df$mu[x])*(1-tau)
      low = qchisq(alpha.FCR/2, 2*m)/2
      up = qchisq(1-alpha.FCR/2, 2*(m+1))/2
      
      cov <- true_mean <= up & (true_mean >= low)
      if(nrej > 0 ){
        CI_length = (up - low)/nrej
      }
      else{
        CI_length = NA
      }
      ret <- c(true_mean,m,sum(cov),length(cov),CI_length)
    }
    else{
      m <- mean(df$g_Y[x])
      true_mean <- mean(df$mu[x])
      sd = sqrt((1+tau^(-2))/nrej)
      up <- m + qnorm(1-alpha.FCR/2,mean=0,sd=sd)
      low <- m - qnorm(1-alpha.FCR/2,mean=0,sd=sd)
      cov <- true_mean <= up & (true_mean >= low)
      ret <- c(true_mean,m,sum(cov),length(cov),mean(up-low))
      
    }

    return(ret)
  })
  
  
  ada_mask_CI <- apply(df$AdaPT.mask.object$s,2,function(s) {
    x <- (df$pvals_mask <= s)  
    nrej <- sum(x)
    if(type == "poisson"){
      m <- sum(df$g_Y[x])
      true_mean <- sum(df$mu[x])*(1-tau)
      low = qchisq(alpha.FCR/2, 2*m)/2
      up = qchisq(1-alpha.FCR/2, 2*(m+1))/2
      
      cov <- true_mean <= up & (true_mean >= low)
      if(nrej > 0 ){
        CI_length = (up - low)/nrej
      }
      else{
        CI_length = NA
      }
      ret <- c(true_mean,m,sum(cov),length(cov),CI_length)
    }
    else{
      m <- mean(df$g_Y[x])
      true_mean <- mean(df$mu[x])
      sd = sqrt((1+tau^(-2))/nrej)
      up <- m + qnorm(1-alpha.FCR/2,mean=0,sd=sd)
      low <- m - qnorm(1-alpha.FCR/2,mean=0,sd=sd)
      cov <- true_mean <= up & (true_mean >= low)
      ret <- c(true_mean,m,sum(cov),length(cov),up-low)
    }
    
    return(ret)
  })
  STAR_mask_CI <- apply(df$STAR.mask.object$mask,2,function(x) {
    nrej <- sum(x)
    if(type == "poisson"){
      m <- sum(df$g_Y[x])
      true_mean <- sum(df$mu[x])*(1-tau)
      low = qchisq(alpha.FCR/2, 2*m)/2
      up = qchisq(1-alpha.FCR/2, 2*(m+1))/2
      
      cov <- true_mean <= up & (true_mean >= low)
      if(nrej > 0 ){
        CI_length = (up - low)/nrej
      }
      else{
        CI_length = NA
      }
      ret <- c(true_mean,m,sum(cov),length(cov),CI_length)
    }
    else{
      m <- mean(df$g_Y[x])
      true_mean <- mean(df$mu[x])
      sd = sqrt((1+tau^(-2))/nrej)
      up <- m + qnorm(1-alpha.FCR/2,mean=0,sd=sd)
      low <- m - qnorm(1-alpha.FCR/2,mean=0,sd=sd)
      cov <- true_mean <= up & (true_mean >= low)
      ret <- c(true_mean,m,sum(cov),length(cov),up-low)
    }
    
    return(ret)
  })
  return(list(BH_full_CI=BH_full_CI,ada_full_CI=ada_full_CI,STAR_full_CI=STAR_full_CI,BH_mask_CI=BH_mask_CI
              ,ada_mask_CI=ada_mask_CI,STAR_mask_CI=STAR_mask_CI))
}





get_summary <-function(filename,alpha.list,mu,tau) {
  load(filename)
  BH.mask.FDP<- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$BH.result.mask$FDP)))
  AdaPT.mask.FDP <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$adapt.result.mask$FDP)))
  STAR.mask.FDP <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$STAR.result.mask$FDP)))
  BH.full.FDP<- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$BH.result.full$FDP)))
  AdaPT.full.FDP <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$adapt.result.full$FDP)))
  STAR.full.FDP <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$STAR.result.full$FDP)))
  df_FDP <- cbind(BH.mask.FDP,AdaPT.mask.FDP,STAR.mask.FDP,BH.full.FDP,AdaPT.full.FDP,STAR.full.FDP)
  rownames(df_FDP) <- alpha.list
  colnames(df_FDP) <- c('BH Blurred','AdaPT Blurred','STAR Blurred','BH Full','AdaPT Full','STAR Full')
  
  df_FDP <- melt(df_FDP)
  df_FDP<- cbind(df_FDP,do.call(rbind,strsplit(as.character(df_FDP$Var2),split=" ")))
  
  
  BH.mask.power<- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$BH.result.mask$power)))
  AdaPT.mask.power <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$adapt.result.mask$power)))
  STAR.mask.power <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$STAR.result.mask$power)))
  BH.full.power<- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$BH.result.full$power)))
  AdaPT.full.power <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$adapt.result.full$power)))
  STAR.full.power <- colMeans(do.call(rbind.data.frame, lapply(results,function(x) x$STAR.result.full$power)))
  df_power <- cbind(BH.mask.power,AdaPT.mask.power,STAR.mask.power,BH.full.power,AdaPT.full.power,STAR.full.power)
  rownames(df_power) <- alpha.list
  colnames(df_power) <- c('BH Blurred','AdaPT Blurred','STAR Blurred','BH Full','AdaPT Full','STAR Full')
  
  
  df_power <- melt(df_power)
  df_power<- cbind(df_power,do.call(rbind,strsplit(as.character(df_power$Var2),split=" ")))
  colnames(df_power) <- c("TargetLevel","ID","Power","Procedure","Splitting")
  
  
  agg_results <- list()
  counter = 0 
  for(i in 1:length(results)){
    if(length(results[[i]]) != 0){
      counter = counter +1
      agg_results[[counter]] <- agg_stats(results[[i]],0.2)
    }
  }
  
  BH.mask.CI <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) x$BH_mask_CI[5,])),na.rm=TRUE)
  AdaPT.mask.CI <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) x$ada_mask_CI[5,])),na.rm=TRUE)
  STAR.mask.CI <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) x$STAR_mask_CI[5,])),na.rm=TRUE)
  BH.full.CI <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) x$BH_full_CI[5,])),na.rm=TRUE)
  AdaPT.full.CI <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) x$ada_full_CI[5,])),na.rm=TRUE)
  STAR.full.CI <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) x$STAR_full_CI[5,])),na.rm=TRUE)
  df_CI_length <- cbind(BH.mask.CI,AdaPT.mask.CI,STAR.mask.CI,BH.full.CI,AdaPT.full.CI,STAR.full.CI)
  rownames(df_CI_length) <- alpha.list
  colnames(df_CI_length) <- c('BH Blurred','AdaPT Blurred','STAR Blurred','BH Full','AdaPT Full','STAR Full')
  
  
  df_CI_length <- melt(df_CI_length)
  df_CI_length<- cbind(df_CI_length,do.call(rbind,strsplit(as.character(df_CI_length$Var2),split=" ")))
  
  BH.mask.MSE <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) (x$BH_mask_CI[1,] -x$BH_mask_CI[2,])**2 )),na.rm=TRUE)
  AdaPT.mask.MSE <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) (x$ada_mask_CI[1,] -x$ada_mask_CI[2,])**2 )),na.rm=TRUE)
  STAR.mask.MSE <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) (x$STAR_mask_CI[1,] -x$STAR_mask_CI[2,])**2 )),na.rm=TRUE)
  BH.full.MSE <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) (x$BH_full_CI[1,] -x$BH_full_CI[2,])**2 )),na.rm=TRUE)
  AdaPT.full.MSE <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) (x$ada_full_CI[1,] -x$ada_full_CI[2,])**2 )),na.rm=TRUE)
  STAR.full.MSE <- colMeans(do.call(rbind.data.frame, lapply(agg_results,function(x) (x$STAR_full_CI[1,] -x$STAR_full_CI[2,])**2 )),na.rm=TRUE)
  df_MSE <- cbind(BH.mask.MSE,AdaPT.mask.MSE,STAR.mask.MSE,BH.full.MSE,AdaPT.full.MSE,STAR.full.MSE)
  rownames(df_MSE) <- alpha.list
  colnames(df_MSE) <- c('BH Blurred','AdaPT Blurred','STAR Blurred','BH Full','AdaPT Full','STAR Full')
  
  df_MSE <- melt(df_MSE)
  df_MSE<- cbind(df_MSE,do.call(rbind,strsplit(as.character(df_MSE$Var2),split=" ")))
  
  
  BH.mask.num_in <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) x$BH_mask_CI[3,])),na.rm=TRUE)
  AdaPT.mask.num_in <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) x$ada_mask_CI[3,])),na.rm=TRUE)
  STAR.mask.num_in <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) x$STAR_mask_CI[3,])),na.rm=TRUE)
  BH.full.num_in <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) x$BH_full_CI[3,])),na.rm=TRUE)
  AdaPT.full.num_in <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) x$ada_full_CI[3,])),na.rm=TRUE)
  STAR.full.num_in <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) x$STAR_full_CI[3,])),na.rm=TRUE)
  df_num_in <- cbind(BH.mask.num_in,AdaPT.mask.num_in,STAR.mask.num_in,BH.full.num_in,AdaPT.full.num_in,STAR.full.num_in)
  rownames(df_num_in) <- alpha.list
  colnames(df_num_in) <- c('BH Blurred','AdaPT Blurred','STAR Blurred','BH Full','AdaPT Full','STAR Full')
  
  
  BH.mask.num <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) !is.na(x$BH_mask_CI[3,])*x$BH_mask_CI[4,])),na.rm=TRUE)
  AdaPT.mask.num <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) !is.na(x$ada_mask_CI[3,])*x$ada_mask_CI[4,])),na.rm=TRUE)
  STAR.mask.num <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) !is.na(x$STAR_mask_CI[3,])*x$STAR_mask_CI[4,])),na.rm=TRUE)
  BH.full.num <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) (!is.na(x$BH_full_CI[3,])*1)*x$BH_full_CI[4,])),na.rm=TRUE)
  AdaPT.full.num <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) (!is.na(x$ada_full_CI[3,])*1)*x$ada_full_CI[4,])),na.rm=TRUE)
  STAR.full.num <- colSums(do.call(rbind.data.frame, lapply(agg_results,function(x) (!is.na(x$STAR_full_CI[3,])*1)*x$STAR_full_CI[4,])),na.rm=TRUE)
  df_num <-cbind(BH.mask.num,AdaPT.mask.num,STAR.mask.num,BH.full.num,AdaPT.full.num,STAR.full.num)
  rownames(df_num) <- alpha.list
  colnames(df_num) <- c('BH Blurred','AdaPT Blurred','STAR Blurred','BH Full','AdaPT Full','STAR Full')
  
  
  df_FCR <- 1- df_num_in / df_num
  rownames(df_FCR) <- alpha.list
  colnames(df_FCR) <- c('BH Blurred','AdaPT Blurred','STAR Blurred','BH Full','AdaPT Full','STAR Full')
  
  df_FCR <- melt(df_FCR)
  df_FCR<- cbind(df_FCR,do.call(rbind,strsplit(as.character(df_FCR$Var2),split=" ")))
  df_FCR <- cbind("FCR",df_FCR)
  df_power <- cbind("Power",df_power)
  df_CI_length <- cbind("CI Length",df_CI_length)
  df_MSE <- cbind("MSE",df_MSE)
  df_FDP <- cbind("FDP",df_FDP)
  
  colnames(df_FDP) <- c("Metric","TargetLevel","ID","Value","Procedure","Splitting")
  colnames(df_power) <- c("Metric","TargetLevel","ID","Value","Procedure","Splitting")
  colnames(df_CI_length) <- c("Metric","TargetLevel","ID","Value","Procedure","Splitting")
  colnames(df_MSE) <- c("Metric","TargetLevel","ID","Value","Procedure","Splitting")
  colnames(df_FCR) <- c("Metric","TargetLevel","ID","Value","Procedure","Splitting")
  res<- cbind(mu,tau,rbind(df_FDP,df_power,df_CI_length,df_MSE,df_FCR))
  
  return(res)
  
}

load("results_PSC/interactive_poisson_0.5_100_50_50_500")
results[[495]]
        
        
alpha.list <- seq(0.01, 0.3, 0.01)
t1 <- get_summary("results_PSC/interactive_poisson_0.9_10_5_50_500",alpha.list,10,0.9)
t2 <- get_summary("results_PSC/interactive_poisson_0.8_10_5_50_500",alpha.list,10,0.8)
t3 <- get_summary("results_PSC/interactive_poisson_0.7_10_5_50_500",alpha.list,10,0.7)
t4 <- get_summary("results_PSC/interactive_poisson_0.6_10_5_50_500",alpha.list,10,0.6)
t5 <- get_summary("results_PSC/interactive_poisson_0.5_10_5_50_500",alpha.list,10,0.5)
t6 <- get_summary("results_PSC/interactive_poisson_0.4_10_5_50_500",alpha.list,10,0.4)
t7 <- get_summary("results_PSC/interactive_poisson_0.3_10_5_50_500",alpha.list,10,0.3)
t8 <- get_summary("results_PSC/interactive_poisson_0.2_10_5_50_500",alpha.list,10,0.2)
t9 <- get_summary("results_PSC/interactive_poisson_0.1_10_5_50_500",alpha.list,10,0.1)

alpha.list <- seq(0.01, 0.3, 0.01)
t1 <- get_summary("results_PSC/interactive_poisson_0.9_100_50_50_500",alpha.list,10,0.9)
t2 <- get_summary("results_PSC/interactive_poisson_0.8_100_50_50_500",alpha.list,10,0.8)
t3 <- get_summary("results_PSC/interactive_poisson_0.7_100_50_50_500",alpha.list,10,0.7)
t4 <- get_summary("results_PSC/interactive_poisson_0.6_100_50_50_500",alpha.list,10,0.6)
t5 <- get_summary("results_PSC/interactive_poisson_0.5_100_50_50_500",alpha.list,10,0.5)
t6 <- get_summary("results_PSC/interactive_poisson_0.4_100_50_50_500",alpha.list,10,0.4)
t7 <- get_summary("results_PSC/interactive_poisson_0.3_100_50_50_500",alpha.list,10,0.3)
t8 <- get_summary("results_PSC/interactive_poisson_0.2_100_50_50_500",alpha.list,10,0.2)
t9 <- get_summary("results_PSC/interactive_poisson_0.1_100_50_50_500",alpha.list,10,0.1)
agg <- rbind(t1,t2,t4,t5,t6,t7,t8,t9)


ds <- agg[agg$Splitting=="Blurred",]
for(i in 1:9/10) {
  t_temp <- t1[t1$Splitting=="Full",]
  t_temp$tau = i
  ds = rbind(ds,t_temp)
}

p1 <- ggplot(ds[ds$Metric == "FDP" & ds$TargetLevel==0.2,],
             aes(x = tau, y = Value, fill = Procedure,linetype=Splitting)) +
  geom_line(aes(color = Procedure,linetype=Splitting), size = 0.8) +
  geom_point(aes(shape = Procedure, color = Procedure), size = 1.5) +
  geom_hline(yintercept=0.2,linetype="dashed",size=1) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  scale_linetype_manual(NULL,values = c('dashed','solid')) + 
  xlab(TeX("$\\tau$"))+
  ylab("False discovery rate") +
  scale_y_continuous(breaks = seq(0, 0.4, len = 5),limits=c(0,0.4))


p2 <- ggplot(ds[ds$Metric == "Power" & ds$TargetLevel==0.2,],
             aes(x = tau, y = Value, fill = Procedure,linetype=Splitting)) +
  geom_line(aes(color = Procedure,linetype=Splitting), size = 0.8) +
  geom_point(aes(shape = Procedure, color = Procedure), size = 1.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  scale_linetype_manual(NULL,values = c('dashed','solid')) + 
  xlab(TeX("$\\tau$"))+
  ylab("Power")


p3 <- ggplot(ds[ds$Metric == "FCR" & ds$TargetLevel==0.2,],
             aes(x = tau, y = Value, fill = Procedure,linetype=Splitting)) +
  geom_line(aes(color = Procedure,linetype=Splitting), size = 0.8) +
  geom_point(aes(shape = Procedure, color = Procedure), size = 1.5) +
  geom_hline(yintercept=0.2,linetype="dashed",size=1) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15), legend.text = element_text(size = 15)) +
  xlab(TeX("$\\tau$")) +
  ylab("False coverage rate / miscoverage rate") +
  scale_linetype_manual(NULL,values = c('dashed','solid'))

p4 <- ggplot(ds[ds$Metric == "CI Length" & ds$TargetLevel==0.2,],
             aes(x = tau, y = Value, fill = Procedure,linetype=Splitting)) +
  geom_line(aes(color = Procedure,linetype=Splitting), size = 0.8) +
  geom_point(aes(shape = Procedure, color = Procedure), size = 1.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab(TeX("$\\tau$")) +
  scale_linetype_manual(NULL,values = c('dashed','solid'))+
  ylab("CI Length")


grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)


alpha.list <- seq(0.01, 0.3, 0.01)
t1 <- get_summary("results_PSC/interactive_normal_0.1_2_0_50_500",alpha.list,2,0.1)
t2 <- get_summary("results_PSC/interactive_normal_0.2_2_0_50_500",alpha.list,2,0.2)
t3 <- get_summary("results_PSC/interactive_normal_0.3_2_0_25_500",alpha.list,2,0.3)
t4 <- get_summary("results_PSC/interactive_normal_0.4_2_0_50_500",alpha.list,2,0.4)
t5 <- get_summary("results_PSC/interactive_normal_0.5_2_0_50_500",alpha.list,2,0.5)
t6 <- get_summary("results_PSC/interactive_normal_0.6_2_0_50_500",alpha.list,2,0.6)
t7 <- get_summary("results_PSC/interactive_normal_0.7_2_0_50_500",alpha.list,2,0.7)
t8 <- get_summary("results_PSC/interactive_normal_0.8_2_0_50_500",alpha.list,2,0.8)
t9 <- get_summary("results_PSC/interactive_normal_0.9_2_0_50_500",alpha.list,2,0.9)



agg <- rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9)


ds <- agg[agg$Splitting=="Blurred",]
for(i in 1:9/10) {
  t_temp <- t1[t1$Splitting=="Full",]
  t_temp$tau = i
  ds = rbind(ds,t_temp)
}




filename = "results/interactive_normal_0.1_2_0_50_5"
load(filename)
t <- results[[1]]
x <- t$x
source('STAR_code/STAR_convex.R')

pdf(file="figures/convex_example_true.pdf")
plot.rej.convex(x, ifelse(t$mu == 2,1,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "", xaxt = "n", yaxt = "n")
dev.off()

pdf(file="figures/convex_example_run.pdf")
par(mfrow=c(2,3))
khat <- max(c(0,which(sort(t$pvals)<=0.2*(1:n)/n)))
alpha <- 0.1 * khat / n
plot.rej.convex(x, ifelse(t$pvals < alpha,1,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "BH (Full)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$pvals < t$AdaPT.full.object$s[,20],2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "AdaPT (Full)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$STAR.full.object$mask[,20],2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "STAR (Full)", xaxt = "n", yaxt = "n")
khat <- max(c(0,which(sort(t$pvals_mask)<=0.1*(1:n)/n)))
alpha <- 0.1 * khat / n
plot.rej.convex(x, ifelse(t$pvals_mask < alpha,2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "BH (Blurred)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$pvals_mask < t$AdaPT.mask.object$s[,20],2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "AdaPT (Blurred)", xaxt = "n", yaxt = "n")
plot.rej.convex(x, ifelse(t$STAR.mask.object$mask[,20],2,0), cex = 0.5,col.bg = "#A9A9A9", col.fg = "#000000",main = "STAR (Blurred)", xaxt = "n", yaxt = "n")
dev.off()



