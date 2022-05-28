library(glmnet)
library(clubSandwich)
library(sandwich)
library(MASS)
library(rootSolve)
library(ggforce)
library(genlasso)
library(splines)
library(cplm)
library(quantreg)
library(latex2exp)
library(parallel)
library(ggplot2)
library(gridExtra)

unif_CI = function(n, knots, Y, X, alpha,deg=1){
  k = length(knots) + deg
  if (length(knots) == 0) {
    bs_X = cbind(rep(1, n), 1:n)
  } else if (k + 1 < n) {
    basis = tp(1:n, knots = knots, degree=deg, k = k)
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
  # cvu <- critval(kappa, alpha = alpha, rdf = n - k - 1)
  # cvu <- critval(kappa, alpha = alpha)
  v = n - k - 1
  cvu <- uniroot(function(x) {kappa*(1 + x^2/v)^(-v/2)/pi + 2*(1 - pt(x, df = v)) - alpha},
                 c(0, max(2*kappa/alpha/pi, 10)))$root
  
  true_trend = lm(X ~ bs_X)$fitted.values
  temp = predict(lm(Y ~ bs_X), se.fit = TRUE, level = 1 - alpha)
  # mean_se[["masking"]] = mean(temp$se.fit)
  width = cvu*temp$se.fit
  # width = cvu*temp$se.fit/temp$residual.scale*sqrt(sum((Y - true_trend)^2)/(n - k - 1))
  return(list(predicted = cbind(temp$fit, temp$fit - width, temp$fit + width, true_trend),
              c = cvu, mean_se = mean(temp$se.fit)))
}


get_trendfilter <- function(Y,pos,degree,alpha){
  n <- length(Y)
  sigma_hat <- sd(Y[-n] - Y[-1])
  
  noise = rnorm(n, sd = sigma_hat)
  g_Y = Y + noise
  h_Y = Y - noise
  
  zz <- file("knots.csv", open = "wt")
  sink(zz)
  select_model = trendfilter(g_Y,pos=pos,ord=degree,verbose=TRUE)
  sink()
  closeAllConnections()
  knots <- read.csv('knots.csv',header=FALSE)
  lambda <- as.numeric(sapply(strsplit(knots[,1],"="), function(x) x[2]))
  additions <- grepl("adding",knots[,2])
  deletions <- grepl("deleting",knots[,2])
  coordinates <- as.numeric(sapply(strsplit(knots[,2],"coordinate "), function(x) x[2]))
  knot_list <- cbind(lambda,additions,deletions,coordinates)
  
  
  cv = cv.trendfilter(select_model)
  knot_subset = knot_list[knot_list[,1] > round(cv$lambda.min,3),]
  coordinates <- aggregate(knot_subset[,2:3],by=list(additions=knot_subset[,4]),FUN=sum)
  knots <- coordinates[which(coordinates[,2]>coordinates[,3]),1]
  
  
  bs_X = bs(1:n, knots = knots, degree = degree); n_x = ncol(bs_X) 
  temp = predict(lm(h_Y ~ bs_X), interval="confidence", se.fit = TRUE, level = 1 - alpha)
  temp2 = predict(lm(h_Y ~ bs_X), interval="prediction", se.fit = TRUE, level = 1 - alpha)
  U_CI <- unif_CI(n = n, knots = knots, Y = h_Y, X = pos, alpha = alpha,deg=degree)$predicted[,2:3]
  
  df <- data.frame(cbind(pos,Y,temp$fit,U_CI,temp2$fit[,2:3]))
  colnames(df) <- c("Position","Actual","Predicted","CI_Low","CI_High","UCI_Low","UCI_High","PI_Low","PI_High")
  return(df)
  
}

quasar <- readRDS("D:/GitHub/masking_regression/trendfiltering-astro-data/fig3top_quasar_spectrum.rds")
galaxy<- readRDS("D:/GitHub/masking_regression/trendfiltering-astro-data/fig3middle_galaxy_spectrum.rds")
stellar <- readRDS("D:/GitHub/masking_regression/trendfiltering-astro-data/fig3bottom_stellar_spectrum.rds")
quasar <- quasar[quasar$wavelength <= 5500 ,]
quasar <- quasar[quasar$wavelength >= 4000 ,]
galaxy <- galaxy[galaxy$wavelength <= 5500 ,]
galaxy <- galaxy[galaxy$wavelength >= 4000 ,]
stellar <- stellar[stellar$wavelength <= 5500 ,]
stellar <- stellar[stellar$wavelength >= 4000 ,]

res_quasar <- get_trendfilter(quasar$flux,quasar$wavelength,1,0.2)
res_quasar_2 <- get_trendfilter(quasar$flux,quasar$wavelength,2,0.2)
res_quasar_3 <- get_trendfilter(quasar$flux,quasar$wavelength,3,0.2)
res_galaxy <- get_trendfilter(galaxy$flux,galaxy$wavelength,1,0.2)
res_galaxy_2 <- get_trendfilter(galaxy$flux,galaxy$wavelength,2,0.2)
res_galaxy_3 <- get_trendfilter(galaxy$flux,galaxy$wavelength,3,0.2)
res_stellar <- get_trendfilter(stellar$flux,stellar$wavelength,1,0.2)
res_stellar_2 <- get_trendfilter(stellar$flux,stellar$wavelength,2,0.2)
res_stellar_3 <- get_trendfilter(stellar$flux,stellar$wavelength,3,0.2)



quasar_linear = ggplot(res_quasar,aes(Position,Actual,linetype="Actual")) +geom_line(size=0.7) +
  geom_line(aes(Position,Predicted,linetype="Fitted"),size=0.7)+
  geom_line(aes(Position,UCI_Low),color="orange",linetype="longdash",size=0.7) +
  geom_line(aes(Position,CI_Low),color="blue",linetype="twodash",size=0.7) + 
  geom_line(aes(Position,CI_High,color="Pointwise CI"),linetype="twodash",size=0.7)+
  geom_line(aes(Position,UCI_High,color="Uniform CI"),linetype="longdash",size=0.7) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "grey"),
    plot.title = element_text(hjust = 0.5,size=20),
    text = element_text(size = 15),
    legend.position="none"
  ) +
  xlab(TeX("$\\lambda(???) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/???) $"))+
  ggtitle(parse(text= paste("Quasar ( ~", expression(z %~~% 2.4),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))

quasar_linear_detailed = ggplot(res_quasar,aes(Position,Actual,linetype="Actual")) +geom_line() +
  geom_line(aes(Position,Predicted,linetype="Fitted"),size=1)+
  geom_line(aes(Position,UCI_Low),color="orange",linetype="longdash",size=1) +
  geom_line(aes(Position,CI_Low),color="blue",linetype="twodash",size=1) + 
  geom_line(aes(Position,CI_High,color="Pointwise CI"),linetype="twodash",size=1)+
  geom_line(aes(Position,UCI_High,color="Uniform CI"),linetype="longdash",size=1) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "grey"),
    plot.title = element_text(hjust = 0.5,size=20),
    text = element_text(size = 15),
    legend.position="none"
  ) +
  xlab(TeX("$\\lambda(???) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/??? ) $"))+
  ggtitle(parse(text= paste("Quasar ( ~", expression(z %~~% 2.4),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))+
  facet_zoom(xlim = c(4700, 4800),ylim=c(1,15),zoom.size=0.5,show.area=TRUE,horizontal=TRUE,split=FALSE)

ggsave(filename = paste(getwd(),"/figures/","quasar_linear",".png",sep=""),
       plot = quasar_linear, width = 20, height = 20)

ggsave(filename = paste(getwd(),"/figures/","quasar_linear_detailed",".png",sep=""),
       plot = quasar_linear_detailed, width = 40, height = 20)





galaxy_quadratic = ggplot(res_galaxy_2,aes(Position,Actual,linetype="Actual")) +geom_line() +
  geom_line(aes(Position,Predicted,linetype="Fitted"),size=1.2)+
  geom_line(aes(Position,UCI_Low),color="orange",linetype="longdash",size=1.2) +
  geom_line(aes(Position,CI_Low),color="blue",linetype="twodash",size=1.2) + 
  geom_line(aes(Position,CI_High,color="Pointwise CI"),linetype="twodash",size=1.2)+
  geom_line(aes(Position,UCI_High,color="Uniform CI"),linetype="longdash",size=1.2) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "grey"),
    plot.title = element_text(hjust = 0.5,size=20),
    text = element_text(size = 15),
    legend.position="none"
  ) +
  xlab(TeX("$\\lambda(???) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/??? ) $"))+
  ggtitle(parse(text= paste("Galaxy ( ~", expression(z %~~% 0.14),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))

galaxy_quadratic_detailed = ggplot(res_galaxy_2,aes(Position,Actual,linetype="Actual")) +geom_line() +
  geom_line(aes(Position,Predicted,linetype="Fitted"),size=1.2)+
  geom_line(aes(Position,UCI_Low),color="orange",linetype="longdash",size=1.2) +
  geom_line(aes(Position,CI_Low),color="blue",linetype="twodash",size=1.2) + 
  geom_line(aes(Position,CI_High,color="Pointwise CI"),linetype="twodash",size=1.2)+
  geom_line(aes(Position,UCI_High,color="Uniform CI"),linetype="longdash",size=1.2) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "grey"),
    plot.title = element_text(hjust = 0.5,size=20),
    text = element_text(size = 15),
    legend.position="none"
  ) +
  xlab(TeX("$\\lambda(???) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/??? ) $"))+
  ggtitle(parse(text= paste("Galaxy ( ~", expression(z %~~% 2.4),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))+
  facet_zoom(xlim = c(4700, 4800),ylim=c(5,9),zoom.size=0.5,show.area=TRUE,horizontal=TRUE,split=FALSE)

ggsave(filename = paste(getwd(),"/figures/","galaxy_quadratic",".png",sep=""),
       plot = galaxy_quadratic, width = 7, height = 7)

ggsave(filename = paste(getwd(),"/figures/","galaxy_quadratic_detailed",".png",sep=""),
       plot = galaxy_quadratic_detailed, width = 7, height = 16)

stellar_cubic =ggplot(res_stellar_2,aes(Position,Actual,linetype="Actual")) +geom_line(size=1) +
  geom_line(aes(Position,Predicted,linetype="Fitted"),size=1)+
  geom_line(aes(Position,UCI_Low),color="orange",linetype="longdash",size=1) +
  geom_line(aes(Position,CI_Low),color="blue",linetype="twodash",size=1) + 
  geom_line(aes(Position,CI_High,color="Pointwise CI"),linetype="twodash",size=1)+
  geom_line(aes(Position,UCI_High,color="Uniform CI"),linetype="longdash",size=1) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "grey"),
    plot.title = element_text(hjust = 0.5,size=20),
    text = element_text(size = 15),
    legend.position="none"
  ) +
  xlab(TeX("$\\lambda(???) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/??? ) $"))+
  ggtitle(parse(text= paste("Star ( ~", expression(z %~~% 0.00),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))


stellar_cubic_detailed = ggplot(res_stellar_2,aes(Position,Actual,linetype="Actual")) +geom_line() +
  geom_line(aes(Position,Predicted,linetype="Fitted"),size=1.2)+
  geom_line(aes(Position,UCI_Low),color="orange",linetype="longdash",size=1.2) +
  geom_line(aes(Position,CI_Low),color="blue",linetype="twodash",size=1.2) + 
  geom_line(aes(Position,CI_High,color="Pointwise CI"),linetype="twodash",size=1.2)+
  geom_line(aes(Position,UCI_High,color="Uniform CI"),linetype="longdash",size=1.2) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "grey"),
    plot.title = element_text(hjust = 0.5,size=20),
    text = element_text(size = 15),
    legend.position="none"
  ) +
  xlab(TeX("$\\lambda(???) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/??? ) $"))+
  ggtitle(parse(text= paste("Star ( ~", expression(z %~~% 0.0),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))+
  facet_zoom(xlim = c(4800, 5000),ylim=c(30,65),zoom.size=0.5,show.area=TRUE,horizontal=TRUE,split=FALSE)

ggsave(filename = paste(getwd(),"/figures/","stellar_cubic",".png",sep=""),
       plot = stellar_cubic, width = 7, height = 7)

ggsave(filename = paste(getwd(),"/figures/","stellar_cubic_detailed",".png",sep=""),
       plot = stellar_cubic_detailed, width = 7, height = 16)
