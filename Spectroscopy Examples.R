###################################################################################################
# Reproduce spectroscopy example in section 6.3
#
# Note that trend filtering package by Collin Politsch is a prerequisite to run these experiments. 
# Please download at: https://capolitsch.github.io/trendfiltering/
###################################################################################################

# Load in required scripts and data
source("regression_code/set_up_regression.R")
quasar <- readRDS("trendfiltering-astro-data/fig3top_quasar_spectrum.rds")
galaxy<- readRDS("trendfiltering-astro-data/fig3middle_galaxy_spectrum.rds")
stellar <- readRDS("trendfiltering-astro-data/fig3bottom_stellar_spectrum.rds")
quasar <- quasar[quasar$wavelength <= 5500 ,]
quasar <- quasar[quasar$wavelength >= 4000 ,]
galaxy <- galaxy[galaxy$wavelength <= 5500 ,]
galaxy <- galaxy[galaxy$wavelength >= 4000 ,]
stellar <- stellar[stellar$wavelength <= 5500 ,]
stellar <- stellar[stellar$wavelength >= 4000 ,]

compute_tf <- function(X,Y,std,type="SURE",onesd_rule=1,deg=1,alpha=0.1){
  n=length(Y)
  noise = rnorm(n, sd = std)
  g_Y = Y + noise
  h_Y = Y - noise
  
  if(type=="SURE"){
    tf = sure_trendfilter(X,g_Y, weights =1/std,k=deg)
    lambda_min = tf$lambda_min
    lambda_1se= tf$lambda_1se
  } else{
    tf = cv_trendfilter(X,g_Y, weights =1/std,k=deg)
    lambda_min = tf$lambda_min[3]
    lambda_1se= tf$lambda_1se[3]
  }
  if(onesd_rule) {
    fit_y = tf$fitted_values[,which(tf$lambda == lambda_1se)]
  }else{
    fit_y = tf$fitted_values[,which(tf$lambda == lambda_min)]
  }
  knots = find_knots(fit_y,X,k=deg)
  
  k= length(knots)+deg
  basis = tp(1:n, knots = knots, degree=deg,k=k)
  bs_X = cbind(rep(1, n), basis$X, basis$Z)
  n_x = ncol(bs_X) 
  P_CI = predict(lm(h_Y ~ bs_X), interval="confidence", se.fit = TRUE, weights=1/std,level = 1 - alpha)
  U_CI = unif_CI(n = n, knots = knots, Y = h_Y, X = 1:n, alpha = alpha,degree=deg)$predicted[,2:3]
  
  df <- data.frame(cbind(X,Y,P_CI$fit,U_CI))
  colnames(df) <- c("Position","Actual","Predicted","CI_Low","CI_High","UCI_Low","UCI_High")
  return(df)
}



quasar_fit = compute_tf(quasar$wavelength,quasar$flux,quasar$flux_std_err,deg=1)
galaxy_fit = compute_tf(galaxy$wavelength,galaxy$flux,galaxy$flux_std_err,deg=2)
stellar_fit = compute_tf(stellar$wavelength,stellar$flux,stellar$flux_std_err,deg=2,type="CV")



quasar_linear_detailed = ggplot(quasar_fit,aes(Position,Actual,linetype="Actual")) +geom_line() +
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
  xlab(TeX("$\\lambda(Å) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/Å ) $"))+
  ggtitle(parse(text= paste("Quasar ( ~", expression(z %~~% 2.4),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))+
  facet_zoom(xlim = c(4700, 4800),ylim=c(1,15),zoom.size=0.5,show.area=TRUE,horizontal=TRUE,split=FALSE)

galaxy_quadratic_detailed = ggplot(galaxy_fit,aes(Position,Actual,linetype="Actual")) +geom_line() +
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
  xlab(TeX("$\\lambda (Å)"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/Å ) $"))+
  ggtitle(parse(text= paste("Galaxy ( ~", expression(z %~~% 2.4),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))+
  facet_zoom(xlim = c(4700, 4800),ylim=c(5,9),zoom.size=0.5,show.area=TRUE,horizontal=TRUE,split=FALSE)



stellar_quadratic_detailed = ggplot(stellar_fit,aes(Position,Actual,linetype="Actual")) +geom_line() +
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
  xlab(TeX("$\\lambda(Å) $"))+
  ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/Å ) $"))+
  ggtitle(parse(text= paste("Star ( ~", expression(z %~~% 0.0),")",sep=""))) + 
  scale_color_manual(NULL,values = c('blue','orange'))+
  scale_linetype_manual(NULL,values = c('solid','dashed'))+
  facet_zoom(xlim = c(4000, 4200),ylim=c(35,100),zoom.size=0.5,show.area=TRUE,horizontal=TRUE,split=FALSE)
