
source("set_up_regression.R")
source('input_regression.R')
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
  library(reshape2)
  library(ggforce)
  library(matrixStats)
  
})

wd = paste(getwd(),"/results/",sep="")
outdir = paste(getwd(),"/figures/",sep="")

################################################################################################################################
################################################################################################################################
# Section 4 Figures: Linear Regression
################################################################################################################################
################################################################################################################################
get_metrics_linear <- function(single_res,x,beta) {
  if(is.character(single_res[[1]][1])){
    error_FCR = NA
    CI_length = NA
    FSR = NA
    power_sign = NA
    power_selected = NA
    precision_selected = NA
    return(c(error_FCR,CI_length,FSR,power_sign,power_selected,precision_selected))
  }
  if (!is.na(single_res$projected[[x]])) {
    error_FCR = mean((single_res$CIs[[x]][,1] > single_res$projected[[x]]) |  (single_res$CIs[[x]][,2] < single_res$projected[[x]]))
    CI_length = mean(single_res$CIs[[x]][,2] - single_res$CIs[[x]][,1])         
    FSR <- mean((single_res$projected[[x]] > 0 & single_res$CIs[[x]][,2] < 0) | (single_res$projected[[x]] <0 & single_res$CIs[[x]][,1] > 0)) 
    power_sign <-mean((single_res$projected[[x]] > 0 & single_res$CIs[[x]][,1] > 0) | (single_res$projected[[x]] <0 & single_res$CIs[[x]][,2] <0))
    power_selected <- sum(abs(beta[single_res$selected[[x]]]) > 0) / (sum(abs(beta)>0))
    precision_selected <- sum(abs(beta[single_res$selected[[x]]]) > 0) / length(single_res$selected[[x]])    
  }
  else{
    error_FCR = 0
    CI_length = NA
    FSR = NA
    power_sign = NA
    power_selected = 0 
    precision_selected = 0
  }
  
  return(c(error_FCR,CI_length,FSR,power_sign,power_selected,precision_selected))
}


# Figure 9 
model = "regression_linear_independent"
p=100
load(file = paste(wd, model,".Rdata", sep = ""))
beta= c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
scale_seq = seq(0, 0.2, length.out = 5)
for(scale in scale_seq){
  agg_masking <- cbind("masking",scale,data.frame(t(sapply(result[[as.character(scale)]], function(x) get_metrics_linear(x,"masking",beta)))))
  agg_split <- cbind("split",scale,data.frame(t(sapply(result[[as.character(scale)]], function(x) get_metrics_linear(x,"split",beta)))))
  agg_full <- cbind("full",scale,data.frame(t(sapply(result[[as.character(scale)]], function(x) get_metrics_linear(x,"full",beta)))))
  
  colnames(agg_masking) <- c("type","scale","error_FCR","CI_length","FSR","power_sign","power_selected","precision_selected")
  colnames(agg_split) <- colnames(agg_masking) 
  colnames(agg_full) <- colnames(agg_masking) 
  agg <- rbind(agg_masking,agg_split,agg_full)
  if(scale == scale_seq[1]){
    res_agg = agg
  }
  else{
    res_agg = rbind(res_agg,agg)
  }
}

df = aggregate(res_agg$error_FCR ~ res_agg$scale + res_agg$type, FUN = mean)
colnames(df) <-c("signal","type","FCR")
FCR_plot <- df %>%
  ggplot( aes(x=signal, y=FCR, group=type, color=type)) +
  geom_line(aes(linetype = type, color = type), size = 0.8) +
  geom_point(aes(shape = type, color = type), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("FCR") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 

df = aggregate(res_agg$CI_length ~ res_agg$scale + res_agg$type, FUN = median)
colnames(df) <-c("signal","type","CI Length")
CI_length_plot <- df %>%
  ggplot( aes(x=signal, y=`CI Length`, group=type, color=type)) +
  geom_line(aes(linetype = type, color = type), size = 0.8) +
  geom_point(aes(shape = type, color = type), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("CI Length") +
  scale_x_continuous(breaks = seq(0, 0.25, 0.05), limits = c(0,0.2)) +
  scale_y_continuous(breaks = seq(0, 0.15, 0.05), limits = c(0,0.15)) 


################################################################################################################################
################################################################################################################################
# Section 5 Figures: Trend Filtering
################################################################################################################################
################################################################################################################################
n = 200
trendfilter_result_regression = function(single_res) {
  methods = names(single_res$CI_bend)
  error_FCR = sapply(methods, function(x) {
    sum(single_res$CI_bend[[x]][,2] > single_res$project_trend[[x]] |
          single_res$CI_bend[[x]][,3] < single_res$project_trend[[x]])/n
  })
  error_typeI = sapply(methods, function(x) {
    any(single_res$CI_bend[[x]][,2] > single_res$project_trend[[x]] |
          single_res$CI_bend[[x]][,3] < single_res$project_trend[[x]])
  })
  CI_length = sapply(single_res$CI_bend, function(x) {
    if(nrow(x) > 0) {mean(x[,3] - x[,2])}
    else NA
  })
  return(list(CI_length = CI_length, error_FCR = error_FCR, error_typeI = error_typeI,
              c = unlist(single_res$c), mean_se = unlist(single_res$mean_se),
              sigma_hat = single_res$sigma_hat, nk_selected = unlist(single_res$nk_selected),
              nk_true = single_res$nk_true
  ))
}


get_metrics <-function(post_result){
  
  typeI = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMeans(sapply(y, function(x) x$error_typeI), na.rm = TRUE)})})
  FCR = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMeans(sapply(y, function(x) x$error_FCR), na.rm = TRUE)})})
  CI_length = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMeans(sapply(y, function(x) x$CI_length), na.rm = TRUE)})})
  CI_length_median = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMedians(sapply(y, function(x) x$CI_length), na.rm = TRUE)})})
  se = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMeans(sapply(y, function(x) x$mean_se), na.rm = TRUE)})})
  se_median = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMedians(sapply(y, function(x) x$mean_se), na.rm = TRUE)})})
  multiplier = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMeans(sapply(y, function(x) x$c), na.rm = TRUE)})})
  multiplier_median = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMedians(sapply(y, function(x) x$c), na.rm = TRUE)})})
  nk_true = sapply(post_result, function(z) {
    sapply(z, function(y) {mean(sapply(y, function(x) x$nk_true), na.rm = TRUE)})})
  nk_selected = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMeans(sapply(y, function(x) x$nk_selected), na.rm = TRUE)})})
  nk_selected_max = sapply(post_result, function(z) {
    sapply(z, function(y) {apply(sapply(y, function(x) x$nk_selected), 1, max)})})
  nk_ratio = sapply(post_result, function(z) {
    sapply(z, function(y) {rowMeans(sapply(y, function(x) x$nk_selected/x$nk_true), na.rm = TRUE)})})
  sigma_hat = sapply(post_result, function(z) {
    sapply(z, function(y) {mean(sapply(y, function(x) x$sigma_hat), na.rm = TRUE)})})
  return(list(typeI=typeI,
              FCR=FCR,
              CI_length=CI_length,
              CI_length_median=CI_length_median,
              se=se,
              se_median=se_median,
              multiplier=multiplier,
              multiplier_median=multiplier_median,
              nk_true=nk_true,
              nk_selected=nk_selected,
              nk_ratio=nk_ratio,
              sigma_hat=sigma_hat))
}

model = "regression_trendfilter_varySigmaProb"
load(file = paste(wd, model,".Rdata", sep = ""))
sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
slope_seq = seq(0.1, 0.5, length.out = 5)

post_result = list()
for (sigma in sigma_seq) {
  post_result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    post_result[[as.character(sigma)]][[as.character(prob)]] =
      lapply(result[[as.character(sigma)]][[as.character(prob)]], function(x)
        trendfilter_result_regression(x))
  }
}
results_trendfilter_varySigmaProb = get_metrics(post_result)


###########################################################
#Figure 13: FCR and type I error control (uniform)
###########################################################
blurring_FCR = results_trendfilter_varySigmaProb$FCR[seq(1,10,2),]
full_FCR = results_trendfilter_varySigmaProb$FCR[seq(2,10,2),]
diff_FCR = blurring_FCR - full_FCR


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_FCR))

plot_blurring_FCR_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_FCR))


plot_full_FCR_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))


blurring_typeI = results_trendfilter_varySigmaProb$typeI[seq(1,10,2),]
full_typeI = results_trendfilter_varySigmaProb$typeI[seq(2,10,2),]


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_typeI))

plot_blurring_typeI_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_typeI))


plot_full_typeI_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))


###########################################################
#Figure 14: CI length as noise and knots vary (uniform)
###########################################################
blurring_CI = results_trendfilter_varySigmaProb$CI_length_median[seq(1,10,2),]
full_CI = results_trendfilter_varySigmaProb$CI_length_median[seq(2,10,2),]
diff_CI = blurring_CI - full_CI


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_CI))

plot_blurring_CI_uniforrm = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_diff") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                # prob = rep(prob_seq[-1], length(sigma_seq)),
                values = as.vector(full_CI))

plot_full_CI_uniforrm = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_diff") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1.0))



df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                # prob = rep(prob_seq[-1], length(sigma_seq)),
                values = as.vector(diff_CI))

plot_diff_CI_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_diff") + xlab("Noise SD") +
  # ylab("Prob of new knots") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))



###########################################################
#Figure 24: Understanding relationship of multiplier to changing confidence interval lengths
###########################################################


blurring_se = results_trendfilter_varySigmaProb$se_median[seq(1,10,2),]
df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_se))


plot_blurring_se_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "Median SE") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 0.2))



blurring_multiplier = results_trendfilter_varySigmaProb$multiplier_median[seq(1,10,2),]
df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_multiplier))


plot_blurring_multiplier_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "Median multiplier") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(2.5, 6.5))



blurring_nk_selected = results_trendfilter_varySigmaProb$nk_selected[seq(1,10,2),]
df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_nk_selected))


plot_blurring_nk_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "Num. Knots Selected") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 200))

###########################################################
#Figure 25: Multiplier trajectory
###########################################################


n=200
alpha = 0.2
cvu = vector(length = 8); i = 1
for (nk in 190:197) {
  knots = sort(sample(2:(n-1), nk))
  k = length(knots) + 1
  if (length(knots) == 0) {
    bs_X = cbind(rep(1, n), 1:n)
  } else if (k + 1 < n) {
    basis = tp(1:n, knots = knots, degree=1, k = k)
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
  cvu[i] <- uniroot(function(x) {kappa*(1 + x^2/v)^(-v/2)/pi + 2*(1 - pt(x, df = v)) - alpha},
                    c(0, max(2*kappa/alpha/pi, 10)))$root
  i = i+1
}


################################################################################################################################
################################################################################################################################
# Section 5 Figures: Trend Filtering (re-do)
################################################################################################################################
################################################################################################################################

#SURE (pointwise)
model = "regression_trendfilter_varySigmaProb_pointwise_SURE"
load(file = paste(wd, model,".Rdata", sep = ""))
sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
prob_seq = prob_seq[-5]
slope_seq = seq(0.1, 0.5, length.out = 5)

post_result = list()
for (sigma in sigma_seq) {
  post_result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    post_result[[as.character(sigma)]][[as.character(prob)]] =
      lapply(result[[as.character(sigma)]][[as.character(prob)]], function(x)
        trendfilter_result_regression(x))
  }
}

regression_trendfilter_varySigmaProb_pointwise_SURE = get_metrics(post_result)

blurring_FCR = regression_trendfilter_varySigmaProb_pointwise_SURE$FCR[seq(1,12,3),]
full_FCR = regression_trendfilter_varySigmaProb_pointwise_SURE$FCR[seq(2,12,3),]
split_FCR = regression_trendfilter_varySigmaProb_pointwise_SURE$FCR[seq(3,12,3),]
diff_FCR = blurring_FCR - full_FCR


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_FCR))

plot_blurring_FCR_pointwise = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_pointwise_blurring_SURE.png",sep=""),plot=plot_blurring_FCR_pointwise)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_FCR))

plot_full_FCR_pointwise  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))

ggsave(paste(outdir,"FCR_pointwise_full_SURE.png",sep=""),plot=plot_full_FCR_pointwise)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_FCR))

plot_split_FCR_pointwise  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_pointwise_split_SURE.png",sep=""),plot=plot_split_FCR_pointwise)



model = "regression_trendfilter_varySigmaProb_uniform_SURE"
load(file = paste(wd, model,".Rdata", sep = ""))
sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
slope_seq = seq(0.1, 0.5, length.out = 5)

post_result = list()
for (sigma in sigma_seq) {
  post_result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    post_result[[as.character(sigma)]][[as.character(prob)]] =
      lapply(result[[as.character(sigma)]][[as.character(prob)]], function(x)
        trendfilter_result_regression(x))
  }
}


regression_trendfilter_varySigmaProb_uniform_SURE = get_metrics(post_result)

blurring_FCR = regression_trendfilter_varySigmaProb_uniform_SURE$FCR[seq(1,15,3),]
full_FCR = regression_trendfilter_varySigmaProb_uniform_SURE$FCR[seq(2,15,3),]
split_FCR = regression_trendfilter_varySigmaProb_uniform_SURE$FCR[seq(3,15,3),]
diff_FCR = blurring_FCR - full_FCR



df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_FCR))

plot_blurring_FCR_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_uniform_blurring_SURE.png",sep=""),plot=plot_blurring_FCR_uniform)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_FCR))

plot_full_FCR_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))

ggsave(paste(outdir,"FCR_uniform_full_SURE.png",sep=""),plot=plot_full_FCR_uniform)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_FCR))

plot_split_FCR_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_uniform_split_SURE.png",sep=""),plot=plot_split_FCR_uniform)


blurring_typeI = regression_trendfilter_varySigmaProb_uniform_SURE$typeI[seq(1,15,3),]
full_typeI = regression_trendfilter_varySigmaProb_uniform_SURE$typeI[seq(2,15,3),]
split_typeI = regression_trendfilter_varySigmaProb_uniform_SURE$typeI[seq(3,15,3),]
diff_typeI = blurring_typeI - full_typeI



df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_typeI))

plot_blurring_typeI_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "typeI") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"typeI_uniform_blurring_SURE.png",sep=""),plot=plot_blurring_typeI_uniform)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_typeI))

plot_full_typeI_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "typeI") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))

ggsave(paste(outdir,"typeI_uniform_full_SURE.png",sep=""),plot=plot_full_typeI_uniform)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_typeI))

plot_split_typeI_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "typeI") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"typeI_uniform_split_SURE.png",sep=""),plot=plot_split_typeI_uniform)



blurring_CI_length_median = regression_trendfilter_varySigmaProb_uniform_SURE$CI_length_median[seq(1,15,3),]
full_CI_length_median = regression_trendfilter_varySigmaProb_uniform_SURE$CI_length_median[seq(2,15,3),]
split_CI_length_median = regression_trendfilter_varySigmaProb_uniform_SURE$CI_length_median[seq(3,15,3),]
diff_CI_length_median = blurring_CI_length_median - full_CI_length_median



df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_CI_length_median))

plot_blurring_CI_length_median_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_length_median") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0,0.3))
ggsave(paste(outdir,"CI_length_median_uniform_blurring_SURE.png",sep=""),plot=plot_blurring_CI_length_median_uniform)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_CI_length_median))

plot_full_CI_length_median_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_length_median") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0,0.3))

ggsave(paste(outdir,"CI_length_median_uniform_full_SURE.png",sep=""),plot=plot_full_CI_length_median_uniform)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_CI_length_median))

plot_split_CI_length_median_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_length_median") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0,0.3))
ggsave(paste(outdir,"CI_length_median_uniform_split_SURE.png",sep=""),plot=plot_split_CI_length_median_uniform)



#CV (pointwise)
model = "regression_trendfilter_varySigmaProb_pointwise_CV"
load(file = paste(wd, model,".Rdata", sep = ""))
sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
prob_seq = prob_seq[-5]
slope_seq = seq(0.1, 0.5, length.out = 5)

post_result = list()
for (sigma in sigma_seq) {
  post_result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    post_result[[as.character(sigma)]][[as.character(prob)]] =
      lapply(result[[as.character(sigma)]][[as.character(prob)]], function(x)
        trendfilter_result_regression(x))
  }
}

regression_trendfilter_varySigmaProb_pointwise_CV = get_metrics(post_result)

blurring_FCR = regression_trendfilter_varySigmaProb_pointwise_CV$FCR[seq(1,12,3),]
full_FCR = regression_trendfilter_varySigmaProb_pointwise_CV$FCR[seq(2,12,3),]
split_FCR = regression_trendfilter_varySigmaProb_pointwise_CV$FCR[seq(3,12,3),]
diff_FCR = blurring_FCR - full_FCR


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_FCR))

plot_blurring_FCR_pointwise = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_pointwise_blurring_CV.png",sep=""),plot=plot_blurring_FCR_pointwise)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_FCR))

plot_full_FCR_pointwise  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))

ggsave(paste(outdir,"FCR_pointwise_full_CV.png",sep=""),plot=plot_full_FCR_pointwise)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_FCR))

plot_split_FCR_pointwise  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_pointwise_split_CV.png",sep=""),plot=plot_split_FCR_pointwise)



model = "regression_trendfilter_varySigmaProb_uniform_CV"
load(file = paste(wd, model,".Rdata", sep = ""))
sigma_seq = seq(0.05, 0.2, length.out = 4)
prob_seq = seq(0.81, 0.99, length.out = 5)
slope_seq = seq(0.1, 0.5, length.out = 5)

post_result = list()
for (sigma in sigma_seq) {
  post_result[[as.character(sigma)]] = list()
  for (prob in prob_seq) { 
    post_result[[as.character(sigma)]][[as.character(prob)]] =
      lapply(result[[as.character(sigma)]][[as.character(prob)]], function(x)
        trendfilter_result_regression(x))
  }
}


regression_trendfilter_varySigmaProb_uniform_CV = get_metrics(post_result)

blurring_FCR = regression_trendfilter_varySigmaProb_uniform_CV$FCR[seq(1,15,3),]
full_FCR = regression_trendfilter_varySigmaProb_uniform_CV$FCR[seq(2,15,3),]
split_FCR = regression_trendfilter_varySigmaProb_uniform_CV$FCR[seq(3,15,3),]
diff_FCR = blurring_FCR - full_FCR



df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_FCR))

plot_blurring_FCR_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_uniform_blurring_CV.png",sep=""),plot=plot_blurring_FCR_uniform)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_FCR))

plot_full_FCR_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))

ggsave(paste(outdir,"FCR_uniform_full_CV.png",sep=""),plot=plot_full_FCR_uniform)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_FCR))

plot_split_FCR_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "FCR") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"FCR_uniform_split_CV.png",sep=""),plot=plot_split_FCR_uniform)


blurring_typeI = regression_trendfilter_varySigmaProb_uniform_CV$typeI[seq(1,15,3),]
full_typeI = regression_trendfilter_varySigmaProb_uniform_CV$typeI[seq(2,15,3),]
split_typeI = regression_trendfilter_varySigmaProb_uniform_CV$typeI[seq(3,15,3),]
diff_typeI = blurring_typeI - full_typeI



df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_typeI))

plot_blurring_typeI_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "typeI") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"typeI_uniform_blurring_CV.png",sep=""),plot=plot_blurring_typeI_uniform)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_typeI))

plot_full_typeI_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "typeI") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))

ggsave(paste(outdir,"typeI_uniform_full_CV.png",sep=""),plot=plot_full_typeI_uniform)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_typeI))

plot_split_typeI_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "typeI") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
ggsave(paste(outdir,"typeI_uniform_split_CV.png",sep=""),plot=plot_split_typeI_uniform)



blurring_CI_length_median = regression_trendfilter_varySigmaProb_uniform_CV$CI_length_median[seq(1,15,3),]
full_CI_length_median = regression_trendfilter_varySigmaProb_uniform_CV$CI_length_median[seq(2,15,3),]
split_CI_length_median = regression_trendfilter_varySigmaProb_uniform_CV$CI_length_median[seq(3,15,3),]
diff_CI_length_median = blurring_CI_length_median - full_CI_length_median



df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(blurring_CI_length_median))

plot_blurring_CI_length_median_uniform = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_length_median") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0,0.3))
ggsave(paste(outdir,"CI_length_median_uniform_blurring_CV.png",sep=""),plot=plot_blurring_CI_length_median_uniform)

df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(full_CI_length_median))

plot_full_CI_length_median_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_length_median") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0,0.3))

ggsave(paste(outdir,"CI_length_median_uniform_full_CV.png",sep=""),plot=plot_full_CI_length_median_uniform)


df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = 1- rep(prob_seq, length(sigma_seq)),
                values = as.vector(split_CI_length_median))

plot_split_CI_length_median_uniform  = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_length_median") + xlab("Noise SD") +
  ylab("Probability of new knots") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0,0.3))
ggsave(paste(outdir,"CI_length_median_uniform_split_CV.png",sep=""),plot=plot_split_CI_length_median_uniform)




################################################################################################################################
################################################################################################################################
# Section A.1 Figures: Poisson Regression
################################################################################################################################
################################################################################################################################


###########################################################
#Figure 16: Example Confidence Intervals
###########################################################



###########################################################
#Figure 17: Simulations varying signal strength
###########################################################



###########################################################
#Figure 18: Simulations varying correlation parameter
###########################################################\


model = "regression_poisson_dependent_varyRho"
load(file = paste(wd, model,".Rdata", sep = ""))


################################################################################################################################
################################################################################################################################
# Section A.2 Figures: Logistic Regression
################################################################################################################################
################################################################################################################################


process_result = function(single_res, mu_0 = 0) {
  power = sapply(single_res$rejections, function(x) {sum(single_res$mu > mu_0 & x)/sum(single_res$mu > mu_0)})
  error_FDR = sapply(single_res$rejections, function(x) {sum(single_res$mu == mu_0 & x)/max(sum(x), 1)})
  error_FCR = sapply(single_res$CIs, function(x) {
    sum(x[,1] > single_res$mu | x[,2] < single_res$mu, na.rm = TRUE)/max(sum(!is.na(x[,1])), 1)
  })
  s = lapply(single_res$CIs, function(x) {
    s = rep(0, nrow(x))
    s = 1*(x[,1] > mu_0) + (-1)*(x[,2] < mu_0)
    return(s)
  })
  error_FSR = sapply(s, function(x) {
    sum((single_res$mu == mu_0 & x != 0) | (single_res$mu > mu_0 & x != 1), na.rm = TRUE)/max(sum(!is.na(x)), 1)
  })
  error_FPR = sapply(s, function(x) {
    sum(single_res$mu == mu_0 & x == 1, na.rm = TRUE)/max(sum(!is.na(x)), 1)
  })
  rejected_FNR = sapply(s, function(x) {
    sum(single_res$mu > mu_0 & x == -1, na.rm = TRUE)/max(sum(!is.na(x) & single_res$mu > mu_0), 1)
  })
  rejected_FZR = sapply(s, function(x) {
    sum(single_res$mu > mu_0 & x == 0, na.rm = TRUE)/max(sum(!is.na(x) & single_res$mu > mu_0), 1)
  })
  true_FNR = sapply(s, function(x) {
    sum(single_res$mu > mu_0 & x == -1, na.rm = TRUE)/max(sum(single_res$mu > mu_0), 1)
  })
  true_FZR = sapply(s, function(x) {
    sum(single_res$mu > mu_0 & x == 0, na.rm = TRUE)/max(sum(single_res$mu > mu_0), 1)
  })
  type_II_sign = sapply(s, function(x) {
    sum(single_res$mu > mu_0 & (is.na(x) | x != 1))/max(sum(single_res$mu > mu_0), 1)
  })
  type_II_coverage = sapply(single_res$CIs, function(x) {
    sum(single_res$mu > mu_0 & (is.na(x[,1]) | x[,1] > single_res$mu | x[,2] < single_res$mu))/max(sum(single_res$mu > mu_0), 1)
  })
  
  CI_length = sapply(single_res$CIs, function(x) {
    if(nrow(x) > 0) {mean(x[,2] - x[,1], na.rm = TRUE)}
    else NA
  })
  return(list(power = power, CI_length = CI_length,
              error_FDR = error_FDR, error_FCR = error_FCR,
              error_FSR = error_FSR, error_FPR = error_FPR,
              rejected_FNR = rejected_FNR, rejected_FZR = rejected_FZR,
              true_FNR = true_FNR, true_FZR = true_FZR,
              type_II_sign = type_II_sign, type_II_coverage = type_II_coverage
  ))
}

process_result_regression = function(single_res, beta, scale) {
  error_FCR = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      0
    } else {
      sum(single_res$CIs[[x]][,1] > single_res$projected[[x]] |
            single_res$CIs[[x]][,2] < single_res$projected[[x]],
          na.rm = TRUE)/max(length(single_res$selected[[x]]), 1)
    }
  })
  if (length(single_res$selected[["full"]]) == 0) {
    FCR_sim = 0
  } else {
    FCR_sim = sum(single_res$nCIs[["full"]][,1] > single_res$projected[["sim"]] |
                    single_res$nCIs[["full"]][,2] < single_res$projected[["sim"]],
                  na.rm = TRUE)/max(length(single_res$selected[["full"]]), 1)
  }
  
  beta_sqdist = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      0
    } else {
      sqrt(sum((single_res$projected[[x]] - beta[single_res$selected[[x]]]*scale)^2)/
             length(single_res$selected[[x]]))
    }
  })
  CI_length = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      NA
    } else {
      mean(single_res$CIs[[x]][,2] - single_res$CIs[[x]][,1],na.rm = TRUE)
    }
  })
  CI_length_scaled = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      NA
    } else {
      mean((single_res$CIs[[x]][,2] - single_res$CIs[[x]][,1])/
             (single_res$nCIs[[x]][,2] - single_res$nCIs[[x]][,1]), na.rm = TRUE)
    }
  })
  power_sign = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      0
    } else {
      sum(single_res$projected[[x]] > 0 & single_res$CIs[[x]][,1] > 0, na.rm = TRUE)/
        max(sum(single_res$projected[[x]] > 0, na.rm = TRUE), 1)
    }
  })
  FSR = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      0
    } else {
      sum(single_res$projected[[x]] > 0 & single_res$CIs[[x]][,2] < 0, na.rm = TRUE)/
        max(sum((single_res$CIs[[x]][,1] > 0 | single_res$CIs[[x]][,2] < 0), na.rm = TRUE), 1)
    }
  })
  precision_selected = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      0
    } else {
      sum(scale*beta[single_res$selected[[x]]] > 0)/length(single_res$selected[[x]])
    }
  })
  power_selected = sapply(methods, function(x) {
    if (length(single_res$selected[[x]]) == 0) {
      NA
    } else {
      sum(scale*beta[single_res$selected[[x]]] > 0)/max(sum(scale*beta[single_res$selected[[x]]] > 0),1)
    }
  })
  ll_ratio = sapply(methods, function(x) {
    if(is.null(single_res$ll_ratio[[x]])) {
      NA
    } else {
      single_res$ll_ratio[[x]]
    }
  })
  return(list(precision_selected = precision_selected, power_selected = power_selected,
              power_sign = power_sign,
              beta_sqdist = beta_sqdist,
              FSR = FSR, error_FCR = error_FCR, FCR_sim = FCR_sim,
              ratio_B = mean(single_res$ratio_B), ratio_M = mean(single_res$ratio_M),
              ratio_fi = mean(single_res$ratio_fi), ll_ratio =ll_ratio,
              CI_length = CI_length, CI_length_scaled = CI_length_scaled
  ))
}



###########################################################
#Figure 19: Example Confidence Intervals
###########################################################

model = "regression_logistic_varyscale"
load(file = paste(getwd(),"/results2/", model,".Rdata", sep = ""))
single_res = result[[0.5]][[1]]

lb = single_res$CIs[[x]][,1]; lb[is.na(lb)] = 0
ub = single_res$CIs[[x]][,2]; ub[is.na(ub)] = 0
projected_beta = single_res$projected[[x]]
plot(1:100, beta*scale, ylab = "coefficients", cex = 0.3, xlab = "",
     ylim=range(c(-0.5, 0.5)))
points(single_res$selected[[x]], projected_beta, col = "blue", pch = 4)
arrows(single_res$selected[[x]], lb, single_res$selected[[x]], ub,
       length=0.05, angle=90, code=3)

ind = as.vector(single_res$CIs[[x]][,1] > projected_beta |
                  single_res$CIs[[x]][,2] < projected_beta)
ind[is.na(ind)] = FALSE
arrows(single_res$selected[[x]][which(ind)], lb[ind],
       single_res$selected[[x]][which(ind)], ub[ind], 
       length=0.05, angle=90, code=3, col = "red")
sum(ind)/length(single_res$selected[[x]])
###########################################################
#Figure 20: Simulations varying signal strength
###########################################################

model = "regression_logistic_varyscale"
load(file = paste(getwd(),"/results2/", model,".Rdata", sep = ""))

beta = c(1, 0, rep(1,20), rep(0, 100 - 31), rep(2,9))
beta= c(0,beta)
methods = c("masking","full","split")
scale_seq = seq(0, 0.5, length.out = 6)
post_result = list()
for (scale in scale_seq) {
  post_result[[as.character(scale)]] = lapply(result[[as.character(scale)]], function(x)
    process_result_regression(x, beta = beta, scale = scale))
}


for (scale in scale_seq) {
  for(iter in 1:length(result[[as.character(scale)]])){
    CI = result[[as.character(scale)]][[iter]]$CIs$masking_update
    nCI = result[[as.character(scale)]][[iter]]$nCIs$masking_update
    CI_length = mean(CI[,2] - CI[,1])
    CI_length_scaled = mean((CI[,2] - CI[,1])/(nCI[,2] - nCI[,1]))
    true_beta = beta[result[[as.character(scale)]][[iter]]$selected$masking]*scale
    
    FCR = sum((CI[,1] > true_beta) | (CI[,2] < true_beta))/max(length(true_beta),1)
    FSR = sum((CI[,1] <0 & true_beta > 0) | (CI[,2] > 0 & true_beta < 0))/max(length(which(true_beta != 0)),1)
    power_sign = sum(true_beta >0  & CI[,1] > 0) / max(length(which(true_beta != 0)),1)
    
    post_result[[as.character(scale)]][[iter]]$CI_length = c(post_result[[as.character(scale)]][[iter]]$CI_length,"masking_update" = CI_length)
    post_result[[as.character(scale)]][[iter]]$CI_length_scaled = c(post_result[[as.character(scale)]][[iter]]$CI_length_scaled,"masking_update" = CI_length_scaled)
    post_result[[as.character(scale)]][[iter]]$error_FCR  = c(post_result[[as.character(scale)]][[iter]]$error_FCR,"masking_update" = FCR)
    post_result[[as.character(scale)]][[iter]]$FSR  = c(post_result[[as.character(scale)]][[iter]]$FSR,"masking_update" = FSR)
    post_result[[as.character(scale)]][[iter]]$power_sign  = c(post_result[[as.character(scale)]][[iter]]$power_sign,"masking_update" = power_sign)
    
  }
}


CI_length = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$CI_length), na.rm = TRUE)})
CI_length_scaled = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$CI_length_scaled), na.rm = TRUE)})
FCR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$error_FCR), na.rm = TRUE)})
power_sign = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$power_sign), na.rm = TRUE)})
FSR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$FSR), na.rm = TRUE)})
power_selected = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$power_selected), na.rm = TRUE)})
precision_selected = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$precision_selected), na.rm = TRUE)})


FCR_plot <- melt(t(FCR)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("FCR") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 


CI_plot <- melt(t(CI_length)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("CI length") +
  scale_y_continuous(breaks = seq(0, 1.5, 0.3), limits = c(0,1.5)) 


FSR_plot <- melt(t(FSR)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("FSR") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))


power_sign_plot <- melt(t(power_sign)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("Power Sign") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))


power_selected_plot <- melt(t(power_selected)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("Power Selected") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))



precision_selected_plot <- melt(t(precision_selected)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("Precision Selected") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))



###########################################################
#Figure 21: Simulations varying blurring parameter
###########################################################

model = "regression_logistic_varyprob"
load(file = paste(dirname, model,".Rdata", sep = ""))

beta = c(1, 0, rep(1,20), rep(0, 100 - 31), rep(2,9))
beta= c(0,beta)
methods = c("masking","full","split")
prob_seq = seq(0.05, 0.45, 0.05)
scale = 0.5
post_result = list()

for (prob in prob_seq) {
  post_result[[as.character(prob)]] = lapply(result[[as.character(prob)]], function(x)
    process_result_regression(x, beta = beta, scale = 0.5))
}

for (prob in prob_seq) {
  for(iter in 1:length(result[[as.character(prob)]])){
    CI = result[[as.character(prob)]][[iter]]$CIs$masking_update
    nCI = result[[as.character(prob)]][[iter]]$nCIs$masking_update
    CI_length = mean(CI[,2] - CI[,1])
    CI_length_probd = mean((CI[,2] - CI[,1])/(nCI[,2] - nCI[,1]))
    true_beta = beta[result[[as.character(prob)]][[iter]]$selected$masking]*scale
    
    FCR = sum((CI[,1] > true_beta) | (CI[,2] < true_beta))/max(length(true_beta),1)
    FSR = sum((CI[,1] <0 & true_beta > 0) | (CI[,2] > 0 & true_beta < 0))/max(length(which(true_beta != 0)),1)
    power_sign = sum(true_beta >0  & CI[,1] > 0) / max(length(which(true_beta != 0)),1)
    
    post_result[[as.character(prob)]][[iter]]$CI_length = c(post_result[[as.character(prob)]][[iter]]$CI_length,"masking_update" = CI_length)
    post_result[[as.character(prob)]][[iter]]$CI_length_probd = c(post_result[[as.character(prob)]][[iter]]$CI_length_probd,"masking_update" = CI_length_probd)
    post_result[[as.character(prob)]][[iter]]$error_FCR  = c(post_result[[as.character(prob)]][[iter]]$error_FCR,"masking_update" = FCR)
    post_result[[as.character(prob)]][[iter]]$FSR  = c(post_result[[as.character(prob)]][[iter]]$FSR,"masking_update" = FSR)
    post_result[[as.character(prob)]][[iter]]$power_sign  = c(post_result[[as.character(prob)]][[iter]]$power_sign,"masking_update" = power_sign)
    
  }
}

CI_length = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$CI_length), na.rm = TRUE)})
CI_length_scaled = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$CI_length_scaled), na.rm = TRUE)})
FCR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$error_FCR), na.rm = TRUE)})
power_sign = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$power_sign), na.rm = TRUE)})
FSR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$FSR), na.rm = TRUE)})
power_selected = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$power_selected), na.rm = TRUE)})
precision_selected = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$precision_selected), na.rm = TRUE)})


FCR_plot <- melt(t(FCR)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Flip probability") +
  ylab("FCR") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 


CI_plot <- melt(t(CI_length)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Flip probability") +
  ylab("CI length") +
  scale_y_continuous(breaks = seq(0, 1.5, 0.3), limits = c(0,1.5)) 


FSR_plot <- melt(t(FSR)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +sapply(result[[as.character(influ)]], function(x) get_metrics_linear(x,"masking",beta))
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Flip probability") +
  ylab("FSR") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))


power_sign_plot <- melt(t(power_sign)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("Power Sign") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))


power_selected_plot <- melt(t(power_selected)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("Power Selected") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))



precision_selected_plot <- melt(t(precision_selected)) %>%
  ggplot( aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_line(aes(linetype = Var2, color = Var2), size = 0.8) +
  geom_point(aes(shape = Var2, color = Var2), size = 2.5) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("signal") +
  ylab("Precision Selected") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))

