
source("regression_code/set_up_regression.R")

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
  if (!is.na(single_res$selected[[x]][1])) {
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



get_graphs_linear <- function(wd,model,label,vary_seq,xlab){
  load(file = paste(wd, model,".Rdata", sep = ""))
  for(par in vary_seq){
    agg_masking <- cbind("masking",par,data.frame(t(sapply(result[[as.character(par)]], function(x) get_metrics_linear(x,"masking",beta)))))
    agg_split <- cbind("split",par,data.frame(t(sapply(result[[as.character(par)]], function(x) get_metrics_linear(x,"split",beta)))))
    agg_full <- cbind("full",par,data.frame(t(sapply(result[[as.character(par)]], function(x) get_metrics_linear(x,"full",beta)))))
    colnames(agg_masking) <- c("type","scale","error_FCR","CI_length","FSR","power_sign","power_selected","precision_selected")
    colnames(agg_split) <- colnames(agg_masking) 
    colnames(agg_full) <- colnames(agg_masking) 
    agg <- rbind(agg_masking,agg_split,agg_full)
    if(par == vary_seq[1]){
      res_agg = agg
    }
    else{
      res_agg = rbind(res_agg,agg)
    }
  }
  res_agg$type = as.character(res_agg$type)
  res_agg[res_agg$type == "masking",'type'] = "blur"
  df = aggregate(res_agg$error_FCR ~ res_agg$scale + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","FCR")
  FCR_plot <- df %>%
    ggplot( aes(x=signal, y=FCR, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 25),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(label) +
    ylab("FCR") +
    geom_hline(yintercept=0.2)+
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  df = aggregate(res_agg$CI_length ~ res_agg$scale + res_agg$type, FUN = median)
  colnames(df) <-c("signal","type","CI Length")
  CI_length_plot <- df %>%
    ggplot( aes(x=signal, y=`CI Length`, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 25),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(label) +
    ylab("CI Length (given nonzero selection)") 
  
  
  df = aggregate(res_agg$FSR~ res_agg$scale + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","FSR")
  FSR_plot <- df %>%
    ggplot( aes(x=signal, y=FSR, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 25),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(label) +
    ylab("FSR") +
    geom_hline(yintercept=0.2)+
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  df = aggregate(res_agg$power_sign ~ res_agg$scale + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","power_sign")
  power_sign_plot <- df %>%
    ggplot( aes(x=signal, y=power_sign, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 25),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(label) +
    ylab("Power Sign") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  
  df = aggregate(res_agg$power_selected ~ res_agg$scale + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","power_selected")
  power_selected_plot <- df %>%
    ggplot( aes(x=signal, y=power_selected, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 25),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(label) +
    ylab("Power Selected") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  df = aggregate(res_agg$precision_selected ~ res_agg$scale + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","precision_selected")
  precision_selected_plot <- df %>%
    ggplot( aes(x=signal, y=precision_selected, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 25),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(label) +
    ylab("Precision Selected") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  ggsave(paste(outdir,model, "_FCR.pdf",sep=""),FCR_plot)
  ggsave(paste(outdir, model, "_CI_length.pdf",sep=""),CI_length_plot)
  ggsave(paste(outdir, model, "_FSR.pdf",sep=""),FSR_plot)
  ggsave(paste(outdir, model, "_power_sign.pdf",sep=""),power_sign_plot)
  ggsave(paste(outdir, model, "_power_selected.pdf",sep=""),power_selected_plot)
  ggsave(paste(outdir, model, "_precision_selected.pdf",sep=""),precision_selected_plot)
  return(list(FCR=FCR_plot,
              CI_length=CI_length_plot,
              FSR=FSR_plot,
              power_sign=power_sign_plot,
              power_selected=power_selected_plot,
              precision_selected=precision_selected_plot))
}


get_metrics <- function(single_res,x,type="fahrmeir") {
  if(class(single_res[1]) != "list"){
    return(c(NA,NA,NA,NA,NA,NA))
  }
  CI = single_res$CIs[[x]][[type]]
  selected = which(!is.na(CI[,1]))
  beta =  single_res$beta
  if(sum(is.nan(CI[,1])) > 0){
    return(c(NA,NA,NA,NA,NA,NA))
  }
  if(length(selected) > 0){
    projected = matrix(NA, nrow = nrow(CI), ncol = 1)
    projected[selected,] = single_res$projected[[x]]
    
    error_FCR = mean(!((CI[,1] < projected) & (CI[,2] > projected)),na.rm=TRUE)
    CI_length = mean(CI[,2] - CI[,1],na.rm=TRUE)
    FSR = mean((projected > 0 & CI[,2] < 0) |  (projected < 0 & CI[,1] >0),na.rm=TRUE)  
    power_sign = mean((projected >0 & CI[,1] > 0) | (projected <0 & CI[,2] <0),na.rm=TRUE)
    power_selected = sum(abs(beta[selected] > 0)) / (sum(abs(beta)>0))
    precision_selected = sum(abs(beta[selected]) > 0) / length(selected)   
  }else {
    error_FCR = NA
    CI_length = NA
    FSR = NA
    power_sign = NA
    power_selected = 0 
    precision_selected = 0
  }
  return(c(error_FCR,CI_length,FSR,power_sign,power_selected,precision_selected))
}


get_graphs <- function(wd,model,label,vary_seq,type,xlab){
  load(file = paste(wd, model,".Rdata", sep = ""))
  
  for(par in vary_seq){
    print(par)
    agg_masking <- cbind("fission",par,data.frame(t(sapply(result[[as.character(par)]], function(x) get_metrics(x,"masking",type=type)))))
    agg_split <- cbind("split",par,data.frame(t(sapply(result[[as.character(par)]], function(x) get_metrics(x,"split",type=type)))))
    agg_full <- cbind("full",par,data.frame(t(sapply(result[[as.character(par)]], function(x) get_metrics(x,"full",type=type)))))
    
    
    colnames(agg_masking) <- c("type",label,"error_FCR","CI_length","FSR","power_sign","power_selected","precision_selected")
    colnames(agg_split) <- colnames(agg_masking) 
    colnames(agg_full) <- colnames(agg_masking) 
    agg <- rbind(agg_masking,agg_split,agg_full)
    if(par ==vary_seq[1]){
      res_agg = agg
    }
    else{
      res_agg = rbind(res_agg,agg)
    }
  }
  
  
  df = aggregate(res_agg$error_FCR ~ res_agg[,2] + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","FCR")
  FCR_plot <- df %>%
    ggplot( aes(x=signal, y=FCR, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"), legend.position= "none",
          text = element_text(size = 15), legend.text = element_text(size = 15)) +
    xlab(xlab) +
    ylab("FCR") +
    geom_hline(yintercept=0.2)+
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  df = aggregate(res_agg$CI_length ~ res_agg[,2] + res_agg$type, FUN = median)
  colnames(df) <-c("signal","type","CI Length")
  CI_length_plot <- df %>%
    ggplot( aes(x=signal, y=`CI Length`, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"), 
          text = element_text(size = 15),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(xlab) +
    ylab("CI Length (given nonzero selection)") 
  
  df = aggregate(res_agg$FSR~ res_agg[,2] + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","FSR")
  FSR_plot <- df %>%
    ggplot( aes(x=signal, y=FSR, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(xlab) +
    ylab("FSR") +
    geom_hline(yintercept=0.2)+
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  df = aggregate(res_agg$power_sign ~ res_agg[,2] + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","power_sign")
  power_sign_plot <- df %>%
    ggplot( aes(x=signal, y=power_sign, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(xlab) +
    ylab("Power Sign") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  df = aggregate(res_agg$power_selected ~ res_agg[,2] + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","power_selected")
  power_selected_plot <- df %>%
    ggplot( aes(x=signal, y=power_selected, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),
          text = element_text(size = 15),
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab(xlab) +
    ylab("Power Selected") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  df = aggregate(res_agg$precision_selected ~ res_agg[,2] + res_agg$type, FUN = mean)
  colnames(df) <-c("signal","type","precision_selected")
  precision_selected_plot <- df %>%
    ggplot( aes(x=signal, y=precision_selected, group=type, color=type)) +
    geom_line(aes(linetype = type, color = type), size = 1.5) +
    geom_point(aes(shape = type, color = type), size = 3) +
    theme(legend.title = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "grey"),legend.position = "none",
          text = element_text(size = 15),
          legend.text = element_text(size = 15)) +
    xlab(xlab) +
    ylab("Precision Selected") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 
  
  
  ggsave(paste(outdir, type, "_", model, "_FCR.pdf",sep=""),FCR_plot)
  ggsave(paste(outdir, type, "_", model, "_CI_length.pdf",sep=""),CI_length_plot)
  ggsave(paste(outdir, type, "_", model, "_FSR.pdf",sep=""),FSR_plot)
  ggsave(paste(outdir, type, "_", model, "_power_sign.pdf",sep=""),power_sign_plot)
  ggsave(paste(outdir, type, "_", model, "_power_selected.pdf",sep=""),power_selected_plot)
  ggsave(paste(outdir, type, "_", model, "_precision_selected.pdf",sep=""),precision_selected_plot)
  return(list(FCR=FCR_plot,
              CI_length=CI_length_plot,
              FSR=FSR_plot,
              power_sign=power_sign_plot,
              power_selected=power_selected_plot,
              precision_selected=precision_selected_plot))
}



get_graphs_linear(wd,"regression_linear_influential","scale",seq(2, 6, length.out = 5),"Leverage Parameter")
get_graphs_linear(wd,"regression_linear_independent","scale",seq(0, 0.2, length.out = 5),"Scale")
get_graphs_linear(wd,"regression_linear_dependent","rho",seq(-0.5, 0.5, length.out = 5),"Rho")
get_graphs(wd,"regression_poisson_influential","scale",vary_seq = seq(2, 6, length.out = 5),type='fahrmeir',xlab="Leverage Parameter")
get_graphs(wd,"regression_poisson_influential","scale",vary_seq = seq(2, 6, length.out = 5),type='sandwich',xlab="Leverage Parameter")
get_graphs(wd,"regression_poisson_varyscale","scale",vary_seq = seq(0, 0.5, length.out = 5),type='fahrmeir',xlab="Scale")
get_graphs(wd,"regression_poisson_varyscale","scale",vary_seq = seq(0, 0.5, length.out = 5),type='sandwich',xlab="Scale")
get_graphs(wd,"regression_poisson_varyRho","rho",vary_seq = seq(-0.5, 0.5, length.out = 5),type='fahrmeir',xlab="Correlation (rho)")
get_graphs(wd,"regression_poisson_varyRho","rho",vary_seq = seq(-0.5, 0.5, length.out = 5),type='sandwich',xlab="Correlation (rho)")
get_graphs(wd,"regression_logistic_varyscale","scale",vary_seq = seq(0, 0.5, length.out = 6),type='fahrmeir',xlab="Scale")
get_graphs(wd,"regression_logistic_varyscale","scale",vary_seq = seq(0, 0.5, length.out = 6),type='sandwich',xlab="Scale")
get_graphs(wd,"regression_logistic_varyprob","prob",vary_seq = seq(0.05, 0.45, 0.05),type='fahrmeir',xlab="Masking Probability")
get_graphs(wd,"regression_logistic_varyprob","prob",vary_seq = seq(0.05, 0.45, 0.05),type='sandwich',xlab="Masking Probability")




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
# Section 5 Figures: Trend Filtering 
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

