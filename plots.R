source("setup.R")
library(ggpubr)
library(ggforce)
rowSes = function(m) {apply(m, 1, function(x) sd(x, na.rm =  TRUE))/sqrt(rowSums(m, na.rm = TRUE))} #standard error
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

# mu = rep(0, p)
# Sigma = toeplitz(rho^(0:(p/5-1)))
# Sigma = bdiag(lapply(1:5, function(x) Sigma))
process_result_regression = function(single_res, beta, scale) {
    error_FCR = sapply(methods, function(x) {
      if (is.na(single_res$selected[[x]])) {
        0
      } else {
        # if(any(single_res$projected[[x]] != beta[single_res$selected[[x]]]*scale)) {
        #   print(scale)
        #   stop("not correct project")
        # }
        sum(single_res$CIs[[x]][,1] > single_res$projected[[x]] |
              single_res$CIs[[x]][,2] < single_res$projected[[x]],
            na.rm = TRUE)/max(length(single_res$selected[[x]]), 1)
      }
    })
    # if (length(single_res$selected[["full"]]) == 0) {
    #   FCR_sim = 0
    # } else {
    #   FCR_sim = sum(single_res$nCIs[["full"]][,1] > single_res$projected[["sim"]] |
    #                   single_res$nCIs[["full"]][,2] < single_res$projected[["sim"]],
    #                 na.rm = TRUE)/max(length(single_res$selected[["full"]]), 1)
    # }
    
    # beta_sqdist = sapply(methods, function(x) {
    #   if (length(single_res$selected[[x]]) == 0) {
    #     0
    #   } else {
    #     sqrt(sum((single_res$projected[[x]] - beta[single_res$selected[[x]]]*scale)^2)/
    #            length(single_res$selected[[x]]))
    #   }
    # })
    # beta_lkdist = sapply(methods, function(x) {
    #   if (length(single_res$selected[[x]]) == 0) {
    #     0
    #   } else {
    #     sum(single_res$projected[[x]] - beta[selected[["x"]]]*scale)
    #   }
    # })

    CI_length = sapply(methods, function(x) {
      if (is.na(single_res$selected[[x]])) {
        NA
      } else {
        mean(single_res$CIs[[x]][,2] - single_res$CIs[[x]][,1], na.rm = TRUE)
      }
    })
    # CI_length_scaled = sapply(methods, function(x) {
    #   if (length(single_res$selected[[x]]) == 0) {
    #     NA
    #   } else {
    #     # if (x == "split") {
    #     #   X = cbind(rep(1, n - length(single_res$split_ind)),
    #     #             single_res$X[-single_res$split_ind,single_res$selected[[x]] + 1])
    #     #   exp_temp = as.vector(exp(X %*% fit_beta[[x]]))
    #     #   bread_mat = solve(t(X) %*% diag(exp_temp/(1 + exp_temp)^2) %*% X) #~bread_temp/n
    #     #   meat_mat = t(X) %*% diag((single_res$Y[-single_res$split_ind] - 1/(1 + 1/exp_temp))^2) %*% X
    #     #   cov_mat = bread_mat %*% meat_mat %*% bread_mat
    #     # } else {
    #     #   X = cbind(rep(1, n), single_res$X[,single_res$selected[[x]] + 1])
    #     #   exp_temp = as.vector(exp(X %*% fit_beta[[x]]))
    #     #   bread_mat = solve(t(X) %*% diag(exp_temp/(1 + exp_temp)^2) %*% X) #~bread_temp/n
    #     #   meat_mat = t(X) %*% diag((single_res$Y - 1/(1 + 1/exp_temp))^2) %*% X
    #     #   cov_mat = bread_mat %*% meat_mat %*% bread_mat
    #     # }
    #     # mean((single_res$CIs[[x]][-1,2] - single_res$CIs[[x]][-1,1])/
    #     #        sqrt(diag(cov_mat)[-1]), na.rm = TRUE)
    #     mean((single_res$CIs[[x]][,2] - single_res$CIs[[x]][,1])/
    #            (single_res$nCIs[[x]][,2] - single_res$nCIs[[x]][,1]), na.rm = TRUE)
    #     # mean((single_res$CIs[[x]][-1,2] - single_res$CIs[[x]][-1,1])/
    #     #        single_res$cov_vec[[x]][-1], na.rm = TRUE)
    #   }
    # })
    
    power_sign = sapply(methods, function(x) {
      if (is.na(single_res$selected[[x]])) {
        0
      } else {
        (sum(single_res$projected[[x]] > 0 & single_res$CIs[[x]][,1] > 0, na.rm = TRUE) +
          sum(single_res$projected[[x]] < 0 & single_res$CIs[[x]][,2] < 0, na.rm = TRUE)) /
          max(sum(single_res$projected[[x]] != 0, na.rm = TRUE), 1)
      }
    })
    
    # FPR_sign = sapply(methods, function(x) {
    #   if (length(single_res$selected[[x]]) == 0) {
    #     0
    #   } else {
    #     sum(fit_beta[[x]][-1] == 0 & single_res$CIs[[x]][-1,1] > 0, na.rm = TRUE)/
    #       max(sum(single_res$CIs[[x]][-1,1] > 0, na.rm = TRUE), 1)
    #   }
    # })
    FSR = sapply(methods, function(x) {
      if (is.na(single_res$selected[[x]])) {
        0
      } else {
        (sum(single_res$projected[[x]] > 0 & single_res$CIs[[x]][,2] < 0, na.rm = TRUE) + 
          sum(single_res$projected[[x]] < 0 & single_res$CIs[[x]][,1] > 0, na.rm = TRUE))/
          max(sum((single_res$CIs[[x]][,1] > 0 | single_res$CIs[[x]][,2] < 0), na.rm = TRUE), 1)
      }
    })
    # error_posbias = sapply(methods, function(x) {
    #   sum(fit_beta[[x]] > 0 & single_res$CIs[[x]][,1] > fit_beta[[x]],
    #       na.rm = TRUE)/max(sum(!is.na(single_res$CIs[[x]][,1])), 1)
    # })
    # error_negbias = sapply(methods, function(x) {
    #   sum(fit_beta[[x]] > 0 & single_res$CIs[[x]][,2] < fit_beta[[x]],
    #       na.rm = TRUE)/max(sum(!is.na(single_res$CIs[[x]][,1])), 1)
    # })
    # n_signs = sapply(methods, function(x) {
    #   sum((single_res$CIs[[x]][,1] > 0 | single_res$CIs[[x]][,2] < 0), na.rm = TRUE)})
    precision_selected = sapply(methods, function(x) {
      if (is.na(single_res$selected[[x]])) {
        0
      } else {
        sum(scale*beta[single_res$selected[[x]]] != 0)/length(single_res$selected[[x]])
      }
    })
    power_selected = sapply(methods, function(x) {
      if (is.na(single_res$selected[[x]])) {
        NA
      } else {
        sum(scale*beta[single_res$selected[[x]]] != 0)/max(sum(scale*beta != 0),1)
      }
    })
    # ll_ratio = sapply(methods, function(x) {
    #   if(is.null(single_res$ll_ratio[[x]])) {
    #     NA
    #   } else {
    #     single_res$ll_ratio[[x]]
    #   }
    # })
  if (0) {
    fit_beta[[method]] = beta
    error_FCR = sapply(single_res, function(x) {
      sum(x[,1] > beta | x[,2] < beta,
          na.rm = TRUE)/max(sum(!is.na(x[,1])), 1)
    })
    CI_length = sapply(single_res, function(x) {
      if(nrow(x) > 0) {mean(x[,2] - x[,1], na.rm = TRUE)}
      else NA
    })
    power_sign = sapply(single_res, function(x) {
      sum(beta > 0 & x[,1] > 0, na.rm = TRUE)/
        sum(beta > 0, na.rm = TRUE)
    })
    FPR_sign = sapply(single_res, function(x) {
      sum(beta == 0 & x[,1] > 0, na.rm = TRUE)/
        max(sum(x[,1] > 0, na.rm = TRUE), 1)
    })
    FSR = sapply(single_res, function(x) {
      sum(beta > 0 & x[,2] < 0, na.rm = TRUE)/
        max(sum((x[,1] > 0 | x[,2] < 0), na.rm = TRUE), 1)
    })
    error_posbias = sapply(single_res, function(x) {
      sum(beta > 0 & x[,1] > beta,
          na.rm = TRUE)/max(sum(!is.na(x[,1])), 1)
    })
    error_negbias = sapply(single_res, function(x) {
      sum(beta > 0 & x[,2] < beta,
          na.rm = TRUE)/max(sum(!is.na(x[,1])), 1)
    })
    n_signs = sapply(single_res, function(x) {
      sum((x[,1] > 0 | x[,2] < 0), na.rm = TRUE)})
    precision_selected = sapply(single_res, function(x) {
      sum(beta > 0 & !is.na(x[,1]))/max(sum(!is.na(x[,1])), 1)})
    power_selected = sapply(single_res, function(x) {sum(beta > 0 & !is.na(x[,1]))/
        sum(beta > 0, na.rm = TRUE)})
  }
  return(list(precision_selected = precision_selected, power_selected = power_selected,
              power_sign = power_sign, #FPR_sign = FPR_sign,
              # beta_sqdist = beta_sqdist,
              FSR = FSR, error_FCR = error_FCR, # FCR_sim = FCR_sim,#n_signs = n_signs,
              # ratio_B = mean(single_res$ratio_B), ratio_M = mean(single_res$ratio_M),
              # ratio_fi = mean(single_res$ratio_fi), ll_ratio =ll_ratio,
              #error_posbias = error_posbias, error_negbias = error_negbias,
              CI_length = CI_length#, CI_length_scaled = CI_length_scaled
              ))
}


process_result_independent = function(single_res, beta, scale) {
  error_FCR = sapply(methods, function(x) {
    if (sum(!is.na(single_res$CIs[[x]])) == 0) {
      0
    } else {
      sum(single_res$CIs[[x]][,1] > beta*scale |
            single_res$CIs[[x]][,2] < beta*scale,
          na.rm = TRUE)/max(sum(!is.na(single_res$CIs[[x]][,1])), 1)
    }
  })
  
  CI_length = sapply(methods, function(x) {
    if (sum(!is.na(single_res$CIs[[x]])) == 0) {
      NA
    } else {
      mean(single_res$CIs[[x]][,2] - single_res$CIs[[x]][,1], na.rm = TRUE)
    }
  })
  
  power_sign = sapply(methods, function(x) {
    if (sum(!is.na(single_res$CIs[[x]])) == 0) {
      0
    } else {
      sum(beta*scale > 0 & single_res$CIs[[x]][,1] > 0, na.rm = TRUE)/
        max(sum(beta*scale > 0, na.rm = TRUE), 1)
    }
  })
  
  FSR = sapply(methods, function(x) {
    if (sum(!is.na(single_res$CIs[[x]])) == 0) {
      0
    } else {
      sum(beta*scale > 0 & single_res$CIs[[x]][,2] < 0, na.rm = TRUE)/
        max(sum((single_res$CIs[[x]][,1] > 0 | single_res$CIs[[x]][,2] < 0), na.rm = TRUE), 1)
    }
  })
  
  precision_selected = sapply(methods, function(x) {
    if (sum(!is.na(single_res$CIs[[x]])) == 0) {
      0
    } else {
      sum(scale*beta[!is.na(single_res$CIs[[x]][,1])] > 0)/length(sum(!is.na(single_res$CIs[[x]])))
    }
  })
  
  power_selected = sapply(methods, function(x) {
    if (sum(!is.na(single_res$CIs[[x]])) == 0) {
      NA
    } else {
      sum(scale*beta[!is.na(single_res$CIs[[x]][,1])] > 0)/max(sum(scale*beta[!is.na(single_res$CIs[[x]][,1])] > 0),1)
    }
  })

  return(list(precision_selected = precision_selected, power_selected = power_selected,
              power_sign = power_sign, #FPR_sign = FPR_sign,
              # beta_sqdist = beta_sqdist,
              FSR = FSR, error_FCR = error_FCR, # FCR_sim = FCR_sim,#n_signs = n_signs,
              # ratio_B = mean(single_res$ratio_B), ratio_M = mean(single_res$ratio_M),
              # ratio_fi = mean(single_res$ratio_fi), ll_ratio =ll_ratio,
              #error_posbias = error_posbias, error_negbias = error_negbias,
              CI_length = CI_length#, CI_length_scaled = CI_length_scaled
  ))
}


######################################################
##################### Figure 7 #######################
model = "regression_linear_dependent"
load(file = paste(getwd(),"/results/", model,".Rdata", sep = ""))

scale_seq = seq(0, 0.2, length.out = 5)
post_result = list()
for (scale in scale_seq) {
  post_result[[as.character(scale)]] = lapply(result[[as.character(scale)]], function(x)
    process_result_regression(x, beta = beta, scale = scale))
}


# post_result = list()
# for (scale in scale_seq) {
#   post_result[[as.character(scale)]] = lapply(result[[as.character(scale)]], function(x)
#     process_result_independent(x, beta = beta, scale = scale))
# }

# for(i in 1:100) {
#   temp = process_result_regression(
#     result[[1]][[i]], beta = beta, scale = scale, type = "binomial", prob = prob)
# }
# 

# post_result = list()
# for (vary in vary_seq) {
#   post_result[[as.character(vary)]] = lapply(result[[as.character(vary)]], function(x)
#     trendfilter_result_regression(x))
# }
# 
# 
# 
# CI_length = sapply(result, function(y) {mean(sapply(y, function(x) {
#   (x$CIs$test[2] - x$CIs$test[1])
# }), na.rm = TRUE)})


# for (scale in scale_seq) {
#   print(mean(sapply(result[[as.character(scale)]], function(x) {
#     (x$CIs$test[2] - x$CIs$test[1])*sqrt(sum(x$X^2*exp(x$X*scale)/(1 + exp(x$X*scale))^2))
#   }), na.rm = TRUE))
# }


power_selected = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$power_selected), na.rm = TRUE)})
precision_selected = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$precision_selected), na.rm = TRUE)})
#sd_power = sapply(post_result, function(y) {rowSes(sapply(y, function(x) x$power))})
CI_length = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$CI_length), na.rm = TRUE)})
# CI_length_scaled = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$CI_length_scaled), na.rm = TRUE)})
FCR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$error_FCR), na.rm = TRUE)})
power_sign = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$power_sign), na.rm = TRUE)})
# FPR_sign = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$FPR_sign), na.rm = TRUE)})
FSR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$FSR), na.rm = TRUE)})
# FCR_sim = sapply(post_result, function(y) {mean(sapply(y, function(x) x$FCR_sim), na.rm = TRUE)})
# ratio_B = sapply(post_result, function(y) {mean(sapply(y, function(x) x$ratio_B), na.rm = TRUE)})
# ratio_M = sapply(post_result, function(y) {mean(sapply(y, function(x) x$ratio_M), na.rm = TRUE)})
# ratio_fi = sapply(post_result, function(y) {mean(sapply(y, function(x) x$ratio_fi), na.rm = TRUE)})
# beta_sqdist = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$beta_sqdist), na.rm = TRUE)})
# ll_ratio = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$ll_ratio), na.rm = TRUE)})
# n_sign = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$n_signs), na.rm = TRUE)})
# error_posbias = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$error_posbias), na.rm = TRUE)})
# error_negbias = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$error_negbias), na.rm = TRUE)})
# se_posbias = sapply(post_result, function(y) {rowSes(sapply(y, function(x) x$error_posbias))})
#FSR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$error_FSR), na.rm = TRUE)})
#type_II = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$type_II_sign), na.rm = TRUE)})
#type_II_coverage = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$type_II_coverage), na.rm = TRUE)})

# rejected_FNR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$rejected_FNR),
#                                                          na.rm = TRUE)})
# rejected_FZR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$rejected_FZR),
#                                                          na.rm = TRUE)})
# FNR_sign = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$true_FNR),
#                                                          na.rm = TRUE)})
# FZR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$true_FZR),
#                                                          na.rm = TRUE)})
# 
# FPR = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$error_FPR),
#                                                          na.rm = TRUE)})

mu_seq = as.numeric(colnames(FCR))  #+ 1; 
# a_seq = c(0.33, 1/2, 1, 2, 3)
# a_seq = c(0.5, seq(0.2, 0.8, by = 0.2))
# a_seq = 1:5
#temp_name = c("BY", paste("i-FCR-", a_seq, sep = ""))
temp_name = rownames(FCR);
# temp_name[3:6] = paste("mask-", a_seq, sep = "")
# temp_name[4:8] = paste("mask-", a_seq, sep = "")
# temp_name[seq(3,12, by = 2)] = paste("mask-", a_seq, sep = "")
# temp_name[seq(4,12, by = 2)] = paste("mask(v2)-", a_seq, sep = "")
# temp_name[seq(3,10, by = 2)] = paste("mask-", a_seq, sep = "")
# temp_name[3:7] = paste("mask-", a_seq, sep = "")
# temp_name[3] = c("split-0.5")
# temp_name[seq(3,11, by = 2)] = paste("mask-", a_seq, sep = "")
legend_name = factor(temp_name, levels = temp_name)
# exclude_methods = temp_name[-seq(5,14, by = 2)]
# exclude_methods = temp_name[-c(1,2,9)]
# exclude_methods = temp_name[seq(4,10, by = 2)]
# exclude_methods = temp_name[2]
# exclude_methods = temp_name[4:7]
exclude_methods = c()

# CI_length_slice =  CI_length_median_uniform[1:2,] #matrix(CI_length_median_uniform[,1], nrow = 2)
mode = "precision_selected"; result_mat = get(mode)
# if (mode %in% c("power_selected", "precision_selected", "power_sign")) {
#   result_mat[,1] = 0
# }
df_power = data.frame(mu_seq = rep(mu_seq, each = nrow(result_mat)),#*sqrt(n_unit),
                      power = as.vector(result_mat),
                      grp = rep(legend_name,
                                ncol(result_mat)))
p = ggplot(data = subset(df_power, !(grp %in% exclude_methods)),
           aes(x = mu_seq, y = power, group = grp, fill = grp)) +
  geom_line(aes(linetype = grp, color = grp), size = 0.8) +
  geom_point(aes(shape = grp, color = grp), size = 2.5) +
  # geom_hline(yintercept = alpha) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  # xlab(expression(paste("correlation parameter ", rho))) + 
  # xlab("Prob of new knots") +
  xlab("signal") + ylab(mode) +
  # ylab("CI length (given nonzero selection)") +
  # ylab("Median CI length") +
  # ylab("Mean standard error") + #(given nonzero rejections)
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  # scale_x_continuous(breaks = c(0.2, 0.3, 0.4, 0.6, 0.7, 0.8), limits = c(0.2,0.8)) +
  # scale_y_continuous(limits = c(0,ceiling(max(result_mat)))) +
  guides(fill=guide_legend(byrow=TRUE))
plot(p)
ggsave(filename = paste(dirname(getwd()),"/figures/",model,"_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", #width = 4.5, height = 5
       width = 4, height = 3.6)

# leg <- get_legend(p); p2 = as_ggplot(leg)
# plot(p2)
# ggsave(filename = paste(dirname(getwd()),"/figures/legend_FCR_uniform.png",sep=""),
#        plot = p2, width = 1, height = 3)
# 
# for (scale in seq(0, 2, length.out = 5)) {
#   dat = generate_logistic(n = 1000, p = 100, beta = beta*scale)
#   exp_gy1 = dat$exp_y/(dat$exp_y + (1 - dat$exp_y)*(prob/(1 - prob)))
#   exp_gy0 = dat$exp_y/(dat$exp_y + (1 - dat$exp_y)*((1 - prob)/prob))
#   logit_g1 = log(exp_gy1/(1 - exp_gy1))
#   logit_g0 = log(exp_gy0/(1 - exp_gy0))
#   logit_y = log(dat$exp_y/(1 - dat$exp_y))
#   print(summary(logit_g0-logit_y))
#   print(summary(logit_g1-logit_y))
# }
# 
# 
# mode = "FCR_sim"; result_mat = get(mode)
# df_single = data.frame(mu_seq = as.numeric(names(result_mat)),#*sqrt(n_unit),
#                       power = result_mat)
# p = ggplot(data = df_single,
#            aes(x = mu_seq, y = power)) +
#   geom_line( size = 0.8, color = "blue") +
#   geom_point( size = 2.5, color = "blue") +
#   geom_hline(yintercept = alpha) +
#   theme(legend.title = element_blank(),
#         panel.background = element_rect(fill = "white", colour = "black"),
#         panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
#         panel.grid.minor = element_line(colour = "grey"),
#         text = element_text(size = 15),
#         legend.position = "none", legend.text = element_text(size = 15)) +
#   xlab("signal") +
#   ylab("ratio of inverse of Fisher") +
#   ylab("FCR_split_2n") +
#   scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0,1.2)) +
#   guides(fill=guide_legend(byrow=TRUE))
# plot(p)
# ggsave(filename = paste(dirname(getwd()),"/figures/",model,"_", mode, ".pdf", sep = ""),
#        plot = p, device = "pdf", #width = 4.5, height = 5
#        width = 4, height = 3.6)
# 
# 
x = "split"
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

plot(beta, ylab = "coefficients", cex = 0.3,
     ylim=range(c(-0.5, 2)))

###########################################################
###################### trendfilter ########################
###########################################################
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

model = "regression_trendfilter_varySigmaProb"
load(file = paste(getwd(),"/results2/", model,".Rdata", sep = ""))
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
# c = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$c), na.rm = TRUE)})
# mean_se = sapply(post_result, function(y) {rowMeans(sapply(y, function(x) x$mean_se), na.rm = TRUE)})
# post_result = lapply(result, function(y) {lapply(y, function(x)
#   process_result_regression(x, beta = beta))})

result_mat = CI_length_median[seq(1,10,2),] - CI_length_median[seq(2,10,2),]
mat = FCR[seq(1,10,2),]
mat = CI_length_uniform[seq(1,10,2),] - CI_length[seq(1,10,2),] 
mode = "CI_length_median";  method = "masking"
if(method == "masking") {
  result_mat = get(mode)[seq(1,10,2),]
} else {
  result_mat = get(mode)[seq(2,10,2),]
}
df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = rep(prob_seq, length(sigma_seq)),
                # prob = rep(prob_seq[-1], length(sigma_seq)),
                values = as.vector(result_mat))
p = ggplot(df, aes(sigma, prob, fill = values)) + geom_tile() + 
  theme_bw() + theme(legend.position = "right") + labs(fill = "CI_diff") + xlab("Noise SD") +
  # ylab("Prob of new knots") +
  ylab("Slope range") +
  scale_fill_gradient(high="violetred1", low= "slategray1", limits = c(0, 1))
plot(p)
ggsave(filename = paste(dirname(getwd()),"/figures/", model, "_", mode,"_",method,".png",sep=""),
       plot = p, width = 5.5, height = 5.5
       #width = 4, height = 5
       )
df = data.frame(sigma = rep(sigma_seq, each = length(prob_seq)),
                prob = rep(prob_seq, length(sigma_seq)),
                typeI = as.vector(typeI[seq(1,10,2),]))

leg <- get_legend(p); p2 = as_ggplot(leg)
plot(p2)
ggsave(filename = paste(dirname(getwd()),"/figures/legend_sigmahat.png",sep=""),
       plot = p2, width = 2, height = 3)

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

