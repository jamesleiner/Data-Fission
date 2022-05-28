##############################################################################
## Generic STAR with Conditional One-Group Model
## Lihua Lei, Aaditya Ramdas & William Fithian,
##   "STAR: A general interactive framework for FDR control under structural constraints"
## Available from https://arxiv.org/pdf/1710.02776.pdf
##############################################################################

library("latex2exp")

generic.STAR <- function(
    pvals, covar, 
    fun.list = create.fun("SeqStep", pstar = 0.5),
    alpha.list = seq(0.05, 0.3, 0.01),
    type = c("naive", "model-assist"),
    score.fun = NULL, score0 = NULL,
    score.params = list(),
    find.candid.fun = NULL,
    find.candid.params = list(),
    update.mask.fun = NULL,
    update.mask.params = list(),
    plot.fun = NULL,
    plot.params = list(),
    num.steps.update.score = 5,
    prop.high = 0.05,
    prop.low = 0.01,
    plot.quiet = TRUE,
    print.quiet = FALSE,
    gif.root.title = NULL, delay = 50, loop = 1,
    width = 1440, height = 720){

    type <- type[1]
    n <- length(pvals)
    m <- length(alpha.list)
    mask <- rep(TRUE, n)
    fdp <- fdp.hat.STAR(pvals, mask, fun.list)

    score.params <- c(list(pvals = pvals,
                           covar = covar,
                           fun.list = fun.list),
                      score.params)
    find.candid.params <- c(list(covar = covar),
                            find.candid.params)
    plot.params <- c(list(covar = covar),
                     plot.params)

    if (type == "naive"){
        score <- pmin(pvals, fun.list[["g"]](pvals))
    } else {
        score <- do.call(score.fun,
                         c(list(mask = mask, score0 = score0),
                           score.params))
    }

    alphaind <- m
    alpha <- alpha.list[alphaind]
    step <- 0
    mask.return <- matrix(FALSE, n, m)
    fdp.return <- rep(0, m)
    score.return <- matrix(0, n, m)    
    num_rej.return <- rep(0, m)
    
    if (!is.null(gif.root.title)){
        gif.title <- paste0(gif.root.title, ".gif")
    }
    
    if (!plot.quiet | !is.null(gif.root.title)){
        title <- TeX(paste0("Step 0: #rej. = ", n, 
                                ", $\\widehat{FDP}$ = ",
                            round(fdp, 2)))
        if (!is.null(gif.root.title)){
            filename <- paste0(gif.root.title, "_0.png")
            png(filename, width = width, height = height)
            do.call(plot.fun,
                    c(list(score = score,
                           mask = mask,
                           main = title),
                      plot.params))
            dev.off()
        } else {
            do.call(plot.fun,
                    c(list(score = score,
                           mask = mask,
                           main = title),
                      plot.params))
            Sys.sleep(3)
        }
    }
    
    while (alphaind > 0 | fdp > alpha){
        step <- step + 1
        prop <- ifelse(alphaind < m, prop.low, prop.high)

        candid <- do.call(find.candid.fun,
                          c(list(mask = mask,
                                 prop = prop),
                            find.candid.params))
        mask <- do.call(update.mask.fun,
                        c(list(candid = candid,
                               score = score,
                               mask = mask,
                               prop = prop),
                               update.mask.params))

        fdp <- fdp.hat.STAR(pvals, mask, fun.list)
        R <- sum(mask)

        if (step %% num.steps.update.score == 0 &&
            type == "model-assist"){
            score <- do.call(score.fun,
                             c(list(mask = mask,
                                    score0 = score),
                               score.params))
        }
        
        if (!plot.quiet | !is.null(gif.root.title)){
            title <- TeX(paste0("Step ",
                                step, ": #rej. = ", R, 
                                ", $\\widehat{FDP}$ = ",
                                round(fdp, 2)))
            if (!is.null(gif.root.title)){
                filename <- paste0(gif.root.title, "_", step,
                                   ".png")
                png(filename, width = width, height = height)
                do.call(plot.fun,
                        c(list(score = score,
                               mask = mask,
                               main = title),
                          plot.params))
                dev.off()
            } else {
                do.call(plot.fun,
                        c(list(score = score,
                               mask = mask,
                               main = title),
                          plot.params))
                Sys.sleep(3)
            }
        }

        if (!print.quiet){
            print(paste0("Step ", step, ": FDP ", fdp, ", Number of Rej. ", R))
        }

        if (fdp <= alpha) {
            temp.alphaind <- min(which(alpha.list >= fdp))            
            mask.return[, temp.alphaind:alphaind] <- mask
            fdp.return[temp.alphaind:alphaind] <- fdp            
            score.return[, temp.alphaind:alphaind] <- score
            num_rej.return[temp.alphaind:alphaind] <- R
            if (alphaind <= 1 || all(alpha.list >= fdp)) break
            alphaind <- max(which(alpha.list < fdp))
            if (alphaind >= 1){
                alpha <- alpha.list[alphaind]
            } else {
                alpha <- 0
            }
        }

        if (fun.list$Const / (1 + R) > alpha) {
            break
        }
    }

    if (!is.null(gif.root.title)){
        system(paste0("convert -delay ", delay, " $(for i in $(seq 0 1 ", step, "); do echo ", gif.root.title, "_${i}.png; done) -loop ", loop, " ", gif.title))
        unlink(paste0(gif.root.title, "*.png"))
    }

    return(list(num_rej = cummax(num_rej.return),
                fdp = fdp.return, 
                score = score.return,
                mask = mask.return))
}


fdp.hat.STAR <- function(pvals, mask, 
                        fun.list = create.fun("SeqStep")){
    R <- sum(mask)
    C <- fun.list$Const
    h <- fun.list$h
    fdp <- (C + sum(h(pvals[mask]))) / (1 + R) 
    return(fdp)
}



create.fun <- function(type = c("SeqStep", "ForwardStop", "HingeExp"),
                       pstar = 0.5){
    type <- type[1]
    if (type == "SeqStep"){
        Const <- 1 / (1 - pstar)
        h <- function(pvals){Const * (pvals > pstar)}
        g <- function(pvals){pstar / (1 - pstar) * (1 - pvals)}
        sinv <- function(pvals){1 - (1 - pstar) / pstar * pvals}
        s.deriv <- function(pvals){
            rep(-pstar / (1 - pstar), length(pvals))
        }
    } else if (type == "ForwardStop"){
        Const <- -log(0.01)
        h <- function(pvals){-log(1-pvals)}
        psstar <- 1 - exp(exp(-Const) - 1)
        mu <- 1 - exp(-Const)
        p.thresh <- 1 - exp(-Const)
        tdp.thresh <- exp(-Const) * (1 - Const - exp(-Const))
        g <- function(pvals){
            sapply(pvals, function(p){
                if (p < psstar){
                    val <- p
                } else if (p < p.thresh){
                    target <- exp(-Const) * p + (1 - p) * log(1 - p)
                    val <- uniroot(function(x){
                        exp(-Const) * x + (1 - x) * log(1 - x) - target
                    }, c(0, psstar))$root
                } else {
                    target <- (1 - p) * (1 - Const - exp(-Const))
                    val <- uniroot(function(x){
                        exp(-Const) * x + (1 - x) * log(1 - x) - target
                    }, c(0, psstar))$root
                }
                return(val)
            })
        }
        sinv <- function(tdpvals){
            sapply(tdpvals, function(tdp){
                target <- exp(-Const) * tdp + (1 - tdp) * log(1 - tdp)
                if (tdp < tdp.thresh){
                    val <- 1 - target / (1 - Const - exp(-Const))
                } else {
                    val <- uniroot(function(x){
                        exp(-Const) * x + (1 - x) * log(1 - x) - target
                    }, c(psstar, 1-1e-6))$root
                }
                return(val)
            })
        }
        s.deriv <- function(pvals){
            vals <- (pmin(Const, h(pvals)) - mu) /
                (pmin(Const, h(g(pvals))) - mu)
            vals <- pmin(vals, -1)
            return(vals)
        }
    } else if (type == "HingeExp"){
        Const <- -log(0.01)/(1-pstar)
        h <- function(pvals){
            ifelse(pvals>pstar, log((1-pstar)/(1-pvals))/(1-pstar), 0)
        }
        psstar <- 1-(1-pstar)*exp((1-pstar)*(exp(-Const*(1-pstar))-1))
        mu <- 1-exp(-Const*(1-pstar))
        p.thresh <- 1-(1-pstar)*exp(-Const*(1-pstar))
        tdp.thresh <- (1-pstar)*exp(-Const*(1-pstar))*(1-Const-exp(-Const*(1-pstar)))
        HC.left <- function(p){
            tmp <- (1-p)/(1-pstar)
            ifelse(p<pstar, -mu*p, exp(-Const*(1-pstar))*p+tmp*log(tmp)-pstar*tmp)
        }
        HC.right <- function(p){
            (1-p)*(1-Const-exp(-Const*(1-pstar)))
        }
        g <- function(pvals){
            sapply(pvals, function(p){
                if (p < psstar){
                    val <- p
                } else if (p < p.thresh){
                    target <- HC.left(p)
                    val <- uniroot(function(x){
                        HC.left(x) - target
                    }, c(0, psstar))$root
                } else {
                    target <- HC.right(p)
                    val <- uniroot(function(x){
                        HC.left(x) - target
                    }, c(0, psstar))$root
                }
                return(val)
            })
        }
        sinv <- function(tdpvals){
            sapply(tdpvals, function(tdp){
                target <- HC.left(tdp)
                if (tdp < tdp.thresh){
                    val <- 1-target/(1-Const-exp(-Const*(1-pstar)))
                } else {
                    val <- uniroot(function(x){
                        HC.left(x) - target
                    }, c(psstar, 1-1e-6))$root
                }
                return(val)
            })
        }
        s.deriv <- function(pvals){
            vals <- (pmin(Const, h(pvals)) - mu) /
                (pmin(Const, h(g(pvals))) - mu)
            vals <- pmin(vals, -1)
            return(vals)
        }
    }

    return(list(h = h, g = g, sinv = sinv, s.deriv = s.deriv,
                Const = Const))
}
