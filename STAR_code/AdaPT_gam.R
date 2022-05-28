##############################################################################
## These functions add generalized additive models to AdaPT:
## See the paper by Lihua Lei & William Fithian,
##   "AdaPT: An interactive procedure for multiple testing with side information"
## Available from http://arxiv.org/abs/1609.06035.
##############################################################################

library("mgcv")

## Function to initialize the GAM.
## Inputs: 
##   x: covariates
##   p: p-values
##   pmin: s(x)
##   pmax: 1 - s(x)
##   cov.formula: the formula for the covariates.
gamexp.censor <- function(x, p, pmin, pmax, cov.formula){
    n <- length(p)
    z <- ifelse(p < pmin | p > pmax, 1, (pmin + 1 - pmax) / (pmin - pmax))
    df <- data.frame(x, z)
    d <- ncol(x)
    if (is.null(colnames(x))){
        names(df)[1:d] <- paste0("cor", 1:d)
    } else {
        names(df)[1:d] <- colnames(x)
    }
    names(df)[d + 1] <- "ifcensor"
    formula <- as.formula(paste("ifcensor ~", cov.formula))
    mod <- gam(formula, data = df, family = gaussian())
    pix <- pmax(predict(mod, type = "response"), 0)
    imputed.p <- sapply(1:n, function(i){
        if (p[i] < pmin[i] | p[i] > pmax[i]){
            temp <- runif(1) < pix[i]
            return(ifelse(temp, runif(1, 0, pmin[i]), runif(1, 0, 1)))
        } else {
            return(p[i])
        }        
    })
    imputed.p[imputed.p == 1] <- 0.9999
    df$ifcensor <- -log(imputed.p)
    names(df)[d + 1] <- "imputed.logpvals"
    formula2 <- as.formula(paste("imputed.logpvals ~", cov.formula))  
    mod2 <- gam(formula2, data = df, family = Gamma())
    mux <- predict(mod2, type = "response")
    return(list(mux = mux, pix = pix))
}

## Function to update pi and mu by GAMs
## Inputs: 
##   x: covariates
##   p: p-values
##   pmin: s(x)
##   pmax: 1 - s(x)
##   pix0: pi(x) at step 0. Initialized by gamexp.censor if not set.
##   mux0: mu(x) at step 0. Initialized by gamexp.censor if not set.
##   cov.formula: the formula for the covariates.
##   num.steps: maximum number of iterations for EM algorithm.
gam.mixem.censor <- function(x, p, pmin, pmax, 
                             pix0 = NULL, mux0 = NULL,
                             cov.formula,                         
                             num.steps = 10){
    n <- length(p)
    if (is.null(pix0) || is.null(mux0)){
        init <- gamexp.censor(x, p, pmin, pmax, cov.formula)
        if (is.null(pix0)){
            pix0 <- init$pix
        }
        if (is.null(mux0)){
            mux0 <- pmax(init$mux, 1)
        }
    } 
    pix <- pix0
    mux <- mux0
    df <- data.frame(x, rep(0, n))
    d <- ncol(x)
    if (is.null(colnames(x))){
        names(df)[1:d] <- paste0("cor", 1:d)
    } else {
        names(df)[1:d] <- colnames(x)
    }
    converge <- FALSE
    for (i in 1: num.steps){
        z <- ifelse(p < pmin | p > pmax, 
                    1 / (1 + 2 * (1 - pix) / pix * mux / (p^(1 / mux - 1) + (1 - p)^(1 / mux - 1))),
                    1 / (1 + (1 - pix) / pix * mux * p^(1 - 1 / mux)))
        z[z==0] <- 0.00001
        p.new <- ifelse(p < pmin | p > pmax,
                        exp((log(p) * p^(1 / mux - 1) + log(1 - p) * (1 - p)^(1 / mux - 1)) / (p^(1 / mux - 1) + (1 - p)^(1 / mux - 1))),
                        p)
        p.new[p.new == 0] <- 0.00001
        df[d + 1] <- z
        names(df)[d + 1] <- "ifcensor"
        formula <- as.formula(paste("ifcensor ~", cov.formula))
        mod <- gam(formula, data = df, family = gaussian())
        pix.new <- pmax(predict(mod, type = "response"), 0)
        df[d + 1] <- -log(p.new)
        names(df)[d + 1] <- "imputed.logpvals"
        formula2 <- as.formula(paste("imputed.logpvals ~", cov.formula))  
        mod2 <- gam(formula2, data = df, family = Gamma(), weights = z)
        mux.new <- predict(mod2, type = "response")
        if (max(abs(pix.new - pix)) < 0.01) {
            converge <- TRUE
            break
        }
        pix <- pmax(pmin(pix.new, 1), 0)
        mux <- pmin(pmax(mux.new, 1), quantile(mux.new, 0.999))
    }
    pix <- pmax(pmin(pix.new, 1), 0)
    mux <- pmin(pmax(mux.new, 1), quantile(mux.new, 0.999))
    return(list(mux = mux, pix = pix, converge = converge))
}

## Function to update s(x) using GAMs on two-groups model.
## Inputs: 
##   x: covariates
##   pvals: p-values
##   s: s(x)
##   delta: tolerance parameter (discussed in section 4.2)
##   pix0: pi(x) at step 0. Initialized as in Appendix A.2 if not set.
##   mux0: mu(x) at step 0. Initialized as in Appendix A.2 if not set.
##   num.steps: maximum number of iterations for EM algorithm.
##   if.update: update pix and mux if TRUE.
##   ...: other args passed to knot.gen and glm.mixem.censor
find.snew.mix.gam <- function(x, pvals, s, delta, cov.formula,
                              pix0 = NULL, mux0 = NULL,
                              num.steps = 10,
                              if.update = TRUE,
                              ...){
    n <- length(s)
    if (if.update){
        temp <- gam.mixem.censor(x, pvals, s, 1-s, pix0, mux0,
                                 cov.formula, num.steps)
        pix <- temp$pix
        mux <- temp$mux
    } else {
        pix <- pix0
        mux <- mux0
    }
    lower.c <- 0
    upper.c <- 1
    level.curve.frac <- function(c, mux, s, pix){
        if (c == 0){return(0)}
        pix[pix == 0] <- 0.00001
        s.new <- ((1-pix)/pix*mux*(1-c)/c)^(mux / (1 - mux))
        s.new[mux <= 1] <- 0
        s.new <- pmin(s.new, s)
        sum(s.new[s.new > 0] * (1 - pix[s.new > 0])) /
            sum(s[s.new > 0] * (1 - pix[s.new > 0]))
    }
    lower.frac <- level.curve.frac(lower.c, mux, s, pix)
    upper.frac <- level.curve.frac(upper.c, mux, s, pix)
    while (lower.frac < 1 - delta - 0.01 & upper.frac > 1 - delta + 0.01) {
        mid.c <- (lower.c + upper.c) / 2
        mid.frac <- level.curve.frac(mid.c, mux, s, pix)
        if (mid.frac < 1 - delta) {
            lower.c <- mid.c
            lower.frac <- mid.frac
        } else {
            upper.c <- mid.c
            upper.frac <- mid.frac
        }
    }
    c <- ifelse(lower.frac > 1 - delta - 0.01, lower.c, upper.c)
    s.new <- pmin(((1-pix)/pix*mux*(1-c)/c)^(mux / (1 - mux)), s)
    s.new[pix <= 0.05] <- 0
    s.new <- ifelse(s.new < (1 - delta) * s, (1 - delta) * s, s.new)    
    return(list(s.new = s.new, pix = pix, mux = mux))
}

## Updated version of AdaPT 
## Inputs: 
##   x: covariates
##   pvals: p-values
##   s0: s(x) at step 0
##   delta.high: tolerance parameter when FDPhat > max(q.list)
##   delta.low: tolerance parameter when FDPhat <= max(q.list)
##   q.list: a list of FDR levels (discussed in section 4.3)
##   disp.quiet: display the FDPhat and number of rejections
##   plot.quiet: display the plot of s(x), \hat{\pi}(x) and \hat{\mu}(x)
##   num.steps.update: number of steps between two updates of pi and mu
##   find.snew.fun: the function to update s(x)
##   ...: other args passed to knot.gen and find.snew.fun
AdaPT <- function(x, pvals, s0 = rep(0.45, length(pvals)), 
                  delta.high = 0.3, delta.low = 0.1,
                  q.list = seq(0.05, 0.3, 0.01), 
                  disp.quiet = FALSE, plot.quiet = FALSE,
                  num.steps.update = 5,
                  find.snew.fun = find.snew.mix.gam,
                  R.max = 0,
                  R.buffer.size = 3, force.delta = 0.1,
                  ...){
    n <- length(pvals)
    m <- length(q.list)
    s <- s0
    fdp <- fdp.hat(pvals, s)    
    qind <- m
    q <- q.list[qind]
    step <- 0
    s.return <- matrix(0, n, m)
    pi.return <- matrix(0, n, m)
    mu.return <- matrix(0, n, m)
    fdp.return <- rep(0, m)
    num_rej.return <- rep(0, m)
    R.buffer <- rep(NA, R.buffer.size)
    R.buffer[R.buffer.size] <- Inf
    while (fdp > q | qind > 0){
        step <- step + 1
        if (step == 1){
            pix <- NULL
            mux <- NULL
        }
        if.update <- step %% num.steps.update == 1 ||
            num.steps.update == 1
        delta <- ifelse(qind < m, delta.low, delta.high)
        temp <- find.snew.fun(x, pvals, s, delta,
                              pix0 = pix,
                              mux0 = mux,
                              if.update = if.update,
                              cov.formula = cov.formula)
        s.new <- temp$s.new
        pix <- temp$pix
        mux <- temp$mux
        if (!plot.quiet){
            par(mfrow = c(3, 1))            
            plot(s, type = 'l')
            lines(s.new, col = "red")
            plot(pix, type = 'l')
            plot(1 / (1 + mux), type = 'l')
            abline(h = 1, col = "blue")
        }
        s <- s.new
        R <- sum(pvals <= s)        
        while (all(R.buffer == R)){
            s <- s * (1 - force.delta)
            R <- sum(pvals <= s)
        }
        fdp <- fdp.hat(pvals, s)
        R.buffer[1:(R.buffer.size-1)] <- R.buffer[2:R.buffer.size]
        R.buffer[R.buffer.size] <- R
        if (R <= R.max) {
            break
        }
        if (fdp <= q) {
            temp.qind <- min(which(q.list >= fdp))
            s.return[, temp.qind:qind] <- s
            pi.return[, temp.qind:qind] <- pix
            mu.return[, temp.qind:qind] <- mux
            num_rej.return[temp.qind:qind] <- R
            fdp.return[temp.qind:qind] <- fdp
            qind <- max(which(q.list < fdp))
            if (qind <= 1) break
            q <- q.list[qind]
        }
        if (!disp.quiet){
            print(paste0("Step ", step, ": FDP ", fdp, ", Number of Rej. ", R))
        }
    }
    return(list(num_rej = cummax(num_rej.return), fdp = fdp.return, s = s.return,
                pi = pi.return, mu = mu.return))
}


## Function to update s(x) using GAMs on normal models with
## s(x) \in {0, 0.5}.
## Inputs: 
##   x: covariates
##   pvals: p-values
##   s: s(x) \in {0, 0.5}
##   delta: tolerance parameter
##   mux0: mu(x) at step 0. Initialized as in Appendix A.2 if not set.
##   num.steps: maximum number of iterations for EM algorithm.
##   if.update: update pix and mux if TRUE.
##   ...: other args passed to knot.gen and glm.mixem.censor
find.snew.normal.gam <- function(x, pvals, s, delta, cov.formula,
                                 pix0 = NULL, mux0 = NULL,
                                 num.steps = 10,
                                 if.update = TRUE){
    n <- length(s)
    pstar <- 0.5
    mask <- pvals <= s | pvals >= 1 - s
    if (if.update){
        mux <- ISS.gam.em.censor(x, pvals, mask, pstar,
                                 cov.formula, mux0,
                                 num.steps)
    } else {
        mux <- mux0
    }
    thresh <- quantile(mux[mask], delta)
    s.new <- ifelse(mask & (mux > thresh), rep(0.5, n), rep(0, n))
    return(list(s.new = s.new, pix = rep(NA, n), mux = mux))
}

