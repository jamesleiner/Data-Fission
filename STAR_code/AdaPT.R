##############################################################################
## These functions implement AdaPT from the paper:
## Lihua Lei & William Fithian,
##   "AdaPT: An interactive procedure for multiple testing with side information"
## Available from http://arxiv.org/abs/1609.06035
##############################################################################

## Function to calculate FDPhat (1 + At) / (Rt v 1).
## Inputs:
##    pvals: p-values
##    s: s(x)
fdp.hat <- function(pvals, s) {
    A <- sum(pvals > 1 - s)
    R <- sum(pvals <= s)
    return((1 + A) / max(R, 1))
}

## Function to calculate the number of rejections (p-values below s(x)).
## Inputs:
##    pvals: p-values
##    s: s(x)
num.rej <- function(pvals, s) {
    sum(pvals <= s)
}

## Function to generate the knots for deriving phi(x) using natural splines in algorithm 2
## when x is a 1-dimensional position variable. 
## Inputs: 
##   n: number of hypotheses
##   good.prior: the fraction of hypotheses thought as promising, same as 
##               gamma in section 4.3.
##   num.good.knots: number of knots in the first gamma-fraction of hypotheses
##   num.bad.knots: number of knots in the last (1-gamma)-fraction of hypotheses
knots.gen <- function(n, good.prior = 0.5,
                      num.good.knots = 3, num.bad.knots = 5){
    ind.thresh <- floor(n * good.prior)
    knots.good <- floor(seq(ceiling(ind.thresh / num.good.knots),
                            ind.thresh, length.out = num.good.knots))
    knots.bad <- floor(seq(ind.thresh, n,
                           length.out = num.bad.knots + 2))
    knots <- c(knots.good, tail(head(knots.bad, -1), -1))
    return(knots)
}

## Function to initialize Algorithm 2 (discussed in Appendix A.2)
## Inputs: 
##   x: covariates
##   p: p-values
##   pmin: s(x)
##   pmax: 1 - s(x)
##   knots: a sequence of knots to be used for natural spline. Set as the 
##          default of knot.gen if knots=NULL
##   ...: other args passed to knot.gen
init.em <- function(x, p, pmin, pmax, 
                    knots = NULL, ...){
    n <- length(p)
    z <- ifelse(p < pmin | p > pmax, 1, (pmin + 1 - pmax) / (pmin - pmax))
    if (is.null(knots)){
        knots <- knots.gen(n, ...)
    }
    mod <- lm(z ~ ns(x, knots = knots))
    pix <- pmax(predict(mod, type = "response"), 0)
    imputed.p <- sapply(1:n, function(i){
        if (p[i] < pmin[i] | p[i] > pmax[i]){
            temp <- runif(1) < (1 + pix[i]) / 2
            return(ifelse(temp, min(p[i], 1 - p[i]), max(p[i], 1 - p[i])))
        } else {
            return(p[i])
        }        
    })
    imputed.p[imputed.p == 1] <- 0.9999
    mod2 <- glm(-log(imputed.p) ~ ns(x, knots = knots), family = Gamma())
    mux <- predict(mod2, type = "response")
    return(list(mux = mux, pix = pix))
}

## Function to update pi and mu in Algorithm 2.
## Inputs: 
##   x: covariates
##   p: p-values
##   pmin: s(x)
##   pmax: 1 - s(x)
##   pix0: pi(x) at step 0. Initialized as in Appendix A.2 if not set.
##   mux0: mu(x) at step 0. Initialized as in Appendix A.2 if not set.
##   knots: a sequence of knots to be used for natural spline. Set as the 
##          default of knot.gen if knots=NULL
##   num.steps: maximum number of iterations for EM algorithm.
##   ...: other args passed to knot.gen and init.em
glm.mixem.censor <- function(x, p, pmin, pmax,
                             pix0 = NULL, mux0 = NULL,
                             knots = NULL,
                             num.steps = 5, ...){
    n <- length(p)
    if (is.null(pix0) || is.null(mux0)){
        init <- init.em(x, p, pmin, pmax, knots, ...)
        if (is.null(pix0)){
            pix0 <- init$pix
        }
        if (is.null(mux0)){
            mux0 <- pmax(init$mux, 1)
        }
    } 
    pix <- pix0
    mux <- mux0
    if (is.null(knots)){
        knots <- knots.gen(n, ...)
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
        mod <- suppressWarnings(glm(z ~ ns(x, knots = knots), family = binomial()))
        pix.new <- predict(mod, type = "response")
        mod2 <- glm(-log(p.new) ~ ns(x, knots = knots), family = Gamma(), weights = z)
        mux.new <- predict(mod2, type = "response")
        if (max(abs(pix.new - pix)) < 0.01) {
            converge <- TRUE
            break
        }
        pix <- pmax(pmin(pix.new, 1), 0)
        mux <- pmax(mux.new, 1)
    }
    pix <- pmax(pmin(pix.new, 1), 0)
    mux <- pmax(mux.new, 1)
    return(list(mux = mux, pix = pix, converge = converge))
}

## Function to update s(x) (discussed in section 4.2)
## Inputs: 
##   x: covariates
##   pvals: p-values
##   s: s(x)
##   delta: tolerance parameter (discussed in section 4.2)
##   pix0: pi(x) at step 0. Initialized as in Appendix A.2 if not set.
##   mux0: mu(x) at step 0. Initialized as in Appendix A.2 if not set.
##   ...: other args passed to knot.gen and glm.mixem.censor
find.snew.mix <- function(x, pvals, s, delta,
                          pix0 = NULL, mux0 = NULL, ...){
    n <- length(s)
    temp <- glm.mixem.censor(x, pvals, s, 1-s, pix0, mux0, ...)
    pix <- temp$pix
    mux <- temp$mux
    lower.c <- 0
    upper.c <- 1
    level.curve.frac <- function(c, mux, s, pix){
        if (c == 0){return(0)}
        pix[pix == 0] <- 0.00001
        s.new <- (1/c + (1-pix)/pix*mux*(1-c)/c)^(mux / (1 - mux))
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
    s.new <- pmin((1/c + (1-pix)/pix*mux*(1-c)/c)^(mux / (1 - mux)), s)
    s.new[pix <= 0.05] <- 0
    s.new <- ifelse(s.new < (1 - delta) * s, (1 - delta) * s, s.new)    
    return(list(s.new = s.new, pix = pix, mux = mux))
}

## Main Function (AdaPT)
## Inputs: 
##   x: covariates
##   pvals: p-values
##   s0: s(x) at step 0
##   delta: tolerance parameter (discussed in section 4.2)
##   q.list: a list of FDR levels (discussed in section 4.3)
##   quiet: display the information of intermediate steps, including 
##          FDPhat, number of rejections, the plot of s(x), \hat{\pi}(x) 
##          and \hat{\mu}(x).
##   find.snew.fun: the function to update s(x)
##   ...: other args passed to knot.gen and find.snew.fun
library('splines')
AdaPT <- function(x, pvals, s0 = rep(0.45, length(pvals)), 
                  delta = 0.1, 
                  q.list = seq(0.05, 0.3, 0.01), 
                  quiet = FALSE,
                  find.snew.fun = find.snew.mix,
                  R.max = 0,
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
    par(mfrow = c(3, 1))
    while (fdp > q | qind > 0){
        step <- step + 1
        if (step == 1){
            pix <- NULL
            mux <- NULL
        }
        temp <- find.snew.fun(x, pvals, s, delta, pix, mux, ...)
        s.new <- temp$s.new
        pix <- temp$pix
        mux <- temp$mux
        s <- s.new
        fdp <- fdp.hat(pvals, s)
        R <- sum(pvals <= s)
        if (R <= R.max) {
            break
        }
        if (fdp <= q) {
            tmp <- which(q.list < fdp)
            if (length(tmp) == 0) break
            new.qind <- max(tmp)            
            for (j in (new.qind + 1):qind){
                s.return[, j] <- s
                pi.return[, j] <- pix
                mu.return[, j] <- mux
                num_rej.return[j] <- R
                fdp.return[j] <- fdp
            }
            if (qind == 1) break
            qind <- new.qind
            q <- q.list[qind]
        }
        if (!quiet){
            print(paste0("Step ", step, ": FDP ", fdp, ", Number of Rej. ", R))
        }
    }
    return(list(num_rej = num_rej.return,
                fdp = fdp.return, s = s.return,
                pi = pi.return, mu = mu.return))
}
