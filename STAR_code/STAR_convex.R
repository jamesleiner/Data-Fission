##############################################################################
## STAR for Convex Region Detection
##############################################################################

source("STAR_code/generic_STAR.R")
library("mgcv")

#### Find_Candidate_Fun
dir.cut <- function(cloud, angle, prop){
    x <- cloud[, 1]
    y <- cloud[, 2]
    n <- length(x)
    if (min(abs(angle - c(0.5, 1, 1.5, 2))) > 1e-4){
        dist <- (y - tanpi(angle) * x) * sinpi(angle)
    } else if (abs(angle - 0.5) < 1e-4) {
        dist <- -x
    } else if (abs(angle - 1.5) < 1e-4) {
        dist <- x
    } else if (abs(angle - 1) < 1e-4) {
        dist <- -y
    } else if (abs(angle - 2) < 1e-4) {
        dist <- y
    } 
    thresh <- quantile(dist, probs = 1 - prop)
    reveal <- (dist >= thresh)
    return(reveal)
}

find.corners <- function(covar, mask, prop, 
                         num.angle = 100,
                         min.cut = 1){
    if (!class(covar) %in% c("data.frame", "matrix")){
        covar <- matrix(covar, ncol = 1)
    }
    n <- nrow(covar)
    basis <- 2 / num.angle
    prop <- max(prop, min.cut/length(mask))
    candids <- lapply(basis*(1:num.angle), function(angle){
        temp <- rep(FALSE, n)
        temp[mask] <- dir.cut(covar[mask, ], angle, prop)
        temp
    })
    return(candids)
}

#### Update_Mask_Fun
## Naive Update
convex.mask.update <- function(candid, score, mask, prop,
                               dir = c("max", "min")){
    dir <- dir[1]
    candid.vals <- sapply(candid, function(set){
        mean(score[set])
    })
    if (dir == "max"){
        reveal.inds <- candid[[which.max(candid.vals)]]
    } else {
        reveal.inds <- candid[[which.min(candid.vals)]]
    }
    mask[reveal.inds] <- FALSE
    return(mask)
}

#### Score_Fun
STAR.gam.em.censor <- function(covar, pvals, mask,
                               cov.formula,
                               fun.list = create.fun("SeqStep"),
                               score0 = NULL,
                               num.steps = 5){
    n <- length(pvals)
    h <- fun.list$h
    g <- fun.list$g
    sinv <- fun.list$sinv
    s.deriv <- fun.list$s.deriv
    Const <- fun.list$Const
    init.tdpvals <- pmin(pvals, g(pvals))
    ## init.ref.tdpvals <- pmax(pvals, g(pvals))
    ## This is a mistake pointed by Heejong Bong from CMU
    init.ref.tdpvals <- sinv(init.tdpvals)
    if (is.null(score0)){
        mux0 <- mean(-log(init.tdpvals))
        mux <- rep(mux0, n)        
    } else {
        mux <- score0
    }
    df <- data.frame(covar, rep(0, n))
    d <- ncol(covar)
    if (is.null(colnames(covar))){
        names(df)[1:d] <- paste0("cor", 1:d)
    } else {
        names(df)[1:d] <- colnames(covar)
    }
    names(df)[d+1] <- "imputed_logpvals"
    converge <- FALSE
    tdpvals <- ifelse(mask, init.tdpvals, pvals)
    ref.tdpvals <- ifelse(mask, init.ref.tdpvals, pvals)
    for (i in 1:num.steps){
        temp <- -1 / s.deriv(ref.tdpvals)
        imputed.logpvals <- ifelse(
            mask,
            (tdpvals^(1/mux-1)*(-log(tdpvals))+ref.tdpvals^(1/mux-1)*temp*(-log(ref.tdpvals))) / (tdpvals^(1/mux-1)+ref.tdpvals^(1/mux-1)*temp),
            -log(tdpvals))

        imputed.logpvals[imputed.logpvals <= 1e-6] <- 1e-6
        df[d + 1] <- imputed.logpvals
        formula <- as.formula(paste("imputed_logpvals ~", cov.formula))  
        mod <- gam(formula, data = df, family = Gamma())
        mux.new <- predict(mod, type = "response")
        
        if (max(abs(mux.new - mux)) < 0.01) {
            converge <- TRUE
            break
        }
        mux <- pmax(mux.new, 1)
    }
    return(mux)
}


#### Plot_Fun
plot.rej.convex <- function(x, mask, main, cex,
                            col.bg = "#FFB6C1",
                            col.fg = "#800000",
                            ...){
    par(...)
    color <- ifelse(mask, col.fg, col.bg)
    plot(x[, 1], x[, 2], type = "n", xlab = "", ylab = "",
         main = main)
    points(x[, 1], x[, 2], col = color, cex = cex)
}

plot.strength.convex <- function(x, strength, main, cex,
                                 col.fg = "#000080",
                                 col.bg = "#ADD8E6",
                                 ...){
    par(...)    
    rbPal <- colorRampPalette(c(col.bg, col.fg))
    color <- rbPal(100)[as.numeric(cut(strength, breaks=100))]
    plot(x[, 1], x[, 2], type = "n", xlab = "", ylab = "",
         main = main)
    points(x[, 1], x[, 2], col = color, cex = cex)
}
    
convex.plot <- function(covar, mask, score,
                        main.rej, main.strength,
                        add.plot = NULL,
                        col.bg.rej = "#FFB6C1",
                        col.fg.rej = "#800000",
                        col.bg.strength = "#000080",
                        col.fg.strength = "#ADD8E6",
                        mar = c(2, 5, 5, 2),
                        cex.main = 2,
                        cex = 3,
                        ...){
    par(mar = mar, cex.main = cex.main)
    if (!is.null(add.plot)){
        par(mfrow = c(1, 3))
        eval(add.plot)
    } else {
        par(mfrow = c(1, 2))
    }
    plot.rej.convex(covar, mask, main = main.rej,
                    cex = cex,
                    col.bg = col.bg.rej,
                    col.fg = col.fg.rej,
                    ...)
    plot.strength.convex(covar, score, main = main.strength,
                         cex = cex,
                         col.bg = col.bg.strength,
                         col.fg = col.fg.strength,
                         ...)
}

####
STAR.convex <- function(
    pvals, x, 
    fun.list = create.fun("SeqStep", pstar = 0.5),
    alpha.list = seq(0.01, 0.3, 0.01),
    type = c("naive", "model-assist"),
    ...){

    type <- type[1]
    find.candid.fun <- find.corners
    plot.fun <- convex.plot
    update.mask.fun <- convex.mask.update
    if (type == "naive"){
        score.fun <- NULL
    } else if (type == "model-assist"){
        score.fun <- STAR.gam.em.censor
    }

    generic.STAR(pvals = pvals, covar = x, type = type,
                fun.list = fun.list, alpha.list = alpha.list,
                score.fun = score.fun,
                find.candid.fun = find.candid.fun,
                update.mask.fun = update.mask.fun,
                plot.fun = plot.fun,
                ...)
}
