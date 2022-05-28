## BH procedure (Benjamini & Hochberg, 1995)
BH <- function(pvals, alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    rej <- pvals <= alpha
    return(rej)
}
    
summary.BH <- function(pvals, H0,
                       alpha.list = seq(0.01, 0.3, 0.01)){
    n <- length(pvals)
    nfrej <- sapply(alpha.list, function(alpha){
        khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
        alpha <- alpha * khat / n
        sum(pvals[!H0] < alpha, na.rm = TRUE)
    })
    ntrej <- sapply(alpha.list, function(alpha){
        khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
        alpha <- alpha * khat / n
        sum(pvals[H0] < alpha, na.rm = TRUE)
    })
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

## Storey's BH Procedure (Storey et al. 2004)
summary.storey <- function(pvals, H0, thr = 0.5,
                           alpha.list = seq(0.01, 0.3, 0.01)){
    est_proportion_nulls=min(1,mean(pvals>thr)/(1-thr))
    pvals[pvals>thr] = Inf
    n <- length(pvals)
    nfrej <- sapply(alpha.list, function(alpha){
        khat <- max(c(0,which(sort(pvals)<=alpha/est_proportion_nulls*(1:n)/n)))
        alpha <- alpha * khat / n
        sum(pvals[!H0] < alpha, na.rm = TRUE)
    })
    ntrej <- sapply(alpha.list, function(alpha){
        khat <- max(c(0,which(sort(pvals)<=alpha/est_proportion_nulls*(1:n)/n)))
        alpha <- alpha * khat / n
        sum(pvals[H0] < alpha, na.rm = TRUE)
    })
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

summary.STAR <- function(STAR.obj, H0){
    nfrej <- apply(STAR.obj$mask, 2, function(x){
        sum(x[!H0], na.rm = TRUE)
    })
    ntrej <- apply(STAR.obj$mask, 2, function(x){
        sum(x[H0], na.rm = TRUE)
    })
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0),1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

summary.AdaPT <- function(adapt, H0, pvals){
    nfrej <- apply(adapt$s, 2, function(s){
        tmp <- (pvals <= s)
        sum(tmp[!H0], na.rm = TRUE)
    })
    ntrej <- apply(adapt$s, 2, function(s){
        tmp <- (pvals <= s)        
        sum(tmp[H0], na.rm = TRUE)
    })
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0),1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

summary.hyptree <- function(pvals, tree.el, H0,
                            alpha.list = seq(0.01, 0.3, 0.01)){
    nfrej <- sapply(alpha.list, function(alpha){
        alpha <- alpha / 2
        hyp.tree <- hFDR.adjust(pvals, tree.el, alpha)
        sum(hyp.tree@p.vals[!H0, 2] < alpha, na.rm = TRUE)
    })
    ntrej <- sapply(alpha.list, function(alpha){
        alpha <- alpha / 2
        hyp.tree <- hFDR.adjust(pvals, tree.el, alpha)
        sum(hyp.tree@p.vals[H0, 2] < alpha, na.rm = TRUE)
    })
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0), 1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

summary.lg <- function(pvals, tree, H0,
                       alpha.list = seq(0.01, 0.3, 0.01),
                       type = 1){
    rejs <- lynch.guo(pvals, tree, alpha.list = alpha.list,
                      type = type)
    nrej <- apply(rejs, 2, sum)
    FDP <- apply(rejs, 2, function(rej){
        sum(rej & (!H0)) / pmax(sum(rej), 1)
    })
    power <- apply(rejs, 2, function(rej){
        sum(rej & H0) / max(sum(H0), 1)
    })
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

summary.Lynch <- function(Lynch.DAG, H0){
    nfrej <- apply(Lynch.DAG$rej, 2, function(x){
        sum(x[!H0], na.rm = TRUE)
    })
    ntrej <- apply(Lynch.DAG$rej, 2, function(x){
        sum(x[H0], na.rm = TRUE)
    })
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0),1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}

summary.cr <- function(CR.DAG, H0){
    nfrej <- apply(CR.DAG$rej, 2, function(x){
        sum(x[!H0], na.rm = TRUE)
    })
    ntrej <- apply(CR.DAG$rej, 2, function(x){
        sum(x[H0], na.rm = TRUE)
    })
    nrej <- nfrej + ntrej
    FDP <- nfrej / pmax(nrej, 1)
    power <- ntrej / max(sum(H0),1)
    df <- data.frame(nrej = nrej, FDP = FDP, power = power)
    return(df)
}
