## Generate p-values
pvals.gen <- function(n, mu, rho, type){
    if (rho >= 0){
        z1 <- rnorm(n)
        z2 <- rnorm(1)
        z <- z1 * sqrt(1 - rho^2) + z2 * rho + mu
    } else {
        if (type == 1){
            new.rho <- rho / n
            Sigma <- diag(rep(1 - new.rho, n)) + new.rho
            z <- mvrnorm(1, mu, Sigma)
        } else if (type == 2){
            m <- ceiling(n / 2)
            ind1 <- sample(n, m)
            ind2 <- (1:n)[-ind1]
            Sigma <- diag(rep(0, n))
            Sigma[ind1, ind1] <- abs(rho)
            Sigma[ind2, ind2] <- abs(rho)
            Sigma[ind1, ind2] <- -abs(rho)
            Sigma[ind2, ind1] <- -abs(rho)
            diag(Sigma) <- 1            
            z <- mvrnorm(1, mu, Sigma)
        } else if (type == 3){
            Sigma <- diag(rep(0, n))
            for (i in 1:n){
                Sigma[i, ] <- rho * (1 - 2 * abs(rho))^abs(1:n - i)
            }
            diag(Sigma) <- 1
            z <- mvrnorm(1, mu, Sigma)
        }
    }
    pvals <- 1 - pnorm(z)
    return(pvals)
}

DAG.gen <- function(m, type,
                    nlayers = 4, nlinks = 3){
    if (type == "layers"){
        edges <- lapply(1:((nlayers-1)*m), function(node){
            level.node <- ceiling(node / m)
            candids <- (level.node*m+1):((level.node+1)*m)
            children <- sample(candids, nlinks)
            data.frame(from = rep(node, nlinks), to = children)
        })
        edges <- Reduce(rbind, edges)
        graph_from_data_frame(edges)
    } else if (type == "reg_layers"){
        edges <- lapply(1:((nlayers-1)*m), function(node){
            level.node <- ceiling(node / m)
            candids <- (node - m * level.node + 1:nlinks) %% m + 1
            children <- m * level.node + candids
            data.frame(from = rep(node, nlinks), to = children)
        })
        edges <- Reduce(rbind, edges)
        graph_from_data_frame(edges)
    } else if (type == "nreg") {
        nnodes <- m
        nlevels <- length(m)
        nodes <- lapply(1:nlevels, function(i){
            if (i == 1){
                1:nnodes[i]
            } else {
                (1+sum(nnodes[1:(i-1)])):sum(nnodes[1:i])
            }
        })
        edges <- lapply(1:(nlevels-1), function(level){
            parents <- nodes[[level]]
            children <- nodes[[level + 1]]
            if (length(parents) >= length(children)){
                to <- as.numeric(
                    sapply(1:length(parents),
                           function(k){
                               ind <- (k:(k+nlinks-1))%%length(children)+1
                               children[ind]
                           }))
                from <- rep(parents, each = nlinks)
                return(data.frame(from = from, to = to))
            } else {
                from <- as.numeric(
                    sapply(1:length(children),
                           function(k){
                               ind <- (k:(k+nlinks-1))%%length(parents)+1
                               parents[ind]
                           }))
                to <- rep(children, each = nlinks)
                return(data.frame(from = from, to = to))
            }
        })
        edges <- Reduce(rbind, edges)
        graph_from_data_frame(edges)
    }
}

DAG.strong.null <- function(DAG, roots, max.num = NULL){
    old.null.set <- roots
    candids <- neighborhood(DAG, order = 1,
                            nodes = roots,
                            mode = "out",
                            mindist = 1)
    candids <- unique(unlist(candids))
    while (length(candids) > 0){
        valid <- sapply(candids, function(node){
            parents <- neighborhood(DAG, 1, nodes = node,
                                    mode = "in", mindist = 1)[[1]]
            all(parents %in% old.null.set)
        })
        new.null.set <- candids[valid]
        old.null.set <- c(old.null.set, new.null.set)
        candids <- neighborhood(DAG, order = 1,
                                nodes = new.null.set,
                                mode = "out",
                                mindist = 1)
        candids <- unique(unlist(candids))
    }
    m <- ifelse(is.null(max.num), length(old.null.set),
                min(max.num, length(old.null.set)))
    return(old.null.set[1:m])
}
