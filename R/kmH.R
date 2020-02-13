kmns.soln <- function(k, x,nstart =  k * prod(dim(x))) {
    if (k == 1)
        kmeans(x, centers = 1) else 
                                   kmeans(x, centers = k, iter.max = 100, nstart = nstart)
}

sc.trW <-  function(k, res, x) (res[[k]]$tot.withinss * k^(2/dim(x)[2]))

dif <- function(k, res, x) (sc.trW(k = k - 1, res = res, x = x) - sc.trW(k = k, res = res, x = x))

KL.meas <- function(k, res, x, nscat) {
    kl <- ifelse(min(res[[k]]$size) <= nscat, -Inf, ifelse(k == 1, NA, abs(dif(k = k, res = res, x = x)/dif(k = k + 1, res = res, x = x))))
    cat("k = ", k, " WSS = ", res[[k]]$tot.withinss," Krzanowski and Lai criterion = ", kl, "\n")
    kl
}

kmns.KL <- function(x, maxclus, scat = 0.001, verbose = F, desired.ncores = 2, nstart =   prod(dim(x))) {
    
    ## x = dataset (matrix of observations: each row is an observation vector)
    ## maxclus = the maximum number of clusters.

    ## identify scatter
    ##
    km.sol <- kmns.soln(k = maxclus + 1, x = x, nstart = nstart)
    nscat <- max(1, round(scat * dim(x)[1]))

    x.red <- x[unlist(lapply(km.sol$cluster, FUN =
                                                 function(x) (any(x==setdiff(1:(maxclus+1), 
                                                 (1:(maxclus+1))[km.sol$size <= nscat]))))), ]
    x.scat <- x[unlist(lapply(km.sol$cluster, FUN =
                                                  function(x) (any(x==(1:(maxclus+1))[km.sol$size <= nscat])))), ]
    id.scat <- (1:nrow(x))[unlist(lapply(km.sol$cluster, FUN =
                                                  function(x) (any(x==(1:(maxclus+1))[km.sol$size <= nscat]))))]

    ## remove these groups from the clustering and downgrade maxclus

    maxclus <- sum(km.sol$size > nscat) - 1
    desired.ncores <- min(detectCores(),desired.ncores)
    registerDoMC(desired.ncores)

    kmns.results <- foreach(i = 1:(maxclus+1),.export=c("kmns.soln")) %dopar% { 
        if (verbose) {
            cat("Beginning k-means for k = ", i, "\n")
        }
        kmns.soln(i, x = x.red,nstart = nstart)
    }
    kmns.KL <- lapply(X = 1:maxclus, FUN = KL.meas, res = kmns.results, x = x.red, nscat = nscat)
    list(kmns.results = kmns.results, KL = kmns.KL, x.red = x.red, x.scat = x.scat, maxclus = maxclus, id.scat = id.scat)
}

prob.dist <- function(i, j, means, wss, size) {
    p <- ncol(means)
    w1 <- wss/(size-1)/p
    mu.sq.diff <- apply((means[i,] - means[j,])^2, 1, sum)
    sig2i <- w1[i]
    sig2j <- w1[j]
    ncpar <-  sig2i * mu.sq.diff / (sig2i - sig2j)^2
    qpar <- sig2j *  mu.sq.diff / (sig2i - sig2j)^2
    ncpar[is.nan(ncpar)] <- 1 ## placeholder
    qpar[is.nan(qpar)] <- 1 ## placeholder
    ##
    ## the chisquared cdf is not calculated well for ncp >= 1e5, there are lots of errors, 
    ## and the calculation also slows to a crawl. So for ncp >= 1e5, we will use the 
    ## normal approximation as per Muirhead, pages 22-24 and problem 1.8 with reference 
    ## Muirhead, R. (2005) Aspects of Multivariate Statistical Theory (2nd Edition). Wiley. ISBN 0-471-76985-1
    ##
    ## ifelse, however goes through the whole calculation, so we will replicate
    ## ncpar and set the values of ncpar > = 1e5 to be zero in 1 
    ## and ncpar can be used for the normal approximation.
    ##

    ##
    ## Actually, it appears that the issue also holds for ncpar >= 1e2.
    ## 

    nc <- 1e2
    
    ncpar1 <- ncpar
    ncpar1[ncpar >= nc] <- 0
        
    pv <- ifelse(ncpar < nc, pchisq(q = qpar, df = p, ncp = ncpar1, lower.tail = TRUE) * (sig2i > sig2j) + pchisq(q = qpar, df = p, ncp = ncpar1, lower.tail = FALSE) * (sig2i < sig2j),
                 ## this is the part now with the normal approximation
                 pnorm(q = qpar, mean = p + ncpar, sd = sqrt(2 * (p + 2 * ncpar)), lower.tail = TRUE) * (sig2i > sig2j) + pnorm(q = qpar, mean = p + ncpar, sd = sqrt(2 * (p + 2 * ncpar)), lower.tail = FALSE) * (sig2i < sig2j))              
    + pnorm(q = 0, mean = mu.sq.diff/sig2i, sd = 2 * sqrt(mu.sq.diff/sig2i)) * (sig2i == sig2j) * (mu.sq.diff != 0)
    pv
}		 


kmH.hc <- function(KL.res, L, K0, Kstar = NULL) { 

    Ko <- K0
    Mu <- KL.results$kmns.results[[Ko]]$centers
    W <- KL.results$kmns.results[[Ko]]$withinss
    N <- KL.results$kmns.results[[Ko]]$size
    cl <- KL.results$kmns.results[[Ko]]$cluster

    cl.mat <- rbind(cl)
    dstar <- NULL

    while (Ko >= max(c(2, Kstar))) {
        Ydist <- outer(1:Ko, 1:Ko, FUN = prob.dist, means = Mu, wss = W, size = N)
        pdistmat <- 1 - (Ydist + t(Ydist))/2
        ## this is the separation index where higher values indicate higher-separation and lower values indicate low separation.
        ##
        ## Note that the diagonals should be zero but we are keeping it 1 to simpllify the finding of the minimum. 
        
        dstar <- c(dstar, min(pdistmat))
        min.row.col <- which(pdistmat == min(pdistmat), arr.ind = TRUE)[1, ]
        
        ##
        ## find the row and column corresponding to the lowest value
        ##
        
        mu.updtd <- (N[min.row.col[2]] * Mu[min.row.col[2], ] + N[min.row.col[1]] * Mu[min.row.col[1], ]) / (N[min.row.col[2]] + N[min.row.col[1]])
        ## 
        ## update the Mu of the merged group (note that the row index of min.row.col is always larger than the second (column) index
        
        W[min.row.col[2]] <- W[min.row.col[2]] + N[min.row.col[2]] * sum((Mu[min.row.col[2], ] - mu.updtd))^2 + W[min.row.col[1]] + N[min.row.col[1]] * sum((Mu[min.row.col[1], ] - mu.updtd)^2) ## check

        N[min.row.col[2]] <- N[min.row.col[1]] + N[min.row.col[2]]
        
        Mu[min.row.col[2],] <- mu.updtd
        
        ## Now replace the vacated cluster ids with those from the largest group id (Ko).

        Mu[min.row.col[1],] <- Mu[Ko, ]
        W[min.row.col[1]] <- W[Ko]
        N[min.row.col[1]] <- N[Ko]
        
        ## Then truncate the vectors to be Ko-1

        Mu <- Mu[-Ko, ]
        length(W) <- Ko-1
        length(N) <- Ko-1
        
        ##
        ## merge the larger (row id) with smallest value
        ##
        
        cl[cl == min.row.col[2]] <- min.row.col[1]
        cl[cl == Ko] <- min.row.col[2]
        cl.mat <- rbind(cl.mat, cl)
        
        Ko <- Ko - 1
    }
    cp <- -diff(dstar)

    if (Kstar) 
    return(list(dstar = dstar, class = cl.mat, cp = cp))
}

kmH.hc.single <- function(KL.res, L, K0, Kstar = NULL, verbose = T) {    
    Mu <- KL.res$kmns.results[[K0]]$centers
    W <- KL.res$kmns.results[[K0]]$withinss
    N <- KL.res$kmns.results[[K0]]$size
    cl <- KL.res$kmns.results[[K0]]$cluster

    if (verbose)
        cat("K0 = ", K0, "Kstar = ", Kstar, "\n")
    
    if (is.null(Kstar)) {
        cl.mat <- NULL
        if (K0 == 2) {
            Ydist <- outer(1:K0, 1:K0, FUN = prob.dist, means = Mu, wss = W, size = N)
            pdistmat <- 1 - (Ydist + t(Ydist))/2
            dstar <- min(pdistmat)
            for (i in 1:L) 
            cl.mat <- rbind(cl.mat, cl)
            cp <- 1-dstar
        } else {
            Ydist <- outer(1:K0, 1:K0, FUN = prob.dist, means = Mu, wss = W, size = N)
            pdistmat <- 1 - (Ydist + t(Ydist))/2
            ## this is the separation index where higher values indicate higher-separation and lower values indicate low separation.
            ##
            diag(pdistmat) <- 0
            pd <- hclust(d = as.dist(pdistmat), method = "single")
            dstar <- pd$height
        cp <- diff(c(dstar, 1))
            which.cp <- order(diff(c(dstar, 1)), decreasing = T)[1:L]
            if (verbose) {
                print(dstar)
                print(cp)
                print(which.cp)
            }
            for (k in 1:min(K0-1,L)) {
                cl.tmp <- cl
                mrgs <- cutree(pd, h = dstar[which.cp[k]] - 0.001 * min(cp))
                for (j in 1:K0)
                    cl.tmp[cl == j] <- mrgs[j]
                cl.mat <- rbind(cl.mat, cl.tmp)
            }
            if (K0 == L) 
        	cl.mat <- rbind(cl.mat, rep(1, length(cl.tmp)))
        }
    } else {
        Ydist <- outer(1:K0, 1:K0, FUN = prob.dist, means = Mu, wss = W, size = N)
        pdistmat <- 1 - (Ydist + t(Ydist))/2
        ## this is the separation index where higher values indicate higher-separation and lower values indicate low separation.
        ##
        diag(pdistmat) <- 0
        pd <- hclust(d = as.dist(pdistmat), method = "single")
        mrgs <- cutree(pd, k = Kstar)
        if (verbose)
            print(mrgs)
        cl.tmp <- cl
        for (j in 1:K0)
            cl.tmp[cl == j] <- mrgs[j]
        cl.mat <- t(cl.tmp)
        cp <- NULL
        dstar <- NULL
    }
    return(list(cp = cp, cluster = cl.mat, dstar = dstar))
}

kmH.main <- function(M, L, KL.res, verbose = T, Kstar = NULL, 
                     minK0 = min(3*Kstar + 1, length(KL.res$KL) - M - 1)) {
    
    KL.res$KL[1:minK0] <- -Inf

    if (!is.null(Kstar))
        L <- 1

    K0 <- order(unlist(KL.res$KL), decreasing = T)[1:M]
    
    for (i in 1:M) {
        kmh <- kmH.hc.single(KL.res, L = L, K0 = K0[i], Kstar = Kstar, verbose = verbose)
	if (i == 1) {
            part.mat <- array(data = kmh$cluster, dim = c(dim(kmh$cluster), 1))
        }  else
            part.mat <- abind(part.mat, kmh$cluster)
    }

    W <- diag(L*M)
    

    
    ind.mat <- NULL
    m <- 1
    for (i in 1:L)
        for (j in 1:M) {
            n <- 1
            for (k in 1:L) 
                for (l in 1:M) {
                    if (m != n) 
                        W[m, n] <- RandIndex(id1 = part.mat[i, , j], id2 = part.mat[k, , l])$AR
                    n <- n + 1
                    ind.mat <- rbind(ind.mat, c(i, j))
                }
            m <- m + 1
        }

    Wi <- which.max(apply(X = W, MARGIN = 1, FUN = mean))
    if (verbose) {
        print(apply(X = W, MARGIN = 1, FUN = mean))
        print(Wi)
    }
    opt.part <- part.mat[ind.mat[Wi, 1], , ind.mat[Wi, 2]]
    return(list(opt.cluster = opt.part, all.partitions = part.mat))
}



##
## kmns.results: if provided, list with the following components:
##   kmns.results: a list with results of candidate kmeans competitors (assumed to be from 1 through maxclus) and ontain
##     KL: the Krzanaowski and Lai values for each K in the kmns.results.
##     x.red: the dataset minus any scatter that may have been identified. 
##     x.scat: the scatter points (not used in calculations)
##     scatter.indices: the indices for the observations identified as scatter. (not used in calculations)
##     id.scat: the id's of the observations identified as scatter
##     maxclus: the maximum number of clusters considered: same as length(kmns.results$KL)
## B: number of  samples (see paper)
## M: number of K_0s tried (see paper, default seems reasonable)
## L: number of K_*s tried (see paper, default 3)
##
## hmap.pdf.file: file to store the pdf of the heatmap (NULL, skipped by default)
## Value:
##  list object with the following components:
##  kmns.kmH: contains the kmeans results, the KL and the initial partition
##  kstar.distn: contains the replicated kstars
## final.partition: final partitioning

kmH <- function(x, kmns.results = NULL, nstart =  prod(dim(x)),B = 100, 
                M = min(round(0.1*sqrt(prod(dim(x)))), 10), L = 3, 
                verbose = F, maxclus = max(sqrt(nrow(x)),50), hmap.pdf.file = NULL, desired.ncores = 2,...) {

    if (is.null(kmns.results)) {
        KL.results <- kmns.KL(x = x, maxclus = maxclus, verbose = verbose)
        kmns.results <- KL.results
    }
    
    xx.kmH <- kmH.main(KL.res=kmns.results, minK0 = 1, M = M, L = L, verbose = verbose)
    zz <- list(kmeans.results = kmns.results, kmH.res = xx.kmH)


    
    a <- xx.kmH$all.partitions 
    n<-dim(a)[2]
    xx <- matrix(apply(a[,rep(1:n,n),]==a[,rep(1:n,each=n),],2,sum),nrow=n)/prod(dim(a)[-2])

    if (!is.null(hmap.pdf.file)) {
  hmap(kmH.soln = xx.kmH,...)
    }

    mu.xx <- 1 - mean(as.dist(1-xx))
    sd.xx <- sd(as.dist(1-xx))
    
    link.meth <- ifelse(sd.xx/mu.xx <= 0.975, ifelse((mu.xx < 0.26) | (mu.xx > 0.33), "single","complete"), "complete")

    if (verbose) 
        cat("mean and sd ",  mu.xx, sd.xx, "\n")
    
    desired.ncores <- min(detectCores(),desired.ncores)
    registerDoMC(desired.ncores)
    
    kstar <- NULL

    samp.size <- min(max(nrow(x)*0.25, 500), nrow(x))
    kstars <- foreach(i = 1:B) %dopar% { 
        if (nrow(zz$kmeans.results$x.red) <= samp.size) 
            a <- xx.kmH$all.partitions else 
                                           a <- xx.kmH$all.partitions[,sample(1:dim(xx.kmH$all.partitions)[2], size = samp.size),]
        n<-dim(a)[2]
        xx <- matrix(apply(a[,rep(1:n,n),]==a[,rep(1:n,each=n),],2,sum),nrow=n)/prod(dim(a)[-2])
        xx.hc <- hclust(as.dist(1-xx), method = link.meth)
        kstar <- max(cutree(xx.hc, h = 0.5))
        if (verbose)
            cat("Done i = ", i, kstar, "\n")
        kstar
    }
    if (verbose)
        print(table(unlist(kstars)))
    optk <- round(median(unlist(kstars)))

    fin.kmH <- kmH.main(KL.res=zz$kmeans.results, Kstar = max(2, optk), M = M, L = 1, verbose = verbose)$opt.cluster
    if (verbose)
        print(link.meth)
   
    list(kmns.kmH = zz, kstar.distn = unlist(kstars), final.partition = fin.kmH)
}

