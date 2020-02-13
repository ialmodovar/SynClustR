##******************************************************
##
## @file: run-kmeans.R
## 
## Find homogenoues groups using a parallelized version of k-means
## Estimate the number of groups using 
## Jump method of Sugar and James (2003) 
## Krzanwoski and Lai (1985) 
##
## Author: 
## Ranjan Maitra
## Iowa State University
## Department of Statistics 
## Ames, IA 50010
## email: maitra@iastate.edu
##
## Israel Almodovar-Rivera, PhD                       
## University of Puerto Rico                          
## Medical Science Campus                             
## Graduate School of Public Health                   
## Department of Biostatistics and Epidemiology       
## San Juan, PR, USA 00953
## email: israel.almodovar@upr.edu
##******************************************************

kmns.all.soln <- function(k, x, x.hc, nstarts, h.c = TRUE) {
  
  if (k == 1){
    kmeans(x, centers = 1)
  } else {
    if (h.c) {
      cl <- cutree(x.hc, k = k)
      x.df <- data.frame(x = x, cl = as.factor(cl))
      km.center <- as.matrix(aggregate(formula = . ~ cl, FUN  = mean, data = x.df)[,-1])
      km1 <- kmeans(x, centers = km.center, iter.max = 100)
      km2 <- kmeans(x, centers = k, iter.max = 100, nstart =  nstarts * k)
      if (km1$tot.withins < km2$tot.withinss)
        km1 else
          km2
    } else {
      kmeans(x, centers = k, iter.max = 100, nstart =  nstarts * k)
    }
  }
}



kmeans.all <- function(x, maxclus, nstarts = prod(dim(x)),desired.ncores = 2, h.c = TRUE) {
    ## x = dataset (matrix of observations: each row is an observation vector)
    ## maxclus = the maximum number of clusters.
    ## h.c = if hierarchical clustering should also be done
  

    distortion <-  function(k, res, x) (res[[k]]$tot.withinss/prod(dim(x)))

    jump <- function(k, dstrtn, x, Y = ncol(x)/2) (ifelse(k == 1, dstrtn[[1]]^(-Y), dstrtn[[k]]^(-Y) - dstrtn[[k-1]]^(-Y)))

    sc.trW <-  function(k, res, x) (res[[k]]$tot.withinss * k^(2/dim(x)[2]))

    dif <- function(k, res, x) (sc.trW(k = k - 1, res = res, x = x) - sc.trW(k = k, res = res, x = x))

    KL <- function(k, res, x) (ifelse(k == 1, NA, abs(dif(k = k, res = res, x = x)/dif(k = k + 1, res = res, x = x))))
    
    x.hc <- hclust(dist(x), method = "ward.D2")
    
    desired.ncores <- min(detectCores(),desired.ncores)
    
    cl <- makeCluster(desired.ncores)
    clusterExport(cl, list("kmns.all.soln"))
    kmns.results <- parLapply(cl,1:maxclus, function(i){
      kmns.all.soln(k = i,x=x, x.hc = x.hc, nstarts = nstarts, h.c =h.c)
    })
    stopCluster(cl)
    

    kmns.dstrtn <- sapply(X = 1:maxclus, FUN = distortion, res = kmns.results, x = x)
    jump.stat <- sapply(1:maxclus, FUN = jump, dstrtn = kmns.dstrtn, x = x)
    kl.stat <- sapply(1:(maxclus-1),FUN = KL, res = kmns.results,x = x)
    
    list(kmns.results = kmns.results, distortions = kmns.dstrtn, jump.stat = jump.stat,
         kl.stat = kl.stat, Khat.jump = which.max(jump.stat),Khat.KL = which.max(kl.stat))

}

##

