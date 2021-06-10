#' Overlap Matrix KDE
#' 
#' @description Compute pairwise overlap matrix using smooth estimation of the distribution function
#' 
#' @param X dataset of size n x p
#' @param Means clusters centroids (means)
#' @param ids cluster memberships
#' @param b smoothing parameter
#' @param kernel choice of the kernel
#' @param inv.roles inverse roles of the gamma kernel
#' @param desired.ncores number of desired cores
#' 
#' @examples 
#' 
#' data("iris", package = "datasets")
#' id <- as.integer(iris[, 5])
#' #estimate mixture parameters
#' Mu <- t(sapply(1:K, function(k){ colMeans(iris[id == k, -5]) }))
#' overlapKDE(X = iris[,-5], Means = Mu,ids = id)
#' 
#' @export
#' 
overlapKDE <- function(X,Means, ids, b = NULL,kernel = "RIG",inv.roles = FALSE,desired.ncores=2)
{
  
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  K <- max(ids)
  
  kernel <- match.arg(kernel,choices = c("gamma","RIG","gaussian"))
  
  residuals.norm <- norm.res(X = X, Means = Means, ids = ids,desired.ncores = desired.ncores)
  psi <- residuals.norm$Psi
  psi.pseudo <- residuals.norm$PseudoPsi
  
  if(is.null(b)){
    wss <- sum(psi * psi)
    b <- gsl_bw_mise_cdf((psi*psi/(wss/(p*(n-K)))))   
  }
  Fhat.Psi <- t(apply(psi.pseudo,1,function(z){kcdf(x = psi,b = b,kernel=kernel,xgrid=z,inv.roles=inv.roles)$Fhat}))
  
  Omega.lk <- array(0,dim=c(K,K))
  ## Return \hat Omega _{l|k}
  for(k in 1:K){
    for(l in k:K){
      Omega.lk[l,k] <-  1-mean(Fhat.Psi[ids==k,l])
    }
  }
  diag(Omega.lk) <- 1
  Omega <- (t(Omega.lk)+Omega.lk)/2
  list(Omega = Omega,
       OmegaConditional = Omega.lk,
       MaxOverlapKDE = max(Omega[!lower.tri(Omega,diag = TRUE)]),
       MeanOverlapKDE = mean(Omega[!lower.tri(Omega,diag = TRUE)]),
       GenOverlapKDE = generalized.overlap(Omega))
}



