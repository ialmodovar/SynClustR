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
#' @export
#' 
overlapKDE <- function(X,Means, ids, b = NULL,kernel = "RIG",inv.roles = FALSE,desired.ncores=2)
{
  
  n <- nrow(X)
  p <- ncol(X)
  K <- max(ids)
  
  kernel <- match.arg(kernel,choices = c("gamma","RIG","gaussian"))
  
  residuals.norm <- norm.res(X = X, Means = Means, ids = ids,desired.ncores = desired.ncores)
  psi <- residuals.norm$Psi
  psi.all <- residuals.norm$PsiAll
  psi.pseudo <- residuals.norm$PseudoPsi
  
  if(is.null(b)){
    wss <- psi * psi
    b <- gsl_bw_mise_cdf((psi*psi/(wss/(p*(n-K)))))   
  }
  Fhat.Psi <- t(apply(psi.all,1,function(z){kcdf(x = psi,b = b,kernel=kernel,xgrid=z,inv.roles=inv.roles)$Fhat}))
  
  Omega.lk <- array(0,dim=c(K,K))
  for(k in 1:K){
    for(l in k:K){
      Omega.lk[k,l] <- Omega.lk[l,k] <- 1-mean(Fhat.Psi[ids==k,l])+1-mean(Fhat.Psi[ids==l,k])
    }
  }
  diag(Omega.lk) <- 1
  list(Omega = Omega.lk,
       MaxOverlapKDE = max(Omega.lk[!lower.tri(Omega.lk,diag = TRUE)]),
       MeanOverlapKDE = mean(Omega.lk[!lower.tri(Omega.lk,diag = TRUE)]),
       GenOverlapKDE = generalized.overlap(Omega.lk))
}

