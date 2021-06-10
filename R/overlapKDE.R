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
#' @return The following is return
#' 
#' \begin{itemize}
#' \item Omega- Overlap matrix of size K x K.
#' \item OmegaConditional- Overlap matrix of missclassification probabilities of 
#' \item MaxOverlap- Maximum overlap of Omega
#' \item MeanOverlap- Average overlap of Omega
#' \item GenOverlap- Generalized overlap of Omega. See Maitra (2010)
#' \end{itemize}
#' 
#' @references
#' 
#' Almodóvar-Rivera, I., & Maitra, R. (2020). Kernel-estimated Nonparametric Overlap-Based Syncytial Clustering. J. Mach. Learn. Res., 21, 122-1.
#'  
#' @examples 
#' 
#' data("iris", package = "datasets")
#' id <- as.integer(iris[, 5])
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
  
  desired.ncores <- min(detectCores(),desired.ncores)
  residuals.norm <- norm.res(X = X, Means = Means, 
                             ids = ids,desired.ncores = desired.ncores)
  psi <- residuals.norm$Psi
  psi.pseudo <- residuals.norm$PseudoPsi
  
  if(is.null(b)){
    wss <- sum(psi * psi)
    b <- gsl_bw_mise_cdf((psi*psi/(wss/(p*(n-K)))))   
  }
  Fhat.Psi <- t(apply(psi.pseudo,1,function(z){kcdf(x = psi,b = b,kernel=kernel,xgrid=z,inv.roles=inv.roles)$Fhat}))
  
  Omega.lk <- array(0,dim=c(K,K))
  ## Compute \hat Omega _{l|k}
  for(k in 1:K){
    for(l in k:K){
      Omega.lk[l,k] <-  1-mean(Fhat.Psi[ids==k,l])
    }
  }
  Omega <- (t(Omega.lk)+Omega.lk) ## compute overlap matrix 
diag(Omega.lk) <- diag(Omega) <- 1
  list(Omega = Omega,
       OmegaConditional = Omega.lk,
       MaxOverlap = max(Omega[!lower.tri(Omega,diag = TRUE)]),
       MeanOverlap = mean(Omega[!lower.tri(Omega,diag = TRUE)]),
       GenOverlap = generalized.overlap(Omega))
}



