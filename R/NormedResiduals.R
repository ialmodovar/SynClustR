
##*************************************************************************************************************************************
#' normed.residuals
#' 
#' @title normed.residuals
#' 
#' @description Compute Normed and Pseudo Normed Residuals
#' 
#' @param X dataset of size n x p
#' @param Means Centroids 
#' @param ids Cluster membership
#' @param desired.ncores Desired number of cores to be used. Default is 2, however the function will determine min(detectCores(),desired.ncores)
#' @return A list with the following
#' \begin{itemize}
#' \item Psi- Compute actual normed residuals for each observation, || X_i -mu_k ||
#' \item PseudoPsi - Return a matrix of size n x K, where K is the maximum number of groups that contains the distance of each observation with the corresponding mean
#' \item Eps- Return a matrix of size n x p that contains the residuals, i.e., e_i = X_i-mu_k
#' \end{itemize}
#' @examples 
#' @export
##***********************************************************************************************************************************


norm.res <- function (X, Means, ids, desired.ncores = max(availableCores(),2)) 
{
  if ((ncol(X) != ncol(Means)) | (length(ids) != nrow(X)) | 
      (nrow(Means) != max(ids))) {
    stop("Sizes do not match \n")
  }
  X <- as.matrix(X)
  Means <- as.matrix(Means)
  n <- nrow(X)
  p <- ncol(X)
  K <- max(ids)
  ## compute residuals e = X_i - \mu_k
  Eps <- lapply(1:n, function(z) {
    find.nan <- which(!is.na(X[z, ]))
    X[z, find.nan] - Means[ids[z], find.nan]
  })
  ## compute normed residuals psi = || X_i - \mu_k||
  psi <- sapply(Eps, norm, type = "2")
  names(psi) <- NULL
  ## adjusted weight for the presence of missing dimensions
  wgt <- apply(X, 1, function(z) sum(!is.na(z)))/p
  psi <- psi * wgt
  if (p > 1) {
    Diff.Psi <- t(sapply(1:n, function(i) {
      find.nan <- which(!is.na(X[i, ]))
      if(length(find.nan) > 1){
        apply(t(sapply(1:K, function(k) {
          X[i, find.nan] - Means[k, find.nan]
        })), 1, norm, type = "2") * wgt[i]
      } else{
        apply(t(sapply(1:K, function(k) {
          X[i, find.nan] - Means[k, find.nan]
        })),2,norm,type="2") *wgt[i]
      }
    }))
  } else {
    Diff.Psi <- sqrt(t(sapply(1:n,function(i){
      sapply(1:K,function(k){
        eps <- X[i,]-Means[k,]
        eps*eps
      })
    })))
  }
  nks <- as.numeric(table(ids))
  ## compute pseudo normed residuals as defined 
  ## in Almodovar-Rivera and Maitra (2020)
  cl <- makeCluster(desired.ncores)
  clusterExport(cl, list("norm", "sum"))
  pseudo.psi <- t(parSapply(cl, 1:n, FUN = function(i) {
    k <- ids[i]
    Means2 <- Means
    nk <- nks[k]
    find.nan <- which(!is.na(X[i, ]))
    ## if dimension avaible  is one
    Means2[k, find.nan] <- (nk * Means[k, find.nan] - X[i, find.nan])/(nk - 1)
    sapply(1:K, function(l) {
      nl <- nks[l]
      Means2[l, find.nan] <- (nl * Means[l, find.nan] + 
                                X[i, find.nan])/(nl + 1)
      Psi.rev <- norm(X[i, find.nan] - Means2[l, find.nan], 
                      type = "2") * wgt[i]
      Psi.rev
    })
  }))
  stopCluster(cl)
  list(Psi = psi, E = Eps, PseudoPsi = Diff.Psi)
}
