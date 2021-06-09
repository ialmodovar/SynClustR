##************************************************************
#' Generalized Overlap
#' 
#' @description Compute the generalized overlap of Maitra (2010)
#' 
#' @param overlap.mat symmetric matrix of size K x K
#' 
#' @export
#' 
##***************************************************************

generalized.overlap <- function(overlap.mat) {
  p <- nrow(overlap.mat)
  if (is.null(p)) 1 else {
    x <- eigen(overlap.mat, symmetric = TRUE, only.values = TRUE)
    (x$values[1] - 1)/(p - 1)
  }
}
