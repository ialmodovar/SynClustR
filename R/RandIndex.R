##**********************************************************************************
##
## @file:RandIndex.R
## 
## Call the C code for adjusted Rand
## Index. Code belong to MixSim.
##
##
##  Volodymyr Melnykov, Wei-Chen Chen, Ranjan Maitra (2012). MixSim: An R Package for
## Simulating Data to Study Performance of Clustering Algorithms. Journal of
## Statistical Software, 51(12), 1-25. URL http://www.jstatsoft.org/v51/i12/. 
##**********************************************************************************

RandIndex <- function(id1, id2){

	if (length(id1) != length(id2)) stop("Lengths of partitioning vectors do not match...\n")

	n <- length(id1)
	

	A <- as.factor(id1)
	B <- as.factor(id2)
	K1 <- nlevels(A)
	K2 <- nlevels(B)

	for (i in 1:nlevels(A)){
		ind <- A == levels(A)[i]
		id1[ind] <- i - 1
	}
	for (i in 1:nlevels(B)){
		ind <- B == levels(B)[i]
		id2[ind] <- i - 1
	}


	Q <- .C("runAdjRand", n = as.integer(n), K1 = as.integer(K1), K2 = as.integer(K2), id1 = as.integer(id1), id2 = as.integer(id2), Rand = as.double(0), aRand = as.double(0), F = as.double(0),PACKAGE="SynClustR")

       	list(R = Q$Rand, AR = Q$aRand, F = Q$F, M = n * (n - 1) * (1 - Q$Rand))
}
