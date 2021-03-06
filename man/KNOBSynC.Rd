\name{KNOBSynC}
\alias{KNOBSynC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel-estimated Nonparametric Overlap-Based Syncytial Clustering
}
\description{
Merge clustering components using smooth estimation of the missclassification probabilities
}
\usage{
KNOBSynC(x, kmns.results=NULL, Kmax = NULL,EstK=NULL,
kernel = "RIG", kappa = NULL, b = NULL,desired.ncores=2,
trueids=NULL,inv.roles=FALSE,
ret.steps = FALSE,verbose=FALSE,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector or matrix}
  \item{kmns.results}{Default is NULL. If k-means results is provided will  use it, if not will run k-means until Kmax and estimate the number  of groups.}
  \item{Kmax}{Default is NULL.If not provided will choose Kmax =
    max{50,sqrt(n)}}
  \item{EstK}{If need to estimate the number of groups will use the Jump method of Sugar and James (2003). If EstK = "KL" will estimate the number of groups using Krznaswoki and Lai (1985).}
  \item{kernel}{Choice of the kernel to be use in the smooth estimation. Default is Reciprocal Inverse Gaussian (RIG), other choices are Gamma kernel, Gamma kernel with inverse roles and Gaussian kernel.}
  \item{kappa}{Merging parameter. This defined as kappa * gen.overlap. See Almodovar-Rivera and Maitra (2018).}
  \item{b}{Smoothing parameter (bandwith) to be use in the smooth estimation of
    the distribution function. If not value is provided will use a
    estimate of the MISE}
  \item{desired.ncores}{Set the desired number of processor cores to be
    use. Default is min(desired.ncores, maximum number of cores).}
   \item{trueids}{If provided will compute the adjusted Rand Index at each step}
   \item{inv.roles}{Inverse role when using a gamma kernel estimator,
     Jeon and Kim (2013). Default is FALSE using Chen (2000).}
   \item{ret.steps}{Return steps calculation for each iteration. Default
   is FALSE}
 \item{verbose}{FALSE}{Show progress in the console.}
 \item{...}{Arguments to be use if k-means.all() need to be call.}
}

\references{
Almodovar-Rivera, I., & Maitra, R. (2018). Kernel-estimated Nonparametric Overlap-Based Syncytial Clustering. arXiv preprint arXiv:1805.09505.

}
\author{
Israel Almodovar-Rivera and Ranjan Maitra.
}

 
\examples{
\donttest{
set.seed(787)
## Example 1: Merging will occur. Best clustering solution when generalized overlap = 
data(Bullseye)
oo <- KNOBSynC(x = Bullseye[,-3],verbose = TRUE)
Bullseye$IdsKmeans <- oo$Ids
Bullseye$IdsKNOBSynC <- oo$IdsMerge

par(mfrow=c(1,3))
with(Bullseye,plot(x = x,y = y, col=Ids,main="True"))
with(Bullseye,plot(x = x,y = y, col=IdsKmeans,main="k-means"))
with(Bullseye,plot(x = x,y = y, col=IdsKNOBSynC,main="KNOB-SynC"))

## Example 2 Merging will not occur since generalized-overlap approx max-overlap 
data(Spherical7)
oo <- KNOBSynC(x = Spherical7[,-3],verbose = TRUE)
Spherical7$IdsKmeans <- oo$Ids
Spherical7$IdsKNOBSynC <- oo$IdsMerge

par(mfrow=c(1,3))
with(Spherical7,plot(x = x,y = y, col=Ids,main="True"))
with(Spherical7,plot(x = x,y = y, col=IdsKmeans,main="k-means"))
with(Spherical7,plot(x = x,y = y, col=IdsKNOBSynC,main="KNOB-SynC"))

## Example 3 a more difficult data  take some time to finish

data(SSS)
oo <- KNOBSynC(x = SSS[,-3],verbose = TRUE,nstart = 100)
SSS$IdsKmeans <- oo$Ids
SSS$IdsKNOBSynC <- oo$IdsMerge

par(mfrow=c(1,3))
with(SSS,plot(x = x,y = y, col=Ids,main="True"))
with(SSS,plot(x = x,y = y, col=IdsKmeans,main="k-means"))
with(SSS,plot(x = x,y = y, col=IdsKNOBSynC,main="KNOB-SynC"))
}

}

