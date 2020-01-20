\name{KmH}
\alias{KmH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merging K-means with hierarchical clustering for identifying general-shaped groups
}
\description{
A hybrid non-parametric clustering approach that amalgamates two methods
(K-means and Hierarchical) to identify general-shaped clusters and that can apply to larger dataset
}
\usage{
kmH <- function(x, kmns.results = NULL, nstart =  prod(dim(x)),B = 100, 
                M = min(round(0.1*sqrt(prod(dim(x)))), 10), L = 3, 
                verbose = FALSE, maxclus = max(sqrt(nrow(x)),50), hmap.pdf.file = NULL, desired.ncores = 2,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector or matrix}
  \item{kmns.results}{Default is NULL. If k-means results is provided will use it, if not will run k-means until maxclus and estimate the
    number of groups.}
  \item{nstart}{Number of initialization to be use in k-means}
  \item{B}{number of  samples (see paper)}
\item{M}{number of K_0s tried (see paper, default seems reasonable)}
\item{L}{number of K_*s tried (see paper, default 3)}
\item{verbose}{Show progress. Default is FALSE.}
\item{maxclus}{Default is NULL.If not provided will choose Kmax =  max{50,sqrt(n)}}
\item{hmap.pdf.file}{file to store the pdf of the heatmap (NULL, skipped by default)}
\item{desired.ncores}{Number of desired cores. By default is min(detectCores(),desired.ncores=2)}
}
\details{
}
\references{
Almodóvar-Rivera, I., & Maitra, R. (2018). Kernel-estimated Nonparametric Overlap-Based Syncytial Clustering. arXiv preprint arXiv:1805.09505.
Peterson, A. D., Ghosh, A. P., & Maitra, R. (2018). Merging K‐means with
hierarchical clustering for identifying general‐shaped groups. Stat, 7(1), e172.

}
\author{
Israel Almodóvar-Rivera and Ranjan Maitra.
}


\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
 
\examples{
\dontrun{
set.seed(787)
data("Bullseye")
kk <- kmH(x = Bullseye[-3],maxclus = sqrt(nrow(Bullseye)),verbose = TRUE)
Khat <- which.max(unlist(kk$kmns.kmH$kmeans.results$KL))
## extract k-means solution
km <- kk$kmns.kmH$kmeans.results$kmns.results[[Khat]]
Bullseye$IdsKmeans <- km$cluster
Bullseye$IdsKmH<- kk$final.partition

par(mfrow=c(1,3))
with(Bullseye,plot(x = x,y = y, col=Ids,main="True"))
with(Bullseye,plot(x = x,y = y, col=IdsKmeans,main="k-means"))
with(Bullseye,plot(x = x,y = y, col=IdsKmH,main="K-mH"))
}

}