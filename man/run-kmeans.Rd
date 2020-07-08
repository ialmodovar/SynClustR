\name{kmeans.all}
\alias{kmeans.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run parallezied k-means algorithm
}
\description{
Run k-means for a range of cluster and estimate the number of groups using the Jump Method and Krzanowski and Lai.
}
\usage{
kmeans.all(x, maxclus, nstarts,desired.ncores = 2, h.c = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).}
\item{maxclus}{Maximum number of clusters in the given range. It will estimate the number of groups from [1,maxclus].}
\item{nstarts}{How many random sets should be chosen?}
\item{desired.ncores}{Number of cores set to be used by the
  user. Default is 2.}
\item{h.c}{If use hierarchical clustering as a initialization}
}
\details{
The data given by x are clustered by the k-means method, which aims to partition the points into k groups such that the sum of squares from points to the assigned cluster centres is minimized. At the minimum, all cluster centres are at the mean of their Voronoi sets (the set of data points which are nearest to the cluster centre).
}
\references{
Forgy, E. W. (1965). Cluster analysis of multivariate data: efficiency vs interpretability of classifications. Biometrics, 21, 768–769.

Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-means clustering algorithm. Applied Statistics, 28, 100–108. doi: 10.2307/2346830.

Krzanowski, W. J., & Lai, Y. T. (1988). A criterion for determining the number of groups in a data set using sum-of-squares clustering. Biometrics, 23-34.

Lloyd, S. P. (1957, 1982). Least squares quantization in PCM. Technical Note, Bell Laboratories. Published in 1982 in IEEE Transactions on Information Theory, 28, 128–137.

MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, eds L. M. Le Cam & J. Neyman, 1, pp. 281–297. Berkeley, CA: University of California Press.

Sugar, C. A., & James, G. M. (2003). Finding the number of clusters in a dataset: An information-theoretic approach. Journal of the American Statistical Association, 98(463), 750-763.
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
\donttest{
set.seed(787)
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
## Run kmeans until 10 groups with 200 initializations
kk <- kmeans.all(x = x,maxclus = 10, nstart = 200)
## Jump approach to estimate number of groups
jump <- unlist(kk$jump.stat)
Khat.jump <- which.max(jump)

## Krzanowski and Lai approach to estimate number of groups
kl <- unlist(kk$kl.stat)
Khat.kl <- which.max(kl)

cbind(Khat.jump, Khat.kl)
}
}
