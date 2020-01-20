\name{RandIndex}
\alias{RandIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rand's Index}
\description{
Computes Rand, adjusted Rand, Fowlkes and Mallows, and Merkin indices.
}
\usage{
RandIndex(id1, id2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id1}{ first partitioning vector. }
  \item{id2}{ second partitioning vector. }
}

\value{
  \item{R}{Rand's index.}
  \item{AR}{adjusted Rand's index.}
  \item{F}{Fowlkes and Mallows index.}
  \item{M}{Mirkin metric.}
}
\references{
Rand, W.M. (1971) ``Objective criteria for the evaluation of clustering methods'', Journal of the American Statistical Association, 66:336, 846-850.

Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Meila, M. (2006) ``Comparing clusterings - an information based distance'', Journal of Multivariate Analysis, 98, 873-895.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }
\examples{
%\dontrun{
id1 <- c(rep(1, 50), rep(2,100))
id2 <- rep(1:3, each = 50)
RandIndex(id1, id2)
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
