\name{kcdf}
\alias{kcdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel estimation of the cumulative distribution function
}
\description{
Perform kernel estimation of the cumulative distribution function.
}
\usage{
 kcdf(x, b = NULL, kernel = "RIG", xgrid = NULL, m = 100, inv.roles = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{A vector (univariate). Default is use of asymmetric kernel defined in the positive real line. These kernel are Reciprocal Inverse Gaussian (RIG), Gamma and Gamma with inverse roles. If the vector contain at least a negative number will use Gaussian kernel.}
\item{b}{Smoothing parameter (Bandwidth). Default is NULL and will estimate the bandwith using a rule-of-thumb that minimize the mean integrated squared error.}
\item{kernel}{Choice of the kernel to be use in the smooth estimation. Default is Reciprocal Inverse Gaussian (RIG), other choices are Gamma kernel, Gamma kernel with inverse roles and Gaussian kernel.} 
\item{xgrid}{Range of grid values that will be evaluated.}
\item{m}{Maximum number of grid.} 
\item{inv.roles}{Only for Gamma kernel. Default is FALSE, which mean kde will be perform using Chen (2000), if TRUE will use Jeon and Kim (2013) with inverse roles.}
}
\details{
Default is use of asymmetric kernel defined in the positive real line. These kernel are Reciprocal Inverse Gaussian (RIG), Gamma and Gamma with inverse roles. If the vector contain at least a negative number will use Gaussian kernel
}
\references{
Chen, S. X. (2000). Probability density function estimation using gamma kernels. Annals of the Institute of Statistical Mathematics, 52(3), 471-480.

Jeon, Y., & Kim, J. H. (2013). A gamma kernel density estimation for insurance loss data. Insurance: Mathematics and Economics, 53(3), 569-579.

Scaillet, O. (2004). Density estimation using inverse and reciprocal inverse Gaussian kernels. Nonparametric statistics, 16(1-2), 217-226.
}
\author{
Israel Almodóvar-Rivera.
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
## gamma(5/2,2) = chisquared with 5 degrees
x <- rgamma(n = 500,shape= 5/2,scale = 2)
xg <- sort(x)
## reciprocal inverse gaussian
F.rig <- kcdf(x, kernel = "RIG",xgrid=xg)
## Gamma kernel of Chen (2000)
F.gamma1 <- kcdf(x, kernel = "gamma",xgrid=xg)
## Gamma kernel of Jeon and Kim (2013)
F.gamma2 <- kcdf(x, kernel = "gamma",inv.roles=TRUE,xgrid=xg)
## Gaussian kernel
F.gaussian <- kcdf(x,kernel="gaussian",xgrid=xg)
## True cdf
F.true <- pgamma(sort(x), shape = 5/2,scale= 2)
## sort data only for display purposes
plot(xg,F.true,col="black",type="l",lwd=2)
lines(xg,F.rig$Fhat,col="red",lwd=2)
lines(xg,F.gamma1$Fhat,col="blue",lwd = 2)
lines(xg,F.gamma2$Fhat,col="orange",lwd = 2)
lines(xg,F.gaussian$Fhat,col="green",lwd=2)
legend("bottomright",col=c("black","red","blue","orange","green"),lty=rep(1,5),lwd=rep(2,5),
legend=c("True","RIG kernel","Gamma kernel","Gamma kernel with inverse roles","Gaussian kernel"))
}
}

