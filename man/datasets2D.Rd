\name{Dataset}
\alias{7-Spherical}
\alias{Aggregation}
\alias{Bananas-Arcs}
\alias{Banana-Clump}
\alias{Bullseye}
\alias{Bullseye-Cigarette}
\alias{Compound}
\alias{Half-ringed-clusters}
\alias{Pathbased}
\alias{SCX-Bananas}
\alias{Spiral}
\alias{SSS}
\alias{XXXX}
\docType{data}
\title{2D Toy datasets}
\description{
  Toy datasets to test and demonstrate \pkg{SynClustR}.
}

\format{

    \code{Spherical7} is a 2D dimensional dataset of size 500 x 3. The third
    column contains the true memberships.
    
  \code{Aggregation}, 2D dimensional data set of size 788 x 3. The third
  column contains the true memberships.

  \code{Bananas-Arcs}, 2D dimensional data set of size 4515 x 3. The third
  column contains the true memberships.

  \code{Bananas-Sphere}, 2D dimensional data set of size 3015 x 3. The third
  column contains the true memberships.

  \code{Bananas-Clump}, 2D dimensional data set of size 200 x 3. The third
  column contains the true memberships.
  
  \code{Bullseye} is a 2D dimensional dataset of size 400 x 3. The third
  column contains the true memberships.

  \code{Bullseye-Cigarette} is a 2D dimensional dataset of size 3025 x 3. The third
  column contains the true memberships.

  \code{Compound} is a 2D dimensional dataset of size 788 x 3. The third
  column contains the true memberships.

  \code{D31} is a 2D dimensional dataset of size 3100 x 3. The third
  column contains the true memberships.

  \code{Flame} is a 2D dimensional dataset of size 240 x 3. The third
  column contains the true memberships.

  \code{Half.ringed.clusters} is a 2D dimensional dataset of size 373 x 3. The third
  column contains the true memberships.

  \code{Pathbased} is a 2D dimensional dataset of size 300 x 3. The third
  column contains the true memberships.

  \code{R15} is a 2D dimensional dataset of size 600 x 3. The third
    column contains the true memberships.
    
  \code{SCX-Bananas} is a 2D dimensional dataset of size 3420 x 3. The third
  column contains the true memberships.
  
  \code{Spiral} is a 2D dimensional dataset of size 312 x 3. The third
  column contains the true memberships.
  
  \code{SSS} is a 2D dimensional dataset of size 5015 x 3. The third
  column contains the true memberships.

    \code{XXXX} is a 2D dimensional dataset of size 415 x 3. The third
  column contains the true memberships.
}
\details{
    \code{Spherical7} is a 2D dimensional dataset of size 500 x 3. The third
    column contains the true memberships.
    
  \code{Aggregation}, 2D dimensional data set of size 788 x 3. The third
  column contains the true memberships.

  \code{Bananas-Arcs}, 2D dimensional data set of size 4515 x 3. The third
  column contains the true memberships.

  \code{Bananas-Sphere}, 2D dimensional data set of size 3015 x 3. The third
  column contains the true memberships.

  \code{Bananas-Clump}, 2D dimensional data set of size 200 x 3. The third
  column contains the true memberships.
  
  \code{Bullseye} is a 2D dimensional dataset of size 400 x 3. The third
  column contains the true memberships.

  \code{Bullseye-Cigarette} is a 2D dimensional dataset of size 3025 x 3. The third
  column contains the true memberships.

  \code{Compound} is a 2D dimensional dataset of size 788 x 3. The third
  column contains the true memberships.

  \code{D31} is a 2D dimensional dataset of size 3100 x 3. The third
  column contains the true memberships.

  \code{Flame} is a 2D dimensional dataset of size 240 x 3. The third
  column contains the true memberships.

  \code{Half.ringed.clusters} is a 2D dimensional dataset of size 373 x 3. The third
  column contains the true memberships.

  \code{Pathbased} is a 2D dimensional dataset of size 300 x 3. The third
  column contains the true memberships.

  \code{R15} is a 2D dimensional dataset of size 600 x 3. The third
    column contains the true memberships.
    
  \code{SCX-Bananas} is a 2D dimensional dataset of size 3420 x 3. The third
  column contains the true memberships.
  
  \code{Spiral} is a 2D dimensional dataset of size 312 x 3. The third
  column contains the true memberships.
  
  \code{SSS} is a 2D dimensional dataset of size 5015 x 3. The third
  column contains the true memberships.

    \code{XXXX} is a 2D dimensional dataset of size 415 x 3. The third
  column contains the true memberships.


}

\author{
Israel A. Almodovar-Rivera and Ranjan Maitra
}
\references{
Almodóvar-Rivera, I., & Maitra, R. (2018). Kernel-estimated
Nonparametric Overlap-Based Syncytial Clustering. arXiv preprint
arXiv:1805.09505.

H. Chang and D.-Y. Yeung (2008). “Robust path-based spectral clustering.” Pattern Recognition, 41(1):191 –
203. ISSN 0031-3203.

L. Fu and E. Medico, FLAME, a novel fuzzy clustering method for the analysis of DNA microarray data. BMC bioinformatics, 2007. 8(1): p. 3. 

A. Gionis, H. Mannila, and P. Tsaparas, Clustering aggregation. ACM Transactions on Knowledge Discovery from Data (TKDD), 2007. 1(1): p. 1-30.

A. Jain and M. Law, Data clustering: A user's dilemma. Lecture Notes in Computer Science, 2005. 3776: p. 1-10. 

Peterson, A. D., Ghosh, A. P., & Maitra, R. (2018). Merging K‐means with
hierarchical clustering for identifying general‐shaped groups. Stat, 7(1), e172.

W. Stuetzle and R. Nugent, “A generalized single linkage method for estimating the cluster tree of a density,” Journal of Computational
and Graphical Statistics, 2010.

C. J. Veenman, M. J. T. Reinders, and E. Backer (2002). “A maximum
variance cluster algorithm.” IEEE Transactions on Pattern Analysis and
Machine Intelligence, 24 (9):1273–1280.

C.T. Zahn, Graph-theoretical methods for detecting and describing gestalt clusters. IEEE Transactions on Computers, 1971. 100(1): p. 68-86. 

}


\keyword{datasets}

\examples{

##
par(mfrow=c(3,4),mar=rep(1.75,4))
data("Spherical7")
with(Spherical7,plot(x,y,col=Ids,main="7-Spherical",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("Bananas-Arcs")
with(Bananas.Arcs,plot(x,y,col=Ids,main="Bananas-Arcs",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("Banana-Clump")
with(Banana.Clump,plot(x,y,col=Ids,main="Banana-Clump",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("Bullseye")
with(Bullseye,plot(x,y,col=Ids,main="Bullseye",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("Bullseye-Cigarette")
with(Bullseye.Cigarette,plot(x,y,col=Ids,main="Bullseye-Cigarette",xlab = "",ylab = "",
xlim=c( 0.6575804,11.7114793),axes = FALSE,frame.plot = TRUE))

data("Compound")
with(Compound,plot(x,y,col=Ids,main="Compound",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))


data("Half-ringed-clusters")
with(Half.ringed.clusters,plot(x,y,col=Ids,main="Half-ringed-clusters",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("Path-based")
with(Pathbased,plot(x,y,col=Ids,main="Path-based",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("SCX-Bananas")
with(SCX.Bananas,plot(x,y,col=Ids,main="SCX-Bananas",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("Spiral")
with(Spiral,plot(x,y,col=Ids,main="Spiral",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("SSS")
with(SSS,plot(x,y,col=Ids,main="SSS",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))

data("XXXX")
with(XXXX,plot(x,y,col=Ids,main="XXXX",
xlab = "",ylab = "",axes = FALSE,frame.plot = TRUE))
}