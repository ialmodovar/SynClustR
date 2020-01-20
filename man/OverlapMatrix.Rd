\name{overlap.map}
\alias{overlap.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display Pairwise Overlap Measures
}
\description{
Display the map of pairwise overlap measures of Maitra and Melnykov (JCGS, 2012)
}
\usage{
overlap.map(overlap.mat, map.col = c("#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB",
"#41B6C4","#1D91C0","#225EA8","#253494","#081D58"), linescol = "magenta", 
             map.range = NULL, lab.col = 1, lab.cex = 0.95, map.cex = 0.85, legend.cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{overlap.mat}{matrix of total pairwise overlaps (necessarily symmetric)}
\item{map.col}{colormap for the mapping. Default is Green to Blue palette from RColorBrewer}
\item{linescol}{color for the lines drawing the squares}
\item{map.range}{range of the overlap map (default: minimum and maximum of lower triangle of the matrix)}
\item{lab.col}{color of the labels (same as nrow(matrix) if provided)}
\item{lab.cex}{character size of the label}
\item{map.cex}{character size of the overlap values laid on the map}
\item{legend.cex}{character size of the legend text (does not work always)}
}
\details{
Provides map of overlap values for each group of mixture model.
}
\references{
Maitra, R. (2010). A re-defined and generalized percent-overlap-of-activation measure for studies of fMRI reproducibility and its use in identifying outlier activation maps. Neuroimage, 50(1), 124-135.
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

}

