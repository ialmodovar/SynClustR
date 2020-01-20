##***************************************************************************
##
## @file:overlap.mat.R
##
## Display the map of pairwise overlap measures of Maitra and Melnykov
## (JCGS, 2012)
##
## overlap.mat = matrix of total pairwise overlaps (necessarily symmetric)
## map.col = colormap for the mapping
## linescol = color for the lines drawing the squares
## map.range = range of the overlap map (default: minimum and maximum of lower
##               triangle of the matrix)
## lab.col = color of the labels (same as nrow(matrix) if provided)
## lab.cex = character size of the label
## map.cex = character size of the overlap values laid on the map
## legend.cex = character size of the legend text (does not work always)
##
## provides map of overlap values for each group of mixture model
##
## written Ranjan Maitra, Ames, IA 50011-1210, June 28, 2009
##
## modified Ranjan Maitra, Ames, IA 50011-1090, December 23, 2016.
## modification to bring in specifications for label color, labels, maps and
## legend character size
##
## Last Modified by
## Israel Almodovar-Rivera, PhD                       
## University of Puerto Rico                          
## Medical Science Campus                             
## Graduate School of Public Health                   
## Department of Biostatistics and Epidemiology       
## San Juan, PR, USA 00953
## email: israel.almodovar@upr.edu
##*********************************************************************************

overlap.map <- function(overlap.mat, map.col = c("#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"), linescol = "magenta", map.range = NULL, lab.col = 1, lab.cex = 0.95, map.cex = 0.85, legend.cex = 1)
{
    oxmat <- overlap.mat
    oxmat[lower.tri(oxmat)] <- NA
    diag(oxmat) <- NA
    p <- ncol(oxmat)
    newox <- oxmat[-p, -1]
    newox <- cbind(rbind(NA, newox), NA)[, p:1]

#    layout(matrix(c(rep(1, 4*p^2), rep(2, 2*p), rep(3, 2*p), rep(4, 2*p), rep(4, 2*p)), nrow = 2*p, ncol = 2*p + 4 ))
     layout(matrix(c(rep(1, p^2), rep(2, p), rep(3, p), rep(4, p)), nrow = p, ncol = p + 3))
#      layout(matrix(c(rep(1, p^2), rep(2, p), rep(3, p), rep(4, p),rep(4,p)), nrow = p, ncol = p + 4))
     
    
    par(mar = c(0.1,0.1,0.75,0.1))
    if (is.null(map.range)) map.range <- range(newox, na.rm = T)  
    image(x = 1:p, y = 1:p, z = newox, axes = F, xlab = "", ylab = "",col = map.col, zlim = map.range)
    text(y = 2:p, x = rep(1, p-1), labels = p:2, cex = lab.cex, col = lab.col[p:2])
    text(x = 2:p, y = rep(1, (p-1)), labels = 1:(p-1), cex = lab.cex, col = lab.col[1:(p-1)])

    for(i in 1:p) {
        for(j in i:p) {
            lines(x = c(i+0.5, i+0.5), y = c(p-j+1,p-j)+1.5, col = linescol, lwd = 0.5)
            lines(y = c(i+0.5, i+0.5), x = c(p-j+1,p-j)+1.5, col = linescol, lwd = 0.5)
        }
    }
    for(i in 2:p) text(x=1:p, y = i, labels = round(newox[,i], 2),
                       col = ifelse(newox[,i] < median(map.range), "black", "white"), cex = map.cex)

    frame()
    
    savepar <- par(cex=0.75, lwd=0.25, mar = c(1, 0.5, 1, 2),
                   xaxs="i", yaxs="i")
    plot.new()
    length.col <- length(map.col) + 1
    ra <- seq(from = map.range[1], to = map.range[2], length=length.col)
    plot.window(xlim=c(0,0.1), ylim= c(map.range[1], map.range[2]))
    rect(0, ra[-length.col], 1, ra[-1], col = map.col, border = NULL)
    axis(4, at = ra, labels = round(ra, digits = 3), las = 1, line = NA, cex.axis = legend.cex)
    rect(0, 0, 1, ra[length.col], col = NULL)
    frame()
}

