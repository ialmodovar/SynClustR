##***************************************************************************************
## 
## @file: Hmap.R
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## Authors:
##
## Israel Almodovar-Rivera, PhD                       Ranjan Maitra, PhD
## University of Puerto Rico                          Department of Statistics
## Medical Science Campus                             Iowa State University
## Graduate School of Public Health                   Ames, IA, 50010
## Department of Biostatistics and Epidemiology       email: maitra@iastate.edu
## San Juan, PR, USA 00953
## email: israel.almodovar@upr.edu
##
## Copyright December 2019
##***************************************************************************************


hmap <- function(kmH.soln,mycol = c("#FFFFFF","#FFF7FB","#ECE2F0","#D0D1E6","#A6BDDB","#67A9CF","#3690C0","#02818A","#016C59","#014636"),
                 fname="myhmap", h = 4.15, w = 4.775,save.plot=FALSE){
  a <- kmH.soln$all.partitions 
  n<-dim(a)[2]
  xx <- matrix(apply(a[,rep(1:n,n),]==a[,rep(1:n,each=n),],2,sum),nrow=n)/prod(dim(a)[-2])
  
  hU <- heatmap(xx, Rowv = FALSE, symm = TRUE, col = mycol,
                distfun = function(c) as.dist(1 - c),
                hclustfun = function(d) hclust(d, method = "single"),
                keep.dendro = FALSE)
  
  ##
  ## sample the matrix for smaller figures
  ##
  
  mat <- cbind( matrix( 1 , ncol = 4, nrow = 4), rep(2, 4))
  layout(mat, widths = c(rep(83,4), 50), heights = rep(83, 4), respect = F)
  
  par(mar = rep(0.1, 4))
  
  if (ncol(xx) <= 512)
    ind <- 1:ncol(xx) else
      ind <- sort(sample(1:ncol(xx), size = 512))
  
  image(1:length(ind), 1:length(ind), (xx[hU$rowInd, hU$rowInd])[ind, ind], xlab = "", ylab = "", col = mycol, zlim = c(0, 1))
  
  savepar <- par(cex=0.7, lwd=0.25, mar = c(1, 0.25, 1, 3), xaxs="i", yaxs="i")
  
  zp <- c(0, 1)
  plot.new()
  plot.window(xlim=c(0,0.1), ylim= zp)
  rect(0, seq(zp[1], zp[2], length=11)[-11],
       1, seq(zp[1], zp[2], length=11)[-1],
       col = mycol)
  
  axis(4, at = seq(zp[1], zp[2], length = 11), labels = round(seq(zp[1], zp[2], length = 11), digits = 10), las = 1, line = NA)
  rect(0, seq(zp[1], zp[2], length=11)[-11],
       1, seq(zp[1], zp[2], length=11)[-1],
       col = mycol)
  if(save.plot){
  dev.copy2pdf(file = fname, height = h, width = w)
  dev.off()
  } 
}
