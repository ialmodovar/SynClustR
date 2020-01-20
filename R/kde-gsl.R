##***************************************************************************************
## 
## @file: kde-gsl.R
##
##  Perform kernel estimation for density and distribution function.
##  Kernels in use
##  1) Gamma kernel estimator of Chen (2000).
##  2) Gamma kernel estimator with inverse roles of Kim (2013).
##  3) Reciprocal Inverse Gaussian of Scaillet 2004.
##  4) Classical gaussian kernel.
##
##  Require GSL 
##  All the bandwidths in here are chosen based on the minimum of the
##  mean integrated squared error (MISE). 
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
## Copyright February 2019
##***************************************************************************************

##************************
## Rule-of-thumb bandwidth
## using gamma kernel
## estimate for cdf
##************************

gsl_bw_mise_cdf <- function(x){
  n <- as.integer(length(x))
  ret <- .C("bandwidth_mise_cdf",x = as.double((x)),n = n,bw=as.double(0),PACKAGE="SynClustR")
  ret$bw
}

##************************
## Rule-of-thumb bandwidth
## using gaussian kernel
## estimate for cdf
##************************


gsl_bw_mise_gaussian <- function(x){
  n <- as.integer(length(x))
  ret <- .C("bandwidth_mise_gaussian",x = as.double(x),n = n,bw=as.double(0),PACKAGE="SynClustR")
  ret$bw
}

##************************
## Rule-of-thumb bandwidth
## using gamma kernel
## estimate for pdf
##************************

gsl_bw_mise_pdf <- function(x){
  n <- as.integer(length(x))
  ret <- .C("bandwidth_mise_gamma_pdf",x = as.double((x)),n = n,bw=as.double(0),PACKAGE="SynClustR")
  ret$bw
}

##************************
## Rule-of-thumb bandwidth
## using gaussian kernel
## estimate for pdf
##************************

gsl_bw_mise_gaussian_pdf <- function(x){
  n <- as.integer(length(x))
  ret <- .C("bandwidth_mise_pdf_gaussian",x = as.double((x)),n = n,bw=as.double(0),PACKAGE="SynClustR")
  ret$bw
}

##****************************
## Density estimation using
## gamma kernel. If inverse
## roles is true will perform
## kernel estimate of Kim 2013.
## Default is Chen 2000.
##******************************

gamma_kpdf <- function(x,b = NULL,xgrid=NULL,m = 100,inv.roles=FALSE){
  if(is.null(b)){
    b <- gsl_bw_mise_pdf(x)
  }
  n <- length(x)
  if(is.null(xgrid)){
    xgrid <- seq(from=min(x),to = max(x),length.out = m)
  }
  m <- length(xgrid)
  fhat <- vector(mode = "double",length = m)
  
  if(!inv.roles){
    ff <- .C("gamma_kpdf",x = as.double(x),n=as.integer(n),
             xgrid= as.double(xgrid),
             m = as.integer(m),
             bw=as.double(b),fhat=fhat,PACKAGE="SynClustR")
  } else {
    ff <- .C("gamma_kpdf_inv",x = as.double(x),n=as.integer(n),
             xgrid= as.double(xgrid),
             m = as.integer(m),
             bw=as.double(b),fhat=fhat,PACKAGE="SynClustR")
  }
 list(fhat = ff$fhat, bw=b)
}

##********************************************
## Smooth estimation of the distribution
## function using gamma kernel. If inverse
## roles is true will perform
## kernel estimate of Kim 2013.
## Default is Chen 2000.
##*********************************************

gamma_kcdf <- function(x,b = NULL,xgrid=NULL,m = 100, inv.roles=FALSE){
  if(is.null(b)){
    b <- gsl_bw_mise_cdf(x)
  }
  
  n <- length(x)
  if(is.null(xgrid)){
    xgrid <- seq(from=min(x),to = max(x),length.out = m)
  }
  m <- length(xgrid)
  Fb <- vector(mode = "double",length = m)
  if(!inv.roles){
    ff <- .C("gamma_kcdf",x = as.double(x),
             n=as.integer(n), xgrid = as.double(xgrid),
             m = as.integer(m),
             bw=as.double(b),Fhat=Fb,PACKAGE="SynClustR")
  } else{
    ff <- .C("gamma_kcdf_inv",x = as.double(x),n=as.integer(n),
             xgrid = as.double(xgrid),
             m = as.integer(m),
             bw=as.double(b),Fhat=Fb,PACKAGE="SynClustR")
  }
 list(Fhat = ff$Fhat, bw=b)
}

##****************************
## Smooth estimation of the distribution
## function using reciprocal
## inverse gaussian kernel.
## See Scaillet 2004.
##******************************

RIG_kcdf <- function(x,b = NULL,xgrid=NULL, m = 100){
  if(is.null(b)){
    b <- gsl_bw_mise_cdf(x)
  }
  n <- length(x)
  if(is.null(xgrid)){
    xgrid <- seq(from=min(x),to = max(x),length.out = m)
  }
  m <- length(xgrid)
  Fb <- vector(mode = "double",length = m)
  ff <- .C("RIG_kcdf",x = as.double(x),n=as.integer(n),
           xgrid=as.double(xgrid), 
           m = as.integer(m),
           bw=as.double(b),Fhat=Fb,PACKAGE="SynClustR")
  list(Fhat = ff$Fhat, bw=b)
}

##***************************************
## Smooth estimation of the distribution
## function using gaussian kernel.
##***************************************

gaussian_kcdf <- function(x,b = NULL,xgrid=NULL, m = 100){
  if(is.null(b)){
    b <- gsl_bw_mise_gaussian(x)
  }
  n <- length(x)
  if(is.null(xgrid)){
    xgrid <- seq(from=min(x),to = max(x),length.out = m)
  }
  m <- length(xgrid)
  Fb <- vector(mode = "double",length = m)
  ff <- .C("gaussian_kcdf",x = as.double(x),
           n=as.integer(n),xgrid = as.double(xgrid),
           m = as.integer(m),
           bw=as.double(b),Fhat=Fb,PACKAGE="SynClustR")
  list(Fhat = ff$Fhat, bw=b)
}

##*************************************
## Wrapper to perform kernel estimates
##*************************************

kcdf <- function(x,b = NULL,kernel= "RIG",xgrid=NULL,m = 100, inv.roles=FALSE){
  kernel <- match.arg(kernel,choices=c("RIG","gamma","gaussian"))
  if(any(x <= 0)){
    kernel <- "gaussian"
  } 
  
  if(kernel == "RIG"){
    FF <- RIG_kcdf(x = x,b = b,xgrid = xgrid,m = m)  
  } 
  if(kernel == "gamma"){
    FF <- gamma_kcdf(x = x, b = b, xgrid = xgrid, m = m, inv.roles = inv.roles)
  }
  if(kernel == "gaussian"){
    FF <-gaussian_kcdf(x = x, b = b, xgrid = xgrid,m = m)
  }
 list(Fhat=FF$Fhat,b = FF$bw)
}
