##***************************************************************************************
## 
## @file: KNOB-SynC.R
##
## Merge cluster using Kernel nonparametric overlap based syncital (KNOB-SynC).
## kernel estimator choices are RIG (default), gamma (original and inverse roles), 
## and gaussian. See Almodovar-Maitra (2018) for details.
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
## Copyright May 2018
##***************************************************************************************


which.max.matrix <- function(mat) (which(x = mat == max(mat), arr.ind=T))

is.integer0 <- function(x){is.integer(x) && length(x) == 0L}

##****************************************
## calculate generalized overlap
##****************************************

generalized.overlap <- function(overlap.mat) {
  p <- nrow(overlap.mat)
  if (is.null(p)) 1 else {
    x <- eigen(overlap.mat, symmetric = TRUE, only.values = TRUE)
    (x$values[1] - 1)/(p - 1)
  }
}


norm.res <- function(X, Means, ids){
  ##******************************
  ##
  ## compute normed residual
  ##
  ##******************************
  if((ncol(X)!=ncol(Means)) | (length(ids) !=  nrow(X)) | (nrow(Means) != max(ids))){
    stop("Sizes do not match \n")
  }
  X <- as.matrix(X)
  Means <- as.matrix(Means)
  n <- nrow(X)
  p <- ncol(X)
  K <- max(ids)
  
  Eps <- lapply(1:n,function(z) {
    find.nan <- which(!is.na(X[z,])); 
    X[z,find.nan]-Means[ids[z],find.nan]}
  )
  
  psi <- sapply(Eps,norm,type = "2")
  
    names(psi) <- NULL


   na.total <-  apply(X,1,function(z) sum(!is.na(z)))/p
     psi  <- psi * na.total
  
  ##******************************
  ## compute psi and pseudo psi
  ##******************************
  if(p > 1){
    
    Diff.Psi <- t(sapply(1:n,function(i) {
      apply(t(sapply(1:K, function(k) {
        find.nan <- which(!is.na(X[i,])); 
        X[i,find.nan]-Means[k,find.nan]})),1,norm,type="2") * sum(!is.na(X[i,]))/p
    }
      ))
    
  } else{
    Diff.Psi <- t(sqrt(apply(X,1,function(z) {
      find.nan <- which(!is.na(z));
      eps <- z[find.nan]-Means[,find.nan]
      ##eps * eps
       eps * eps * sum(!is.na(z))/p
    })))
  }
  
  
  ##**************************
  ## Compute Pseudo residuals
  ##**************************
  nks <- as.numeric(table(ids))
  pseudo.psi<-  t(sapply(1:n,function(i){
    k <- ids[i]
    Means2 <- Means
    nk <- nks[k]
    find.nan <- which(!is.na(X[i,]));
    Means2[k,find.nan] <- (nk * Means[k,find.nan] - X[i,find.nan])/(nk-1)
    sapply(1:K, function(l){
      nl <- nks[l]
      Means2[l,find.nan] <- (nl * Means[l,find.nan] + X[i,find.nan])/(nl+1)
      Psi.rev<- norm(X[i,find.nan]-Means2[l,find.nan],type = "2") * sum(!is.na(X[i,]))/p
      Psi.rev
    })
  }))
  
  list(Psi = psi, PsiAll=Diff.Psi, E = Eps,PseudoPsi=pseudo.psi)
}

##KNOBSynC
KNOBSynC <- function(x, kmns.results=NULL, Kmax = NULL,EstK="jump",kernel = "RIG",
                     kappa = NULL, b = NULL, trueids=NULL,inv.roles=FALSE, 
                    ret.steps = FALSE,verbose=FALSE,...)
{
  ##*************************************************************************************************************************************
  ## X: dataset of size n x p
  ## kmns.results: kmeans results 
  ## nstart: number of initializations to be use along k-means. Only need it if kmns.results are not provided
  ## kernel: kernel estimator to be use, choices are RIG (default), gamma with/without inverse roles and gaussian  
  ## b: Bandwidth to be use, if null a estimated on the MISE will be use
  ## trueids: if provided Adjusted Rand Index will be return for each step
  ## inv.roles: if TRUE will use gamma kernel of Kim (2013), default is Chen 2000
  ## kappa: multiple the generalized overlap. Complexity of merging Higher values of kappa less clusters are merge
  ## display.overlap.mat= Display overlap matrix Maitra 2010
  ## verbose: if TRUE will print each step  
  ##***********************************************************************************************************************************
  X <- as.matrix(x)
  n <- nrow(X);p <- ncol(X)
  if(verbose){
    cat(paste("#",paste(rep("=",100),collapse=""),"#\n",sep=""))
    cat(" \t \t KNOB-SynC: Kernel-estimated Nonparametric Overlap-Based Syncytial Clustering  \n\n")  
  }
  if((is.null(kmns.results))) {
    ##EstK <- ifelse(is.null(EstK),ifelse(p^2 <= n,"jump","KL" ),match.arg(EstK,choices = c("jump","KL")))
    EstK <- match.arg(EstK,choices = c("jump","KL"))
    if(is.null(Kmax)){
      ##  Kmax <- ifelse(is.null(Kmax),ifelse(n < 50, round(sqrt(n),digits = 0),max(round(sqrt(n),digits = 0),50)),ifelse(Kmax> n,round(sqrt(n),digits = 0),Kmax))
      Kmax <- ifelse(n < 50, round(sqrt(n),digits = 0),max(round(sqrt(n),digits = 0),50))
    }
    if(verbose){
      cat("Printing set-up:\n")
      cat(sprintf("Performing k-means range of number of groups using %s method [1,%d]. \n",EstK,Kmax))
    }
    
    kmeans.results <- kmeans.all(x = X,maxclus = Kmax,...)
    if(EstK=="jump"){
      K <- which.max(unlist(kmeans.results$jump.stat))
      kmns.results <- kmeans.results$kmns.results[[K]]
      Means <- kmns.results$centers
      ids <- kmns.results$cluster
      wss <- kmns.results$tot.withinss
    }  else{
      K <- which.max(unlist(kmeans.results$kl.stat))
      kmns.results <- kmeans.results$kmns.results[[K]]
      Means <- kmns.results$centers
      ids <- kmns.results$cluster     
      wss <- kmns.results$tot.withinss
    }
  }  else{
    if(verbose){
      cat("k-means solution was provided.  \n")
    }
    Means <- kmns.results$centers
    ids <- kmns.results$cluster
    wss <- kmns.results$tot.withinss
    K <- max(unique(ids))
  }
  
  kernel <- match.arg(kernel,choices = c("gamma","RIG","gaussian"))
  if(verbose){
    if(!is.null(b)){
      cat("Kernel: ", kernel, "\nBandwidth:", b,"\n")
      cat("Initial number of groups:", K, "\n")
    }else{
      cat("Kernel: ", kernel, "\nBandwidth will be estimated based on ruled-of-thumb for MISE. \n")
      cat("Initial number of groups:", K, "\n")
    }
  }
  if(is.null(kappa)){
    if(verbose){
      cat("Kappa was not provided following, 1,2,3,4,5 \n")
    }
    Kappa <- c(1,2,3,4,5)
  } else{
    if(verbose){
      cat(sprintf("Kappa was provided merging clusters will occur %s times generalized overlap \n", as.character(kappa)))
    }
    Kappa <- length(kappa)
    ##Kappa <- kappa
  }
  ##************************************************
  ## compute Psi = || X_i - \mu_ik|| and pseudo Psi
  ##************************************************
  if(verbose){
    cat("Computing normed residuals and pseudo-normed residuals. \n")
  }
  residuals.norm <- norm.res(X = X, Means = Means, ids = ids)
  
  psi <- residuals.norm$Psi 
  pseudo.psi <- residuals.norm$PsiAll
  psi.rev <- residuals.norm$PseudoPsi
  
  ## minimum generazid overlap
  min.gen.overlap <- 1/prod(dim(X))
  ##*************************************
  ##
  ## Choice of the kernel to estimate
  ## the missclasification probabilities
  ##
  ##*************************************
  if(verbose){
    cat("Estimating the missclasification overlaps using kernel estimation.\n")
  }
  if(is.null(b)){
    b <- gsl_bw_mise_cdf((psi*psi/(wss/(p*(n-K)))))
  }
  Fhat.Psi <- t(apply(pseudo.psi,1,function(z){kcdf(x = psi,b = b,kernel=kernel,xgrid=z,inv.roles=inv.roles)$Fhat}))
  
  ##*********************************
  ##
  ## compute \hat{omega}_{kl}
  ##
  ##*********************************
  
  iter <- 1
  gen.overlap <- rep(0,K)  ## store compute generalized overlap
  max.overlap <- rep(0,K)  ## store compute max overlap
  mean.overlap <- rep(0,K)  ## store compute bar overlap
  Kstep <- rep(0,K)
  
  Omega.lk <- array(0,dim=c(K,K))
  for(k in 1:(K-1)){
    for(l in (k+1):K){
      Omega.lk[k,l] <- 1-mean(Fhat.Psi[ids==k,l]) ## Omega_{l|k}
    }
  }
  
  omega.mat <- Omega.lk
  Omega <- (t(omega.mat) + (omega.mat))
  diag(Omega) <- 1.00
  ##***************************************************
  ##
  ## this asymmetric kernel does not integrate to 1
  ## therefore some values can be larger than 1.
  ## compute generalized overlap of Maitra 2010
  ##
  ##***************************************************
  
  if(!is.null(trueids)){
    ar <- rep(0,K)
    ar[iter] <- RandIndex(id1=trueids,id2=ids)$AR
  }

  idsMerge <- ids
  Kmerge <- K
  Kstep[iter] <- Kmerge
  Omega.lk.merge <- Omega.lk
  Omega.initial <- Omega
  Omega.initial.merge <- Omega.lk.merge
  GG <- list()
  gen.overlap[iter] <- generalized.overlap(overlap.mat = Omega)
  max.overlap[iter] <- max(Omega.lk.merge[!lower.tri(Omega.lk.merge,diag = TRUE)])
  mean.overlap[iter] <- mean(Omega.lk.merge[!lower.tri(Omega.lk.merge,diag = TRUE)])
  ##************************************************************************************
  ## If the generalized overlap for the k-means clustering is very low, there's not need to merge them.
  ## Forcing them to merge will create unstable clustering.
  ## \check{\omega}/\genoverlap >= 5 is the same as Infinity see Almodovar and Maitra 2018
  ##**********************************************************************************
  if(max.overlap[iter]/gen.overlap[iter] <= 4){
    if(verbose){
      cat(sprintf("|Iteration |  K  |  Max Overlap  | Mean Overlap |  Generalized Overlap \t| \n"))
      cat(sprintf("| %d \t   | %d  |  %.6f \t |  %.6f \t|  %.6f \t        |\n",iter,Kmerge,max.overlap[iter],mean.overlap[iter],gen.overlap[iter]))
    }
    if(!is.null(trueids)){
      knob.sync.kappa <- list(KmeansSoln = kmns.results,OmegaMapKmns = Omega.lk, OmegaMapKNS=Omega.lk.merge,Ids = ids, IdsMerge = idsMerge, 
                              Groups = GG, MaxOverlap=max.overlap[1:iter],MeanOverlap = mean.overlap[1:iter], GenOverlap = gen.overlap[1:(iter)],Psi=psi, Fhat = Fhat.Psi, AR= ar[1:iter],bw = b)
    } else{
      knob.sync.kappa <- list(KmeansSoln = kmns.results,OmegaMapKmns = Omega.lk, OmegaMapKNS=Omega.lk.merge,Ids = ids, IdsMerge = idsMerge, 
                              Groups = GG, MaxOverlap=max.overlap[1:iter],MeanOverlap = mean.overlap[1:iter], GenOverlap = gen.overlap[1:(iter)],Psi=psi, Fhat = Fhat.Psi,bw = b)
    }
    Kmerge <- max(unique(knob.sync.kappa$IdsMerge))
    if(verbose){
      if(!is.null(trueids)){
        ar <- knob.sync.kappa$AR[length(knob.sync.kappa$AR)]
        cat("C =", Kmerge,"\nAdjusted Rand Index =", ar, "\n") 
      } else{
        cat("C =", Kmerge,"\n") 
      }
      cat(paste("#",paste(rep("=",100),collapse=""),"#\n",sep=""))
    }
    knob.sync.kappa
    
  } else {
    
    knob.sync.kappa <- list()
    for(kappa in Kappa){
      if(verbose){
        cat(sprintf("Kappa = %s,\n", as.character(kappa)))
      }
      iter <- 1
      idsMerge <- ids
      Kmerge <- K
      Kstep[iter] <- Kmerge
      GG <- list()
      gen.overlap[iter] <- generalized.overlap(overlap.mat = Omega.initial)
      max.overlap[iter] <- max(Omega.initial.merge[!lower.tri(Omega.initial.merge,diag = TRUE)])
      mean.overlap[iter] <- mean(Omega.initial.merge[!lower.tri(Omega.initial.merge,diag = TRUE)])
      if(verbose){
        cat(sprintf("|Iteration |  K  |  Max Overlap  | Mean Overlap |  Generalized Overlap \t| \n"))
        cat(sprintf("| %d \t   | %d  |  %.6f \t |  %.6f \t|  %.6f \t        |\n",iter,Kmerge,max.overlap[iter],mean.overlap[iter],gen.overlap[iter]))
      }
      
      index <- which(Omega.lk > min(Omega.lk[Omega.lk  > kappa*gen.overlap[iter]]),arr.ind = TRUE)
      if(length(index) > 0){
        index.order <- index[order(index[,1]),]
        if(!is.null(dim(index.order))){
          idstobeMerge <- unique(index.order[,1])
        } else{
          index.order <- t(as.matrix(index.order))
          idstobeMerge <- unique(index.order[,1])
        }
        
        index <- index.order[,2]
        names(index) <- paste(index.order[,1])
        
        idstobeMerge <- (as.numeric(names(index)))
        ids.unique <- unique(index)
        Kp <- length(ids.unique)
        K.ids <- length(index)
        group.merge <- NULL
        idsMerge[idsMerge == index[K.ids]] <- idstobeMerge[K.ids]
        ggMerge <- c(index[K.ids],idstobeMerge[K.ids])
        group.merge <- c(group.merge,ggMerge)
        GG[[1]] <- group.merge
        
        for(i in (K.ids-1):1){
          
          ggMerge <- c(index[i],idstobeMerge[i])
          
          index.GG <- which(unlist(lapply(lapply(GG,function(z) ggMerge %in% z),function(w) any(w==TRUE))))
          
          if(is.integer0(index.GG)){
            group.merge <- NULL
            ng <- length(GG)
            GG[[ng+1]] <- ggMerge
            group.merge <- c(group.merge,ggMerge)
            group.merge <- unique(group.merge)
            GG[[ng+1]] <- group.merge
          } else{
            if(length(index.GG) > 1){
              np <- length(index.GG)
              sort.index.GG <- sort(index.GG)
              min.index.GG <- min(sort.index.GG)
              for(j in sort.index.GG[sort.index.GG>min.index.GG]){
                GG[[min.index.GG]] <- unique(c(GG[[min.index.GG]],GG[[j]]))
                GG <- GG[-j]
                if(length(GG)==1){
                  break;
                }
              }
            } else{
              group.merge <- GG[[index.GG]]
              group.merge <- c(group.merge,ggMerge)
              group.merge <- unique(group.merge)
              GG[[index.GG]] <- group.merge
            }
          }
        }
        
        ## unique members
        unique.g <- unique(ids)
        uni.g <- unique.g[!(unique.g %in% unlist(GG))]
        if(!is.integer0(uni.g)){
          N.GG <- length(GG)
          for(i in 1:length(uni.g)){
            GG[[N.GG+i]] <- uni.g[i]
          }
        }
        
        Kp <- length(GG)
        
        for(i in 1:Kp){
          ids.of.GG <- GG[[i]]
          nk.gg <- length(ids.of.GG)
          for(j in 1:nk.gg){
            idsMerge[idsMerge==ids.of.GG[j]] <- K+i # to avoid any conflict with ids
          }
        }
        
        idsMerge[idsMerge > K] <- idsMerge[idsMerge > K] - K
        Kmerge <- length(names(table(idsMerge)))
        
        idstobeMerge <- as.numeric(names(table(idsMerge)))
        for(i in 1:Kmerge){
          idsMerge[idsMerge==idstobeMerge[i]] <- i
        }
        ##*********************************
        ## Compute omega_{C_k C_l}
        ##********************************
        Omega.lk.merge <- array(0,dim=c(Kmerge,Kmerge))
        
        for(k in 1:(Kmerge-1)){
          for(l in (k+1):(Kmerge)){
            ids.of.GG <- GG[[l]]
            nk <- length(ids.of.GG)
            max.omega.gg <- rep(0,nk)
            for(j in 1:nk){
              max.omega.gg[j] <- (1-mean(apply(as.matrix(Fhat.Psi[ids == ids.of.GG[j],GG[[k]]]),1,min)))^nk
            }
            Omega.lk.merge[k,l] <- max(max.omega.gg)
          }
        }
        
        iter <- iter+1
        Kstep[iter] <- Kmerge
        omega.mat <- Omega.lk.merge
        Omega <- (t(omega.mat) + (omega.mat))
        diag(Omega) <- 1.00
        
        gen.overlap[iter] <- generalized.overlap(overlap.mat = Omega)
        max.overlap[iter] <- max(Omega.lk.merge[!lower.tri(Omega.lk.merge,diag = TRUE)])
        mean.overlap[iter] <- mean(Omega.lk.merge[!lower.tri(Omega.lk.merge,diag = TRUE)])
        if(verbose){
          cat(sprintf("| %d \t   | %d  |  %.6f \t |  %.6f \t|  %.6f \t        |\n",iter,Kmerge,max.overlap[iter],mean.overlap[iter],gen.overlap[iter]))
        }
        if(!is.null(trueids)){
          ar[iter] <- RandIndex(id1=trueids,id2=idsMerge)$AR
        }
        while((max.overlap[iter] > kappa*gen.overlap[iter]) & (gen.overlap[iter] > min.gen.overlap) & (gen.overlap[iter] <= gen.overlap[iter-1])) {
          
          ##***********************************************************
          ##
          ## If merging occurs, then generalized overlap will start
          ## decreasing. If we merge cluster with the highest overlap.
          ##
          ##***********************************************************
          min.merge <- Omega.lk.merge[Omega.lk.merge > kappa*gen.overlap[iter]]
          if((length(min.merge) == 1)&(length(min.merge) != 0) ){
            index <- which(Omega.lk.merge == min.merge,arr.ind = TRUE)
          } else{
            index <- which(Omega.lk.merge > min(min.merge),arr.ind = TRUE)
          }
          index.order <- index[order(index[,1]),]
          if(!is.null(dim(index.order))){
            idstobeMerge <- unique(index.order[,1])
          } else{
            index.order <- t(as.matrix(index.order))
            idstobeMerge <- unique(index.order[,1])
          }
          
          index <- index.order[,2]
          names(index) <- paste(index.order[,1])
          idstobeMerge <- (as.numeric(names(index)))
          ids.unique <- unique(index)
          Kp <- length(ids.unique)
          K.ids <- length(index)
          
          if(length(index)> 0){
            max.ggMerge <- rep(0,K.ids)
            for(i in (K.ids):1){
              ggMerge <- c(index[i],idstobeMerge[i])
              min.ggMerge <- min(ggMerge)
              max.ggMerge[i] <- max(ggMerge)
              GG[[min.ggMerge]] <- unique(c(GG[[min.ggMerge]],GG[[max.ggMerge[i]]]))
              GG[[max.ggMerge[i]]] <- NA
              idsMerge[idsMerge == max.ggMerge[i]] <- min.ggMerge
            }
            
            GG <- GG[-max.ggMerge]
            if(any(unlist(lapply(GG,is.na)))){
              GG <- lapply(GG, function(x) x[!is.na(x)])
            }
            Kmerge <- length(names(table(idsMerge)))
            idstobeMerge <- as.numeric(names(table(idsMerge)))
            
            for(i in 1:Kmerge){
              idsMerge[idsMerge==idstobeMerge[i]] <- i
            }
            
            
            ##**************************************
            ##
            ## generalized overlap is defined as
            ## omega = (\lambda_{(1)}-1)/(K-1)
            ##
            ## Compute overlap measure
            ## \hat{\omega}_{C_k C_l}
            ## Almodovar and Maitra 2018
            ##
            ##*********************************
            Omega.lk.merge <- array(0,dim=c(Kmerge,Kmerge))
            for(k in 1:(Kmerge-1)){
              for(l in (k+1):(Kmerge)){
                ids.of.GG <- GG[[l]]
                nk <- length(ids.of.GG)
                max.omega.gg <- rep(0,nk)
                for(j in 1:nk){
                  max.omega.gg[j] <- (1-mean(apply(as.matrix(Fhat.Psi[ids == ids.of.GG[j],GG[[k]]]),1,min)))^nk
                }
                Omega.lk.merge[k,l] <- max(max.omega.gg)
              }
            }
            
            iter <- iter+1
            Kstep[iter] <- Kmerge
            omega.mat <- Omega.lk.merge
            Omega <- (t(omega.mat) + (omega.mat))
            diag(Omega) <- 1
            gen.overlap[iter] <- generalized.overlap(overlap.mat = Omega) ## compute generalized overlap
            max.overlap[iter] <- max(Omega.lk.merge[!lower.tri(Omega.lk.merge,diag = TRUE)]) ## compute maximum overlap
            mean.overlap[iter] <- mean(Omega.lk.merge[!lower.tri(Omega.lk.merge,diag = TRUE)]) ## compute mean overlap
            
            if(verbose){
              cat(sprintf("| %d \t   | %d  |  %.6f \t |  %.6f \t|  %.6f \t        |\n",iter,Kmerge,max.overlap[iter],mean.overlap[iter],gen.overlap[iter]))
              
            }
            if(!is.null(trueids)){
              ar[iter] <- RandIndex(id1=trueids,id2=idsMerge)$AR
            }
          } else{
            break;
          }
        }
      }
      
      if(!is.null(trueids)){
        ##********************************
        ## Return with Adjusted Rand Index
          ##********************************
          Omega.lk2 <- Omega.lk
          Omega.lk2 <- Omega.lk2+ t(Omega.lk2)
          diag(Omega.lk2) <- 1

          Omega.lk2.merge <- Omega.lk.merge
          Omega.lk2.merge <- Omega.lk2.merge+ t(Omega.lk2.merge)
          diag(Omega.lk2.merge) <- 1

          
        knob.sync.kappa[[kappa]] <- list(KmeansSoln = kmns.results,OmegaMapKmns = Omega.lk2, OmegaMapKNS=Omega.lk2.merge,Ids = ids, IdsMerge = idsMerge, 
                                         Groups = GG, MaxOverlap=max.overlap[1:iter],MeanOverlap = mean.overlap[1:iter], GenOverlap = gen.overlap[1:(iter)],
                                         Psi=psi, Fhat = Fhat.Psi, AR= ar[1:iter],bw = b)
      } else{
                    Omega.lk2 <- Omega.lk
          Omega.lk2 <- Omega.lk2+ t(Omega.lk2)
          diag(Omega.lk2) <- 1

          Omega.lk2.merge <- Omega.lk.merge
          Omega.lk2.merge <- Omega.lk2.merge+ t(Omega.lk2.merge)
          diag(Omega.lk2.merge) <- 1

        knob.sync.kappa[[kappa]] <- list(KmeansSoln = kmns.results,OmegaMapKmns = Omega.lk2,OmegaMapKNS=Omega.lk2.merge,Ids = ids, IdsMerge = idsMerge, 
                                         Groups = GG, MaxOverlap=max.overlap[1:iter],MeanOverlap = mean.overlap[1:iter],GenOverlap = gen.overlap[1:(iter)],
                                         Psi=psi, Fhat = Fhat.Psi, bw = b)
      }
    }
    
    if(length(Kappa) > 1){
      kappa.min <- which.min(sapply(knob.sync.kappa,function(z) z$GenOverlap[length(z$GenOverlap)]))
    } else{
      kappa.min <- which(sapply(knob.sync.kappa,is.null)==FALSE)
    }
    Kmerge <- max(unique(knob.sync.kappa[[kappa.min]]$IdsMerge))
    if(verbose){
      if(!is.null(trueids)){
        ar <- knob.sync.kappa[[kappa.min]]$AR[length(knob.sync.kappa[[kappa.min]]$AR)]
        cat("C =", Kmerge,"\nAdjusted Rand Index =", ar, "\n") 
      } else{
        cat("C =", Kmerge,"\n") 
      }
      cat(paste("#",paste(rep("=",100),collapse=""),"#\n",sep=""))
    }
      if(!ret.steps){
          knob.sync.kappa[[kappa.min]]
      }
      else{
       knob.sync.kappa
      }
  }
}
