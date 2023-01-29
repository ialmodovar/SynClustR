##*************************************************************************************************************************************
#' KNOBSynC
#' 
#' @title Kernel-estimated Nonparametric Overlap-Based Syncytial Clustering
#' 
#' @description Merge clustering components using smooth estimation of the missclassification probabilities
#' 
#' @param X dataset of size n x p
#' @param kmns.results clustering solution. Default is NULL, if provided a class kmeans can be input. This can be a list with centers and cluster ids  
#' @param min.gen.overlap Minimum desired generalized overlap of Maitra (2010). Algorithm will stop when this overlap has been reached. Default is 0.00001. 
#' @param kappa  Syncytial clustering parameter. As kappa increase fewer clusters will be merge. For more details see Almodovar-Rivera and Maitra (2020).
#' @param Kmax Maximum value of the range of the number of cluster. Default is NULL, if n >= 50, Kmax = max{50, sqrt(n)} otherwise Kmax=sqrt(n)
#' @param EstK If \code{kmns.results} is NULL, which methodology will be used to estimated the number of groups, options are "jump" of Sugar and James (2003) and "KL" of Krzanoswki and Lai (1985)
#' @param kernel kernel estimator to be use to estimate missclassification probabilites. Default is reciprocal inverse gaussian (RIG), other choices see \code{kcdf}.
#' @param b Smoothing parameter to be use in the estimation of the probabilities. Default is NULL, bandwidth to be use is the one that minimize the mean integrated squared error (MISE).
#' @param inv.roles if TRUE will use gamma kernel of Kim (2013), default is Chen 2000
#' @param desired.ncores Desired number of cores to be used. Default is 2, however the function will determine min(detectCores(),desired.ncores)
#' @param ret.steps If TRUE will return all the information at each step of the algorithm
#' @param verbose If TRUE will print each step  
#' @return A list with the following
#' \begin{itemize}
#' \item KmeansSoln - Return the clustering solution based on k-means. If \code{kmns.results} was NULL, then the clustering solution will be based on EstK. If provided by the user it will be the same as the input. 
#' \item OmegaMapKmns - Overlap matrix based on the k-means solution. 
#' \item OmegaMapKNS - Overlap matrix based on KNOB-SynC solution. 
#' \item Ids - Cluster membership based on k-means.
#' \item IdsMerge - Cluster membership based on KNOB-SynC.
#' \item Groups - A list containing the Groups that were merge together.
#' \item MaxOverlap - Maximum pairwise overlap for the KNOB-SynC Overlap matrix.
#' \item MeanOverlap - Average pairwise overlap for the KNOB-SynC Overlap matrix.
#' \item GenOverlap - Generalized pairwise overlap for the KNOB-SynC Overlap matrix.
#' \item Psi - Normed residuals for the k-means solution
#' \item Fhat - Estimation of the missclassification probability for each observation
#' \item bw - Numeric value of the smoothing parameter. If the default was NULL, then the use estimated value is return, otherwise it returns the user input.
#' \end{itemize}
#' @examples 
#' \dontrun{
#' set.seed(787)
#' ## Example 1
#' data(Bullseye)
#' oo <- KNOBSynC(x = Bullseye[,-3],verbose = TRUE)
#' Bullseye$IdsKmeans <- oo$Ids
#' Bullseye$IdsKNOBSynC <- oo$IdsMerge
#' par(mfrow=c(1,3))
#' with(Bullseye,plot(x = x,y = y, col=Ids,main="True"))
#' with(Bullseye,plot(x = x,y = y, col=IdsKmeans,main="k-means"))
#' with(Bullseye,plot(x = x,y = y, col=IdsKNOBSynC,main="KNOB-SynC"))
#' ## Example 2  
#' data(Spherical7)
#' oo <- KNOBSynC(x = Spherical7[,-3],verbose = TRUE)
#' Spherical7$IdsKmeans <- oo$Ids
#' Spherical7$IdsKNOBSynC <- oo$IdsMerge
#' par(mfrow=c(1,3))
#' with(Spherical7,plot(x = x,y = y, col=Ids,main="True"))
#' with(Spherical7,plot(x = x,y = y, col=IdsKmeans,main="k-means"))
#' with(Spherical7,plot(x = x,y = y, col=IdsKNOBSynC,main="KNOB-SynC"))
#' }
#' @export
##***********************************************************************************************************************************

which.max.matrix <- function(mat) (which(x = mat == max(mat), arr.ind=T))

is.integer0 <- function(x){is.integer(x) && length(x) == 0L}

KNOBSynC <- function(x, kmns.results=NULL, min.gen.overlap = 1e-5,kappa = NULL, 
                      Kmax = NULL,EstK=NULL,kernel = "RIG", b = NULL,inv.roles=FALSE, 
                      desired.ncores=2, 
                     ret.steps = FALSE,verbose=FALSE,...)
{
  X <- as.matrix(x)
  n <- nrow(X);p <- ncol(X)
  if(verbose){
    cat(paste("#",paste(rep("=",100),collapse=""),"#\n",sep=""))
    cat(" \t \t Kernel-estimated Nonparametric Overlap-Based Syncytial Clustering  \n\n")  
  }
  if((is.null(kmns.results))) {
    EstK <- ifelse(is.null(EstK),ifelse(p^2 <= n,"jump","KL" ),match.arg(EstK,choices = c("jump","KL")))
    if(is.null(Kmax)){
      Kmax <- ifelse(n < 50, round(sqrt(n),digits = 0),max(round(sqrt(n),digits = 0),50))
    }
    if(verbose){
      cat("Printing set-up:\n")
      cat(sprintf("Performing k-means range of number of groups using %s method [1,%d]. \n",EstK,Kmax))
    }
    
    kmeans.results <- kmeans.all(x = X,maxclus = Kmax, desired.ncores = desired.ncores,...)
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
      cat("k-means solution was provided.\n")
    }
    Means <- kmns.results$centers
    ids <- kmns.results$cluster
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
    ##    Kappa <- length(kappa)
    Kappa <- kappa
  }
  ##************************************************
  ## compute Psi = || X_i - \mu_ik|| and pseudo Psi
  ##************************************************
  if(verbose){
    cat("Computing normed residuals and pseudo-normed residuals. \n")
  }
  desired.ncores <- max(availableCores(),desired.ncores)
  residuals.norm <- norm.res(X = X, Means = Means, ids = ids,desired.ncores=desired.ncores)
  
  psi <- residuals.norm$Psi 
  pseudo.psi <- residuals.norm$PseudoPsi
  
  ## minimum generazid overlap
  ##  min.gen.overlap <- 1/prod(dim(X))
  ##*************************************
  ##
  ## Choice of the kernel to estimate
  ## the missclasification probabilities
  ##
  ##*************************************
  if(verbose){
    cat("Estimating the missclasification probabilities using kernel estimation.\n")
  }
  if(is.null(b)){
    b <- gsl_bw_mise_cdf((psi*psi/(sum(psi*psi)/(p*(n-K)))))
  }
  
  cl <- makeCluster(desired.ncores,...)
  clusterExport(cl, list("kcdf"))
  Fhat.Psi <- t(parApply(cl = cl,X = pseudo.psi,MARGIN = 1,FUN = function(z){kcdf(x = psi,b = b,kernel=kernel,
                                                                                  xgrid=z,inv.roles=inv.roles)$Fhat}))
  stopCluster(cl)
  
  
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
  Groups.step <- list()
  Ids.step <- list()  
  Omega.lk <- array(0,dim=c(K,K))
  for(k in 1:(K-1)){
    for(l in (k+1):K){
      Omega.lk[k,l] <- 1-mean(Fhat.Psi[ids==k,l]) ## Omega_{l|k}
    }
  }

  # for(k in 1:K){
  #   for(l in 1:K){
  #     Omega.lk[k,l] <- 1-mean(Fhat.Psi[ids==k,l]) ## Omega_{l|k}
  #   }
  # }
  
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
  
  idsMerge <- ids
  Kmerge <- K
  Kstep[iter] <- Kmerge
  Omega.lk.merge <- Omega.lk
  Omega.initial <- Omega
  Omega.initial.merge <- Omega.lk.merge
  GG <- list()
  gen.overlap[iter] <- generalized.overlap(overlap.mat = Omega)
  max.overlap[iter] <- max(Omega[!lower.tri(Omega,diag = TRUE)])
  mean.overlap[iter] <- mean(Omega[!lower.tri(Omega,diag = TRUE)])
  Groups.step[[iter]] <- GG
  Ids.step[[iter]] <- idsMerge 
  ##************************************************************************************
  ## If the generalized overlap for the k-means clustering is very low, there's not need to merge them.
  ## Forcing them to merge will create unstable clustering.
  ## \check{\omega}/\genoverlap >= 5 is the same as Infinity see Almodovar and Maitra 2018
  ##**********************************************************************************
  if(max.overlap[iter] <= 4* gen.overlap[iter]){
    if(verbose){
      cat(sprintf("|Iteration |  K  |  Max Overlap  | Mean Overlap |  Generalized Overlap \t| \n"))
      cat(sprintf("| %d \t   | %d  |  %.6f \t |  %.6f \t|  %.6f \t        |\n",iter,Kmerge,max.overlap[iter],mean.overlap[iter],gen.overlap[iter]))
    }
    knob.sync.kappa <- list(KmeansSoln = kmns.results,OmegaMapKmns = Omega.lk, OmegaMapKNS=Omega.lk.merge,Ids = ids, IdsMerge = idsMerge, StepGroups = Groups.step,
                              Groups = GG, MaxOverlap=max.overlap[1:iter],MeanOverlap = mean.overlap[1:iter], GenOverlap = gen.overlap[1:(iter)],Psi=psi, Fhat = Fhat.Psi,bw = b)
    
    Kmerge <- max(unique(knob.sync.kappa$IdsMerge))
    if(verbose){
      cat("C =", Kmerge,"\n") 
      cat(paste("#",paste(rep("=",100),collapse=""),"#\n",sep=""))
    }
    return(knob.sync.kappa)
    
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
      Groups.step <- list()
      Groups.step[[iter]] <- unique(idsMerge)
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
        Groups.step[[iter]] <- GG
        Ids.step[[iter]] <- idsMerge 
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
            if(Kmerge==1){
              ids.of.GG <- GG[[1]]
              nk <- length(ids.of.GG)
              max.omega.gg <- rep(0,nk)
              
              Omega.lk.merge[1,1] <- 1
              iter <- iter+1
              Kstep[iter] <- Kmerge
              gen.overlap[iter] <- 1 ## compute generalized overlap
              max.overlap[iter] <- 1 ## compute maximum overlap
              mean.overlap[iter] <- 1 ## compute mean overlap            
            } else{
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
            }
            Groups.step[[iter]] <- GG
            Ids.step[[iter]] <- idsMerge 
            if(verbose){
              cat(sprintf("| %d \t   | %d  |  %.6f \t |  %.6f \t|  %.6f \t        |\n",iter,Kmerge,max.overlap[iter],mean.overlap[iter],gen.overlap[iter]))
              
            }
          } else{
            break;
          }
        }
      }
      
        Omega.lk2 <- Omega.lk
        Omega.lk2 <- Omega.lk2+ t(Omega.lk2)
        diag(Omega.lk2) <- 1
        
        Omega.lk2.merge <- Omega.lk.merge
        Omega.lk2.merge <- Omega.lk2.merge+ t(Omega.lk2.merge)
        diag(Omega.lk2.merge) <- 1
        if(ret.steps){
          names(Groups.step) <- paste("Step",1:iter,sep="")
          names(Ids.step) <- paste("Step",1:iter,sep="")
          knob.sync.kappa[[kappa]] <- list(KmeansSoln = kmns.results,OmegaMapKmns = Omega.lk2,
                                           OmegaMapKNS=Omega.lk2.merge,Ids = ids, IdsMerge = idsMerge, 
                                           GroupsStep = Groups.step,  IdsStep = Ids.step,
                                           Groups = GG, MaxOverlap=max.overlap[1:iter],
                                           MeanOverlap = mean.overlap[1:iter],
                                           GenOverlap = gen.overlap[1:(iter)],
                                           Psi=psi, Fhat = Fhat.Psi, bw = b)
        } else{
          knob.sync.kappa[[kappa]] <- list(KmeansSoln = kmns.results,OmegaMapKmns = Omega.lk2,
                                           OmegaMapKNS=Omega.lk2.merge,Ids = ids, IdsMerge = idsMerge, 
                                           Groups = GG, MaxOverlap=max.overlap[1:iter],
                                           MeanOverlap = mean.overlap[1:iter],
                                           GenOverlap = gen.overlap[1:(iter)],
                                           Psi=psi, Fhat = Fhat.Psi, bw = b)
        }
      }
    ##
    ## choose best kappa solution by
    ## finding the minimum generalized overlap
    ##
    if(length(Kappa) > 1){
      kappa.min <- which.min(sapply(knob.sync.kappa,function(z) z$GenOverlap[length(z$GenOverlap)]))
    } else{
      kappa.min <- which(sapply(knob.sync.kappa,is.null)==FALSE)
    }
    Kmerge <- max(unique(knob.sync.kappa[[kappa.min]]$IdsMerge))
    if(verbose){
        cat("C =", Kmerge,"\n") 
      }
      cat(paste("#",paste(rep("=",100),collapse=""),"#\n",sep=""))
    }
    knob.sync.kappa[[kappa.min]]
  }
