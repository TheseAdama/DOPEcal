MultiNeighbor <- function(X, G){
  nx <- nrow(X)
  L <- nrow(G)
  Pop <- array(0, dim = c(L, nx+1, ncol(X)))
  for(l in 1:L){
    Pop[l, , ] <- rbind(X, G[l,])
  }
  
  return(Pop)
}

LHSDset <- function(Dset, Dx, N){
  xx = LHSD(N, Dx)
  L = nrow(Dset)
  G = c()
  for(l in 1:L){
    G <- rbind(G, Dbind(xx, Dset[l,]))
  }
  
  return(G)
}

# Proposal distributon for discrete variable
ProbaMat <- function(N){
  
  if(N<=4) message("N must be >4 !")
  Mat = matrix(0,nrow = N, ncol=N)
  for(i in 3:(N-2)){
    Mat[i,i-2] <-  1/7
    Mat[i,i-1] <- 2/7
    Mat[i,i] <- 1/7
    Mat[i,i+1] <- 2/7
    Mat[i,i+2] <- 1/7
  }
  
  Mat[1,1] <- 1/4
  Mat[1,2] <- 1/2
  Mat[1,3] <- 1/4
  
  Mat[2,1] <- 1/3 
  Mat[2,2] <- 1/6
  Mat[2,3] <- 1/3
  Mat[2,4] <- 1/6 
  
  Mat[N,N-2] <- 1/4
  Mat[N,N-1] <- 1/2
  Mat[N,N] <- 1/4
  
  Mat[N-1,N-3] <- 1/6
  Mat[N-1,N-2] <- 1/3
  Mat[N-1,N-1] <- 1/6
  Mat[N-1,N] <- 1/3
  
  return(Mat)
}

Neighbordisc <- function(xx, Dset, Dx, ProbMat){
  
  # Continuous variables
  x  = xx[1:nrow(Dx)]
  xv = Neighbor(x, Dx, type="uniform")
  
  # Discret variables
  xd = xx[(nrow(Dx)+1):(nrow(Dx)+ncol(Dset))]
  ind = match(xd, Dset)[1]
  xdv = Dset[sample.int(nrow(Dset), 1, prob = ProbMat[ind,]), ]
  
  # Combine the two neighbor
  xxv <- cbind(matrix(xv, nrow=1), matrix(xdv, nrow=1))
  
  return(as.vector(xxv))
}



# Recuit simulé : cas entrée discrete ou categorielle
SAOptimdisc <- function(xinit, Dx, Dset, f, maxiter,
                        c=0.99, Tinit, Paccept="metropolis", schema="geo"){
  
  xinit <- as.matrix(xinit)
  n <- nrow(xinit)
  d <- ncol(xinit)
  
  ProbMat = ProbaMat(nrow(Dset))
  
  Yinit <- f(xinit)
  temp <- Tinit #Initial temp
  
  Xstar <- xinit
  Ystar <- Yinit
  k <- 1
  while(k <= n*maxiter){
    i = k%%n
    if(i==0){i=n}
    X = Xstar
    X[i,] <- Neighbordisc(xx=X[i,], Dset, Dx, ProbMat)
    Y <- f(X)
    if (is.na(Y)){
      next
    }
    Delta <- Ystar - Y
    
    if(Paccept=="metropolis"){
      alpha = min(1, exp(-Delta/temp))
    }else{
      alpha = exp(-Delta/temp)/(1+exp(-Delta/temp))
    }
    
    u <- runif(1)
    if(alpha > u){
      Xstar <- X
      Ystar <- Y
    }
    
    if(k%%n==0){ 
      if(schema=="geo"){
        temp <- c*temp
      }
      if(schema=="log"){
        temp <- temp/log(1+k)
      }
    }
    
    # Calcul du pourcentage d'évolution
    rate <- (k /(n*maxiter))* 100
    
    # Affichage du pourcentage
    cat(sprintf("Progression Optimisation: %.2f%%\r", rate))
    
    # Forcez l'affichage dans la console
    flush.console()
    
    k <- k + 1
  }
  return(Xstar)
}

# Algo glouton : cas discret
ForwardOptimDiscret <- function(f, nD, tomax=TRUE, Dset, Dx, N, echo=TRUE,...){
  
  # Verification des paramètres
  if(!is.function(f)) stop("f must be a function to optimize !")
  if (!is.matrix(Dx)) stop("The domain Dx must be a matrix !")
  if (!is.numeric(nD) || nD < 1) stop("An integer greater than 0 for nD must be entered.")
  if (!is.logical(tomax)) stop("The parameter tomax must be TRUE or FALSE.")
  if (!is.numeric(N) || N < 1) stop("N must be a numeric value greater than 0.")
  if (!is.matrix(Dset)) stop("The domain Dset must be a matrix !")
  
  if(is.null(N)) {
    if(echo) message("Setting the number of candidat N to 10**3")
    N <- 1000
  }
  
  G <- LHSDset(Dset, Dx, N)
  nc = nrow(Dx)+ncol(Dset)
  Xopt <- matrix(ncol=nc, nrow=1)
  for(i in 1:nD){
    
    if(i==1){
      P <- array(0, dim = c(nrow(G), 1, nc))
      P[c(1:nrow(G)), , ] <- G
    }else{
      P <- MultiNeighbor(Xopt, G)
    }
    
    fP = apply(P, 1, function(x){return(f(matrix(x, ncol=nc, byrow = FALSE)))})
    
    if(tomax==TRUE){
      Xopt <- P[which.max(fP), , ]
      Copt <- max(fP)
      G = G[-which.max(fP),]
    }else{
      Xopt <- P[which.min(fP), , ]
      Copt <- min(fP)
      G = G[-which.min(fP),]
    }
    Xopt = matrix(Xopt, ncol=nc, byrow = FALSE)
  }
  
  return(list(Xopt=Xopt, Copt = Copt))
}
