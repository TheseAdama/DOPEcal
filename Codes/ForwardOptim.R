Xaugment <- function(X, Dx, N){
  
  G <- LHSD(N,Dx)
  Pop <- array(0, dim = c(N,nrow(X)+1, nrow(Dx)))
  for(i in 1:N){
    Pop[i, , ] <- rbind(X, G[i,])
  }
  
  return(Pop)
}

# Algorithme d'optimisation glouton
ForwardOptim <- function(f, Dx, nD, N=NULL, tomax=TRUE, echo=TRUE) {
  
  # Verification des paramètres
  if(!is.function(f)) stop("f must be a function to optimize !")
  if (!is.matrix(Dx)) stop("The domain Dx must be a matrix !")
  if (!is.numeric(nD) || nD < 1) stop("An integer greater than 0 for nD must be entered.")
  if (!is.logical(tomax)) stop("The parameter tomax must be TRUE or FALSE.")
  if (!is.numeric(N) || N < 1) stop("N must be a numeric value greater than 0.")
  
  if(is.null(N)) {
    if(echo) message("Setting the number of candidat N to 10**3")
    N <- 1000
  }
  
  # Boucle for
  Xopt <- matrix(ncol=nrow(Dx), nrow=1)
  for(i in 1:nD){
    
    # Candidats
    if(i==1){
      P <- array(0, dim = c(N, 1, nrow(Dx)))
      P[c(1:N), , ] <- LHSD(N, Dx)
    }else{
      P <- Xaugment(Xopt, Dx, N)
    }
    
    fP = apply(P, 1, function(x){return(f(matrix(x, ncol=nrow(Dx), byrow = FALSE)))})
    
    # Best Candidat
    if(tomax==TRUE){
      Xopt <- P[which.max(fP), , ]
      Copt <- max(fP)
    }else{
      Xopt <- P[which.min(fP), , ]
      Copt <- min(fP)
    }
    Xopt = matrix(Xopt,ncol=nrow(Dx),byrow = FALSE)
    
    # Affichage du pourcentage
    cat(sprintf("Progression iteration : %d\r", i))
    
  }
  
  return(list(Xopt=Xopt, Copt = Copt))
}

