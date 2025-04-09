# Recuit simulé 
SAOptim <- function(Xinit, f, Dx, maxiter, c=0.99, Tinit=0.1, Paccept="metropolis", schema="geo"){
  
  n <- nrow(Xinit)
  d <- ncol(Xinit)
  
  Yinit <- f(Xinit)
  temp <- Tinit
  
  Xstar <- Xinit
  Ystar <- Yinit
  k <- 1
  
  while(k <= maxiter){
    i = k%%n+1
    X = Xstar
    X[i,] <- Neighbor(X[i,], Dx, type="uniform")
    Y <- f(X)
    if (is.na(Y)){next}
    Delta <- Ystar - Y
    
    if(Paccept=="metropolis"){
      alpha = min(1,exp(-Delta/temp))
    }else{
      alpha = 1/(1+exp(Delta/temp))
    }
    
    u <- runif(1)
    if(alpha > u){
      Xstar <- X
      Ystar <- Y
    }
    
    if(schema=="geo"){temp <- c*temp}
    if(schema=="log"){temp <- temp/log(1+k)}
    
    # Calcul du pourcentage d'évolution
    rate <- (k /maxiter)* 100
    
    # Affichage du pourcentage
    cat(sprintf("Progression Optimisation: %.2f%%\r", rate))
    
    # Affichage dans la console
    flush.console()
    
    k <- k + 1
  }
  cat("\n")
  
  return(Xstar)
}
