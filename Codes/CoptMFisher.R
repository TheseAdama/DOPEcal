CoptMFisher <- function(X, model, sigeps, Dx, Dtheta, dftheta, typeCopt="Det", L=1000, type='SK'){
  
  if(!is.matrix(X)) stop("X doit etre une matrice ! ")
  THETA = Unif(L, Dtheta)
  C <- rep(0, L)
  for(l in 1:L){
    R = predict(model, newdata=Dbind(X,THETA[l,]), type=type, 
                cov.compute = TRUE, se.compute = FALSE, checkNames = FALSE)
    Covmat <- R$cov + diag(sigeps, nrow(X))
    Mat <- Fishermatrix(X, THETA[l,], dftheta, Covmat=Covmat)
    if (typeCopt=="Det")      C[l] <- det(Mat)
    if (typeCopt=="Tr")       C[l] <- sum(diag(Mat))
    if (typeCopt=="DetInv")   C[l] <- det(solve(Mat))
    if (typeCopt=="TrInv")    C[l] <- sum(diag(solve(Mat)))
  }
  
  Cmean <- mean(C)
  return(Cmean)
}
