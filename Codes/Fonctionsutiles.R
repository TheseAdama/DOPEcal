LHSD<- function(N,D){
  
  G <- lhs::maximinLHS(N,nrow(D))
  for(i in 1:nrow(D)){
    G[,i] <- (D[i,2]-D[i,1])*G[,i] + D[i,1]
  }
  G = as.matrix(G)
  return(G)
}

Unif <- function(N, D){
  Mat <- matrix(0, nrow = N, ncol = nrow(D))
  for(i in 1:nrow(D)){
    Mat[,i] <- stats::runif(N, D[i,1], D[i,2])
  }
  return(Mat)
}

Dbind <- function(D,x){
  mat <- cbind(D, matrix(rep(x,nrow(D)), ncol=length(x), byrow=TRUE))
  return(mat)
}


Fishermatrix <- function(X, theta0, dftheta, Covmat){
  
  J = matrix(0, ncol=length(theta0), nrow=nrow(X))
  for(i in 1:nrow(X)){
    J[i,] = dftheta(X[i,], theta0)
  }
  R = t(J)%*%solve(Covmat)%*%J
  return(R)
}

Poidpost <- function(THETA, Ysim, model,  X, sigeps, Dtheta, dprior, type){
  d = ncol(X)
  n = nrow(X)
  p = ncol(THETA)
  K <- nrow(THETA)
  w <- rep(0, K)
  for(k in 1:K){
    R = predict(model, newdata=Dbind(X, THETA[k,]), type, cov.compute = TRUE, 
                se.compute = FALSE, checkNames = FALSE)
    E = Ysim - R$mean
    VV = R$cov + diag(sigeps, n)
    pzero = dprior(THETA[k,])
    w[k] <- (det(VV)^{-1/2})*exp(-0.5*(t(E)%*%solve(VV)%*%E))*pzero
  }
  w <- as.vector(w)
  w <- w/sum(w)
  
  return(w)
}


# Generer des voisins
Neighbor <- function(x, Dx, type="uniform",...){
  x <- matrix(x, nrow=1)
  xv = matrix(0, ncol = nrow(Dx) , nrow=1)
  
  if(type=="uniform"){
    for(k in 1:nrow(Dx)){
      m = x[1,k] - (x[1,k]-Dx[k,1])/3
      M = x[1,k] + (Dx[k,2]-x[1,k])/3
      xv[1,k] <- stats::runif(1, min = m, max = M)
    }
  }
  
  if(type=="gaussian"){
    for(k in 1:nrow(Dx)){
      sig = sqrt(0.25*(Dx[k,2]-Dx[k,1]))
      xv[1,k] <- rtruncnorm(1, a=Dx[k,1], b=Dx[k,2], mean=x[1,k], sd=sig)
    }
  }
  
  return(xv)
}

