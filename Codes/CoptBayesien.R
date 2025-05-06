# Criteres bayesiens basés sur la densite a posteriori
CoptBayes <- function(X, model, Dx, Dtheta, sigeps,  dprior, rprior, 
                      L=100, K=100, typeCopt='KL', type='SK',...){
  
  if(!is.matrix(X)) stop("X doit etre une matrice ! ")
  
  p <- nrow(Dtheta)
  d <- nrow(Dx)
  n <- nrow(X)
  
  TT = rprior(L)
  THETA = Unif(K, Dtheta)
  Eps <- t(SimDesign::rmvnorm(K, mean=rep(0,n), sigma=diag(1,n)))
  
  CC <- rep(0,L)
  for(l in 1:L){
    R = predict(model, newdata=Dbind(X,TT[l,]), type=type, 
                cov.compute = TRUE, se.compute = FALSE, checkNames = FALSE)
    Mpost = R$mean
    Cpost = t(chol(R$cov + diag(sigeps, n)))
    RR <- rep(0, K)
    for(k in 1:K){
      Ysim = Mpost +  Cpost%*%matrix(Eps[, k], ncol = 1) 
      W <- Poidpost(THETA, Ysim, model,  X, sigeps, Dtheta, dprior, type)
      
      if (typeCopt=="SOV"){
        thetaD =  t(W)%*%THETA
        DD <- THETA - matrix(rep(thetaD, K), nrow=K, ncol=p)
        AA <- as.vector(W)*as.vector(rowSums(DD^2))
      }
      
      if (typeCopt=="MSE"){
        DD <- THETA - matrix(rep(TT[l,],K), nrow=K, ncol=p)
        AA <- as.vector(W)*as.vector(rowSums(DD^2))
      }
      
      if (typeCopt=="KL"){
        PP <- apply(THETA, 1, function(theta){return(dprior(theta))})
        PP = PP/sum(PP)
        AA <- W*log(W+1e-16) - W*log(PP+1e-16)
      }
      RR[k] <- sum(AA)
    }
    CC[l] <- sum(RR)
  }
  R = mean(CC)
  return(R)
}

# Divergence de KL entre pi_0 et pi_post (cas gaussien)
KLprior_post <- function(muprior, Sigmaprior, mupost, Sigmapost) {
  p <- ncol(Sigmaprior)
  A <- log(1e-16+ det(Sigmapost) / det(Sigmaprior))
  B <- sum(diag(solve(Sigmapost) %*% Sigmaprior))
  diff <- matrix(mupost - muprior, ncol = 1)
  C <- t(diff) %*% solve(Sigmapost) %*% diff
  R <- 0.5 * (A - p + B + C)
  return(as.numeric(R))
}

# Criteres d'optimalités bayesienne (code de calcul lineaire)
CoptBayesLin <- function(X, model, Dx, Dtheta,  thetaprior, Sigmaprior, sigeps,
                         L=100, typeCopt='KL', type='SK',...){
  p <- nrow(Dtheta)
  d <- nrow(Dx)
  n <- nrow(X)
  alpha = 1e-5
  R = predict(model, newdata=Dbind(X,thetaprior+alpha), type=type, 
              cov.compute = FALSE, se.compute = FALSE, checkNames = FALSE)
  M1 = R$mean
  R = predict(model, newdata=Dbind(X,thetaprior-alpha), type=type, 
              cov.compute = FALSE, se.compute = FALSE, checkNames = FALSE)
  M2 = R$mean
  H = (M1-M2)/(2*alpha)
  
  SigmapriorInv = solve(Sigmaprior)
  Minf = t(H)%*%H
  if(typeCopt=="Det"){R = det((1/sigeps)*Minf + SigmapriorInv)}
  
  if(typeCopt=="Tr"){R = Tr((1/sigeps)*Minf + SigmapriorInv)}
  
  if(typeCopt=="KL"){
    KLvec = rep(0, L)
    TT = t(SimDesign::rmvnorm(L, mean=thetaprior, sigma=Sigmaprior))
    Eps <- t(SimDesign::rmvnorm(L, mean=rep(0,n), sigma=diag(1,n)))
    for(l in 1:L){
      R = predict(model, newdata=Dbind(X,TT[l,]), type=type, 
                  cov.compute = TRUE, se.compute = FALSE, checkNames = FALSE)
      Mpost = R$mean
      Cpost = t(chol(R$cov + diag(sigeps, n)))
      Ysim = Mpost +  Cpost%*%matrix(Eps[, l], ncol = 1) 
      Sigmapost = solve((1/sigeps)*Minf + SigmapriorInv)
      mupost = Sigmapost%*%((1/sigeps)*t(H)%*%Ysim + SigmapriorInv%*%muprior)
      KLvec[l] = KLprior_post(muprior, Sigmaprior, mupost, Sigmapost)
    }
    R = mean(KLvec)
  }
  
  return(R)
}
