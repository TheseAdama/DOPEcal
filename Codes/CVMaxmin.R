CVMm <- function(X, model, Dx, Dtheta, L=100, type="SK", alpha=0.5){
  n = nrow(X)
  d = nrow(Dx)
  p = nrow(Dtheta)
  THETA = Unif(L, Dtheta)
  MM = matrix(0, ncol=n, nrow=L)
  for(i in 1:L){
    R =  predict(model,newdata = Dbind(X,THETA[i,]),type = type, cov.compute = FALSE, checkNames = "FALSE")
    MM[i,] = R$mean
  }
  CV = mean((MM-colMeans(MM))**2)^{alpha}
  dMm = min(dist(X, method = "euclidean"))^{1-alpha}
  o = CV*dMm
  return(o)
}
# 
# fcv <- function(X){
#   o = CVMm(X, model, Dx, Dtheta, L=30, type="SK", alpha=0.6)
#   return(o)
# }
# 
# R = ForwardOptim(f=fcv, Dx, nD=20, N=100, tomax=TRUE, echo=TRUE)
# 
# image(x, x, zz, col = viridis(1000), xlab = "x1", ylab = "x2", 
#       main ="code")
# points(R$Xopt, col='red', pch=17, lwd=2)
