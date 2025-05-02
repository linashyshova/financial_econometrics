llik_CT_sDVECH <- function(par,x){
  
  a <- exp(par[1])/(1+exp(par[1]))
  b <- exp(par[2])/(1+exp(par[2]))
  
  d <- dim(x)
  n <- d[1]
  
  VECHt <- matrix(0,nrow=n,ncol=3)
  llik <- 0
  
  C <- cov(x)
  VECHt[1,] <- c(C[1,1],C[1,2],C[2,2])
  
  for(t in 2:n){
    
    VECHt[t,1] <- C[1,1]*(1-a-b)+b*VECHt[t-1,1]+a*x[t-1,1]^2
    VECHt[t,3] <- C[2,2]*(1-a-b)+b*VECHt[t-1,3]+a*x[t-1,2]^2
    VECHt[t,2] <- C[1,2]*(1-a-b)+b*VECHt[t-1,2]+a*x[t-1,1]*x[t-1,2]
    
    SIGMAt <- cbind(c(VECHt[t,1],VECHt[t,2]),c(VECHt[t,2],VECHt[t,3]))
    
    llik <- llik-0.5*(log(det(SIGMAt))+x[t,]%*%solve(SIGMAt)%*%t(t(x[t,])))/n
  }
  
  return(llik)
  
}