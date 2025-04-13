# This file contains the average log-likelihood function for GARCH(1,2). 
# The function below takes the data, labeled "x", and the parameter vector, 
# labeled "par", as input and gives the average log-likelihood, labeled "llik", 
# as output. 


llik_fun_GARCH_1_2 <- function(par,x){
  
  n <- length(x)
  
  #set paramter values from the input par using link functions for restrictions
  omega <- exp(par[1])                   #exp() to ensure omega>0
  alpha1 <- exp(par[2])/(1+exp(par[2]))   #logistic()=exp()/(exp()) for 0<alpha_1<1
  alpha2 <- exp(par[3])/(1+exp(par[3]))   #logistic()=exp()/(exp()) for 0<alpha_2<1
  beta <- exp(par[4])/(1+exp(par[4]))    #logistic()=exp()/(exp()) for 0<beta<1
  # cat("| omega =", omega, "| alpha1 =", alpha1, "| alpha2 =", alpha2, "| beta =", beta, "\n")Ç
  
  if (alpha1 + alpha2 + beta >= 1) return(-2000000)
  
  ## Filter Volatility
  sig2 <- rep(0,n)
  sig2[1:2] <- var(x) #initialize volatility at unconditional variance
  
  for(t in 3:n){
    sig2[t] <- omega + alpha1*x[t-1]^2 + alpha2*x[t-2]^2 + beta*sig2[t-1]
    # cat("t =", t, "|  sig2[t] =",  sig2[t],"\n")
  }
  
  ## Calculate Log Likelihood Values
  
  #construct sequence of log-likelihood contributions
  l <- -(1/2)*log(2*pi) - (1/2)*log(sig2) - (1/2)*x^2/sig2
  
  llik <- sum(l)  # obtain the average log-likelihood
  cat("llik =", llik, "| omega =", omega, "| alpha1 =", alpha1, "| alpha2 =", alpha2, "| beta =", beta, "\n")
  return(llik) # return the average log-likelihood as output
}