llik_OD_regression <- function(y, x, par) {
  n <- length(x)
  omega <- par[1]
  phi <- exp(par[2])/(1+exp(par[2]))
  alpha <- exp(par[3])
  s2 <- exp(par[4])
  
  beta <- rep(0, n)
  beta[1] <- omega/(1 - phi)
  for (t in 2:n) {
    beta[t] <- omega + phi * beta[t-1] + alpha * (y[t-1] - beta[t-1]*x[t-1]) * x[t-1]
  }
  
  l <- -(1/2)*log(s2) - (1/2)*(y - beta*x)^2 / s2
  llik <- sum(l)
  return(llik)
}