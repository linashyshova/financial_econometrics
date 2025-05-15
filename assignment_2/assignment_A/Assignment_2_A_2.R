# Plot filtered volatility of a GARCH(1,1) model for the the log-returns of the S&P 500
# in order to compare it with filtered volatility for SV

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(moments)
library(ggplot2)
library(tseries)
library(scales)
source("assignment_2/assignment_A/llik_fun_GARCH_pq.R") 

########################
############ 1. Obtain S&P500 weekly log-returns
########################

#The file market.txt contains weekly prices of the S&P500

sep_prices <- scan("assignment_2/data/market.txt")
x <- diff(log(sep_prices))

# Choose the initial parameter values for the Newton Raphson optimization
a <- 0.1 # initial value for alpha1
b <- 0.6  # initial value for beta
omega <- var(x)*(1-a-b) # initial value for omega

# Transform intitial values using the inverse of the link funtions
par_ini <- c(
  log(omega),
  log(a / (1 - a)),
  log(b / (1 - b))
)

# Optimize log-likelihood function
est <- optim(
  par = par_ini,
  fn = function(par) -llik_fun_GARCH_pq(x, par=par, p = 1, q = 1),
  method = "BFGS",
  control = list(maxit = 100000)
)

# Obtain parameter estimates
omega_hat <- exp(est$par[1])
alpha_hat <- exp(est$par[2])/(1+exp(est$par[2]))
beta_hat <- exp(est$par[3])/(1+exp(est$par[3]))

theta_hat <- c(omega_hat,alpha_hat,beta_hat)
cat("The parameter estimates are:")
round(theta_hat,6)

# Display the log-likelihood. 
cat("The log-likelihood value is:")
-est$value*length(x)

# Display the exit flag to see if convergence of the algorithm has been attained
cat("Exit flag:")
est$convergence # zero indicates succesfull optimization

# Calculate and plot filtered conditional variance
n <- length(x)
sigma2 <- rep(0,n)
sigma2[1:2] <- var(x)

for(t in 3:n){
  sigma2[t] = omega_hat + alpha_hat*x[t-1]^2 + beta_hat*sigma2[t-1]
}

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(x,type="l", main="Time series",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(sigma2,type="l",col=2, main="Filtered sigmat GARCH(1,1)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")