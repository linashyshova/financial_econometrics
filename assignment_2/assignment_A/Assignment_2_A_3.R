### Estimate an SV model with R with a different set of auxiliary statistics

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

library(moments)

### This file contains the code to obtain simulated moments from an SV model
# The input are parameter values and the error vectors et and vt
# For more details see the Lecture Notes.

sim_m_SV <- function(e,par){
  
  omega <- par[1]      
  beta <- exp(par[2])/(1+exp(par[2]))
  sig2f <- exp(par[3]) 
  
  H <- length(e[,1])
  
  epsilon <- e[,1] 
  eta <- sqrt(sig2f)*e[,2]
  
  x <- rep(0,H) 
  f <- rep(0,H) 
  
  f[1] <- omega/(1-beta)
  x[1] = exp(f[1]/2) * epsilon[1]
  
  for(t in 2:H){
    f[t] <- omega + beta * f[t-1] + eta[t]  # state equation
    x[t] <- exp(f[t]/2) * epsilon[t]       # observation equation
  }
  
  xa <- abs(x) 
  acv_15 <- acf(xa, lag.max = 15, type = "covariance", plot = F)$acf
  
  output <- c(mean(xa), acv_15)
  return(output)
}

#This R code contains the function that needs to be minimized,
#with respect to ft at each time period to obtain an approximation for
#the filtered volatility of the SV model. 
#For a more detailed explanation see the Lecture Notes.

filter_SV <- function(yt,ft,ft1,theta){
  
  omega <- theta[1]   
  beta <- theta[2]   
  sig2f <- theta[3]  
  
  output <- yt^2*exp(-ft)+3*ft+(ft-omega-beta*ft1)^2/sig2f
  return(output)
}

########################
############ 1. Obtain S&P500 weekly log-returns
########################

#The file market.txt contains weekly prices of the S&P500

sep_prices <- scan("assignment_2/data/market.txt")
x <- diff(log(sep_prices))

########################
############ 2. Obtain sample moments
########################

n <- length(x)
xa <- abs(x)
acv_15 <- acf(xa, lag.max = 15, type = "covariance", plot = F)$acf
sample_m <- c(mean(xa), acv_15)

########################
############ 3. Generate errors for the simulations from the SV model
########################

set.seed(123)
H <- 30*n
epsilon <- rnorm(H)  
eta <- rnorm(H)  
e <- cbind(epsilon,eta)

########################
############ 4. Choose initial parameter values for the optimization
########################

#choose the initail parameter values for the numerical optimization

b <- 0.90
sig2f <- 0.1
omega <- log(var(x))*(1-b)

par_ini <- c(omega,log(b/(1-b)),log(sig2f))

########################
############ 5. Obtain parameter estimates
########################

# 3a. Minimize criterion function

est <- optim(
  par=par_ini,
  fn=function(par) mean((sim_m_SV(e,par)-sample_m)^2),
  method = "BFGS",
  control = list(maxit = 100000)
)

# Display the exit flag to see if convergence of the algorithm has been attained
cat("Exit flag:")
est$convergence

# 3a. Obtain parameter extimate using the link functions

omega_hat <- est$par[1]
beta_hat <- exp(est$par[2])/(1+exp(est$par[2]))
sig2f_hat <- exp(est$par[3])

theta_hat <- c(omega_hat,beta_hat,sig2f_hat)
cat("The parameter estimates are:")
round(theta_hat,4)


########################
############ 6. Obtain filtered volatility
########################

f <- rep(0,n)
f[1] <- log(var(x))

#Run a for loop minizing the function filter_SV() at each time period to
#obtain the filtered estimate of ft 

for(t in 2:n){ # start recursion from t=2 to t=T
  ft_ini <- f[t-1]
  f_est <- optim(par=ft_ini ,fn= function(ft) filter_SV(x[t],ft,f[t-1],theta_hat), method = "BFGS")
  f[t] <- f_est$par
}

########################
############ 7. Plot series and estimated volatility
########################

par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(x,type="l", main="Time series",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(exp(f),type="l",col=2, main="Filtered sigmat SV (different auxiliary statistics)",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
