# This R file contains code for the estimation of the 
# bivariate sDVECH(1,1) with covariance targeting
# Stock returns of HSI and IBEX are considered

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################
install.packages("yfR")
library(yfR)
source("assignment_1/assignment_B/llik_CT_sDVECH.R") 

########################
############ 1. Read data
########################

hsi <- read.table("assignment_1/data/HSI.txt")
p_hsi <- hsi[[1]]

ibex <- read.table("assignment_1/data/IBEX.txt")
p_ibex = ibex[[1]]

# obtain log-returns for HSI
r_hsi <- diff(log(p_hsi))*100
# obtain log-returns for IBEX
r_ibex <- diff(log(p_ibex))*100

#combine the two series of returns in a single matrix
x <- cbind(r_hsi,r_ibex)

#obtain the sample size
n <- length(r_hsi)

########################
############ 2. plot prices and log-returns
########################

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(p_hsi,type="l",main = "Prices HSI",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(r_hsi,type="l",main = "log-returns HSI",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(p_ibex,type="l",main = "Prices IBEX",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(r_ibex,type="l",main = "Log-returns IBEX",ylab="",xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ 3. Initial paramter value for optimization
########################

alpha_ini <- 0.2 
beta_ini <- 0.6  
par_ini <- c(log(alpha_ini/(1-alpha_ini)),log(beta_ini/(1-beta_ini)))

########################
############ 4. Optimize the log-likelihood function llik_CT_sDVECH()
########################

est <- optim(par=par_ini,fn=function(par)-llik_CT_sDVECH(par,x), method = "BFGS")

########################
############ 5. Display estimation results 
########################

# parameter values
(a_hat <- exp(est$par[1])/(1+exp(est$par[1])))
(b_hat <- exp(est$par[2])/(1+exp(est$par[2])))

cat("log likelihood value:")
-est$value*n


########################
############ 6. Obtain the estimated conditional covariance matrix
########################

## Create nx3 matrix
VECHt <- matrix(0,nrow=n,ncol=3)

# set the initial value of the conditional variance equal to the sample
# covariance

C <- cov(x)
VECHt[1,] <- c(C[1,1],C[1,2],C[2,2])

for(t in 2:n){
  VECHt[t,1] <- C[1,1]*(1-a_hat-b_hat)+b_hat*VECHt[t-1,1]+a_hat*x[t-1,1]^2
  VECHt[t,3] <- C[2,2]*(1-a_hat-b_hat)+b_hat*VECHt[t-1,3]+a_hat*x[t-1,2]^2
  VECHt[t,2] <- C[1,2]*(1-a_hat-b_hat)+b_hat*VECHt[t-1,2]+a_hat*x[t-1,1]*x[t-1,2]
}

########################
############ 7. Plot the estimated conditional variances, covariance and correlation.
########################

sd1t <- sqrt(VECHt[,1])
sd2t <- sqrt(VECHt[,3])
corrt <- VECHt[,2]/(sd1t*sd2t)

var_hsi <- VECHt[,1] 
var_ibex <- VECHt[,3]
cov_hsi_ibex <- VECHt[,2]

par(mfrow=c(2,2),mar=c(4.1,4.1,1.1,2.1))
plot(var_hsi, type="l", main = "Conditional variance HSI", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(var_ibex, type="l", main = "Conditional variance IBEX", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(cov_hsi_ibex, type="l", main = "Conditional covariance", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
plot(corrt, type="l", main = "Conditional correlation", ylab="", xlab="")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################
############ 8. Obtain the conditional variance and the α-VaR at 1% level for the
############    portfolio of the bank (30% in HSI and 70% in IBEX)
########################

# Portfolio weights: 30% HSI, 70% IBEX
w <- c(0.3, 0.7)

port_var <- numeric(n)
VaR_1pct <- numeric(n)

# Loop over time to compute conditional variance and VaR
for (t in 1:n) {
  Ht <- matrix(c(VECHt[t,1], VECHt[t,2], VECHt[t,2], VECHt[t,3]), nrow=2)
  port_var[t] <- t(w) %*% Ht %*% w
  VaR_1pct[t] <- -2.326 * sqrt(port_var[t])  # Normal 1% quantile ≈ -2.326
}

# Plot results
par(mfrow=c(2,1),mar=c(4.1,4.1,1.1,2.1))
plot(sqrt(port_var), type="l", main="Conditional portfolio volatility", ylab="", xlab="")
grid()
plot(VaR_1pct, type="l", main="1% VaR for the portfolio", ylab="", xlab="")
grid()
