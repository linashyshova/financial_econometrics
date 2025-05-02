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
library(ggplot2)
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
############ 2. Plot prices and log-returns
########################

# Create date vectors
start_date <- as.Date("1997-01-01")
date_prices <- seq(from = start_date, by = "week", length.out = length(p_hsi))
date_returns <- seq(from = start_date + 7, by = "week", length.out = length(r_hsi))

# Function to plot with years on x axis
plot_series <- function(dates, values, title, ylab) {
  df <- data.frame(date = dates, value = values)
  ggplot(df, aes(x = date, y = value)) +
    geom_line(color = "blue") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = "Year", y = ylab, title = title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

plot_series(date_prices, p_hsi,  "HSI Index Price - Weekly", "HSI Price")
plot_series(date_returns, r_hsi, "HSI Log-Returns - Weekly", "HSI Return")
plot_series(date_prices, p_ibex, "IBEX Index Price - Weekly", "IBEX Price")
plot_series(date_returns, r_ibex,"IBEX Log-Returns - Weekly", "IBEX Return")

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

plot_series(date_returns, var_hsi,  "Conditional variance HSI", "HSI Variance")
plot_series(date_returns, var_ibex, "Conditional variance IBEX", "IBEX Variance")
plot_series(date_returns, cov_hsi_ibex, "Conditional covariance", "Conditional covariance")
plot_series(date_returns, corrt,"Conditional correlation", "Conditional correlation")

########################
############ 8. Obtain the conditional variance and the α-VaR at 1% level for the
############    portfolio of the bank (30% in HSI and 70% in IBEX)
########################

# Portfolio weights: 30% HSI, 70% IBEX
w1 <- 0.3
w2 <- 0.7

# Initialize vectors to store results
port_var <- numeric(n)     # Array for conditional portfolio variance
VaR <- numeric(n)     # Array for 1% VaR
z_value = qnorm(0.99)

for (t in 1:n) {
  var_hsi  <- VECHt[t, 1]  
  covar    <- VECHt[t, 2]   
  var_ibex <- VECHt[t, 3] 
  port_var[t] <- w1^2 * var_hsi + w2^2 * var_ibex + 2 * w1 * w2 * covar
  VaR[t] <- z_value * sqrt(port_var[t])
}

# Plot portfolio conditional volatility and VaR
plot_series(date_returns, port_var, "Conditional portfolio variance", "Conditional variance")
plot_series(date_returns, VaR, "1% Value-at-Risk of the portfolio", "Loss %")

########################
############ 9. Forecast portfolio volatility for 52 weeks ahead
########################

C_vech <- c(C[1,1], C[1,2], C[2,2])

# Initialize forecast matrix
sigma_forecast <- matrix(0, nrow = 52, ncol = 3)
sigma_forecast[1, ] <- C_vech * (1 - a_hat - b_hat) + b_hat * VECHt[n, ]  # h = 1

# Forecast h = 2 to 52
for (h in 2:52) {
  sigma_forecast[h, ] <- C_vech * (1 - a_hat - b_hat) + b_hat * sigma_forecast[h - 1, ]
}

# compute portfolio variance forecast: w' H w
w <- c(0.3, 0.7)
port_var_forecast <- numeric(52)
for (h in 1:52) {
  sigma_mat <- matrix(c(sigma_forecast[h,1], sigma_forecast[h,2], sigma_forecast[h,2], sigma_forecast[h,3]), nrow = 2)
  port_var_forecast[h] <- t(w) %*% sigma_mat %*% w
}

port_vol_forecast <- sqrt(port_var_forecast)

df_vol <- data.frame(
  week = 1:52,
  volatility = sqrt(port_var_forecast)  # convert variance to volatility
)
plot(port_vol_forecast, type = "l", 
     main = "Forecasted Portfolio Volatility (Next 52 Weeks)", 
     xlab = "Forecast horizon (weeks)", ylab = "Volatility")
grid()


