
# Clear workspace
rm(list=ls())

# Load necessary packages and files
library(ggplot2)
library(tseries)
library(scales)
library(xtable)
library(ggfortify)


# --------------------------------Question 1------------------------------------
# Read the data
msft_data <- scan("MSFT.txt")
bac_data <- scan("BAC.txt")
xom_data <- scan("XOM.txt")
sp500_data <- scan("market.txt")

# Compute the log-returns
msft_log_ret <- diff(log(msft_data))
bac_log_ret <- diff(log(bac_data))
xom_log_ret <- diff(log(xom_data))
sp500_log_ret <- diff(log(sp500_data))


log_returns <- data.frame(
  MSFT = msft_log_ret,
  BAC = bac_log_ret,
  XOM = xom_log_ret,
  MARKET = sp500_log_ret
)
log_returns <- na.omit(log_returns)
# MSFT
ggplot(log_returns, aes(x = 1:nrow(log_returns), y = MSFT)) +
  geom_line(color = "blue") +
  labs(title = "Log-Returns of MSFT", x = "Time", y = "Log-Return") +
  theme_minimal()
#ggsave("msft_log_return.png", width = 8, height = 5, dpi = 300)
# BAC
ggplot(log_returns, aes(x = 1:nrow(log_returns), y = BAC)) +
  geom_line(color = "darkred") +
  labs(title = "Log-Returns of BAC", x = "Time", y = "Log-Return") +
  theme_minimal()
#ggsave("bac_log_return.png", width = 8, height = 5, dpi = 300)
# XOM
ggplot(log_returns, aes(x = 1:nrow(log_returns), y = XOM)) +
  geom_line(color = "darkgreen") +
  labs(title = "Log-Returns of XOM", x = "Time", y = "Log-Return") +
  theme_minimal()
#ggsave("xom_log_return.png", width = 8, height = 5, dpi = 300)
# Market
ggplot(log_returns, aes(x = 1:nrow(log_returns), y = MARKET)) +
  geom_line(color = "purple") +
  labs(title = "Log-Returns of Market", x = "Time", y = "Log-Return") +
  theme_minimal()
#ggsave("market_log_return.png", width = 8, height = 5, dpi = 300)


# Estimate beta of CAPM model for each asset
model1 <- lm(MSFT ~ MARKET, data = log_returns)
model2 <- lm(BAC ~ MARKET, data = log_returns)
model3 <- lm(XOM ~ MARKET, data = log_returns)

cat("Beta MSFT:", coef(model1)[2], "\n")
cat("Beta BAC:", coef(model2)[2], "\n")
cat("Beta XOM:", coef(model3)[2], "\n")


# --------------------------------Question 3------------------------------------
source("llik_OD_regression.R")

estimate_OD_regression <- function(xt, yt){
  # Initialization
  alpha <- 0.1/sd(xt*yt)
  phi <- 0.8
  omega <- (cov(xt,yt)/var(xt))*(1-phi)
  sig2 <- var(yt)
  
  par_init <- c(omega, log(phi/(1-phi)), log(alpha), log(sig2))
  
  # Optimization
  est <- optim(par = par_init,
               fn = function(par) -llik_OD_regression(yt, xt, par),
               method = "BFGS",
               control = list(maxit = 10000))
  
  # Obtain estimates
  omega_hat <- est$par[1]
  phi_hat <- exp(est$par[2])/(1+exp(est$par[2]))
  alpha_hat <- exp(est$par[3])
  sigma2_hat <- exp(est$par[4])

  theta_hat <- c(omega_hat, phi_hat, alpha_hat, sigma2_hat)
  cat("The parameter estimates are:", round(theta_hat,4))
  
  # Filtered beta
  n <- length(xt)
  beta <- rep(0,n)
  beta[1] <- omega_hat/(1-phi_hat)

  for(t in 2:n){
    beta[t] <- omega_hat + phi_hat*beta[t-1] + alpha_hat*(yt[t-1] - beta[t-1]*xt[t-1])*xt[t-1];
  }
  
  # Plot estimated beta
  plot(beta, type = "l", col = "royalblue", main = "βt", ylab = "", xlab = "")
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  
  return(list(theta_hat = theta_hat, beta = beta))
}


msft_result <- estimate_OD_regression(log_returns$MARKET, log_returns$MSFT)

bac_result <- estimate_OD_regression(log_returns$MARKET, log_returns$BAC)

xom_result <- estimate_OD_regression(log_returns$MARKET, log_returns$XOM)


# --------------------------------Question 4------------------------------------
compute_beta_tplus1 <- function(result, x, y){
  omega <- result$theta_hat[1]
  phi <- result$theta_hat[2]
  alpha <- result$theta_hat[3]
  beta_t <- tail(result$beta, 1)
  x_t <- tail(x, 1)
  y_t <- tail(y, 1)
  
  beta_tplus1 <- omega + phi*beta_t + alpha*(y_t - beta_t*x_t)*x_t
  
  return(beta_tplus1)
}

msft_beta_tplus1 <- compute_beta_tplus1(msft_result, log_returns$MARKET, log_returns$MSFT)
bac_beta_tplus1 <- compute_beta_tplus1(bac_result, log_returns$MARKET, log_returns$BAC)
xom_beta_tplus1 <- compute_beta_tplus1(xom_result, log_returns$MARKET, log_returns$XOM)


next_beta_df <- data.frame(
  Asset = c("MSFT", "BAC", "XOM"),
  Beta_Tplus1 = c(msft_beta_tplus1, bac_beta_tplus1, xom_beta_tplus1)
)

print(next_beta_df)



current_beta_df <- data.frame(
  Asset = c("MSFT", "BAC", "XOM"),
  Beta_T = c(tail(msft_result$beta,1), tail(bac_result$beta,1), tail(xom_result$beta,1))
)

print(current_beta_df)

