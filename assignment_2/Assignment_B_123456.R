
# Clear workspace
rm(list=ls())

# Load necessary packages and files
library(ggplot2)
library(tseries)
library(scales)
library(xtable)
library(ggfortify)
library(rugarch)
library(rmgarch)


# --------------------------------Question 1------------------------------------
# Read the data
msft_data <- scan("MSFT.txt")
bac_data <- scan("BAC.txt")
xom_data <- scan("XOM.txt")
sp500_data <- scan("market.txt")

# Compute the log-returns
msft_log_ret <- 100 * diff(log(msft_data))
bac_log_ret  <- 100 * diff(log(bac_data))
xom_log_ret  <- 100 * diff(log(xom_data))
sp500_log_ret <- 100 * diff(log(sp500_data))


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

# --------------------------------Question 5------------------------------------
get_ccc_beta <- function(market_ret, asset_ret) {
  llik_fun_GARCH <- function(par, x) {
    omega <- exp(par[1])
    alpha <- exp(par[2]) / (1 + exp(par[2]))
    beta <- exp(par[3]) / (1 + exp(par[3]))
    T <- length(x)
    sigma2 <- numeric(T)
    sigma2[1] <- var(x)
    for (t in 2:T) {
      sigma2[t] <- omega + alpha * x[t - 1]^2 + beta * sigma2[t - 1]
    }
    ll <- -0.5 * sum(log(2 * pi) + log(sigma2) + x^2 / sigma2)
    return(-ll)
  }
  
  x <- cbind(market_ret, asset_ret)
  alpha_ini <- 0.2; beta_ini <- 0.6
  omega_ini <- var(x[,1])*(1-alpha_ini-beta_ini)
  par_ini <- c(log(omega_ini), log(alpha_ini/(1-alpha_ini)), log(beta_ini/(1-beta_ini)))
  
  # Market
  est1 <- optim(par=par_ini, fn=function(par) llik_fun_GARCH(par, x[,1]))
  p1 <- est1$par
  omega1 <- exp(p1[1]); alpha1 <- exp(p1[2]) / (1 + exp(p1[2])); beta1 <- exp(p1[3]) / (1 + exp(p1[3]))
  
  # Asset
  est2 <- optim(par=par_ini, fn=function(par) llik_fun_GARCH(par, x[,2]))
  p2 <- est2$par
  omega2 <- exp(p2[1]); alpha2 <- exp(p2[2]) / (1 + exp(p2[2])); beta2 <- exp(p2[3]) / (1 + exp(p2[3]))
  
  n <- nrow(x); s1 <- rep(0,n); s2 <- rep(0,n)
  s1[1] <- var(x[,1]); s2[1] <- var(x[,2])
  for (t in 2:n) {
    s1[t] <- omega1 + alpha1 * x[t-1,1]^2 + beta1 * s1[t-1]
    s2[t] <- omega2 + alpha2 * x[t-1,2]^2 + beta2 * s2[t-1]
  }
  
  e1 <- x[,1] / sqrt(s1); e2 <- x[,2] / sqrt(s2)
  rho <- cor(e1, e2)
  beta_t <- rho * sqrt(s2) / sqrt(s1)
  return(beta_t)
}

beta_ccc_msft <- get_ccc_beta(log_returns$MARKET, log_returns$MSFT)
beta_ccc_bac  <- get_ccc_beta(log_returns$MARKET, log_returns$BAC)
beta_ccc_xom  <- get_ccc_beta(log_returns$MARKET, log_returns$XOM)

beta_dynamic_msft <- msft_result$beta
beta_dynamic_bac  <- bac_result$beta
beta_dynamic_xom  <- xom_result$beta

plot(beta_ccc_msft, type = 'l', col = "blue", lwd = 2,
     main = "MSFT: CCC-GARCH vs Q3 Beta_t",
     ylab = expression(beta[t]), xlab = "Time")
lines(beta_dynamic_msft, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("CCC-GARCH", "Beta_t Q3"), col = c("blue", "red"),
       lty = c(1,2), lwd = 2, bty = "n")

plot(beta_ccc_bac, type = 'l', col = "blue", lwd = 2,
     main = "BAC: CCC-GARCH vs Q3 Beta_t",
     ylab = expression(beta[t]), xlab = "Time")
lines(beta_dynamic_bac, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("CCC-GARCH", "Beta_t Q3"), col = c("blue", "red"),
       lty = c(1,2), lwd = 2, bty = "n")

plot(beta_ccc_xom, type = 'l', col = "blue", lwd = 2,
     main = "XOM: CCC-GARCH vs Q3 Beta_t",
     ylab = expression(beta[t]), xlab = "Time")
lines(beta_dynamic_xom, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("CCC-GARCH", "Beta_t Q3"), col = c("blue", "red"),
       lty = c(1,2), lwd = 2, bty = "n")


# ------------------ Question 6:-----------------------------------------------
objective_function <- function(par, x, y, e, sample_m) {
  
  # simulation moments
  sim_m <- simulate_summary_moments(e, x, par)
  weights <- 1/(sample_m^2 + 1e-6)
  return(mean(weights*(sim_m - sample_m)^2))
}


simulate_summary_moments <- function(e, x, par) {
  alpha0 <- par[1]
  alpha1 <- exp(par[2])/(1+exp(par[2]))
  sigma_eta <- exp(par[3])
  sigma_eps <- exp(par[4])
  
  n <- length(x)
  H <- nrow(e)
  M <- H/n
  
  eta <- sqrt(sigma_eta) * e[,1]
  eps <- sqrt(sigma_eps) * e[,2]
  
  n_moments <- 1 + 1 + 5 + 5 + 5  
  output <- matrix(0, nrow=M, ncol=n_moments)
  
  for(m in 1:M) {
    f <- numeric(n)
    beta <- numeric(n)
    
    # initializaiton
    f[1] <- alpha0/(1-alpha1)
    
    for(t in 2:n) {
      f[t] <- alpha0 + alpha1*f[t-1] + eta[(m-1)*n + t]
    }
    
    beta <- exp(f)
    y_sim <- beta * x + eps[((m-1)*n + 1):(m*n)]
    
    hb <- cov(y_sim,x)/var(x)
    yr <- y_sim - hb*x
    
    acf_yr <- acf(yr, lag.max=5, plot=FALSE)$acf[-1]
    
    xy <- yr * x
    acf_xy <- acf(xy, lag.max=5, plot=FALSE)$acf[-1]
    
    acf_y2 <- acf(y_sim^2, lag.max=5, plot=FALSE)$acf[-1]
    
    output[m,] <- c(var(yr), hb, acf_yr, acf_xy, acf_y2)
  }
  
  return(colMeans(output))
}

# the estimation function
estimate_indirect <- function(x, y) {
  n <- length(x)
  M <- 20
  H <- M*n
  set.seed(123)
  e <- cbind(rnorm(H), rnorm(H))
  
  hb <- cov(y,x)/var(x)
  yr <- y - hb*x
  xy <- yr * x
  
  acf_yr <- acf(yr, lag.max=5, plot=FALSE)$acf[-1]
  acf_xy <- acf(xy, lag.max=5, plot=FALSE)$acf[-1]
  acf_y2 <- acf(y^2, lag.max=5, plot=FALSE)$acf[-1]
  
  sample_m <- c(var(yr), hb, acf_yr, acf_xy, acf_y2)
  
  hb <- max(cov(y,x)/var(x), 0.1)  
  a1_ini <- 0.9
  a0_ini <- log(hb)*(1-a1_ini)  
  s_eta_ini <- log(0.01)  
  s_eps_ini <- log(var(yr))
  
  par_ini <- c(a0_ini, log(a1_ini/(1-a1_ini)), s_eta_ini, s_eps_ini)
  
  # constrained optimization
  lower_bounds <- c(-5, -5, -10, -10)
  upper_bounds <- c(5, 5, 5, 5)
  
  est <- optim(par=par_ini,
               fn=function(p) objective_function(p,x,y,e,sample_m),
               method="L-BFGS-B",
               lower=lower_bounds,
               upper=upper_bounds,
               control=list(maxit=10000))
  
  alpha0_hat <- est$par[1]
  alpha1_hat <- exp(est$par[2])/(1+exp(est$par[2]))
  sigma_eta_hat <- exp(est$par[3])
  sigma_eps_hat <- exp(est$par[4])
  
  return(c(alpha0_hat, alpha1_hat, sigma_eta_hat, sigma_eps_hat))
}

params_msft <- estimate_indirect(log_returns$MARKET, log_returns$MSFT)
params_bac <- estimate_indirect(log_returns$MARKET, log_returns$BAC)
params_xom <- estimate_indirect(log_returns$MARKET, log_returns$XOM)

param_df <- data.frame(
  Asset = c("MSFT", "BAC", "XOM"),
  alpha0 = round(c(params_msft[1], params_bac[1], params_xom[1]), 4),
  alpha1 = round(c(params_msft[2], params_bac[2], params_xom[2]), 4),
  sigma_eta = round(c(params_msft[3], params_bac[3], params_xom[3]), 6),
  sigma_eps = round(c(params_msft[4], params_bac[4], params_xom[4]), 4)
)

print(param_df)