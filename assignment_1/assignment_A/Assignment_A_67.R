### Empirical Finance assignment 1

rm(list = ls())
graphics.off()
cat("\014")

prices_HSI <- as.numeric(readLines("assignment_1/data/HSI.txt"))


### Q1a.7
log_returns_HSI <- log(prices_HSI[2:length(prices_HSI)]) - log(prices_HSI[1:(length(prices_HSI)-1)])


# GARCH(1,1) log-likelihood function
llik_fun_GARCH <- function(par, log_returns_HSI) {
  n <- length(log_returns_HSI)
  omega <- exp(par[1])
  alpha <- exp(par[2]) / (1 + exp(par[2]))
  beta  <- exp(par[3]) / (1 + exp(par[3]))
  
  sig2 <- rep(0, n)
  sig2[1] <- var(log_returns_HSI)
  
  for (t in 2:n) {
    sig2[t] <- omega + alpha * log_returns_HSI[t - 1]^2 + beta * sig2[t - 1]
  }
  
  l <- -0.5 * log(2 * pi) - 0.5 * log(sig2) - 0.5 * log_returns_HSI^2 / sig2
  llik <- mean(l)
  
  return(llik)
  
}

a <- 0.2
b <- 0.6
init_omega <- var(log_returns_HSI) * (1-a-b)
init_params <- c(log(init_omega), log(a/(1-a)), log(b/(1-b)))
est_GARCH <- optim(init_params, fn=function(par) - llik_fun_GARCH(par, log_returns_HSI), method = "BFGS")

theta_hat <- est_GARCH$par
omega_hat <- exp(theta_hat[1])
alpha_hat <- exp(theta_hat[2]) / (1 + exp(theta_hat[2]))
beta_hat  <- exp(theta_hat[3]) / (1 + exp(theta_hat[3]))


cat("Estimated GARCH(1,1) parameters:\n")
cat("omega =", omega_hat, "\n")
cat("alpha =", alpha_hat, "\n")
cat("beta  =", beta_hat, "\n")
GARCH_llikelihood <- length(log_returns_HSI) * llik_fun_GARCH(theta_hat, log_returns_HSI)
cat("\nLog-likelihood GARCH(1,1):", GARCH_llikelihood, "\n")
# AIC and BIC
AIC_GARCH11 <- -2 * GARCH_llikelihood + 2 * 3 #3 parameters       
BIC_GARCH11   <- -2 * GARCH_llikelihood + log(length(log_returns_HSI)) * 3
cat("\nAIC_GARCH(1,1):", AIC_GARCH11, "\n")
cat("BIC_GARCH(1,1):", BIC_GARCH11, "\n")


# GJR-GARCH(1,1) 
llik_fun_GJR_GARCH <- function(par, log_returns_HSI) {
  n <- length(log_returns_HSI)
  omega <- exp(par[1])
  alpha <- exp(par[2]) / (1 + exp(par[2]))
  beta  <- exp(par[3]) / (1 + exp(par[3]))
  gamma <- exp(par[4]) / (1 + exp(par[4]))
  
  sig2 <- rep(0, n)
  sig2[1] <- var(log_returns_HSI)
  
  for (t in 2:n) {
    indicator <- ifelse(log_returns_HSI[t - 1] < 0, 1, 0)
    sig2[t] <- omega + alpha * log_returns_HSI[t - 1]^2 + beta * sig2[t - 1] + gamma * indicator * log_returns_HSI[t - 1]^2
  }
  
  l <- -0.5 * log(2 * pi) - 0.5 * log(sig2) - 0.5 * log_returns_HSI^2 / sig2
  llik <- mean(l)
  
  return(llik)
}
a <- 0.2
b <- 0.6
g <- 0.05 
init_omega <- var(log_returns_HSI) * (1 -a-b - 0.5*g)
init_params <- c(log(init_omega), log(a/(1-a)), log(b/(1-b)), log(g/(1-g)))

# Optimization
est_GJR <- optim(
  par = init_params,
  fn = function(par) - llik_fun_GJR_GARCH(par, log_returns_HSI),
  method = "BFGS"
)

theta_hat <- est_GJR$par

omega_hat <- exp(theta_hat[1])
alpha_hat <- exp(theta_hat[2]) / (1 + exp(theta_hat[2]))
beta_hat  <- exp(theta_hat[3]) / (1 + exp(theta_hat[3]))
gamma_hat <- exp(theta_hat[4]) / (1 + exp(theta_hat[4]))

cat("Estimated GJR-GARCH(1,1) parameters:\n")
cat("omega =", omega_hat, "\n")
cat("alpha =", alpha_hat, "\n")
cat("beta  =", beta_hat, "\n")
cat("gamma =", gamma_hat, "\n")
GJR_llikelihood <- length(log_returns_HSI) * llik_fun_GJR_GARCH(theta_hat, log_returns_HSI)
cat("\nLog-likelihood:", GJR_llikelihood, "\n")
# AIC and BIC
AIC_GJR <- -2 * GJR_llikelihood + 2 * 4 #4 parameters with gamma       
BIC_GJR   <- -2 * GJR_llikelihood + log(length(log_returns_HSI)) * 4
cat("\nAIC_GJR:", AIC_GJR, "\n")
cat("BIC_GJR:", BIC_GJR, "\n")



