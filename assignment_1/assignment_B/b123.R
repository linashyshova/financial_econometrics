library(ggplot2)
library(scales)

# --- Load Data ---
hsi <- read.table("HSI.txt")
ibex <- read.table("IBEX.txt")

hsi_ret <- diff(log(hsi$V1))
ibex_ret <- diff(log(ibex$V1))
T <- length(hsi_ret)

# --- GARCH(1,1) Log-Likelihood Function ---
garch_loglik <- function(params, data) {
  omega <- abs(params[1])
  alpha <- abs(params[2])
  beta <- abs(params[3])
  T <- length(data)
  sigma2 <- rep(var(data), T)
  ll <- 0
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * data[t-1]^2 + beta * sigma2[t-1]
    ll <- ll - 0.5 * (log(2 * pi) + log(sigma2[t]) + data[t]^2 / sigma2[t])
  }
  return(-ll)
}

# --- Estimate GARCH(1,1) Models ---
opt_hsi <- optim(c(0.01, 0.05, 0.9), garch_loglik, data = hsi_ret)
opt_ibex <- optim(c(0.01, 0.05, 0.9), garch_loglik, data = ibex_ret)
params_hsi <- abs(opt_hsi$par)
params_ibex <- abs(opt_ibex$par)

# --- Filter Conditional Variances ---
get_sigma2 <- function(params, data) {
  omega <- params[1]; alpha <- params[2]; beta <- params[3]
  sigma2 <- rep(var(data), length(data))
  for (t in 2:length(data)) {
    sigma2[t] <- omega + alpha * data[t-1]^2 + beta * sigma2[t-1]
  }
  return(sigma2)
}

sigma2_hsi <- get_sigma2(params_hsi, hsi_ret)
sigma2_ibex <- get_sigma2(params_ibex, ibex_ret)

# --- Estimate CCC Correlation ---
z_hsi <- hsi_ret / sqrt(sigma2_hsi)
z_ibex <- ibex_ret / sqrt(sigma2_ibex)
rho_hat <- cor(z_hsi, z_ibex)

# --- Construct Conditional Covariance Matrices Ht ---
Ht <- array(NA, c(T, 2, 2))
cov_hsi_ibex <- numeric(T)
for (t in 1:T) {
  D_t <- diag(c(sqrt(sigma2_hsi[t]), sqrt(sigma2_ibex[t])))
  R <- matrix(c(1, rho_hat, rho_hat, 1), 2)
  Ht[t,,] <- D_t %*% R %*% D_t
  cov_hsi_ibex[t] <- Ht[t,1,2]
}

# --- Portfolio Variance and VaR (1%) ---
w <- c(0.3, 0.7)
port_var <- sapply(1:T, function(t) t(w) %*% Ht[t,,] %*% w)
VaR_1pct <- -qnorm(0.01) * sqrt(port_var)

# --- Optimal Sharpe-Ratio Portfolio Weights ---
mu <- c(mean(hsi_ret), mean(ibex_ret))
w_opt <- matrix(NA, nrow = T, ncol = 2)
for (t in 1:T) {
  H_inv <- solve(Ht[t,,])
  temp <- H_inv %*% mu
  w_opt[t,] <- temp / sum(temp)
}

# --- Combine Outputs in a Data Frame ---
results <- data.frame(
  sigma2_hsi = sigma2_hsi,
  sigma2_ibex = sigma2_ibex,
  cov_hsi_ibex = cov_hsi_ibex,
  port_var = port_var,
  VaR_1pct = VaR_1pct,
  weight_hsi = w_opt[,1],
  weight_ibex = w_opt[,2]
)

# --- Plotting ---
par(mfrow = c(3,1))
plot(results$sigma2_hsi, type='l', main="Conditional Variance: HSI", ylab="Variance", xlab="Time")
plot(results$sigma2_ibex, type='l', main="Conditional Variance: IBEX", ylab="Variance", xlab="Time")
plot(results$cov_hsi_ibex, type='l', main="Conditional Covariance: HSI & IBEX", ylab="Covariance", xlab="Time")

par(mfrow = c(2,1))
plot(results$port_var, type='l', main="Portfolio Conditional Variance", ylab="Variance", xlab="Time")
plot(results$VaR_1pct, type='l', main="Portfolio 1% Value-at-Risk (VaR)", ylab="VaR", xlab="Time")

matplot(results[, c("weight_hsi", "weight_ibex")], type='l', col=c("blue","red"),
        main="Optimal Portfolio Weights (Sharpe Ratio)", ylab="Weight", xlab="Time")
legend("topright", legend=c("HSI","IBEX"), col=c("blue","red"), lty=1)

# --- Print Weights for T+1 ---
cat("\nOptimal weights at T+1:\n")
print(round(w_opt[T,], 4))