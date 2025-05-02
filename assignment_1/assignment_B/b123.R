library(ggplot2)
library(scales)

# --- Load Data ---
hsi <- read.table("HSI.txt")
ibex <- read.table("IBEX.txt")

hsi_ret <- diff(log(hsi$V1))
ibex_ret <- diff(log(ibex$V1))
T <- length(hsi_ret)

# --- Create Time Index ---
start_date <- as.Date("1997-01-10")  # adjust as needed
date_returns <- seq(from = start_date, by = "week", length.out = T)

# --- Plot Helper Function ---
plot_series <- function(dates, values, ylab_title) {
  df <- data.frame(date = dates, value = values)
  ggplot(df, aes(x = date, y = value)) +
    geom_line(color = "darkblue", linewidth = 0.6) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = "Year", y = ylab_title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

# =====================================================================
#                                1b.1
# =====================================================================

# GARCH(1,1) log-likelihood
garch_loglik <- function(params, data) {
  omega <- abs(params[1])
  alpha <- abs(params[2])
  beta <- abs(params[3])
  sigma2 <- rep(var(data), length(data))
  ll <- 0
  for (t in 2:length(data)) {
    sigma2[t] <- omega + alpha * data[t-1]^2 + beta * sigma2[t-1]
    ll <- ll - 0.5 * (log(2 * pi) + log(sigma2[t]) + data[t]^2 / sigma2[t])
  }
  return(-ll)
}

# Estimate parameters
opt_hsi <- optim(c(0.01, 0.05, 0.9), garch_loglik, data = hsi_ret)
opt_ibex <- optim(c(0.01, 0.05, 0.9), garch_loglik, data = ibex_ret)
params_hsi <- abs(opt_hsi$par)
params_ibex <- abs(opt_ibex$par)

# Filter conditional variances
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

# Estimate CCC correlation
z_hsi <- hsi_ret / sqrt(sigma2_hsi)
z_ibex <- ibex_ret / sqrt(sigma2_ibex)
rho_hat <- cor(z_hsi, z_ibex)

# Construct Ht and covariances
Ht <- array(NA, c(T, 2, 2))
cov_hsi_ibex <- numeric(T)
for (t in 1:T) {
  D_t <- diag(c(sqrt(sigma2_hsi[t]), sqrt(sigma2_ibex[t])))
  R <- matrix(c(1, rho_hat, rho_hat, 1), 2)
  Ht[t,,] <- D_t %*% R %*% D_t
  cov_hsi_ibex[t] <- Ht[t,1,2]
}

# --- Plots 1b.1 ---
plot_series(date_returns, sigma2_hsi, "HSI Conditional Variance")
plot_series(date_returns, sigma2_ibex, "IBEX Conditional Variance")
plot_series(date_returns, cov_hsi_ibex, "Conditional Covariance")

# =====================================================================
#                                1b.2
# =====================================================================

w <- c(0.3, 0.7)
port_var <- sapply(1:T, function(t) t(w) %*% Ht[t,,] %*% w)
VaR_1pct <- -qnorm(0.01) * sqrt(port_var)

# --- Plots 1b.2 ---
plot_series(date_returns, port_var, "Portfolio Conditional Variance")
plot_series(date_returns, VaR_1pct, "1% Value-at-Risk")

# =====================================================================
#                                1b.3
# =====================================================================

mu <- c(mean(hsi_ret), mean(ibex_ret))
w_opt <- matrix(NA, nrow = T, ncol = 2)
for (t in 1:T) {
  H_inv <- solve(Ht[t,,])
  temp <- H_inv %*% mu
  w_opt[t,] <- temp / sum(temp)
}

# --- Plot Optimal Weights ---
df_w <- data.frame(
  date = date_returns,
  HSI = w_opt[,1],
  IBEX = w_opt[,2]
)

library(tidyr)
df_long <- pivot_longer(df_w, cols = c("HSI", "IBEX"), names_to = "Asset", values_to = "Weight")

ggplot(df_long, aes(x = date, y = Weight, color = Asset)) +
  geom_line(linewidth = 0.7) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Year", y = "Optimal Weight") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# --- Print T+1 weights ---
cat("Optimal weights at T+1:")
cat("HSI:", round(w_opt[T, 1], 4))
cat("IBEX:", round(w_opt[T, 2], 4))
