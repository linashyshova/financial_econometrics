# Clear workspace
rm(list=ls())

# Load necessary packages and files
library(ggplot2)
library(tseries)
library(scales)
library(xtable)


### ------Calculate log-likelihood function------
llik_fun_GARCH_pq <- function(data, par, p, q) {
  n <- length(data)
  omega <- exp(par[1])  # ensure omega > 0

  alpha <- exp(par[2:(1+q)]) / (1 + exp(par[2:(1+q)]))
  beta <- exp(par[(2+q):(1+q+p)]) / (1 + exp(par[(2+q):(1+q+p)]))

  if (omega <= 0 || any(alpha < 0) || any(beta < 0)) return(1e10)
  if (sum(alpha) + sum(beta) >= 1) return(1e10)  # stationarity

  sig2 <- rep(0, n)  # initialize with sample variance
  sig2[1:max(p,q)] <- var(data)

  for (t in (max(p, q) + 1):n) {
    arch_terms <- sum(alpha * data[(t - 1):(t - q)]^2)
    garch_terms <- sum(beta * sig2[(t - 1):(t - p)])
    sig2[t] <- omega + arch_terms + garch_terms
  }

  llik <- -0.5 * sum(log(sig2[(max(p,q)+1):n]) + 
                     data[(max(p,q)+1):n]^2 / sig2[(max(p,q)+1):n])

  return(-llik)
}


############################## Question 1 ###############################
# Read data
hsi_data <- read.table("HSI.txt")

# Add column name
names(hsi_data) <- "price"
# Add another column to indicate the date. Let's assume we started collecting data from 01.01.1997
start_date <- as.Date("1997-01-01")
hsi_data$date <- seq(from = start_date, by = "week", length.out = nrow(hsi_data))

# Plot weekly prices
ggplot(hsi_data, aes(x = date, y = price)) +
  geom_line(color = "steelblue", linewidth = 0.7) +
  labs(x = "Year", y = "Price USD", title = "HSI Index Price - Weekly") +
  theme_minimal()
  
# Calculate the log-returns
hsi_data$log_return <- c(NA, 100 * diff(log(hsi_data$price)))

# Plot weekly log-returns
ggplot(hsi_data, aes(x = date, y = log_return)) +
  geom_line(na.rm = TRUE, color = "steelblue", linewidth = 0.7) +
  labs(x = "Year", y = "Log Return (%)", title = "Weekly Log Returns of Stock Price") +
  theme_minimal()

par(mfrow = c(1,2))
# ACF for weekly log returns
acf(hsi_data$log_return, na.action = na.pass, main = "")

# ACF for weekly squared log returns
acf(hsi_data$log_return^2, na.action = na.pass, main = "")




############################## Question 3 ###############################
estimate_garch_mle <- function(data, p, q) {
  n <- length(data)
  num_params <- 1 + p + q

  init_alpha <- rep(0.1 / q, q)
  init_beta <- rep(0.8 / p, p)

  # Scale if needed
  if (sum(init_alpha) + sum(init_beta) >= 1) {
    scale <- 0.95 / (sum(init_alpha) + sum(init_beta))
    init_alpha <- init_alpha * scale
    init_beta <- init_beta * scale
  }

  init_omega <- max(1e-6, var(data) * (1 - sum(init_alpha) - sum(init_beta)))
  par_init <- c(log(init_omega),
                log(init_alpha / (1 - init_alpha)),
                log(init_beta / (1 - init_beta)))

  opt <- optim(par_init, function(par) llik_fun_GARCH_pq(data, par, p, q),
               method = "BFGS", control = list(maxit = 10000))

  omega_hat <- exp(opt$par[1])
  alpha_hat <- exp(opt$par[2:(1+q)]) / (1 + exp(opt$par[2:(1+q)]))
  beta_hat <- exp(opt$par[(2+q):(1+q+p)]) / (1 + exp(opt$par[(2+q):(1+q+p)]))

  sig2_fitted <- rep(0, n)
  sig2_fitted[1:max(p,q)] <- var(data)
  for (t in (max(p, q) + 1):n) {
    sig2_fitted[t] <- omega_hat +
      sum(alpha_hat * data[(t - 1):(t - q)]^2) +
      sum(beta_hat * sig2_fitted[(t - 1):(t - p)])
  }

  loglik <- -opt$value
  aic <- 2 * num_params - 2 * loglik
  bic <- log(n) * num_params - 2 * loglik

  return(list(
    p = p, q = q,
    omega = omega_hat,
    alpha = alpha_hat,
    beta = beta_hat,
    filtered_var_est = sig2_fitted,
    loglik = loglik,
    aic = aic,
    bic = bic,
    convergence = opt$convergence,
    stationary = sum(alpha_hat) + sum(beta_hat) < 1
  ))
}

hsi_data <- na.omit(hsi_data)
y <- hsi_data$log_return
garch_12 <- estimate_garch_mle(y, p = 1, q = 2)

cat("Omega:")
garch_12$omega
cat("Alpha:")
garch_12$alpha
cat("Beta:")
garch_12$beta
cat("Log-likelihood:")
garch_12$loglik


hsi_data$filter_var <- garch_12$filtered_var_est
ggplot(hsi_data, aes(x = date, y = filter_var)) +
  geom_line(color = "steelblue") +
  labs(x = "Year", y = "Filtered conditional variance", title = "HSI Filtered conditional variance") +
  theme_minimal()


############################## Question 4 ###############################
# Homoskedasticity test by plotting the ACF of squared residuals
u <- y/sqrt(hsi_data$filter_var)
acf(u^2, main="")
title("ACF squared residuals", line = 0.3)

# Jarque-Bera test on normality assumption of residuals
jarque.bera.test(u)



############################## Question 5 ###############################
compare_garch_models <- function(data, max_p, max_q) {
  models <- list()
  results <- data.frame()
  
  for (p in 1:max_p) {
    for (q in 1:max_q) {
      model_name <- paste("GARCH(", p, ",", q, ")", sep = "")
      
      model <- estimate_garch_mle(data, p, q)
      models[[model_name]] <- model
      
      results <- rbind(results, data.frame(
        Model = model_name,
        LogLikelihood = model$loglik,
        AIC = model$aic,
        BIC = model$bic,
        Stationary = model$stationary,
        Convergence = model$convergence
      ))
    }
  }
  
  return(list(models = models, results = results))
}

comparison <- compare_garch_models(hsi_data$log_return, max_p = 5, max_q = 5)
print(comparison$results)

# Convert to LaTeX table
latex_table <- xtable(comparison$results, caption = "GARCH(p,q) Model Comparison", label = "tab:garch_models")
#print(latex_table, include.rownames = FALSE, file = "garch_model_results.tex")

# Select the best model based on AIC value
best_model_name <- comparison$results[which.min(comparison$results$AIC), "Model"]
best_model <- comparison$models[[best_model_name]]
cat("Best GARCH model based on AIC:", best_model_name, "\n")

# Plot filtered variances of the best model
hsi_data$filter_var_best <- best_model$filtered_var_est
ggplot(hsi_data, aes(x = date, y = filter_var_best)) +
  geom_line(color = "steelblue") +
  labs(x = "Year", y = "Filtered conditional variance", title = paste("Filtered conditional variance -", best_model_name)) +
  theme_minimal()
