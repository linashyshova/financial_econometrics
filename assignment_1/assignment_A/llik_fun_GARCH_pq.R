# This file contains the log-likelihood function for a GARCH(p,q) model.
# The function takes the data (x), the parameters (par), and the order of the GARCH model (p, q) as inputs
# and returns the average log-likelihood value.

llik_fun_GARCH_pq <- function(x, par, p, q){
    n <- length(x)
    # Extract parameters
    omega <- exp(par[1]) # Ensure omega > 0
    
    alpha <- numeric(q)
    for (i in (1:q)){
        alpha[i] <- exp(par[i + 1]) / (1 + exp(par[i + 1])) # Ensure alpha_i > 0
    }
    beta <- numeric(p)
    for (i in (1:p)){
        beta[i] <- exp(par[i + q + 1]) / (1 + exp(par[i + q + 1])) # Ensure beta_i > 0
    }

    # Check if the sum of alpha and beta is less than 1
    if (sum(alpha) + sum(beta) >= 1) {
        return(NA)
    }
    
    # Filter volatility
    sig2 <- rep(0,n)
    sig2[1:max(p, q)] <- var(x) # Initialize with sample variance
    for (t in (max(p, q) + 1):n) {
        sig2[t] <- omega + sum(alpha * x[(t - 1):(t - q)]^2) + sum(beta * sig2[(t - 1):(t - p)])
    }

    # Calculate log-likelihood
    l <- -(1/2)*log(2*pi) - (1/2)*log(sig2) - (1/2)*x^2/sig2
    llik <- mean(l)
    
    return (llik)
}