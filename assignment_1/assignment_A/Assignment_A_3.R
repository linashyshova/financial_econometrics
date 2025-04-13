##  MAXIMUM LIKELIHOOD ESTIMATE PARAMETERS OF GARCH(1,2)
#
#  Description: 
#  This code snippet shows how to optimize the log likelihood
#  and estimate the parameters of a GARCH(1,2) model by maximum likelihood.
#  Note tha the log-likelihood is contained in the R file "llik_fun_GARCH_1_2.R"

########################
############ 00. Clean workspace
########################

rm(list=ls())    

########################
############ 0. Load packages and/or functions
########################

# The average log-likelihood function for the GARCH(1,2) is llik_fun_GARCH_1_2(). 
# This function is contained in the R file "llik_fun_GARCH_1_2.R". 

source("assignment_1/assignment_A/llik_fun_GARCH_1_2.R") 

########################
############ 1. Load data
########################

df <- read.table("assignment_1/data/HSI.txt")
names(df) <- "price"  
df$log_return <- c(0, 100 * diff(log(df$price)))
df$sq_log_return <- df$log_return^2
x <- df$sq_log_return
########################
############ 2. Set initialization for numerical optimization
########################

## 2a. choose the initail parameter values for the Newton Raphson optimization

a1 <- 0.2 # initial value for alpha1
a2 <- 0.2 # initial value for alpha2
b <- 0.5  # initial value for beta
omega <- var(na.omit(x))*(1-a1-a2-b) # initial value for omega
## 2b. Transform intitial values using the inverse of the link funtions

# In the likelihood function llik_fun_GARCH_1_2() the parameter inputs are 
# transformed using link functions to ensure the parameter restrictions.
# This is done to avoid numerical problems in the optimization.
# The restrictions are the following (see llik_fun_GARCH_1_2): 
# - exp() to ensure omega>0
# - logistic()=exp()/(exp()) for 0<alpha1<1
# - logistic()=exp()/(exp()) for 0<alpha2<1
# - logistic()=exp()/(exp()) for 0<beta<1
# Therefore, we transform the parameter values with the inverse of these functions
# to give the desired initialization:

par_ini <- c(log(omega),log(a1/(1-a1)),log(a2/(1-a2)),log(b/(1-b)))
par_ini

########################
############ 3. Optimize Log Likelihood function
########################

#Note that the numerical optimizer optim() minimizes a function. We want
#to maximize the log-likelihood. Therefore we need to provide a negative
#log-likelihood. Note also that for numerical reasons we optimize the
#average log-likelihood instead of the log-likelihood.

# Optim input:
# (1) Initial parameter (par): par_ini
# (2) Negative average log-likelihood function (fn): - llik_fun_GARCH()
# (3) Numerical algorithm used (method): SANN

# fmincon output:
# (1) estimates: $par
# (2) log likelihood function value at theta_hat: $value
# (3) exit flag indicating convergence: $convergence

est <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH_1_2(par,x), method = "SANN", control = list(maxit = 1000000))
print(est$value)
########################
############ 4. Obtain parameter estimates
########################

# est$par gives the values that maximize the likelihood
# However, these are not the paramter estimates since we have used the link functions
# to guarantee omega>0, 0<alpha<1 and 0<beta<1.
# Therefore, we apply these link functions to obtain the parameter estimates

omega_hat <- exp(est$par[1])
alpha1_hat <- exp(est$par[2])/(1+exp(est$par[2]))
alpha2_hat <- exp(est$par[3])/(1+exp(est$par[3]))
beta_hat <- exp(est$par[4])/(1+exp(est$par[4]))

theta_hat <- c(omega_hat,alpha1_hat,alpha2_hat,beta_hat)
cat("The parameter estimates are:")
round(theta_hat,6)
# display the log-likelihood. 
# Note that we multiply the average log-likelihood by the sample size and
# -1 to get the log-likelihood.

# display the exit flag to see if convergence of the algorithm has been attained

cat("Exit flag:")
est$convergence # zero indicates succesfull optimization

########################
############ 5. Calculate and plot filtered conditional variance
########################

n <- length(x)
sigma2 <- rep(0,n)
sigma2[1:2] <- var(x)

for(t in 3:n){
  sigma2[t] = omega_hat + alpha1_hat*x[t-1]^2 + alpha2_hat*x[t-2]^2 + beta_hat*sigma2[t-1]
}

df$filtered_var <- sigma2

start_date <- as.Date("1997-01-01")
df$date <- seq(from = start_date, by = "week", length.out = nrow(df))

# Plot filtered conditional variance
ggplot(df, aes(x = date, y = filtered_var)) +
  geom_line(color = "blue") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Year", y = "Filtered conditional variance", title = "HSI Filtered conditional variance") +
  theme(axis.text.x = element_text(angle = 90))

