# 0. Set seed

set.seed(123)

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(tidyverse)))

# 2. Data simulation

n <- 50 # Number of observations

p <- 10 # Number of parameters

s <- p - 4 # Number of zero elements

x <- matrix(data = NA, nrow = n, ncol = p) # Matrix of explanatory variables

x[,1] <- rep(1, n) # Intercept column
x[,2] <- sample(x = 0:1, size = n, replace = TRUE, prob = c(0.5, 0.5)) # Dummy variable
x[,3] <- rpois(n = n, lambda = 1) # Poisson-distributed variable

for (j in 4:p) {
  x[,j] <- rnorm(n = n, mean = seq(from = -3, to = 4, by = 1), sd = 2) # Normal-distributed variables
}

beta_true <- c(0.5, 2, -0.5, 1, rep(0, s)) # Mean parameter

sigma2_true <- 1 # Variance parameter

y <- rnorm(n = n, mean = x%*%beta_true, sd = sqrt(sigma2_true)) # Response variable

# 3. Simulation framework to evaluate Bayesian regression models

e <- c(rep(1602, 18), rep(102, 18), rep(3, 18)) # Shape parameter of inverse-gamma distribution

f <- c(rep(1601, 9), rep(16010, 9), rep(101, 9), rep(1010, 9), rep(2, 9), rep(20, 9)) # Scale parameter of inverse-gamma distribution

g <- c(rep(1600, 3), rep(100, 3), rep(1, 3)) # Shape parameter of gamma distribution

h <- c(1600, 160, 16, 100, 10, 1, 1, 0.1, 0.01) # Rate parameter of gamma distribution

hyperparameter <- cbind(e, f, rep(g, 6), rep(h, 6))

# Before executing the following code, sample_beta_lasso, sample_tau, sample_sigma2_lasso, 
# and sample_lambda_lasso functions must be executed

# 4. Gibbs sampling algorithm

Gibbs_lasso <- function(y, x, hyperparameter, n, p, n_skip, n_sams, n_burn, verbose = TRUE)
{
  X <- t(x)%*%x # Compute X^{\top} X
  X_y <- t(x)%*%y # Compute X^{\top} y
  
  # Number of iterations of the Gibbs sampling algorithm
  B <- n_burn + n_sams*n_skip
  ncat <- floor(0.01*B)
  
  # Objects where the samples of beta, and sigma2 will be stored
  BETA <- array(data = NA, dim = c(nrow(hyperparameter), n_sams, p))
  SIGMA <- array(data = NA, dim = c(nrow(hyperparameter), n_sams, 1))
  
  for (j in 1:nrow(hyperparameter)) {
    # Initialize beta, tau, sigma2, and lambda values
    beta <- rep(0, p)
    tau <- rep(1, p)
    sigma2 <- 1
    lambda <- 1
    
    # Gibbs sampling algorithm
    for (i in 1:B) {
      beta <- sample_beta_lasso(sigma2, tau, X, X_y, p) # Update beta
      tau <- sample_tau(beta, lambda, p) # Update tau2
      sigma2 <- sample_sigma2_lasso(beta, y, x, e = hyperparameter[j,1], f = hyperparameter[j,2], n) # Update sigma2
      lambda <- sample_lambda_lasso(tau, g = hyperparameter[j,3], h = hyperparameter[j,4], p) # Update lambda
      
      # Save effective samples
      if (i > n_burn && (i - n_burn) %% n_skip == 0) {
        t <- (i - n_burn) / n_skip 
        BETA[j, t, ] <- beta
        SIGMA[j, t,] <- sigma2
      }
      # Algorithm progress
      if (verbose && i %% ncat == 0)
        cat(sprintf("%.1f%% completado\n", 100*i/B))
    }
  }
  return(list(BETA = BETA, SIGMA = SIGMA))
}

# 4.1 Gibbs sampling algorithm implementation

M3 <- Gibbs_lasso(y, x, hyperparameter, n, p, 
                  n_skip = 10, # Accounting for Markov chain autocorrelation will require systematic sampling
                  n_sams = 10000, # Set the number of effective samples
                  n_burn = 1000) # Set the number of burn-in samples

# 5. Bayesian inference for beta and sigma2

simulation_lasso <- function(model, hyperparameter, p)
{
  # Object where the inference for beta, and sigma2 will be stored
  table <- data_frame("Escenario de simulación" = numeric(), 
                      "Coeficiente" = numeric(), 
                      "Media posterior" = numeric(),
                      "Límite inferior IC" = numeric(), 
                      "Límite superior IC" = numeric())
  
  table_2 <- table
  
  for (i in 1:nrow(hyperparameter)) {
    # Inference for the i-th simulation scenery
    BETA_HAT <- apply(model$BETA[i,,], MARGIN = 2, FUN = mean) # Posterior mean for beta
    CI_BETA <- apply(model$BETA[i,,], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)) # 95% credible interval
    
    SIGMA2_HAT <- mean(model$SIGMA[i,,]) # Posterior mean for sigma2
    CI_SIGMA <- quantile(x = model$SIGMA[i,,], probs = c(0.025, 0.975)) # 95% credible interval)
    
    for (j in 1:p) {
      table[j,] <- cbind(i, j, BETA_HAT[j], CI_BETA[1, j], CI_BETA[2, j])
    }
    
    table[p + 1, ] <- cbind(i, p + 1, SIGMA2_HAT, CI_SIGMA[1], CI_SIGMA[2])
    
    table_2 <- rbind(table_2, table)
  }
  return(table = table_2)
}

simulation_result <- simulation_lasso(M3, hyperparameter, p)
