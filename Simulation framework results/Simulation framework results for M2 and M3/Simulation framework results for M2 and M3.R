# 0. Set seed

set.seed(123)

# 1. Load dataset and necessary libraries

load("~/Trabajo de grado/Simulation framework/Dataset - Simulation framework.RData")

suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(GIGrvg)))

# 2. Auxiliary functions

# Compute Deviance Information Criterion

compute_DIC <- function(model, y, x, beta, sigma2){
  LL_HAT <- sum(dnorm(x = y, mean = x%*%beta, sd = sqrt(sigma2), log = TRUE))
  LL_B <- model$LL
  
  pDIC <- 2*(LL_HAT - mean(LL_B))
  DIC <- -2*LL_HAT + 2*pDIC
  
  return(DIC)
}

# Compute Watanabe-Akaike Information Criterion

compute_WAIC <- function(model, y, x, n){
  LPPD <- 0
  pWAIC <- 0
  
  for (i in 1:n) {
    # LPPD
    TMP <- dnorm(x = y[i], mean = x[i,]%*%t(model$BETA), sd = sqrt(model$SIGMA))
    LPPD <- LPPD + log(mean(TMP))
    # pWAIC
    TMP_2 <- dnorm(x = y[i], mean = x[i,]%*%t(model$BETA), sd = sqrt(model$SIGMA), log = TRUE)
    pWAIC <- pWAIC + 2*(log(mean(TMP)) - mean(TMP_2))
  }
  WAIC <- -2*LPPD + 2*pWAIC
  
  return(WAIC)
}

# Inference for Bayesian parametric Ridge and LASSO 

hat <- function(model, scenery, y, x, n, p, i){
  # Posterior mean for beta and sigma2
  beta <- colMeans(model$BETA)
  sigma2 <- mean(model$SIGMA)
  
  # 95% credible interval
  beta_ic <- apply(model$BETA, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
  sigma2_ic <- quantile(model$SIGMA, probs = c(0.025, 0.975))
  
  # Object where the Bayesian inference for beta and sigma2 will be stored
  inference <- data.frame(`Escenario` = i, `Coeficiente` = 1:(p + 1),
                          `Media posterior` = c(beta, sigma2),
                          `IC LI` = c(beta_ic[1,], sigma2_ic[1]),
                          `IC LS` = c(beta_ic[2,], sigma2_ic[2]))
  return(list(inference = inference, beta = beta, sigma2 = sigma2))
}

# Cross validation

cross_validation <- function(y, x, n, p, 
                             a, b, c, d,
                             e, f, g, h){
  # Create the train and test dataset
  index <- sample(1:n, size = 0.7*n)
  
  y_train <- y[index]; x_train <- x[index,] # Train dataset
  y_test <- y[-index]; x_test <- x[-index,] # Test dataset
  
  # Objects where the mean absolute error, and the mean squared prediction error will be stored
  mape <- NULL
  mspe <- NULL
  
  # Fit Bayesian Ridge and LASSO regression models
  M2 <- Gibbs_ridge(y_train, x_train, a, b, c, d, n = length(y_train), p,
                    n_burn = 1000,
                    n_sams = 10000,
                    n_skip = 10)
    
  M3 <- Gibbs_lasso(y_train, x_train, e, f, g, h, n = length(y_train), p,
                    n_burn = 1000,
                    n_sams = 10000,
                    n_skip = 10)
    
  # Posterior mean for beta
  beta_ridge <- colMeans(M2$BETA)
  beta_lasso <- colMeans(M3$BETA)
    
  # Linear predictor
  y_hat_ridge <- x_test%*%beta_ridge
  y_hat_lasso <- x_test%*%beta_lasso
    
  # Compute mean absolute prediction error
  mape_ridge <- mean(abs(y_test - y_hat_ridge))
  mape_lasso <- mean(abs(y_test - y_hat_lasso))
    
  # Compute mean squared prediction error
  mspe_ridge <- mean((y_test - y_hat_ridge)^2)
  mspe_lasso <- mean((y_test - y_hat_lasso)^2)
    
  mape <- rbind(mape, c(mape_ridge, mape_lasso))
  mspe <- rbind(mspe, c(mspe_ridge, mspe_lasso))
  return(list(mape = mape, mspe = mspe))
}

# 3. Hyperparameter elicitation

a <- 3; e <- 3 # Shape parameter of inverse-gamma distribution

b <- 2; f <- 2 # Scale parameter of inverse-gamma distribution

c <- 1; g <- 1 # Shape parameter of gamma distribution

d <- 1; h <- 1 # Rate parameter of gamma distribution

# 4. Results

simulation_results <- function(data,
                               a, b, c, d,
                               e, f, g, h){
  # Number of sceneries
  scenery <- length(data$Y)
  
  # Objects where the evaluation metrics will be stored
  RIDGE <- NULL
  LASSO <- NULL
  INF_CRI <- NULL
  CRO_VAL <- NULL
  
    for (i in 1:scenery) {
    y <- data$Y[[i]] # Set the response variable
    x <- data$X[[i]] # Set the matrix containing the explanatory variables
    n <- length(y) # Sample size
    
    # The parametric model assumes k = 1
    beta_true <- data$BETA[[i]][1,]
    p <- length(beta_true) # Number of explanatory variables
    
    # Before executing the following code, Bayesian Ridge.r, and Bayesian LASSO
    # must be executed, since it displays the Gibbs sampling algorithm
    
    M2 <- Gibbs_ridge(y, x, a, b, c, d, n, p,
                      n_burn = 1000,
                      n_sams = 10000,
                      n_skip = 10)
    
    M3 <- Gibbs_lasso(y, x, e, f, g, h, n, p,
                      n_burn = 1000,
                      n_sams = 10000,
                      n_skip = 10)
    
    # Inference for Bayesian parametric Ridge and LASSO
    ridge <- hat(M2, Scenery, y, x, n, p, i)
    lasso <- hat(M3, Scenery, y, x, n, p, i)
    
    RIDGE <- rbind(RIDGE, round(ridge$inference, 2))
    LASSO <- rbind(LASSO, round(lasso$inference, 2))
    
    # Compute Information Criterion
    ridge_dic <- compute_DIC(M2, y, x, ridge$beta, ridge$sigma2)
    ridge_waic <- compute_WAIC(M2, y, x, n)
    
    lasso_dic <- compute_DIC(M3, y, x, lasso$beta, lasso$sigma2)
    lasso_waic <- compute_WAIC(M3, y, x, n)
    
    inf_cri <- data.frame(`Escenario` = i,
                          `DIC Ridge` = ridge_dic, `DIC LASSO` = lasso_dic, 
                          `WAIC Ridge` = ridge_waic, `WAIC LASSO` = lasso_waic)
    
    INF_CRI <- rbind(INF_CRI, round(inf_cri, 2))
    
    # Cross validation
    cro_val <- cross_validation(y, x, n, p, a, b, c, d, e, f, g, h)
    
    CRO_VAL <- rbind(CRO_VAL, c(i, round(cro_val$mape, 2), round(cro_val$mspe, 2)))
    }
  return(list(RIDGE = RIDGE, LASSO = LASSO, INF_CRI = INF_CRI, CRO_VAL = CRO_VAL))
}

results <- simulation_results(data, a, b, c, d, e, f, g, h)
