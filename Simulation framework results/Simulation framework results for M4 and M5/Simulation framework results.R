# 0. Set seed

set.seed(123)

suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(GIGrvg)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(mclust)))

# 1. Load dataset

load("~/Trabajo de grado/Simulation framework/Dataset - Simulation framework.RData")

# 2. Auxiliary functions

# Compute Watanabe-Akaike Information Criterion

compute_WAIC_np <- function(model, y, x, n){
  ite <- nrow(model$XI)
  TMP <- matrix(data = NA, nrow = ite, ncol = n)
  TMP_2 <- matrix(data = NA, nrow = ite, ncol = n)

  for (i in 1:n) {
    for (b in 1:ite) {
      xi <- model$XI
      beta <- model$BETA[[b]][xi[b,i], , drop = FALSE]
      sigma2 <- model$SIGMA[[b]][xi[b,i]]
    
      # LPPD
      TMP[b,i] <- dnorm(x = y[i], mean = x[i,]%*%t(beta), sd = sqrt(sigma2))
      
      #pWAIC
      TMP_2[b,i] <- dnorm(x = y[i], mean = x[i,]%*%t(beta), sd = sqrt(sigma2), log = TRUE) 
    }
  }

  LPPD <- sum(log(colMeans(TMP)))
  pWAIC <- 2*sum((log(colMeans(TMP)) - colMeans(TMP_2)))

  WAIC <- -2*LPPD + 2*pWAIC
  
  return(WAIC)
}

# Compute Incidence Matrix

matrix_A <- function(model, k_pos, n){
  ite <- nrow(model$XI)
  A <- matrix(data = 0, nrow = n, ncol = n)
  
  for (b in 1:ite) {
    xi_b <- model$XI[b,]
    for (k in unique(xi_b)) {
      cluster_i <- which(xi_b == k)
      A[cluster_i, cluster_i] <- A[cluster_i, cluster_i] + 1
    }
  }
  A <- A / ite
  diag(A) <- 1
  
  # Convert A to the incidence matrix
  A <- 1 - A
  
  # Apply multidimensional scaling
  coord <- cmdscale(as.dist(A), k = 4)
  
  # Clustering based on co-clustering probabilities
  xi_hat <- Mclust(data = coord, G = 1:k_pos)$classification
  
  return(xi_hat = xi_hat)
}

# Inference on the number of clusters

hat_k <- function(model){
  K <- apply(model$XI, 1, function(x) length(unique(x))) # Compute the number of clusters at each iteration
  K_table <- as.data.frame(table(K)/length(K))
  
  # Posterior number of clusters
  k_pos <- as.numeric(as.character(K_table$K[which.max(K_table$Freq)]))
  
  return(list(K = K_table, k_pos = k_pos))
}

# Inference for Bayesian nonparametric Ridge and LASSO

hat_np <- function(model, K, p){
  ite <- nrow(model$XI)
  permu <- gtools::permutations(n = K, r = K) # Permutations
  
  # Objects where the samples of beta, and sigma2 will be stored
  BETA_corrected <- lapply(1:p, function(x) matrix(nrow = 0, ncol = K))
  SIGMA2_corrected <- matrix(nrow = 0, ncol = K)
  
  # Objects where the posterior mean of beta, and sigma2 will be stored
  beta_pos <- matrix(0, nrow = K, ncol = p)
  sigma2_pos <- numeric(K)
  
  # Posterior distribution of K
  k <- apply(X = model$XI, MARGIN = 1, function(x) length(unique(x)))
  k_tab <- table(factor(x = k, levels = k, labels = k))
  
  for (b in 1:ite) {
    if (length(unique(model$XI[b,])) == K) {
      beta_pos <- beta_pos +
        (model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE] / max(k_tab))
      
      sigma2_pos <- sigma2_pos + (model$SIGMA[[b]] / max(k_tab))
    }
  }
  
  for(b in 1:ite){
    if (length(unique(model$XI[b,])) == K) {
      # Average over the permuted spaces
      beta_current <- model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE]
      sigma2_current <- model$SIGMA[[b]]
      
      # Reorder according to the permutations, and compute the distance of each sample to its posterior mean
      dist <- apply(X = permu, MARGIN = 1, 
                    FUN = function(perm) {
                      permuted_beta <- sum((beta_current[perm,, drop = FALSE] - beta_pos)^2)
                    }
      )
      # Select the optimum permutation
      best_permu <- permu[which.min(dist),]
      for (j in 1:p) {
        BETA_corrected[[j]] <- rbind(BETA_corrected[[j]], beta_current[best_permu, j]) 
      }
      
      # Reorder according to the permutations, and compute the distance of each sample to its posterior mean
      dist2 <- apply(X = permu, MARGIN = 1, 
                     FUN = function(perm) {
                       permuted_sigma2 <- sigma2_current[perm]
                       sum((permuted_sigma2 - sigma2_pos)^2)
                     }
      )
      # Select the optimum permutation
      best_permu2 <- permu[which.min(dist2),]
      SIGMA2_corrected <- rbind(SIGMA2_corrected, sigma2_current[best_permu2])
    }
  }
  
  # Posterior mean, and 95% credible interval
  inf <- matrix(data = NA, nrow = p + 1, ncol = 4*K)
  id <- function(k)(4*(k - 1) + 1):(4*k)
  
  for (j in 1:p) {
    beta_mean <- colMeans(BETA_corrected[[j]])
    beta_median <- apply(X = BETA_corrected[[j]], MARGIN = 2, FUN = median)
    beta_ic <- apply(X = BETA_corrected[[j]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
  
      for (k in 1:K) {
        inf[j, id(k)] <- c(beta_mean[k], beta_median[k], beta_ic[1, k], beta_ic[2, k])
      }
    }
  sigma2_mean <- colMeans(SIGMA2_corrected)
  sigma2_median <- apply(SIGMA2_corrected, MARGIN = 2, FUN = median)
  sigma2_ic <- apply(SIGMA2_corrected, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
  
  for (k in 1:K) {
    inf[p + 1, id(k)] <- c(sigma2_mean[k], sigma2_median[k], sigma2_ic[1, k], sigma2_ic[2, k])
  }
  
  col <- seq(from = 1, to = ncol(inf), by = 4) # Position of the posterior mean
  beta <- inf[c(1:p), col, drop = FALSE]
  sigma2 <- inf[(p + 1), col]
  
  return(list(inference = inf, beta = beta, sigma2 = sigma2, beta_pos = beta_pos, sigma2_pos = sigma2_pos))
}

# Out of sample prediction

out_sample_ridge <- function(y_test, x_test, K, xi_hat, beta, sigma2, alpha,
                             a, b, c, d, p){
  # Object where the cluster assignment will be stored
  xi_test <- numeric(length(y_test))
  
  for (i in 1:length(y_test)) {
    # Object where the probability for the i-th observation of the test data set
    # to be assigned to a cluster will be stored
    log_probs <- numeric(K + 1)
    
    # Vector of explanatory variables, and response variable of the i-th observation
    # of the test data
    x_i <- matrix(data = x_test[i,], nrow = 1, ncol = p, byrow = FALSE)
    y_i <- y_test[i]
    
    for (k in 1:K) {
      n_k <- sum(xi_hat == k) # Number of observations in the k-th cluster
      log_probs[k] <- log(n_k) + dnorm(y_i, mean = x_i%*%beta[,k], sd = sqrt(sigma2[k]), log = TRUE)
    }
    
    # Sample beta from its prior distribution
    lambda_prior <- rgamma(n = 1, shape = c, rate = d)
    beta_prior <- c(mvtnorm::rmvnorm(n = 1, mean = rep(0, p), sigma = (1/lambda_prior)*(diag(x = 1, p))))
    
    # Sample sigma2 from its prior distribution
    sigma2_prior <- 1/rgamma(n = 1, shape = a, rate = b)
    
    # Probability for the i-th observation to be assigned to a new cluster K + 1
    log_probs[K + 1] <- log(alpha) + dnorm(y_i, mean = x_i%*%beta_prior, sd = sqrt(sigma2_prior), log = TRUE)
    
    # Sample the cluster assignment for the i-th observation of the test set
    cluster <- sample(1:(K + 1), size = 1, prob = exp(log_probs - max(log_probs)))
    xi_test[i] <- cluster
  }
  return(xi_test = xi_test)
}

out_sample_lasso <- function(y_test, x_test, K, xi_hat, beta, sigma2, alpha,
                             e, f, g, h, p){
  # Object where the cluster assignment will be stored
  xi_test <- numeric(length(y_test))
  
  for (i in 1:length(y_test)) {
    # Object where the probability for the i-th observation of the test data set
    # to be assigned to a cluster will be stored
    log_probs <- numeric(K + 1)
    
    # Vector of explanatory variables, and response variable of the i-th observation
    # of the test data
    x_i <- matrix(data = x_test[i,], nrow = 1, ncol = p, byrow = FALSE)
    y_i <- y_test[i]
    
    for (k in 1:K) {
      n_k <- sum(xi_hat == k) # Number of observations in the k-th cluster
      log_probs[k] <- log(n_k) + dnorm(y_i, mean = x_i%*%beta[,k], sd = sqrt(sigma2[k]), log = TRUE)
    }
    
    # Sample beta from its prior distribution
    lambda_prior <- rgamma(n = 1, shape = g, rate = h)
    tau_prior <- rexp(n = p, rate = (0.5*lambda_prior))
    beta_prior <- c(mvtnorm::rmvnorm(n = 1, mean = rep(0, p), sigma = (diag(x = tau_prior, p))))
    
    # Sample sigma2 from its prior distribution
    sigma2_prior <- 1/rgamma(n = 1, shape = e, rate = f)
    
    # Probability for the i-th observation to be assigned to a new cluster K + 1
    log_probs[K + 1] <- log(alpha) + dnorm(y_i, mean = x_i%*%beta_prior, sd = sqrt(sigma2_prior), log = TRUE)
    
    # Sample the cluster assignment for the i-th observation of the test set
    cluster <- sample(1:(K + 1), size = 1, prob = exp(log_probs - max(log_probs)))
    xi_test[i] <- cluster
  }
  return(xi_test = xi_test)
}

# Cross validation

cross_validation <- cross_validation <- function(y, x, n, p, 
                                                 a, b, c, d,
                                                 e, f, g, h, l, m,
                                                 n_burn, n_sams, n_skip){
  # Create the train and test dataset
  index <- sample(1:n, size = 0.7*n)
  
  y_train <- y[index]; x_train <- x[index,] # Train dataset
  y_test <- y[-index]; x_test <- x[-index,] # Test dataset
  
  # Objects where the mean absolute error, and the mean squared prediction error will be stored
  mape <- NULL
  mspe <- NULL
  
  n_train <- length(y_train)
    
  # Fit Bayesian Nonparametric Ridge and LASSO regression models
  M4 <- Gibbs_ridgenp(y_train, x_train, n_burn, n_sams, n_skip, 
                      a, b, c, d, l, m, n_train, p, verbose = TRUE)
    
  M5 <- Gibbs_lassonp(y_train, x_train, n_burn, n_sams, n_skip = 10, 
                      e, f, g, h, l, m, n_train, p, verbose = TRUE)
    
  # Posterior number of clusters
  ridge_k <- hat_k(M4); ridge_k_pos <- ridge_k$k_pos
  lasso_k <- hat_k(M5); lasso_k_pos <- lasso_k$k_pos
  
  # Posterior mean for beta and sigma2
  beta_ridge <- hat_np(M4, ridge_k_pos, p)$beta
  sigma2_ridge <- hat_np(M4, ridge_k_pos, p)$sigma2
    
  beta_lasso <- hat_np(M5, lasso_k_pos, p)$beta
  sigma2_lasso <- hat_np(M5, lasso_k_pos, p)$sigma2
  
  # Posterior mean for alpha
  alpha_ridge <- mean(M4$ALPHA)
  alpha_lasso <- mean(M5$ALPHA)
  
  # Cluster assignment for train dataset
  xi_hat_ridge <- matrix_A(M4, ridge_k_pos, n_train)
  xi_hat_lasso <- matrix_A(M5, lasso_k_pos, n_train)
  
  # Cluster assignment for test dataset
  xi_ridge <- out_sample_ridge(y_test, x_test, ridge_k_pos, xi_hat_ridge, 
                               beta_ridge, sigma2_ridge, alpha_ridge, a, b, c, d, p)
  xi_lasso <- out_sample_lasso(y_test, x_test, lasso_k_pos, xi_hat_lasso, 
                               beta_lasso, sigma2_lasso, alpha_lasso, e, f, g, h, p)
  
  # Linear predictor
  y_hat_ridge <- numeric(length(y_test))
  y_hat_lasso <- numeric(length(y_test))
  
  for (i in 1:length(y_test)) {
    if (max(xi_ridge) > max(xi_hat_ridge)) {
      # Since an observation of the test dataset is assigned to a new cluster, 
      # beta is samples from its prior distribution
      lambda <- rgamma(n = 1, shape = c, rate = d)
      beta_new <- c(mvtnorm::rmvnorm(n = 1, mean = rep(0, p), sigma = (1/lambda)*(diag(x = 1, p))))
      
      beta_ridge <- cbind(beta_ridge, beta_new)
    } else{
      beta_ridge <- beta_ridge
    } 
    
    if (max(xi_lasso) > max(xi_hat_lasso)) {
      # Since an observation of the test dataset is assigned to a new cluster, 
      # beta is samples from its prior distribution
      lambda <- rgamma(n = 1, shape = g, rate = h)
      tau <- rexp(n = p, rate = (0.5*lambda))
      beta_new <- c(mvtnorm::rmvnorm(n = 1, mean = rep(0, p), sigma = (diag(x = tau, p))))
      
      beta_lasso <- cbind(beta_lasso, beta_new)
    } else{
      beta_lasso <- beta_lasso
    }
    
    y_hat_ridge[i] <- x_test[i,]%*%beta_ridge[,xi_ridge[i]]
    y_hat_lasso[i] <- x_test[i,]%*%beta_lasso[,xi_lasso[i]]
  }
    
  # Compute mean absolute error
  mape_ridge <- mean(abs(y_test - y_hat_ridge))
  mape_lasso <- mean(abs(y_test - y_hat_lasso))
    
  # Compute mean squared prediction error
  mspe_ridge <- mean((y_test - y_hat_ridge)^2)
  mspe_lasso <- mean((y_test - y_hat_lasso)^2)
    
  mape <- rbind(mape, c(mape_ridge, mape_lasso))
  mspe <- rbind(mspe, c(mspe_ridge, mspe_lasso))
  
  return(list(M4 = M4, M5 = M5, y_test = y_test, y_train = y_train,
              x_test = x_test, x_train = x_train, mape = mape,
              mspe = mspe))
}

# 3. Hyperparameter elicitation

a <- 3; e <- 3 # Shape parameter of inverse-gamma distribution

b <- 2; f <- 2 # Scale parameter of inverse-gamma distribution

c <- 1; g <- 1 # Shape parameter of gamma distribution

d <- 1; h <- 1 # Rate parameter of gamma distribution

l <- 1; m <- 1

# 4. Results

simulation_results <- function(data,
                               a, b, c, d,
                               e, f, g, h, l, m,
                               n_burn, n_sams, n_skip){
  # Number of sceneries
  scenery <- length(data$Y)
  
  # Objects where the results will be stored
  RIDGE <- list("vector", scenery)
  LASSO <- list("vector", scenery)
  NUM_K <- NULL
  HAT_BNR <- list("vector", scenery)
  HAT_BNL <- list("vector", scenery)
  # Objects where the evaluation metrics will be stored
  INF_CRI <- NULL
  CRO_VAL <- NULL
  CRO_VAL_ALL <- list("vector", scenery)
  
  Y <- data$Y
  X <- data$X
  BETA_TRUE <- data$BETA
  
    for (i in 1:24) {
    y <- Y[[i]] # Set the response variable
    x <- X[[i]] # Set the matrix containing the explanatory variables
    n <- as.numeric(length(y)) # Sample size
    
    beta_true <- BETA_TRUE[[i]][1,]
    p <- as.numeric(length(beta_true)) # Number of explanatory variables
    
    # Before executing the following code, Bayesian Nonparametric Ridge.r
    # and Bayesian Nonparametric LASSO.r must be executed,
    # since it displays the Gibbs sampling algorithm
    
    M4 <- Gibbs_ridgenp(y, x, 
                        n_burn, n_sams, n_skip, 
                        a, b, c, d, l, m, n, p, verbose = TRUE)
    
    M5 <- Gibbs_lassonp(y, x, 
                        n_burn, n_sams, n_skip, 
                        e, f, g, h, l, m, n, p, verbose = TRUE)
    
    # Save the Bayesian Nonparametric Ridge and LASSO regression models results
    RIDGE[[i]] <- M4
    LASSO[[i]] <- M5
    
    # Bayesian Inference on the number of clusters
    ridge_k <- hat_k(M4); ridge_k_pos <- ridge_k$k_pos
    lasso_k <- hat_k(M5); lasso_k_pos <- lasso_k$k_pos
    
    ridge_k <- data.frame(Escenario = i, Modelo = "M4",
                          K = ridge_k$K[,1], Frecuencia = ridge_k$K[,2])
    
    lasso_k <- data.frame(Escenario = i, Modelo = "M5",
                          K = lasso_k$K[,1], Frecuencia = lasso_k$K[,2])
    
    
    NUM_K <- rbind(NUM_K, ridge_k, lasso_k)
    
    ridge <- hat_np(M4, ridge_k_pos, p)
    lasso <- hat_np(M5, lasso_k_pos, p)
    
    HAT_BNR[[i]] <- cbind(i, 1:(p + 1), round(ridge$inference, 2))
    HAT_BNL[[i]] <- cbind(i, 1:(p + 1), round(lasso$inference, 2))
    
    # Compute the Watanabe-Akaike Information Criterion
    ridge_waic <- compute_WAIC_np(M4, y, x, n)
    lasso_waic <- compute_WAIC_np(M5, y, x, n)
    
    inf_cri <- data.frame(`Escenario` = i, `WAIC Ridge` = ridge_waic, `WAIC LASSO` = lasso_waic)
    
    INF_CRI <- rbind(INF_CRI, round(inf_cri, 2))
    
    # Cross validation
    cro_val <- cross_validation(y, x, n, p, a, b, c, d, e, f, g, h, l, m,
                                n_burn, n_sams, n_skip)
    
    CRO_VAL <- rbind(CRO_VAL, c(i, cro_val$mape, cro_val$mspe))
    CRO_VAL_ALL[[i]] <- cro_val
    }
  return(list(RIDGE = RIDGE, LASSO = LASSO, HAT_BNR = HAT_BNR, HAT_BNL = HAT_BNL,
              NUM_K = NUM_K, INF_CRI = INF_CRI, CRO_VAL = CRO_VAL,
              CRO_VAL_ALL = CRO_VAL_ALL))
}

results <- simulation_results(data, a, b, c, d, e, f, g, h, l, m,
                              n_burn = 1000, n_sams = 10000, n_skip = 10)