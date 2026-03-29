rm(list=ls()); set.seed(123)

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(ggplot2)))

# 2. Import database

Data <- read_csv("~/Downloads/Database - Case study 2.txt")

# 2.1 Select the response variable and explanatory variables

y <- Data$Wage # Set the response variable

x <- Data %>%
  select(-c(Wage)) %>% # Set the matrix containing the explanatory variables
  scale(center = TRUE, scale = TRUE) %>% # Standardize the explanatory variables
  as.matrix()

x <- cbind(x, 1) # Create the intercept column
  
n <- length(y) # Sample size

p <- ncol(x) # Number of explanatory variables

# 3. Bayesian Nonparametric Ridge regression

# 3.1 Hyperparameter elicitation

x_b <- Data %>%
    dplyr::select(-c(Wage)) %>% # Set the matrix containing the explanatory variables
    as.matrix()

x_b <- cbind(x_b, 1) # Create the intercept column

beta_OLS <- solve(t(x_b)%*%x_b)%*%t(x_b)%*%y

residuals <- y - x_b%*%beta_OLS

sigma2_OLS <- sum(residuals^2)/(n - p)

a <- 3 # Shape parameter of inverse-gamma distribution

b <- a*sigma2_OLS # Scale parameter of inverse-gamma distribution

c <- 1 # Shape parameter of gamma distribution

d <- 1 # Rate parameter of gamma distribution

l <- 1 # Shape parameter of gamma distribution

m <- 1 # Shape parameter of gamma distribution

# Before executing the following code, Bayesian Nonparametric LASSO.r must be executed, since it displays
# the Gibbs sampling algorithm

# 3.2 Gibbs sampling algorithm implementation

M4 <- Gibbs_ridgenp(y, x,
                    n_burn = 1000, # Set the number of burn-in samples
                    n_sams = 10000, # Set the number of effective samples
                    n_skip = 10, # Accounting for Markov chain autocorrelation will require systematic sampling,
                    a, b, c, d, l, m, n, p)

# 4. Bayesian inference

# 4.1 Display the log-likelihood chain

plot(M4$LL, type = "p", pch = ".", cex = 1.1, col = "darkorchid3", xlab = "Iteración", ylab = "Log-verosimilitud", main = "")
abline(h = mean(M4$LL), lwd = 3, col = "darkorchid3")

# 4.2 Inference on the number of clusters

K <- apply(M4$XI, 1, function(x) length(unique(x))) # Compute the number of clusters at each iteration

K_table <- as.data.frame(table(K)/length(K))

# Plot the posterior distribution of K

ggplot(K_table, aes(x = K, y = Freq)) +
  geom_segment(aes(x = K, xend = K, y = 0, yend = Freq),
               color = "cyan2", lwd = 1.5) +
  ylim(c(0, 1)) +
  ylab("Densidad") + xlab("Número de clústeres") +
  theme_bw(base_size = 16) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Compute Incidence Matrix

matrix_A <- function(model, K, n){
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
  xi_hat <- Mclust(data = coord, G = 1:K)$classification
  
  return(xi_hat = xi_hat)
}

# 4.3 Inference for Bayesian Nonparametric Ridge

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
  inf <- matrix(data = NA, nrow = p + 1, ncol = 5*K)
  id <- function(k)(5*(k - 1) + 1):(5*k)
  
  for (j in 1:p) {
    beta_mean <- colMeans(BETA_corrected[[j]])
    beta_median <- apply(BETA_corrected[[j]], MARGIN = 2, FUN = median)
    beta_sd <- apply(BETA_corrected[[j]], MARGIN = 2, FUN = sd)
    beta_ic <- apply(X = BETA_corrected[[j]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
    
    for (k in 1:K) {
      inf[j, id(k)] <- c(beta_mean[k], beta_median[k], beta_sd[k], beta_ic[1, k], beta_ic[2, k])
    }
  }
  sigma2_mean <- colMeans(SIGMA2_corrected)
  sigma2_median <- apply(SIGMA2_corrected, MARGIN = 2, FUN = median)
  sigma2_sd <- apply(SIGMA2_corrected, MARGIN = 2, FUN = sd)
  sigma2_ic <- apply(SIGMA2_corrected, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
  
  for (k in 1:K) {
    inf[p + 1, id(k)] <- c(sigma2_mean[k], sigma2_median[k], sigma2_sd, sigma2_ic[1, k], sigma2_ic[2, k])
  }
  
  col <- seq(from = 1, to = ncol(inf), by = 5) # Position of the posterior mean
  beta <- inf[c(1:p), col]
  sigma2 <- inf[(p + 1), col]
  
  return(list(inference = inf, beta = as.matrix(beta), sigma2 = sigma2))
}

inference <- hat_np(M4, K = 1, p)

# Bayesian inference for lambda

LAMBDA_MEAN <- round(mean(M4$LAMBDA), 4) # Posterior mean

LAMBDA_MEDIAN <- round(median(M4$LAMBDA), 4) # Posterior median

LAMBDA_SD <- round(sd(M4$LAMBDA), 4) # Posterior standard deviation

CI_LAMBDA <- round(quantile(x = M4$LAMBDA, probs = c(0.025, 0.975)), 4) # 95% credible interval

# 4.3.4 Bayesian inference for alpha

ALPHA_MEAN <- round(mean(M4$ALPHA), 4) # Posterior mean

ALPHA_MEDIAN <- round(median(M4$ALPHA), 4) # Posterior median

ALPHA_SD <- round(sd(M4$ALPHA), 4) # Posterior standard deviation

CI_ALPHA <- round(quantile(x = M4$ALPHA, probs = c(0.025, 0.975)), 4) # 95% credible interval

# 4.4 Compute information criterion

# Watanabe-Akaike Information Criterion

ite <- nrow(M4$XI)
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

# 4.5 Cross validation

# Out of sample prediction

out_sample <- function(y_test, x_test, K, xi_hat, beta, sigma2, alpha, 
                       a, b, c, d){
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

# Cross validation

cross_validation <- cross_validation <- function(y, p, 
                                                 a, b, c, d, l, m){
  x <- Data %>%
    dplyr::select(-c(Wage)) %>% # Set the matrix containing the explanatory variables
    as.matrix()
  x <- cbind(x, 1) # Create the intercept column

  # Create the train and test dataset
  index <- sample(1:n, size = 0.7*n)

  y_train <- y[index]; x_train <- x[index,] # Train dataset
  y_test <- y[-index]; x_test <- x[-index,] # Test dataset

  # Hyperparameter elicitation
  beta_OLS <- solve(t(x_train)%*%x_train)%*%t(x_train)%*%y_train
  residuals <- y_train - x_train%*%beta_OLS
  sigma2_OLS <- sum(residuals^2)/(n - p)

  x <- Data %>%
  select(-c(Wage)) %>% # Set the matrix containing the explanatory variables
  scale(center = TRUE, scale = TRUE) %>% # Standardize the explanatory variables
  as.matrix()

  x <- cbind(x, 1) # Create the intercept column

  x_train <- x[index,] # Standardized train dataset
  x_test <- x[-index,] # Standardized test dataset

  # Objects where the mean absolute error, and the mean squared prediction error will be stored
  mape <- NULL
  mspe <- NULL

  n <- length(y_train)

  # Fit Bayesian Nonparametric Ridge regression model
  M4 <- Gibbs_ridgenp(y_train, x_train, 
                      n_burn = 1000, 
                      n_sams = 10000, 
                      n_skip = 10, 
                      a, b = a*sigma2_OLS, c, d, l, m, n, p, verbose = TRUE)
    
  # Inference on the number of clusters
  K <- apply(M4$XI, 1, function(x) length(unique(x)))
  K_table <- as.data.frame(table(K)/length(K))
    
  # Posterior number of clusters
  k_pos <- as.numeric(K_table$K[which.max(K_table$Freq)])
    
  # Posterior mean for beta and sigma2
  beta <- hat_np(M4, k_pos, p)$beta
  sigma2 <- hat_np(M4, k_pos, p)$sigma2
    
  # Posterior mean for alpha
  alpha <- mean(M4$ALPHA)
    
  # Cluster assignment for train dataset
  xi_hat <- matrix_A(M4, k_pos, n)
    
  # Cross validation
  xi <- out_sample(y_test, x_test, k_pos, xi_hat, beta, sigma2, alpha,
                   a, b = a*sigma2_OLS, c, d)
    
  # Linear predictor
  y_hat_ridge <- numeric(length(y_test))
  for (i in 1:length(y_test)) {
    y_hat_ridge[i] <- x_test[i,]%*%beta[,xi[i]]
  }
    
  # Compute mean absolute prediction error and mean squared prediction error
  mape <- mean(abs(y_test - y_hat_ridge))
  mspe <- mean((y_test - y_hat_ridge)^2)

  return(list(mape = mape, mspe = mspe))
}

CV <- cross_validation(y, p, a, b, c, d, l, m)

# 4.6 Bayesian inference for density function

posterior_density_estimate <- function(model, y_seq, x) {
  B <- length(model$BETA)
  n <- ncol(model$XI)
  M <- length(y_seq)
  
  # Store density estimates
  FE <- matrix(NA, nrow = B, ncol = M)
  
  for (b in 1:B) {
    xi <- model$XI[b, ]
    beta <- model$BETA[[b]]
    sigma2 <- model$SIGMA[[b]]
    
    # Unique cluster labels
    cluster_labels <- sort(unique(xi))
    cluster_counts <- as.numeric(table(xi)[as.character(cluster_labels)])
    
    # Cluster probabilities
    weight_b <- cluster_counts / sum(cluster_counts)
    
    for (i in 1:M) {
      FE[b, i] <- sum(weight_b * dnorm(y_seq[i], mean = x%*%t(beta), sd = sqrt(sigma2)))
    }
  }
  
  f_hat <- colMeans(FE)  # Posterior mean density
  f_inf <- apply(FE, 2, quantile, probs = 0.025)  # 2.5%  credible interval
  f_sup <- apply(FE, 2, quantile, probs = 0.975)  # 97.5% credible interval
  
  return(list(f_hat = f_hat, f_inf = f_inf, f_sup = f_sup))
}

# Define a sequence of y values for density estimation
y_seq <- seq(min(y), max(y), length.out = 150)

# Define a grid of x values
x_seq <- colMeans(x)

# Compute posterior density estimate and credible intervals
density_estimate <- posterior_density_estimate(M4, y_seq, x_seq)

f_hat <- density_estimate$f_hat
f_inf <- density_estimate$f_inf
f_sup <- density_estimate$f_sup

# Plot the histogram
hist(x = y, freq = FALSE,
     ylab = "Densidad", main = "",
     col = alpha("grey", 0.3), cex.label = 1.5, cex.axis  = 1.5)
# Overlay the posterior density estimate as a blue line
polygon(c(y_seq, rev(y_seq)), c(f_inf, rev(f_sup)), col = alpha("darkorchid3", 0.3), border = NA)
lines(y_seq, f_hat, lwd = 2, col = "darkorchid3")
