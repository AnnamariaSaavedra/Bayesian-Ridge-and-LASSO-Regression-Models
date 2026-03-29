rm(list=ls()); set.seed(123)

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(coda)))
suppressMessages(suppressWarnings(library(ggplot2)))

# 2. Import database

Data <- read_csv("~/Trabajo de grado/Database - Case study 2.txt")

# 2.1 Select the response variable and explanatory variables

y <- Data$Wage # Set the response variable

x <- Data %>%
  dplyr::select(-c(Wage)) %>% # Set the matrix containing the explanatory variables
  scale(center = TRUE, scale = TRUE) %>% # Standardize the explanatory variables
  as.matrix()

x <- cbind(x, 1) # Create the intercept column
  
n <- length(y) # Sample size

p <- ncol(x) # Number of explanatory variables

# 3. Bayesian Ridge

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

# Before executing the following code, Bayesian Ridge.r must be executed, since it displays
# the Gibbs sampling algorithm

# 3.2 Gibbs sampling algorithm implementation

M2 <- Gibbs_ridge(y, x, a, b, c, d, n, p,
                  n_skip = 10, # Accounting for Markov chain autocorrelation will require systematic sampling,
                  n_sams = 10000, # Set the number of effective samples
                  n_burn = 1000) # Set the number of burn-in samples

# Compute the effective sample size for model parameters

TEM_beta <- coda::effectiveSize(M2$BETA); summary(TEM_beta) # beta

TEM_sigma2 <- coda::effectiveSize(M2$SIGMA); summary(TEM_sigma2) # sigma2

TEM_lambda <- coda::effectiveSize(M2$LAMBDA); summary(TEM_lambda) # lambda

# Compute the Monte Carlo standard error for model parameters

EEMC_beta <- apply(X = M2$BETA, MARGIN = 2, FUN = sd)/sqrt(TEM_beta); round(summary(EEMC_beta), 4) # beta

EEMC_sigma2 <- sd(M2$SIGMA)/sqrt(TEM_sigma2); round(summary(EEMC_sigma2), 4) # sigma2

EEMC_lambda <- sd(M2$LAMBDA)/sqrt(TEM_lambda); round(summary(EEMC_lambda), 4) # lambda

# 4. Display log-likelihood chain

plot(M2$LL, type = "p", pch = ".", cex = 1.1, col = "#00CD66", xlab = "Iteración", ylab = "Log-verosimilitud", main = "",)
abline(h = mean(M2$LL), lwd = 3, col = "#00CD66")

# 5. Bayesian inference

# Bayesian inference for beta

BETA_MEAN <- round(apply(M2$BETA, MARGIN = 2, FUN = mean), 3) # Posterior mean

BETA_MEDIAN <- round(apply(M2$BETA, MARGIN = 2, FUN = median), 3) # Posterior median

BETA_SD <- round(apply(M2$BETA, MARGIN = 2, FUN = sd), 3) # Posterior standard deviation

CI_BETA <- round(apply(M2$BETA, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)), 3) # 95% credible interval

# Bayesian inference for sigma2

SIGMA2_MEAN <- round(mean(M2$SIGMA), 3) # Posterior mean

SIGMA2_MEDIAN <- round(median(M2$SIGMA), 3) # Posterior median

SIGMA2_SD <- round(sd(M2$SIGMA), 3) # Posterior standard deviation

CI_SIGMA <- round(quantile(x = M2$SIGMA, probs = c(0.025, 0.975)), 3) # 95% credible interval

# Bayesian inference for lambda

LAMBDA_MEAN <- round(mean(M2$LAMBDA), 3) # Posterior mean

LAMBDA_MEDIAN <- round(median(M2$LAMBDA), 3) # Posterior median

LAMBDA_SD <- round(sd(M2$LAMBDA), 3) # Posterior standard deviation

CI_LAMBDA <- round(quantile(x = M2$LAMBDA, probs = c(0.025, 0.975)), 4) # 95% credible interval

# 6. Compute information criterion and cross validation

# Deviance Information Criterion

LL_HAT <- sum(dnorm(x = y, mean = x%*%BETA_MEAN, sd = sqrt(SIGMA2_MEAN), log = TRUE))

LL_B <- M2$LL

pDIC <- 2*(LL_HAT - mean(LL_B))

DIC <- -2*LL_HAT + 2*pDIC

# Watanabe-Akaike Information Criterion

LPPD <- 0

pWAIC <- 0

for (i in 1:n) {
  # LPPD
  TMP <- dnorm(x = y[i], mean = x[i,]%*%t(M2$BETA), sd = sqrt(M2$SIGMA))
  LPPD <- LPPD + log(mean(TMP))
  # pWAIC
  TMP_2 <- dnorm(x = y[i], mean = x[i,]%*%t(M2$BETA), sd = sqrt(M2$SIGMA), log = TRUE)
  pWAIC <- pWAIC + 2*(log(mean(TMP)) - mean(TMP_2))
}

WAIC <- -2*LPPD + 2*pWAIC

# Cross validation

cross_validation <- function(y, n, p, a, c, d){
  x <- Data %>%
    dplyr::select(-c(Wage)) %>% # Set the matrix containing the explanatory variables
    as.matrix()

  x <- cbind(x, 1) # Create the intercept column
  
  # Create the train and test dataset
  index <- sample(1:n, size = 0.7*n)
  
  y_train <- y[index]; x_train <- x[index,] # Train dataset
  y_test <- y[-index]; x_test <- x[-index,] # Test dataset
  
  # Objects where the mean absolute error, and the mean squared prediction error will be stored
  mape <- NULL
  mspe <- NULL
  
  n <- length(y_train)
  
  # Hyperparameter elicitation
  beta_OLS <- solve(t(x_train)%*%x_train)%*%t(x_train)%*%y_train
  residuals <- y_train - x_train%*%beta_OLS
  sigma2_OLS <- sum(residuals^2)/(n - p)
  
  x <- Data %>%
    dplyr::select(-c(Wage)) %>% # Set the matrix containing the explanatory variables
    scale(center = TRUE, scale = TRUE) %>% # Standardize the explanatory variables
  
  x <- cbind(x, 1) # Create the intercept column
  
  x_train <- x[index,] # Standardized train dataset
  x_test <- x[-index,] # Standardized test dataset
  
  # Fit Bayesian Ridge regression model
  M2 <- Gibbs_ridge(y_train, x_train, a, b = a*sigma2_0, c, d, n, p,
                    n_skip = 10, # Accounting for Markov chain autocorrelation will require systematic sampling,
                    n_sams = 10000, # Set the number of effective samples
                    n_burn = 1000) # Set the number of burn-in samples
  
  # Posterior mean for beta
  beta <- colMeans(M2$BETA)
  
  # Linear predictor
  y_hat <- x_test%*%beta
  
  # Compute mean absolute prediction error
  mape <- mean(abs(y_test - y_hat))
  
  # Compute mean squared prediction error
  mspe <- mean((y_test - y_hat)^2)
  
  return(list(mape = mape, mspe = mspe))
}

cross_validation_M2 <- cross_validation(y, n, p, a, c, d)

# 7. Monte Carlo samples from the posterior predictive distribution of test statistics

# Create test statistics function

test_stats <- function(x) {
  c(mean = mean(x),
    median = median(x),
    sd = sd(x),
    iqr = diff(quantile(x, c(0.25, 0.75))),
    min = min(x),
    max = max(x)
  )
}

ts_display <- c("Media", "Mediana", "Desviación estándar", "Rango intercuartílico",
                "Mínimo", "Máximo")

ts <- NULL # Object where the test statistics will be stored

# Simulated statistics

for (b in 1:length(M2$SIGMA)) {
  # Samples from the posterior distribution
  beta <- M2$BETA[b, ]
  sigma2 <- M2$SIGMA[b]
  
  # Posterior predictive datasets, each of size n
  y_tilde <- rnorm(n = n, mean = x%*%beta, sd = sqrt(sigma2))
  
  # Samples from the posterior predictive distribution of test statistics
  ts <- rbind(ts, test_stats(y_tilde)) # Compute test statistics
}

ts_hat <- test_stats(y) # Observed test statistics

# Comparison plots between simulated and observed test statistics

par(mfrow = c(3, 2), mar = c(3, 3, 2, 1), mgp = c(1.75, 0.75, 0))
for (j in 1:length(ts_hat)) {
  test_statistics  <- ts[, j]
  test_statistics_hat <- ts_hat[j]
  
  # Plot histogram
  hist(
    x = test_statistics, freq = FALSE, nclass = 30,
    col = "gray90", border = "gray90",
    xlab = ts_display[j], ylab = "Densidad",
    main = ts_display[j]
  )
  
  abline(
    v = quantile(test_statistics, c(0.025, 0.5, 0.975)),
    col = c(4, 2, 4), lty = c(4, 2, 4), lwd = c(2, 1, 2)
  )
  
  abline(v = test_statistics_hat, lwd = 2)
}

# Compute posterior predictive p-value

ppp <- NULL

for (j in 1:length(ts_hat)) {
  ppp[j] <- round(mean(ts[,j] < ts_hat[j]), 4)
}

# 8. Posterior predictive density estimate

posterior_density_estimate <- function(model, 
                                       y_seq, # Define a sequence of y values for density estimation
                                       x # Define a grid of x values
                                       ) {
  B <- length(model$SIGMA)
  M <- length(y_seq)
  
  # Object where density estimates will be stored
  FE <- matrix(NA, nrow = B, ncol = M)
  
  for (b in 1:B) {
    beta_b <- model$BETA[b,]
    sigma2_b <- model$SIGMA[b]
    
    for (i in 1:M) {
      FE[b, i] <- dnorm(y_seq[i], mean = x%*%beta_b, sd = sqrt(sigma2_b))
    }
  }
  
  f_hat <- colMeans(FE)  # Posterior mean density
  f_inf <- apply(FE, 2, quantile, probs = 0.025)  # 2.5% credible interval
  f_sup <- apply(FE, 2, quantile, probs = 0.975)  # 97.5% credible interval
  
  return(list(f_hat = f_hat, f_inf = f_inf, f_sup = f_sup))
}

# Define a sequence of y values for density estimation
y_seq <- seq(min(y), max(y), length.out = 150)

# Define a grid of x values
x_seq <- colMeans(x)

# Compute posterior density estimate and credible intervals
density_estimate <- posterior_density_estimate(M2, y_seq, x_seq)

f_hat <- density_estimate$f_hat
f_inf <- density_estimate$f_inf
f_sup <- density_estimate$f_sup

# Plot the histogram
hist(x = y, freq = FALSE, xlim = c(11, 18), ylim = c(0, 1),
     ylab = "Densidad", main = "",
     col = alpha("grey", 0.3), cex.lab = 1.5, cex.axis  = 1.5)
# Overlay the posterior density estimate as a blue line
polygon(c(y_seq, rev(y_seq)), c(f_inf, rev(f_sup)), col = alpha("#00CD66", 0.3), border = NA)
lines(y_seq, f_hat, lwd = 2, col = "#00CD66")
