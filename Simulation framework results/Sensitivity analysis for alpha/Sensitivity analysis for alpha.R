# 0. Set seed

rm(list=ls()); set.seed(123)

# 1. Load necessary library

suppressMessages(suppressWarnings(library(ggplot2)))

# 2. Simulation framework to evaluate hyperparameter elicitation

n <- 200 # Number of observations

K <- 2 # Number of clusters

p <- 10 # Number of parameters

p_0 <- 4 # Number of non-zero elements

sigma2_sim <- c(1, 3) # Variance parameter

beta_sim <- matrix(data = c(1, 2, 4, 8, rep(0, 6),
                            -1, -3, -9, -18, rep(0, 6)),
                   nrow = K, ncol = p, byrow = TRUE) # Mean parameter

omega <- rep(1/K, K) # Vector of mixing proportions
xi <- sample(x = 1:K, size = n, replace = TRUE, prob = omega) # Vector of cluster assignment

n_k <- table(xi) # Size of cluster k

x <- matrix(data = NA, nrow = n, ncol = p) # Matrix of explanatory variables

x[,1] <- rep(1, n) # Intercept column
x[,2] <- sample(x = 0:1, size = n, replace = TRUE, prob = c(0.5, 0.5)) # Dummy variable
x[,3] <- rbeta(n = n, shape1 = 1, shape2 = 1) # Beta-distributed variable
x[,4] <- rpois(n = n, lambda = 1) # Poisson-distributed variable

for (j in 5:p) {
  x[,j] <- rnorm(n = n, mean = seq(from = -3, to = 3, by = 1), sd = 2) # Normal-distributed variables
}

# Objects where the response variable will be stored
y <- numeric()
f_true <- numeric()

# Sample response variable for each cluster
for (k in 1:K) {
  y[xi == k] <- rnorm(n = n_k[k], mean = x[xi == k,]%*%beta_sim[k,], sd = sqrt(sigma2_sim[k]))
  f_true[xi == k] <- dnorm(x = y[xi == k], mean = x[xi == k,]%*%beta_sim[k,], sd = sqrt(sigma2_sim[k]))
}

# Before executing the following code, load the Bayesian Nonparametric Ridge and LASSO
# regression models results

load("")

# 3. Plot the posterior distribution of K

K <- apply(hyperparameter$XI, 1, function(x) length(unique(x))) # Compute the number of clusters at each iteration

K_table <- as.data.frame(table(K)/length(K))

ggplot(K_table, aes(x = K, y = Freq)) +
  geom_segment(aes(x = K, xend = K, y = 0, yend = Freq),
               color = "chartreuse4", lwd = 1.5) +
  ylim(c(0, 1)) +
  ylab("Densidad") + xlab("Número de clústeres") +
  theme_bw(base_size = 20) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# 4. Bayesian inference for K =

hat_np <- function(model, K, p){
  ite <- nrow(model$XI)
  permu <- gtools::permutations(n = K, r = K) # Permutations
  
  # Objects where the samples of beta, and sigma2 will be stored
  BETA_corrected <- vector("list", p)
  SIGMA2_corrected <- NULL
  
  # Posterior distribution of K
  k <- apply(X = model$XI, MARGIN = 1, function(x) length(unique(x)))
  k_tab <- table(factor(x = k, levels = k, labels = k))
  
  # Objects where the posterior mean of beta, and sigma2 will be stored
  beta_pos <- matrix(0, nrow = K, ncol = p)
  sigma2_pos <- 0
  
  for (b in 1:ite) {
    if (length(unique(model$XI[b,])) == K) {
      beta_pos <- beta_pos +
        model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE] / max(k_tab)
      
      sigma2_pos <- sigma2_pos + model$SIGMA[[b]] / max(k_tab)
    }
  }
  
  # Average over the permuted spaces
  for (b in 1:ite) {
    if (length(table(model$XI[b,])) == K) {
      beta_current <- model$BETA[[b]][sort(unique(model$XI[b,])), , drop = FALSE]
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
      
      sigma2_current <- model$SIGMA[[b]]
      # Reorder according to the permutations, and compute the distance of each sample to its posterior mean
      dist2 <- apply(X = permu, MARGIN = 1, 
                     FUN = function(perm) {
                       permuted_sigma2 <- sigma2_current[perm]
                       sum((permuted_sigma2 - sigma2_pos)^2)
                     }
      )
      # Select the optimum permutation
      best_permu <- permu[which.min(dist2),]
      SIGMA2_corrected <- rbind(SIGMA2_corrected, sigma2_current[best_permu])
    }
  }
  
  # Posterior mean, and 95% credible interval
  inf <- matrix(data = NA, nrow = p + 1, ncol = 3*K)
  id <- function(k)(3*(k - 1) + 1):(3*k)
  
  for (j in 1:p) {
    beta_mean <- colMeans(BETA_corrected[[j]])
    beta_ic <- apply(X = BETA_corrected[[j]], MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
    
    for (k in 1:K) {
      inf[j, id(k)] <- c(beta_mean[k], beta_ic[1, k], beta_ic[2, k])
    }
  }
  sigma2_mean <- colMeans(SIGMA2_corrected)
  sigma2_ic <- apply(SIGMA2_corrected, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
  
  for (k in 1:K) {
    inf[p + 1, id(k)] <- c(sigma2_mean[k], sigma2_ic[1, k], sigma2_ic[2, k])
  }
  return(inference = inf)
}
