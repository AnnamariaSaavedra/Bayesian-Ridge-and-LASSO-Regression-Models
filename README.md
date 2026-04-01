# Bayesian Nonparametric LASSO and Ridge Regression Models

The GitHub repository presents the algorithms for fitting and assessing Bayesian LASSO and Ridge regression models, the design of simulation scenarios, sensitivity analyses on hyperparameter elicitation, and the fitting of the proposed regression models to two datasets. The RStudio codes enable the reader to replicate the results presented in the master's thesis titled: Non-Parametric Bayesian LASSO and Ridge Regression Models: An Approach to Identifying the Macroeconomic Determinants of Growth.

Bayesian nonparametric LASSO and Ridge regression models computational framework is based on the characterization of the Dirichlet Process mixture model, according to the Chinese Restaurant Process, and the Gibbs Sampler algorithm.

-  *Case study 2* folder contains the fitting and assessing results of five regression models over the dataset described in section 6.2. The reader will find the Monte Carlo Markov Chain samples for each Bayesian regression model.
-  *Case study* folder contains the fitting and assessing results of five regression models over the dataset described in section 6.1. The reader will find the Monte Carlo Markov Chain samples for each Bayesian regression model.
-  *Data* folder contains the case study datasets described in chapter 6 and the simulated datasets described in chapter 6 for fitting and assessing the proposed models and performing the sensitivity analyses.
-  *Gibbs Sampler* folder contains the algorithms for fitting and assessing Bayesian LASSO and Ridge regression models and G-prior Bayesian normal regression model.
-  *Simulation framework results* folder contains the sensitivity analyses results for sigma2, lambda, and alpha and the proposed models results for the simulated datasets.
-  *Theoretical framework* folder contains the RStudio codes for displaying the figures presented in chapters 2 and 3.
