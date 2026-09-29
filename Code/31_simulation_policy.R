## Source: additional ressources Hill (2011) paper
## This script is used to generate data based on the IHDP data set with more than
## two treatments. This data is used in Chapter 5.3 of the thesis.

# libraries -----
library(dplyr)

# paths ----
path <- getwd()
path <- paste0(path, "/Code")

# data ----
# load orginal data set
load(paste0(path, "/../sim.data"))

# reduce data set as in Hill (2011)
data <- imp1[!(imp1$treat==1 & imp1$momwhite==0),]
rm(imp1)

# identify covariates
covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")
covs.cat.n=c("sex","twin","b.marr","mom.lths","mom.hs",	"mom.scoll","cig","first","booze","drugs","work.dur","prenatal","ark","ein","har","mia","pen","tex","was")
covs.ols = c(covs.cont.n,covs.cat.n)
p=length(c(covs.cont.n,covs.cat.n))

# create treatment variable
Trt=data$treat

# standardize data
X = data[,covs.ols]
X[,covs.cont.n]=as.data.frame(t((t(X[,covs.cont.n])-unlist(lapply(X[,covs.cont.n],mean)))/sqrt(unlist(lapply(X[covs.cont.n],var)))))

# save dimensions for further simulations
dimx = ncol(X)
Xmat = as.matrix(X)
N = nrow(X)


# simulation ----
simulate <- function(X, seed = 123) {
  
  set.seed(seed)

  # data
  X <- as.data.frame(X)
  n <- nrow(X)
  p <- ncol(X)
  
  
  # control 
  beta_0 <- sample(x = 0:4, size = p + 1, replace = TRUE,
                   prob = c(0.50, 0.20, 0.15, 0.10, 0.05))
  mu_0 <- as.numeric(cbind(1, as.matrix(X)) %*% beta_0)
  
  sigma_0 <- 1.0
  eps_0 <- rnorm(n, mean = 0, sd = sigma_0)
  
  Y_0 <- mu_0 + eps_0
  
  
  # g(X)
  g <- 0.9 * X$bw - 0.7 * X$b.head + 0.25 * X$preterm
  
  
  # CATEs 
  tau_1 <- -1.0 + 3.0 * plogis(2.0 * (g - 0.5))
  tau_2 <- -1.0 + 3.0 * plogis(-2.0 * (g + 0.5))
  
  
  # treatment 1
  mu_1 <- mu_0 + tau_1
  
  sigma_1 <- 0.8 + 1.2 * plogis(1.5 * g)
  sigma_1 <- pmax(sigma_1, 0.1)
  
  eps_1 <- rnorm(n, mean = 0, sd = sigma_1)
  Y_1 <- mu_1 + eps_1
  

  # treatment 2
  mu_2 <- mu_0 + tau_2
  
  sigma_2 <- 0.8
  eps_2 <- rnorm(n, mean = 0, sd = sigma_2)
  
  Y_2 <- mu_2 + eps_2
  
  
  # treatment assignment 
  Trt <- sample(x = 0:2, size = n, replace = TRUE, prob = c(1/3, 1/3, 1/3))
  Y_obs <- ifelse(Trt == 0, Y_0, ifelse(Trt == 1, Y_1, Y_2))
  

  # observed data
  data_obs <- data.frame(id = seq_len(n), X,
                         Trt = factor(Trt, levels = c(0, 1, 2)), Y = Y_obs)
  
  
  # ground truth
  data_truth <- data.frame(id = seq_len(n), mu_0 = mu_0, mu_1 = mu_1, mu_2 = mu_2,
                           Y_0 = Y_0, Y_1 = Y_1, Y_2 = Y_2, tau_1 = tau_1,
                           tau_2 = tau_2, sigma_0 = sigma_0, sigma_1 = sigma_1,
                           sigma_2 = sigma_2, g = g)
  
  return(list(data_obs = data_obs, ground_truth = data_truth))
}

# loop ----
for (iter in 1:100) {
  
  data <- simulate(Xmat, seed = 30092026 + iter)
  
  IHDP <- data[[1]]
  save(IHDP, file = paste0(path, "/../Data/Policy/IHDP_", iter, ".RData"))
  
  truth <- data[[2]]
  save(truth, file = paste0(path, "/Policy/TRUTH_", iter, ".RData"))

}
