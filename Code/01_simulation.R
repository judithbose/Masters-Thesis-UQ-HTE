## This script is used to generate the simulated data used in Chapter 4.2 
## of the thesis.

# path 
path <- getwd()
path <- paste0(path, "/Code")

# function for generating simulation data
generate_data <- function(n, design = "linear", seed) {
  
  set.seed(seed)
  
  # generate covariates (standard normal)
  d <- 20
  X <- matrix(rnorm(n * d, mean = 0, sd = 1), nrow = n, ncol = d)
  colnames(X) <- paste0("X", 1:d)
  
  # randomize treatment
  Trt <- rbinom(n, size = 1, prob = 0.5)
  
  # generate error term (standard normal)
  epsilon <- rnorm(n, mean = 0, sd = 1)

  if (design == "linear") { # linear design

    mu_0 <- 2 * X[, 1] + X[, 2]
    tau <- X[, 1]
    
  } else if (design == "nonlinear") { # non-linear design

    mu_0 <- sin(pi * X[, 1]) + 2 * (X[, 2] - 0.5)^2
    tau <- sin(pi * X[, 1])
    
  } else if (design == "interaction") { # interaction design

    mu_0 <- X[, 1] * X[, 2]
    tau <- X[, 3] * X[, 4]
  }
  
  # compute observed outcome
  Y <- mu_0 + Trt * tau + epsilon
  
  # return data
  data <- data.frame(X, Trt = Trt, Y = Y, true_cate = tau)
  return(data)
}

# create 6 different scenarios
for (iter in 1:100) {
  
  data <- generate_data(seed = 30092026 + iter, n = 600, design = "linear")
  save(data, file = paste0(path, "/../Data/Simulation/data_L600/data_", iter, ".RData"))
  
  data <- generate_data(seed = 30092026 + iter, n = 1200, design = "linear")
  save(data, file = paste0(path, "/../Data/Simulation/data_L1200/data_", iter, ".RData"))
  
  data <- generate_data(seed = 30092026 + iter, n = 600, design = "nonlinear")
  save(data, file = paste0(path, "/../Data/Simulation/data_N600/data_", iter, ".RData"))
  
  data <- generate_data(seed = 30092026 + iter, n = 1200, design = "nonlinear")
  save(data, file = paste0(path, "/../Data/Simulation/data_N1200/data_", iter, ".RData"))
  
  data <- generate_data(seed = 30092026 + iter, n = 600, design = "interaction")
  save(data, file = paste0(path, "/../Data/Simulation/data_I600/data_", iter, ".RData"))
  
  data <- generate_data(seed = 30092026 + iter, n = 1200, design = "interaction")
  save(data, file = paste0(path, "/../Data/Simulation/data_I1200/data_", iter, ".RData"))

}
