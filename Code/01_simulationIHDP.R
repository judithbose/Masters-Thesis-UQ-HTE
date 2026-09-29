## Source: additional ressources Hill (2011) paper
## This script is used to generate the outcomes for the real data used in Chapter 4.2 
## of the thesis. It is based on the code provided by Hill (2011). This has been
## rewritten solely for better readability and converted into functions to make
## it easier to implement.

path <- getwd()
path <- paste0(path, "/Code")

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


# function to simulate setting A see Hill (2011)
simulate_A <- function(sigy = 1, tau = 4, seed, m){
  set.seed(seed)
  
  betaA = sample(c(0:4),dimx+1,replace=TRUE,prob=c(.5,.2,.15,.1,.05))
  yahat = cbind(rep(1, N), Xmat) %*% betaA

  res <- vector("list", m)
  
  for (j in seq_len(m)) {
    set.seed(01042026 + j)
    
    YB0 <- rnorm(N, yahat, sigy)
    YB1 <- rnorm(N, yahat+tau, sigy)
    
    Y <- YB1
    Y[Trt == 0] <- YB0[Trt == 0]
    tauAis=4
    
    B <- cbind(X, Trt, Y, true_cate = tauAis)
    res[[j]] <- B
  }
  
  return(res)
}


# create data for marginal coverage check
for (iter in 1:100) {
  
  IHDP_list <- simulate_A(seed = 30092026 + iter, m = 1)
  
  for (j in seq_along(IHDP_list)) {
    
    data <- IHDP_list[[j]]
    save(data, file = paste0(path, "/../data/IHDP/dataA/data_", iter, "_", j, ".RData"))
  
  }
}

# create data for conditional coverage check
for (iter in 1:5) {
  
  IHDP_list <- simulate_A(seed = 30092026 + iter, m = 200)
  
  for (j in seq_along(IHDP_list)) {
    
    data <- IHDP_list[[j]]
    save(data, file = paste0(path, "/../data/IHDP/dataA_conditional/data_", iter, "_", j, ".RData"))
    
  }
}


# function to simulate setting B see Hill (2011)
simulate_B <- function(sigy = 1, seed, m){
  set.seed(seed)
  
  betaB <- sample(c(0.0, 0.1, 0.2, 0.3, 0.4), dimx + 1, replace = TRUE, 
                  prob = c(0.6, 0.1, 0.1, 0.1, 0.1))
  
  xb <- cbind(rep(1, N), (Xmat + .5)) %*% betaB
  yb0hat <- exp(xb)
  yb1hat <- xb
  
  offset <- mean(yb1hat[Trt == 0] - yb0hat[Trt == 0]) - 4
  yb1hat <- xb - offset
  
  res <- vector("list", m)
  
  for (j in seq_len(m)) {
    set.seed(01042026 + j)
    
    YB0 <- rnorm(N, yb0hat, sigy)
    YB1 <- rnorm(N, yb1hat, sigy)
    
    Y <- YB1
    Y[Trt == 0] <- YB0[Trt == 0]
    
    tauBis <- yb1hat - yb0hat
    tauB <- mean(tauBis)
    tauBs <- mean(YB1[Trt == 0] - YB0[Trt == 0])
    
    B <- cbind(X, Trt, Y, true_cate = tauBis)
    res[[j]] <- B
  }
  
  return(res)
}


# create data for marginal coverage check
for (iter in 1:100) {
  
  IHDP_list <- simulate_B(seed = 30092026 + iter, m = 1)
  
  for (j in seq_along(IHDP_list)) {
    
    data <- IHDP_list[[j]]
    save(data, file = paste0(path, "/../data/IHDP/dataB/data_", iter, "_", j, ".RData"))
    
  }
}


# create data for conditional coverage check
for (iter in 1:5) {
  
  IHDP_list <- simulate_B(seed = 30092026 + iter, m = 200)
  
  for (j in seq_along(IHDP_list)) {
    
    data <- IHDP_list[[j]]
    save(data, file = paste0(path, "/../data/IHDP/dataB_conditional/data_", iter, "_", j, ".RData"))
  
  }
}
