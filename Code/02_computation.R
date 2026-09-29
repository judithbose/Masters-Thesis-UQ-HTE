# Implementation to be able to run it on a cluster
iter <- as.integer(Sys.getenv("PBS_ARRAYID"))

# Libraries ----

# path for library (specify on cluster)
.libPaths(c("~/R/library", .libPaths()))

library(dplyr)
library(tibble) 
library(bartCause)  # BART
library(cfcausal)   # Lei & Candes Method
library(grf)        # Causal Forest
library(ranger)     # Random Forest
library(xgboost)    # XGBoost
library(gbm)        # Boosting for propensity score
library(quantreg)   # Quantile Regression

# Path / Settings ----

# change path based on data you want to use
path <- file.path(getwd(), "Data", "IHDP", "dataA")

# specify when submitting to cluster that we need 4 cores
n_cores <- 4 

# List for collecting results of THIS iteration
all_results <- list()
k_idx <- 1

# Helper: weighted quantile ----
# compute weighted quantiles according to approach of Tibshirani et al. 
weighted_conformal_quantile <- function(scores_calib, w_calib, w_test, alpha) {

  n <- length(scores_calib)

  # sort scores
  ord <- order(scores_calib)
  scores <- scores_calib[ord]
  w <- w_calib[ord]

  # cumulative weights
  cum_w <- cumsum(w)

  # total weight including test point
  total_w <- sum(w) + w_test

  # weighted CDF
  cdf <- cum_w / total_w

  # find first index reaching 1-alpha
  idx <- which(cdf >= (1 - alpha))[1]

  if (is.na(idx)) {
    return(max(scores))
  } else {
    return(scores[idx])
  }
}

# Helper: density-ratio estimation ----
# fit density ratio
# (true ratio is unknown so ratio needs to be estimated, as described in 
# Tibshrani et al. Random Forests will be used)
fit_density_ratio <- function(source_x, target_x, min_prob = 1e-6, 
                              num.trees = 512) {
  
  # create data format
  dom <- factor(c(rep("source", nrow(source_x)), rep("target", nrow(target_x))),
                levels = c("source", "target"))
  dat <- data.frame(domain = dom, rbind(source_x, target_x))
  
  # use random forest to estimate ratio between source and target population
  fit <- ranger(domain ~ ., data = dat, probability = TRUE,
                num.trees = num.trees, respect.unordered.factors = "order",
                num.threads = n_cores)
  
  list(model = fit, n_source = nrow(source_x), n_target = nrow(target_x),
       min_prob = min_prob)
}

# predict density ratio (function is split up for easier use in the following)
predict_density_ratio <- function(dr_fit, new_x) {
  
  # predict based on input model
  pred <- predict(dr_fit$model, data = new_x)$predictions
  p_target <- as.numeric(pred[, "target"])
  
  # protection against 0/1 and numerical outliers
  p_target <- pmin(pmax(p_target, dr_fit$min_prob), 1 - dr_fit$min_prob)
  w <- (p_target / (1 - p_target)) * (dr_fit$n_source / dr_fit$n_target)
  
  # optional: trimming of huge weights
  w <- pmin(w, 10)
  w <- w / mean(w)
  
  return(w)
}

# Helper: nuisance model ----
# fit nuisance models (mu0, mu1 using random forest, xgboost and propensity score)
fit_nuisance_models <- function(df, xvars, yvar, tvar, num.trees = 512,
                                method) {
  
  # data preparation
  df[[tvar]] <- factor(df[[tvar]], levels = c(0, 1))
  d0 <- df[df[[tvar]] == 0, , drop = FALSE]
  d1 <- df[df[[tvar]] == 1, , drop = FALSE]
  
  # model formulation
  fml_pi  <- as.formula(paste(tvar, "~", paste(xvars, collapse = " + ")))
  fml_y  <- as.formula(paste(yvar, "~", paste(xvars, collapse = " + ")))
  
  # random forest
  if (method == "rf"){
    
    # nuisance parameter mu0 and mu1
    model_mu0 <- ranger(fml_y,  data = d0[, c(yvar, xvars)], num.threads = n_cores,
                        num.trees = num.trees, respect.unordered.factors = "order")
    model_mu1 <- ranger(fml_y,  data = d1[, c(yvar, xvars)], num.threads = n_cores,
                        num.trees = num.trees, respect.unordered.factors = "order")
    
  } else if (method == "xgb") { # XGBoost
    
    x0_mat <- model.matrix(fml_y, data = d0)[, -1]
    x1_mat <- model.matrix(fml_y, data = d1)[, -1]
    
    # nuisance parameter mu0 and mu1
    model_mu0 <- xgboost(x = x0_mat, y = d0[[yvar]], objective = "reg:squarederror", nrounds = 300, nthread = n_cores, verbose = 0)
    model_mu1 <- xgboost(x = x1_mat, y = d1[[yvar]], objective = "reg:squarederror", nrounds = 300, nthread = n_cores, verbose = 0)
    
    }  
  
  # data preparation
  df_gbm <- df
  df_gbm[[tvar]] <- as.numeric(as.character(df_gbm[[tvar]]))
  
  # nuisance parameter propensity score
  model_pi <- gbm(formula = fml_pi, data = df_gbm[, c(tvar, xvars)],
                  distribution = "bernoulli", n.trees = 1000, n.cores = 1,
                  interaction.depth = 3, shrinkage = 0.01, cv.folds = 5, verbose = FALSE)
  
  # return nuisance parameters
  list(model_pi = model_pi, model_mu0 = model_mu0, model_mu1 = model_mu1)
}

# Helper: DR pseudo-outcome ----
# compute DR pseudo-outcome
predict_dr_outcome <- function(df, nuisance_models, xvars, yvar, tvar, 
                               min_prop = 0.01, method) {
  
  # data set
  df[[tvar]] <- as.numeric(as.character(df[[tvar]]))
  T_val <- df[[tvar]]
  Y <- df[[yvar]]
  X <- df[, xvars, drop = FALSE]
  
  if (method == "rf") { # Random Forest
    
    # predict nuisance parameters
    mu0_pred <- as.numeric(predict(nuisance_models$model_mu0, data = X)$predictions)
    mu1_pred <- as.numeric(predict(nuisance_models$model_mu1, data = X)$predictions)
    
  } else if (method == "xgb") { # XGBoost
    
    X_mat <- model.matrix(~ . - 1, data = X)
    
    # predict nuisance parameters
    mu0_pred <- predict(nuisance_models$model_mu0, X_mat)
    mu1_pred <- predict(nuisance_models$model_mu1, X_mat)
  }  
  
  # predict propensity score
  best_pi  <- gbm.perf(nuisance_models$model_pi, method = "cv", plot.it = FALSE)
  pi_pred  <- predict(nuisance_models$model_pi, newdata = X, n.trees = best_pi, type = "response")
  
  # avoid values that are too large or too small
  pi_pred <- pmin(pmax(pi_pred, min_prop), 1 - min_prop)
  
  # compute pseudo-outcome as explained in Chapter 2.2
  pseudo <- (mu1_pred - mu0_pred) + (T_val * (Y - mu1_pred) / pi_pred) -
    ((1 - T_val) * (Y - mu0_pred) / (1 - pi_pred))
  
  # clip pseudo value as well
  clip_val <- quantile(abs(pseudo), 0.99)
  pseudo <- pmax(pmin(pseudo, clip_val), -clip_val)
  
  as.numeric(pseudo)
}

# Helper: heuristic variance deflation ----
# compute multiplication factor (deflate the variance)
propensity_multiplier <- function(df, nuisance_models, xvars, tvar, min_prop = 0.01) {
  
  # data set
  df[[tvar]] <- as.numeric(as.character(df[[tvar]]))
  T_val <- df[[tvar]]
  X <- df[, xvars, drop = FALSE]
  
  # use predicted propensity scores
  best_pi <- gbm.perf(nuisance_models$model_pi, method = "cv", plot.it = FALSE)
  pi_hat <- predict(nuisance_models$model_pi, newdata = X, n.trees = best_pi,
                    type = "response")
  pi_hat <- pmin(pmax(pi_hat, min_prop), 1 - min_prop)
  

  # group T = 1: pi_hat
  mean_treated <- mean(pi_hat[T_val == 1], na.rm = TRUE)
  # group T = 0: 1 - pi_hat
  mean_control <- mean(1 - pi_hat[T_val == 0], na.rm = TRUE)
  
  # create vector
  multiplier <- ifelse(T_val == 1, mean_treated, mean_control)
  
  return(as.numeric(multiplier))
}

# Helper: cross-fitting ----
crossfit_dr_scores <- function(train_df, xvars, yvar, tvar, K = 3, 
                               num.trees = 512, min_prop = 0.01, method) {
  
  n <- nrow(train_df)
  pseudo_oof <- rep(NA_real_, n)
  
  # generate random folds 
  folds <- sample(rep(seq_len(K), length.out = n))
  
  # k folds cross-validation
  for (k in seq_len(K)) {
    idx_tr <- which(folds != k)
    idx_va <- which(folds == k)
    
    # fit nuisance models
    nuisance_k <- fit_nuisance_models(train_df[idx_tr, ], xvars = xvars, 
                                      yvar = yvar, tvar = tvar, 
                                      num.trees = num.trees, method = method)
    
    # calculate pseudo outcomes
    pseudo_oof[idx_va] <- predict_dr_outcome(train_df[idx_va, ],  nuisance_k, 
                                             xvars = xvars, yvar = yvar, 
                                             tvar = tvar, min_prop = min_prop, 
                                             method = method)
  }
  
  # full nuisance model 
  nuisance_full <- fit_nuisance_models(train_df, xvars = xvars, yvar = yvar, 
                                       tvar = tvar, num.trees = num.trees,
                                       method = method)
  
  list(pseudo_oof = pseudo_oof, nuisance_full = nuisance_full)
}
# Main: CQR + DR ----
weighted_cqr_dr_cate <- function(method = c("rf", "xgb"), nuisance_df, quantile_df, 
                                 train_df, calib_df, test_df, yvar, tvar, xvars, 
                                 K = 3, num.trees = 512, min_prop = 0.05) {
  
  method <- match.arg(method)
  
  # cross fit scores
  cf <- crossfit_dr_scores(train_df = train_df, xvars = xvars, yvar = yvar, 
                           tvar = tvar, K = K, num.trees = num.trees, 
                           min_prop = min_prop, method = method)
  nuisance_full <- cf$nuisance_full
  
  # pseudo outcomes (nuisance and quantile)
  pseudo_nuis <- cf$pseudo_oof[1 : nrow(nuisance_df)]
  pseudo_quant <- cf$pseudo_oof[(nrow(nuisance_df) + 1):(nrow(nuisance_df) + nrow(quantile_df))]
  
  mult_train <- propensity_multiplier(df = train_df, nuisance_models = nuisance_full,
                                      xvars = xvars, tvar = tvar, min_prop = min_prop)
  mult_nuis <- mult_train[1:nrow(nuisance_df)]
  mult_quant <- mult_train[(nrow(nuisance_df) + 1):(nrow(nuisance_df) + nrow(quantile_df))]
  
  pseudo_nuis <- pseudo_nuis * mult_nuis
  pseudo_quant <- pseudo_quant * mult_quant
  
  
  # compute conformalized quantile regression
  train_cqr_df <- data.frame(pseudo = c(pseudo_nuis, pseudo_quant), train_df[, xvars, drop = FALSE])
  fml_cqr <- as.formula(paste("pseudo ~", paste(xvars, collapse = " + ")))
  fit_qrf <- ranger::ranger(formula = fml_cqr, data = train_cqr_df, min.node.size = 20,
                            quantreg = TRUE, num.trees = 1000, num.threads = n_cores)

  
  # pseudo outcomes (calibration)  
  pseudo_calib <- predict_dr_outcome(calib_df, nuisance_full, xvars = xvars,
                                     yvar = yvar, tvar = tvar, min_prop = min_prop,
                                     method = method)
  mult_calib <- propensity_multiplier(df = calib_df, nuisance_models = nuisance_full,
                                      xvars = xvars, tvar = tvar, min_prop = min_prop)
  pseudo_calib <- pseudo_calib * mult_calib
  calib_x <- calib_df[, xvars, drop = FALSE]
  
  # pseudo outcomes (test)
  pseudo_test <- predict_dr_outcome(test_df, nuisance_full, xvars = xvars,
                                    yvar = yvar, tvar = tvar, min_prop = min_prop,
                                    method = method)
  mult_test <- propensity_multiplier(df = test_df, nuisance_models = nuisance_full,
                                     xvars = xvars, tvar = tvar, min_prop = min_prop)
  pseudo_test <- pseudo_test * mult_test
  test_x  <- test_df[, xvars, drop = FALSE]

  
  # estimate weights (by density ratio)
  dr_fit <- fit_density_ratio(source_x = train_df[, xvars], target_x = test_x,
                              num.trees = num.trees)

  w_calib <- predict_density_ratio(dr_fit, calib_x)
  w_test  <- predict_density_ratio(dr_fit, test_x)

  
  # compute scores and quantile (95%)
  pred_cal_95 <- predict(fit_qrf, data = calib_x, type = "quantiles", 
                         quantiles = c(0.025, 0.975))$predictions
  scores_cal_95 <- pmax(as.numeric(pred_cal_95[, 1]) - pseudo_calib, 
                        pseudo_calib - as.numeric(pred_cal_95[, 2]))

  qhat_95 <- sapply(w_test, function(wt) weighted_conformal_quantile(scores_calib = scores_cal_95, 
                                                                     w_calib = w_calib, w_test = wt, 
                                                                     alpha = 0.05))
  pred_test_95 <- predict(fit_qrf, data = test_x, type = "quantiles", 
                          quantiles = c(0.025, 0.975))$predictions
  
  # 95% interval
  cate_lower_95 <- as.numeric(pred_test_95[, 1]) - qhat_95
  cate_upper_95 <- as.numeric(pred_test_95[, 2]) + qhat_95

  
  # compute scores and quantile (80%)
  pred_cal_80 <- predict(fit_qrf, data = calib_x, type = "quantiles", 
                          quantiles = c(0.10, 0.90))$predictions
  scores_cal_80 <- pmax(as.numeric(pred_cal_80[, 1]) - pseudo_calib, 
                        pseudo_calib - as.numeric(pred_cal_80[, 2]))

  qhat_80 <- sapply(w_test, function(wt) weighted_conformal_quantile(scores_calib = scores_cal_80, 
                                                                     w_calib = w_calib, w_test = wt, 
                                                                     alpha = 0.2))
  pred_test_80 <- predict(fit_qrf, data = test_x, type = "quantiles", 
                           quantiles = c(0.10, 0.90))$predictions
   
  # 80% interval
  cate_lower_80 <- as.numeric(pred_test_80[, 1]) - qhat_80
  cate_upper_80 <- as.numeric(pred_test_80[, 2]) + qhat_80
  
   
  # CATE point estimation
  if (method == "rf") { # Random Forest
    fit_cate <- ranger::ranger(pseudo ~ ., data = train_cqr_df, 
                               num.trees = num.trees, num.threads = n_cores,
                               respect.unordered.factors = "order")
    cate_hat <- as.numeric(predict(fit_cate, data = test_x)$predictions)
    
  } else if (method == "xgb") { # XGBoost
    X_train <- model.matrix(~ . - 1, data = train_df[, xvars, drop = FALSE])
    X_test  <- model.matrix(~ . - 1, data = test_x)
    
    dtrain_cate <- xgb.DMatrix(data = X_train, label = cf$pseudo_oof)
    fit_cate <- xgb.train(params = list(objective = "reg:squarederror"),
                          data = dtrain_cate, nrounds = 300,
                          nthread = n_cores)
    cate_hat <- as.numeric(predict(fit_cate, X_test))
  }  

  list(cate_point = data.frame(cate_hat = cate_hat),
       intervals_95 = data.frame(lower = cate_lower_95, upper = cate_upper_95),
       intervals_80 = data.frame(lower = cate_lower_80, upper = cate_upper_80),
       true_cate = data[test_idx, "true_cate"], cf = cf, pseudo_test = pseudo_test)
}

# Main: CP + DR  ----
weighted_split_dr_cate <- function(method = c("rf", "xgb"), train_df, calib_df, 
                                   test_df, yvar, tvar, xvars, K = 3, num.trees = 512,
                                   min_prop = 0.05, alpha_95 = 0.05, alpha_80 = 0.20,
                                   n_cores = 1) {
  
  method <- match.arg(method)
  
  # DR pseudo-outcomes on training data
  cf <- crossfit_dr_scores(train_df = train_df, xvars = xvars, yvar = yvar,
                           tvar = tvar, K = K, num.trees = num.trees, 
                           min_prop = min_prop, method = method)
  nuisance_full <- cf$nuisance_full
  
  
  # pseudo-outcomes for calibration and test
  pseudo_calib <- predict_dr_outcome(calib_df, nuisance_full, xvars = xvars, 
                                     yvar = yvar, tvar = tvar, min_prop = min_prop,
                                     method = method)
  pseudo_test <- predict_dr_outcome(test_df, nuisance_full, xvars = xvars, 
                                    yvar = yvar, tvar = tvar, min_prop = min_prop, 
                                    method = method)
  calib_x <- calib_df[, xvars, drop = FALSE]
  test_x  <- test_df[, xvars, drop = FALSE]
  
  
  # density ratio weights
  dr_fit <- fit_density_ratio(source_x = train_df[, xvars, drop = FALSE],
                              target_x = test_x, num.trees = num.trees)
  
  w_calib <- predict_density_ratio(dr_fit, calib_x)
  w_test  <- predict_density_ratio(dr_fit, test_x)
  
  
  # CATE point estimation
  if (method == "rf") {
    
    train_cate_df <- data.frame(pseudo = cf$pseudo_oof, train_df[, xvars, drop = FALSE])
    fit_cate <- ranger::ranger(pseudo ~ ., data = train_cate_df, num.threads = n_cores,
                               num.trees = num.trees, respect.unordered.factors = "order")
    
    pred_cal <- as.numeric(predict(fit_cate, data = calib_x)$predictions)
    cate_hat <- as.numeric(predict(fit_cate, data = test_x)$predictions)
    
  } else if (method == "xgb") {
    X_train <- model.matrix(~ . - 1, data = train_df[, xvars, drop = FALSE])
    X_cal   <- model.matrix(~ . - 1, data = calib_x)
    X_test  <- model.matrix(~ . - 1, data = test_x)
    
    dtrain_cate <- xgb.DMatrix(data = X_train, label = cf$pseudo_oof)
    fit_cate <- xgb.train(params = list(objective = "reg:squarederror"),
                          data = dtrain_cate, nrounds = 300, nthread = n_cores,
                          verbose = 0)
    
    pred_cal <- as.numeric(predict(fit_cate, X_cal))
    cate_hat  <- as.numeric(predict(fit_cate, X_test))
  }
  
  # conformal scores: absolute residuals
  scores_cal <- abs(pseudo_calib - pred_cal)
  
  # 95% interval
  qhat_95 <- sapply(w_test, function(wt) weighted_conformal_quantile(
    scores_calib = scores_cal, w_calib = w_calib, w_test = wt, alpha = alpha_95))
  intervals_95 <- data.frame(lower = cate_hat - qhat_95, upper = cate_hat + qhat_95)
  
  # 80% interval
  qhat_80 <- sapply(w_test, function(wt) weighted_conformal_quantile(
    scores_calib = scores_cal, w_calib = w_calib, w_test = wt, alpha = alpha_80))
  intervals_80 <- data.frame(lower = cate_hat - qhat_80, upper = cate_hat + qhat_80)
  
  list(cate_point = data.frame(cate_hat = cate_hat), intervals_95 = intervals_95,
       intervals_80 = intervals_80, pseudo_test = pseudo_test, weights_calib = w_calib, 
       weights_test = w_test, cf = cf)
}

# iteration (conditional coverage) ----
for (people in 1:5){

  set.seed(30092026 + iter)
  
  # load data
  load(paste0(path, "/IHDP_", people, "_", iter,".RData"))
  
  load(paste0(path, "/IHDP_", iter, "_1.RData"))
  #load(paste0(path, "/data_", iter, ".RData"))
  
  data <- IHDP

  # define covariate names
  covariateNames <- setdiff(names(data), c("Trt", "Y", "true_cate"))
  n <- nrow(data)

  # subset data
  n_nuis <- floor(0.25 * n)                # size nuisance data frame
  n_quant <- floor(0.5 * n)                # size quantile data frame
  n_calib <- floor(0.15 * n)               # size calibration data frame
  n_test <- n - n_nuis - n_quant - n_calib # size test data frame

  # shift test and target distribution
  shift_var <- covariateNames[1]
  mean_val <- mean(data[[shift_var]])
  prob_test <- plogis(data[[shift_var]] - mean_val) 

  # index 
  test_idx <- sample(1:n, size = n_test, replace = FALSE, prob = prob_test)
  remaining_idx <- setdiff(1:n, test_idx) 
  idx_rest <- sample(remaining_idx)

  nuis_idx  <- idx_rest[1:n_nuis]                                # nuisance index
  quant_idx <- idx_rest[(n_nuis + 1):(n_nuis + n_quant)]         # quantile index
  calib_idx <- idx_rest[(n_nuis + n_quant + 1):length(idx_rest)] # calibration index
  train_idx <- c(nuis_idx, quant_idx)                            # train index

  nuisance_df <- data[nuis_idx, ]            # nuisance data frame
  quantile_df <- data[quant_idx, ]           # quantile data frame
  train_df <- data[c(nuis_idx, quant_idx), ] # train data frame
  calib_df    <- data[calib_idx, ]           # calibration data frame
  test_df     <- data[test_idx, ]            # test data frame


  # CQR + DR Methods ----
  cqr_methods <- c("rf" = "cqr_dr_rf", "xgb" = "cqr_dr_xgb")
  saved_cf <- NULL 
  xgboost_cf <- NULL

  for (m_code in names(cqr_methods)) {
    
    set.seed(42 + 40 * iter + 2 * 1) 
  
    out <- weighted_cqr_dr_cate(method = m_code, nuisance_df = nuisance_df, 
                                quantile_df = quantile_df, train_df = train_df,
                                calib_df = calib_df, test_df = test_df,
                                xvars = covariateNames, yvar = "Y", 
                                tvar = "Trt", K = 3, min_prop = 0.01)
  
    if (m_code == "rf") saved_cf <- out$cf
    if (m_code == "xgb") xgboost_cf <- out$cf
  
    # 95% interval
    all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95,
                                   method = cqr_methods[[m_code]], person_id = test_idx, 
                                   cate = out$cate_point$cate_hat, 
                                   lower = out$intervals_95$lower,
                                   upper = out$intervals_95$upper,
                                   true_cate = out$true_cate, dr = out$pseudo_test)
    k_idx <- k_idx + 1
  
    # 80% interval
    all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80,
                                   method = cqr_methods[[m_code]], person_id = test_idx, 
                                   cate = out$cate_point$cate_hat, 
                                   lower = out$intervals_80$lower,
                                   upper = out$intervals_80$upper,
                                   true_cate = out$true_cate, dr = out$pseudo_test)
    k_idx <- k_idx + 1
  }

  
  # Split CP + DR ----
  methods_split <- c("rf" = "dr_rf", "xgb" = "dr_xgb")
  cf_objects <- list(rf = saved_cf, xgb = xgboost_cf)

  for (m in names(methods_split)) {
  
    set.seed(42 + 40 * iter + 2 * 1)
  
    out <- weighted_split_dr_cate(method = m, train_df = train_df, 
                                  calib_df = calib_df, test_df = test_df, 
                                  xvars = covariateNames, yvar = "Y", tvar = "Trt",
                                  n_cores = n_cores, K = 10)
  
    # 95%
    all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, 
                                   method = methods_split[[m]],
                                   person_id = test_idx, cate = out$cate_point$cate_hat,
                                   lower = out$intervals_95$lower, 
                                   upper = out$intervals_95$upper,
                                   true_cate = test_df$true_cate, dr = out$pseudo_test)
    k_idx <- k_idx + 1
  
    # 80%
    all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80, 
                                   method = methods_split[[m]],
                                   person_id = test_idx, cate = out$cate_point$cate_hat,
                                   lower = out$intervals_80$lower, 
                                   upper = out$intervals_80$upper,
                                   true_cate = test_df$true_cate, dr = out$pseudo_test)
    k_idx <- k_idx + 1
  }

  
  # Lei & Candes ----
  set.seed(42 + 40 * iter + 2 * 1)
  
  CIfun_std <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                            data$Trt[train_idx], alpha = 0.05, algo = "nest", 
                            exact = TRUE, type = "mean", useCV = FALSE) 
  ci_std <- CIfun_std(data[test_idx, covariateNames]) 

  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, method = "L&C", 
                                 person_id = test_idx, cate = NA, lower = ci_std$lower, 
                                 upper = ci_std$upper, true_cate = data[test_idx, "true_cate"], 
                                 dr = NA) 
  k_idx <- k_idx + 1 

  
  set.seed(42 + 40 * iter + 2 * 1)
  
  CIfun_std <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                            data$Trt[train_idx], alpha = 0.2, algo = "nest",
                            exact = TRUE, type = "mean", useCV = FALSE) 
  ci_std <- CIfun_std(data[test_idx, covariateNames]) 
  
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.8, method = "L&C", 
                                 person_id = test_idx, cate = NA, lower = ci_std$lower, 
                                 upper = ci_std$upper, true_cate = data[test_idx, "true_cate"],
                                 dr = NA) 
  k_idx <- k_idx + 1 


  # Lei & Candes + CQR ----
  set.seed(42 + 40 * iter + 2 * 1)
  
  CIfun_na <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                           data$Trt[train_idx], alpha = 0.05, algo = "nest", 
                           exact = TRUE, type = "CQR", quantiles = c(0.025, 0.975), 
                           useCV = FALSE) 
  ci_na <- CIfun_na(data[test_idx, covariateNames])
  
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, method = "L&C cqr",
                                 person_id = test_idx, cate = NA, lower = ci_na$lower,
                                 upper = ci_na$upper, true_cate = data[test_idx, "true_cate"],
                                 dr = NA) 
  k_idx <- k_idx + 1 


  set.seed(42 + 40 * iter + 2 * 1)
  CIfun_na <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                           data$Trt[train_idx], alpha = 0.2, algo = "nest", 
                           exact = TRUE, type = "CQR", quantiles = c(0.1, 0.9),
                           useCV = FALSE) 
  ci_na <- CIfun_na(data[test_idx, covariateNames])
  
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.8, method = "L&C cqr", 
                                 person_id = test_idx, cate = NA, lower = ci_na$lower,
                                 upper = ci_na$upper, true_cate = data[test_idx, "true_cate"], 
                                 dr = NA) 
  k_idx <- k_idx + 1 

  
  # Causal Forest ----
  
  set.seed(42 + 40 * iter + 2 * 1)
  cf_auto <- causal_forest(data[train_idx, covariateNames], data$Y[train_idx],
                           data$Trt[train_idx], num.threads = n_cores)

  cate_cf <- predict(cf_auto, data[test_idx, covariateNames], estimate.variance = TRUE)
  tau_hat <- cate_cf$predictions
  se_hat  <- sqrt(cate_cf$variance.estimates)

  # 95% Interval
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, 
                                 method = "causal forest", person_id = test_idx, 
                                 cate = tau_hat, lower = tau_hat - 1.96 * se_hat,
                                 upper = tau_hat + 1.96 * se_hat, 
                                 true_cate = data[test_idx, "true_cate"], dr = NA)
  k_idx <- k_idx + 1

  # 80% Interval
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80, 
                                 method = "causal forest", person_id = test_idx, 
                                 cate = tau_hat, lower = tau_hat - 1.282 * se_hat,
                                 upper = tau_hat + 1.282 * se_hat, 
                                 true_cate = data[test_idx, "true_cate"], dr = NA)
  k_idx <- k_idx + 1
  

  # BART ----
  
  set.seed(42 + 40 * iter + 2 * 1)
  fit <- bartc(confounders = data[train_idx, covariateNames], 
               response = data$Y[train_idx], treatment = data$Trt[train_idx],
               keepTrees = TRUE, n.threads = n_cores)

  ite_draws <- predict(fit, newdata = data[test_idx, covariateNames], type = "ite")
  bart_cate <- apply(ite_draws, 2, mean)

  # 95%
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, 
                                 method = "BART", person_id = test_idx, cate = bart_cate, 
                                 lower = apply(ite_draws, 2, quantile, probs = 0.025), 
                                 upper = apply(ite_draws, 2, quantile, probs = 0.975), 
                                 true_cate = data[test_idx, "true_cate"], dr = NA)
  k_idx <- k_idx + 1

  # 80%
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80, 
                                 method = "BART", person_id = test_idx, cate = bart_cate,
                                 lower = apply(ite_draws, 2, quantile, probs = 0.10),
                                 upper = apply(ite_draws, 2, quantile, probs = 0.90),
                                 true_cate = data[test_idx, "true_cate"], dr = NA)
  k_idx <- k_idx + 1


  # save results ----
  results <- bind_rows(all_results)
  NameOut <- paste0(path, "/results_iter_", people, "_", iter, ".RData")
  save(results, file = NameOut)
}

# iteration (marginal coverage) ----
set.seed(30092026 + iter)
  
# load data
load(paste0(path, "/IHDP_", iter, "_1.RData"))
#load(paste0(path, "/data_", iter, ".RData"))
  
data <- IHDP
  
# define covariate names
covariateNames <- setdiff(names(data), c("Trt", "Y", "true_cate"))
n <- nrow(data)
  
# subset data
n_nuis <- floor(0.25 * n)                # size nuisance data frame
n_quant <- floor(0.5 * n)                # size quantile data frame
n_calib <- floor(0.15 * n)               # size calibration data frame
n_test <- n - n_nuis - n_quant - n_calib # size test data frame
  
# shift test and target distribution
shift_var <- covariateNames[1]
mean_val <- mean(data[[shift_var]])
prob_test <- plogis(data[[shift_var]] - mean_val) 
  
# index 
test_idx <- sample(1:n, size = n_test, replace = FALSE, prob = prob_test)
remaining_idx <- setdiff(1:n, test_idx) 
idx_rest <- sample(remaining_idx)
  
nuis_idx  <- idx_rest[1:n_nuis]                                # nuisance index
quant_idx <- idx_rest[(n_nuis + 1):(n_nuis + n_quant)]         # quantile index
calib_idx <- idx_rest[(n_nuis + n_quant + 1):length(idx_rest)] # calibration index
train_idx <- c(nuis_idx, quant_idx)                            # train index
  
nuisance_df <- data[nuis_idx, ]            # nuisance data frame
quantile_df <- data[quant_idx, ]           # quantile data frame
train_df <- data[c(nuis_idx, quant_idx), ] # train data frame
calib_df    <- data[calib_idx, ]           # calibration data frame
test_df     <- data[test_idx, ]            # test data frame
  
  
# CQR + DR Methods ----
cqr_methods <- c("rf" = "cqr_dr_rf", "xgb" = "cqr_dr_xgb")
saved_cf <- NULL 
xgboost_cf <- NULL
  
for (m_code in names(cqr_methods)) {
    
  set.seed(42 + 40 * iter + 2 * 1) 
    
  out <- weighted_cqr_dr_cate(method = m_code, nuisance_df = nuisance_df, 
                              quantile_df = quantile_df, train_df = train_df,
                              calib_df = calib_df, test_df = test_df,
                              xvars = covariateNames, yvar = "Y", 
                              tvar = "Trt", K = 3, min_prop = 0.01)
    
  if (m_code == "rf") saved_cf <- out$cf
  if (m_code == "xgb") xgboost_cf <- out$cf
    
  # 95% interval
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95,
                                 method = cqr_methods[[m_code]], person_id = test_idx, 
                                 cate = out$cate_point$cate_hat, 
                                 lower = out$intervals_95$lower,
                                 upper = out$intervals_95$upper,
                                 true_cate = out$true_cate, dr = out$pseudo_test)
  k_idx <- k_idx + 1
    
  # 80% interval
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80,
                                 method = cqr_methods[[m_code]], person_id = test_idx, 
                                 cate = out$cate_point$cate_hat, 
                                 lower = out$intervals_80$lower,
                                 upper = out$intervals_80$upper,
                                 true_cate = out$true_cate, dr = out$pseudo_test)
  k_idx <- k_idx + 1
}
  
  
# Split CP + DR ----
methods_split <- c("rf" = "dr_rf", "xgb" = "dr_xgb")
cf_objects <- list(rf = saved_cf, xgb = xgboost_cf)
  
for (m in names(methods_split)) {
    
  set.seed(42 + 40 * iter + 2 * 1)
  
  out <- weighted_split_dr_cate(method = m, train_df = train_df, 
                                calib_df = calib_df, test_df = test_df, 
                                xvars = covariateNames, yvar = "Y", tvar = "Trt",
                                n_cores = n_cores, K = 10)
    
  # 95%
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, 
                                 method = methods_split[[m]],
                                 person_id = test_idx, cate = out$cate_point$cate_hat,
                                 lower = out$intervals_95$lower, 
                                 upper = out$intervals_95$upper,
                                 true_cate = test_df$true_cate, dr = out$pseudo_test)
  k_idx <- k_idx + 1
    
  # 80%
  all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80, 
                                 method = methods_split[[m]],
                                 person_id = test_idx, cate = out$cate_point$cate_hat,
                                 lower = out$intervals_80$lower, 
                                 upper = out$intervals_80$upper,
                                 true_cate = test_df$true_cate, dr = out$pseudo_test)
  k_idx <- k_idx + 1
}
  
  
# Lei & Candes ----
set.seed(42 + 40 * iter + 2 * 1)
  
CIfun_std <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                          data$Trt[train_idx], alpha = 0.05, algo = "nest", 
                          exact = TRUE, type = "mean", useCV = FALSE) 
ci_std <- CIfun_std(data[test_idx, covariateNames]) 
  
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, method = "L&C", 
                               person_id = test_idx, cate = NA, lower = ci_std$lower, 
                               upper = ci_std$upper, true_cate = data[test_idx, "true_cate"], 
                               dr = NA) 
k_idx <- k_idx + 1 
  
  
set.seed(42 + 40 * iter + 2 * 1)
  
CIfun_std <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                          data$Trt[train_idx], alpha = 0.2, algo = "nest",
                          exact = TRUE, type = "mean", useCV = FALSE) 
ci_std <- CIfun_std(data[test_idx, covariateNames]) 
  
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.8, method = "L&C", 
                               person_id = test_idx, cate = NA, lower = ci_std$lower, 
                               upper = ci_std$upper, true_cate = data[test_idx, "true_cate"],
                               dr = NA) 
k_idx <- k_idx + 1 
  
  
# Lei & Candes + CQR ----
set.seed(42 + 40 * iter + 2 * 1)
  
CIfun_na <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                         data$Trt[train_idx], alpha = 0.05, algo = "nest", 
                         exact = TRUE, type = "CQR", quantiles = c(0.025, 0.975), 
                         useCV = FALSE) 
ci_na <- CIfun_na(data[test_idx, covariateNames])
  
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, method = "L&C cqr",
                               person_id = test_idx, cate = NA, lower = ci_na$lower,
                               upper = ci_na$upper, true_cate = data[test_idx, "true_cate"],
                               dr = NA) 
k_idx <- k_idx + 1 
  
  
set.seed(42 + 40 * iter + 2 * 1)
CIfun_na <- conformalIte(data[train_idx, covariateNames], data$Y[train_idx], 
                         data$Trt[train_idx], alpha = 0.2, algo = "nest", 
                         exact = TRUE, type = "CQR", quantiles = c(0.1, 0.9),
                         useCV = FALSE) 
ci_na <- CIfun_na(data[test_idx, covariateNames])
  
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.8, method = "L&C cqr", 
                               person_id = test_idx, cate = NA, lower = ci_na$lower,
                               upper = ci_na$upper, true_cate = data[test_idx, "true_cate"], 
                               dr = NA) 
k_idx <- k_idx + 1 
  
  
# Causal Forest ----
  
set.seed(42 + 40 * iter + 2 * 1)
cf_auto <- causal_forest(data[train_idx, covariateNames], data$Y[train_idx],
                           data$Trt[train_idx], num.threads = n_cores)
  
cate_cf <- predict(cf_auto, data[test_idx, covariateNames], estimate.variance = TRUE)
tau_hat <- cate_cf$predictions
se_hat  <- sqrt(cate_cf$variance.estimates)
  
# 95% Interval
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, 
                               method = "causal forest", person_id = test_idx, 
                               cate = tau_hat, lower = tau_hat - 1.96 * se_hat,
                               upper = tau_hat + 1.96 * se_hat, 
                               true_cate = data[test_idx, "true_cate"], dr = NA)
k_idx <- k_idx + 1
  
# 80% Interval
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80, 
                               method = "causal forest", person_id = test_idx, 
                               cate = tau_hat, lower = tau_hat - 1.282 * se_hat,
                               upper = tau_hat + 1.282 * se_hat, 
                               true_cate = data[test_idx, "true_cate"], dr = NA)
k_idx <- k_idx + 1
  
  
# BART ----
  
set.seed(42 + 40 * iter + 2 * 1)
fit <- bartc(confounders = data[train_idx, covariateNames], 
             response = data$Y[train_idx], treatment = data$Trt[train_idx],
             keepTrees = TRUE, n.threads = n_cores)
  
ite_draws <- predict(fit, newdata = data[test_idx, covariateNames], type = "ite")
bart_cate <- apply(ite_draws, 2, mean)
  
# 95%
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.95, 
                               method = "BART", person_id = test_idx, cate = bart_cate, 
                               lower = apply(ite_draws, 2, quantile, probs = 0.025), 
                               upper = apply(ite_draws, 2, quantile, probs = 0.975), 
                               true_cate = data[test_idx, "true_cate"], dr = NA)
k_idx <- k_idx + 1
  
# 80%
all_results[[k_idx]] <- tibble(sim = iter, rep = 1, conf_level = 0.80, 
                               method = "BART", person_id = test_idx, cate = bart_cate,
                               lower = apply(ite_draws, 2, quantile, probs = 0.10),
                               upper = apply(ite_draws, 2, quantile, probs = 0.90),
                               true_cate = data[test_idx, "true_cate"], dr = NA)
k_idx <- k_idx + 1
  
  
# save results ----
results <- bind_rows(all_results)
NameOut <- paste0(path, "/results_iter_", people, "_", iter, ".RData")
save(results, file = NameOut)
