# Implementation to be able to run it on a cluster
iter <- as.integer(Sys.getenv("PBS_ARRAYID"))

# Libraries ----

# path for library (specify on cluster)
.libPaths(c("~/R/library", .libPaths()))

library(dplyr)
library(tibble)  
library(FNN)
library(purrr)

library(bartCause)  # BART
library(cfcausal)   # Lei & Candes Method
library(grf)        # Causal Forest
library(ranger)     # Random Forest
library(xgboost)    # XGBoost
library(gbm)        # Boosting for propensity score
library(quantreg)   # Quantile Regression
library(mice)       # MICE Imputation


# Path / Settings ----

# change path based on data you want to use
path <- file.path(getwd(), "Data", "IHDP_Missing", "dataA_MAR")

# specify when submitting to cluster that we need 4 cores
n_cores <- 4 

# List for collecting results of THIS iteration
all_results <- list()
k_idx <- 1
nsim <- 1

# Helper: to address missing data (for machine learning algorithms) ----

# XGBoost
to_xgb_matrix <- function(X) {
  X <- data.matrix(as.data.frame(X))
  storage.mode(X) <- "double"
  X <- X + 0.0
  X
}

# BART
prepare_missing_indicators <- function(train_df, test_df) {
  cols_with_na <- names(train_df)[colSums(is.na(train_df)) > 0]
  
  for (col in cols_with_na) {
    # add indicators
    train_df[[paste0(col, "_isNA")]] <- as.numeric(is.na(train_df[[col]]))
    test_df[[paste0(col, "_isNA")]]  <- as.numeric(is.na(test_df[[col]]))
    
    # replace NA by -999
    train_df[[col]][is.na(train_df[[col]])] <- -999
    test_df[[col]][is.na(test_df[[col]])]   <- -999
  }
  return(list(train = train_df, test = test_df))
}
# Helper: Zaffran ----
## Helper: Basics (Zaffran) ----
pattern_to_id <- function(m) {strtoi(paste(as.integer(m), collapse = ""), base = 2)}
 
patterns_to_ids <- function(M) {apply(M, 1, pattern_to_id)}

## Helper: Base Model (Zaffran) ----
fit_basemodel <- function(X_train, Y_train, target = "Quantiles",
                          basemodel = c("rf", "xgb"), num.trees = 512,
                          alpha = 0.1, params = list()) {

  basemodel <- match.arg(basemodel)
  X_train <- as.matrix(X_train)

  if (basemodel == "rf") { # Random Forest
    model <- ranger::ranger(x = as.data.frame(X_train), y = Y_train,
                            num.trees = num.trees, quantreg = TRUE)

  } else if (basemodel == "xgb") { # XGBoost

    X_train <- to_xgb_matrix(X_train)
    dtrain  <- xgboost::xgb.DMatrix(data = X_train, label = Y_train)

    q_low <- xgboost::xgb.train(params = list(objective = "reg:quantileerror",
                                              quantile_alpha = alpha / 2,
                                              nthread = n_cores),
                                data = dtrain, nrounds = 300)

    q_high <- xgboost::xgb.train(params = list(objective = "reg:quantileerror",
                                               quantile_alpha = 1 - alpha / 2,
                                               nthread = n_cores),
                                 data = dtrain, nrounds = 300)

    model <- list(q_low = q_low, q_high = q_high)
  }

  list(basemodel = basemodel, alpha = alpha, model = model)
}

predict_basemodel <- function(object, X_new) {
  X_new <- as.matrix(X_new)

  if (object$basemodel == "rf") { # Random Forest

    # predict
    preds <- predict(object$model, data = as.data.frame(X_new), type = "quantiles",
                     quantiles = c(object$alpha / 2, 1 - object$alpha / 2))$predictions
    list(y_inf = preds[, 1], y_sup = preds[, 2])

  } else if (object$basemodel == "xgb") { # XGBoost

    # predict
    X_new <- to_xgb_matrix(X_new)
    dnew  <- xgboost::xgb.DMatrix(data = X_new)

    list(y_inf = as.numeric(predict(object$model$q_low, dnew)),
         y_sup = as.numeric(predict(object$model$q_high, dnew)))
  }
}


## Helper: Conformal Part (Zaffran) ----
quantile_corrected <- function(scores, alpha) {
  n <- length(scores)
  q <- (1 - alpha) * (1 + 1 / n)

  if (q > 1) return(Inf)

  as.numeric(quantile(scores, probs = q, names = FALSE))
}

# function to calibrate the prediction interval
calibrate_pi <- function(fitted_basemodel, imputer, X_cal, M_cal, Y_cal, X_mis_test,
                         features_test, M_test, groups_test, alpha = 0.1) {
  
  # redefine data
  X_cal <- as.matrix(X_cal)
  M_cal <- as.matrix(M_cal)

  X_mis_test <- as.matrix(X_mis_test)
  M_test <- as.matrix(M_test)

  features_test <- as.matrix(features_test)

  patterns <- unique(M_test, MARGIN = 1)
  ids <- apply(patterns, 1, pattern_to_id)

  n_test <- nrow(X_mis_test)
  q_scores_test <- numeric(n_test)

  # calculate smallest number of points to use the exact quantile
  n_min <- ceiling(1 / alpha) - 1

  # for each pattern
  for (idp in seq_len(nrow(patterns))) {

    id_pattern <- ids[idp]
    pattern <- patterns[idp, ]
    idx_test_pat <- which(groups_test == id_pattern)

    ind_exact <- apply(M_cal[, pattern == 0, drop = FALSE] == 0, 1, all)
    n_exact <- sum(ind_exact)

    # indicator if exact version can be used
    use_exact <- (n_exact >= n_min)

    if (use_exact) { ind_subsample <- ind_exact
    } else { ind_subsample <- rep(TRUE, nrow(M_cal)) }

    empty <- (sum(ind_subsample) == 0)

    # compute masking
    X_imp_cal_masking <- X_cal[ind_subsample, , drop = FALSE]
    M_cal_masking <- M_cal[ind_subsample, , drop = FALSE]
    Y_cal_masking <- Y_cal[ind_subsample]

    pattern_ext <- if (ncol(X_imp_cal_masking) > length(pattern)) {
      c(pattern, rep(0, ncol(X_imp_cal_masking) - length(pattern)))
    } else pattern
    
    
    X_imp_cal_masking[, pattern_ext == 1] <- NA
    M_cal_masking[, pattern == 1] <- 1

    if (!empty) {
      X_imp_cal_masking <- transform_imputer(imputer, X_imp_cal_masking)
      features_cal <- cbind(X_imp_cal_masking, M_cal_masking)

      cal_predictions <- predict_basemodel(fitted_basemodel, features_cal)
      scores <- pmax(cal_predictions$y_inf - Y_cal_masking, Y_cal_masking - cal_predictions$y_sup)

      q_scores_test[idx_test_pat] <- quantile_corrected(scores, alpha)
    } else {
      q_scores_test[idx_test_pat] <- Inf
    }
  }
  pred_test <- predict_basemodel(fitted_basemodel, features_test)

  list(y_inf = pred_test$y_inf - q_scores_test, y_sup = pred_test$y_sup + q_scores_test)
}

## Helper: Imputation (Zaffran) ----
fit_imputer <- function(X) {

  X <- as.matrix(X)
  d <- ncol(X)

  # mean imputation
  col_means <- colMeans(X, na.rm = TRUE)
  col_means[is.nan(col_means)] <- 0

  obj <- list(col_means = col_means)
  obj$fill_value <- col_means

  obj
}

transform_imputer <- function(object, Xnew) {

  Xnew <- as.matrix(Xnew)
  d <- ncol(Xnew)

  Ximp <- Xnew
  for (j in seq_len(d)) {
    na_idx <- is.na(Ximp[, j])
    if (any(na_idx)) Ximp[na_idx, j] <- object$fill_value[j]
  }

  return(Ximp)
}

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
  
  # Random Forest
  if (method == "rf"){ # Random Forest
    
    # nuisance parameter mu0 and mu1
    model_mu0 <- ranger(fml_y,  data = d0[, c(yvar, xvars)], num.threads = n_cores,
                        num.trees = num.trees, respect.unordered.factors = "order")
    model_mu1 <- ranger(fml_y,  data = d1[, c(yvar, xvars)], num.threads = n_cores,
                        num.trees = num.trees, respect.unordered.factors = "order")
    
  } else if (method == "xgb") { # XGBoost
    
    x0_mat <- to_xgb_matrix(d0[, xvars, drop = FALSE])
    x1_mat <- to_xgb_matrix(d1[, xvars, drop = FALSE])
    
    dtrain0 <- xgboost::xgb.DMatrix(data = x0_mat, label = d0[[yvar]])
    dtrain1 <- xgboost::xgb.DMatrix(data = x1_mat, label = d1[[yvar]])
    
    # nuisance parameter mu0 and mu1
    model_mu0 <- xgboost::xgb.train(
      params = list(objective = "reg:squarederror", nthread = n_cores),
      data = dtrain0, nrounds = 300, verbose = 0)
    
    model_mu1 <- xgboost::xgb.train(
      params = list(objective = "reg:squarederror", nthread = n_cores),
      data = dtrain1, nrounds = 300, verbose = 0)
  }  
  
  # data prepartion
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
    X_mat <- to_xgb_matrix(X[, xvars, drop = FALSE])
    dnew  <- xgboost::xgb.DMatrix(data = X_mat)
    
    # predict nuisance parameters
    mu0_pred <- predict(nuisance_models$model_mu0, dnew)
    mu1_pred <- predict(nuisance_models$model_mu1, dnew)

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

# Helper: Helper: multiplication factor ----
# compute multiplication factor (to get unbiased results)
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
  if (method == "rf") { # random forest
    fit_cate <- ranger::ranger(pseudo ~ ., data = train_cqr_df, 
                               num.trees = num.trees, num.threads = n_cores,
                               respect.unordered.factors = "order")
    cate_hat <- as.numeric(predict(fit_cate, data = test_x)$predictions)
    
    
  } else if (method == "xgb") { # XGBoost
    
    X_train <- to_xgb_matrix(train_cqr_df[, xvars, drop = FALSE])
    y_train <- train_cqr_df$pseudo
    X_test <- to_xgb_matrix(test_x[, xvars, drop = FALSE])
    
    dtrain_cate <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
    dtest_cate  <- xgboost::xgb.DMatrix(data = X_test)
    
    fit_cate <- xgboost::xgb.train(
      params = list(objective = "reg:squarederror", nthread = n_cores),
      data = dtrain_cate, nrounds = 300, verbose = 0)
    
    cate_hat <- as.numeric(predict(fit_cate, dtest_cate))

  }  
  
  list(cate_point = data.frame(cate_hat = cate_hat),
       intervals_95 = data.frame(lower = cate_lower_95, upper = cate_upper_95),
       intervals_80 = data.frame(lower = cate_lower_80, upper = cate_upper_80),
       true_cate = test_df$true_cate, cf = cf, pseudo_test = pseudo_test,
       cqr_raw = list(scores_cal_95 = scores_cal_95, scores_cal_80 = scores_cal_80,
                      weights_calib = w_calib, weights_test = w_test,
                      pred_test_95 = data.frame(q_lower = as.numeric(pred_test_95[, 1]),
                                                q_upper = as.numeric(pred_test_95[, 2])),
                      pred_test_80 = data.frame(q_lower = as.numeric(pred_test_80[, 1]),
                                                q_upper = as.numeric(pred_test_80[, 2]))))
}

# Main: zaffran ----
zaffran <- function(method = c("rf", "xgb"), imputation = "mean",
                    nuisance_df, quantile_df, train_df, calib_df,
                    test_df, yvar, tvar, xvars, K = 3,
                    num.trees = 512, min_prop = 0.05, alpha){
  
  method <- match.arg(method)
  
  # cross fit DR-outcomes
  cf <- crossfit_dr_scores(train_df = train_df, xvars = xvars, yvar = yvar,
                           tvar = tvar, K = K, num.trees = num.trees,
                           min_prop = min_prop, method = method)
  
  nuisance_full <- cf$nuisance_full
  
  # pseudo-outcomes
  pseudo_nuis <- cf$pseudo_oof[1 : nrow(nuisance_df)]
  pseudo_quant <- cf$pseudo_oof[(nrow(nuisance_df) + 1):(nrow(nuisance_df) + nrow(quantile_df))]
  
  # heuristic multiplication factor
  mult_train <- propensity_multiplier(df = train_df, nuisance_models = nuisance_full,
                                      xvars = xvars, tvar = tvar, min_prop = min_prop)
  mult_nuis <- mult_train[1:nrow(nuisance_df)]
  mult_quant <- mult_train[(nrow(nuisance_df) + 1):(nrow(nuisance_df) + nrow(quantile_df))]
  
  pseudo_nuis <- pseudo_nuis * mult_nuis
  pseudo_quant <- pseudo_quant * mult_quant
  
  
  # predict dr pseudo-outcomes 
  pseudo_calib <- predict_dr_outcome(calib_df, nuisance_full, xvars = xvars,
                                     yvar = yvar, tvar = tvar, min_prop = min_prop,
                                     method = method)
  mult_calib <- propensity_multiplier(df = calib_df, nuisance_models = nuisance_full,
                                      xvars = xvars, tvar = tvar, min_prop = min_prop)
  pseudo_calib <- pseudo_calib * mult_calib
  
  
  pseudo_test <- predict_dr_outcome(test_df, nuisance_full, xvars = xvars,
                                    yvar = yvar, tvar = tvar, min_prop = min_prop,
                                    method = method)
  mult_test <- propensity_multiplier(df = test_df, nuisance_models = nuisance_full,
                                     xvars = xvars, tvar = tvar, min_prop = min_prop)
  pseudo_test <- pseudo_test * mult_test

  
  # imputer
  imputer <- fit_imputer(train_df[ , covariateNames])

  # impute data 
  X_imp_train <- transform_imputer(imputer, train_df[ , covariateNames])
  X_imp_cal   <- transform_imputer(imputer, calib_df[ , covariateNames])
  X_imp_test  <- transform_imputer(imputer, test_df[ , covariateNames])

  
  # mask
  mask <- as.matrix(as.data.frame(lapply(data[covariateNames], function(x) as.integer(is.na(x)))))
  
  features_train <- cbind(as.matrix(X_imp_train), mask[train_idx, , drop = FALSE])
  features_test  <- cbind(as.matrix(X_imp_test),  mask[test_idx, , drop = FALSE])
  
  # fit model
  fitted_model <- fit_basemodel(features_train, c(pseudo_nuis, pseudo_quant),
                                target = "Quantiles", basemodel = method,
                                alpha = alpha)
  
  groups_test <- patterns_to_ids(mask[test_idx, covariateNames])
  
  # get prediction intervals
  res_exact <- calibrate_pi(fitted_basemodel = fitted_model, imputer = imputer,
                            X_cal = calib_df[ , covariateNames], 
                            M_cal = mask[calib_idx, covariateNames],
                            Y_cal = pseudo_calib, 
                            X_mis_test = data[test_idx, covariateNames],
                            features_test = features_test, groups_test = groups_test,
                            M_test = mask[test_idx, covariateNames], alpha = alpha)
  
  return(list(res_exact = res_exact, pseudo_test = pseudo_test))
  
}
# Imputation ----
run_one <- function(data_imp, iter, proportion, imp_id, nuis_idx, quant_idx, 
                    calib_idx, test_idx, covariateNames, n_cores = 4) {
  
  res_list <- list()
  raw_list <- list()
  k <- 1
  
  nuisance_df <- data_imp[nuis_idx, ]
  quantile_df <- data_imp[quant_idx, ]
  train_df    <- data_imp[c(nuis_idx, quant_idx), ]
  calib_df    <- data_imp[calib_idx, ]
  test_df     <- data_imp[test_idx, ]
  
  # cqr_dr
  cqr_methods <- c("rf" = "cqr_dr_rf", "xgb" = "cqr_dr_xgb")
  
  for (m_code in names(cqr_methods)) {
    
    set.seed(42 + 40 * iter + 1000 * imp_id)
    
    out <- weighted_cqr_dr_cate(method = m_code, nuisance_df = nuisance_df,
                                quantile_df = quantile_df, train_df = train_df,
                                calib_df = calib_df, test_df = test_df,
                                xvars = covariateNames, yvar = "Y", tvar = "Trt",
                                K = 3, min_prop = 0.01)
    
    res_list[[k]] <- tibble(sim = iter, proportion = proportion, imp = imp_id,
                            conf_level = 0.95, method = cqr_methods[[m_code]],
                            person_id = test_idx, cate = out$cate_point$cate_hat,
                            lower = out$intervals_95$lower, upper = out$intervals_95$upper,
                            true_cate = test_df$true_cate, dr = out$pseudo_test)
    k <- k + 1
    
    res_list[[k]] <- tibble(sim = iter, proportion = proportion, imp = imp_id,
                            conf_level = 0.80, method = cqr_methods[[m_code]],
                            person_id = test_idx, cate = out$cate_point$cate_hat,
                            lower = out$intervals_80$lower, upper = out$intervals_80$upper,
                            true_cate = test_df$true_cate, dr = out$pseudo_test)
    k <- k + 1
    
    raw_list[[cqr_methods[[m_code]]]] <- out
  }

  list(results = bind_rows(res_list), raw = raw_list)
}

# Pooling ----

# to pool a matrix
pool_matrix <- function(mat, type = c("mean", "median")) {
  type <- match.arg(type)
  
  if (type == "mean") {
    rowMeans(mat, na.rm = TRUE)
  } else {
    apply(mat, 1, median, na.rm = TRUE)
  }
}

get_cqr_names <- function(alpha) {
  if (abs(alpha - 0.05) < 1e-12) {
    list(score_name = "scores_cal_95", pred_name  = "pred_test_95")
  } else if (abs(alpha - 0.20) < 1e-12) {
    list(score_name = "scores_cal_80", pred_name  = "pred_test_80")
  } 
}

## Pooling mean ----
pool_cqr_scores_mean <- function(imp_runs, method_name, alpha, sim, proportion,
                                 test_idx, true_cate) {
  
  # names
  nm <- get_cqr_names(alpha)
  
  # run imputation
  outs <- lapply(imp_runs, function(z) z$raw[[method_name]])
  
  # create matrices
  scores_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$score_name]]))
  
  weights_calib_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw$weights_calib))
  weights_test_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw$weights_test))
  
  qlo_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$pred_name]]$q_lower))
  qhi_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$pred_name]]$q_upper))
  
  cate_mat <- do.call(cbind, lapply(outs, function(o) o$cate_point$cate_hat))
  
  # pooling
  scores_pool <- rowMeans(scores_mat, na.rm = TRUE)
  weights_calib_pool <- rowMeans(weights_calib_mat, na.rm = TRUE)
  weights_test_pool  <- rowMeans(weights_test_mat, na.rm = TRUE)
  
  qlo_pool <- rowMeans(qlo_mat, na.rm = TRUE)
  qhi_pool <- rowMeans(qhi_mat, na.rm = TRUE)
  cate_pool <- rowMeans(cate_mat, na.rm = TRUE)
  
  qhat_pool <- vapply(weights_test_pool, function(wt) {
      weighted_conformal_quantile(scores_calib = scores_pool, alpha = alpha,
                                  w_calib = weights_calib_pool, w_test = wt)},
    numeric(1))
  
  tibble(sim = sim, proportion = proportion, conf_level = 1 - alpha,
         method = paste0(method_name, "_mean"), person_id = test_idx,
         cate = cate_pool, lower = qlo_pool - qhat_pool,
         upper = qhi_pool + qhat_pool, true_cate = true_cate)
}

# Pooling median ----
pool_cqr_scores_median <- function(imp_runs, method_name, alpha, sim,
                                   proportion, test_idx, true_cate) {
  
  # names
  nm <- get_cqr_names(alpha)
  
  # run imputation
  outs <- lapply(imp_runs, function(z) z$raw[[method_name]])
  
  # create matrices
  scores_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$score_name]]))
  
  weights_calib_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw$weights_calib))
  weights_test_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw$weights_test))
  
  qlo_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$pred_name]]$q_lower))
  qhi_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$pred_name]]$q_upper))
  
  cate_mat <- do.call(cbind, lapply(outs, function(o) o$cate_point$cate_hat))
  
  # pooling
  scores_pool <- apply(scores_mat, 1, median, na.rm = TRUE)
  weights_calib_pool <- apply(weights_calib_mat, 1, median, na.rm = TRUE)
  weights_test_pool  <- apply(weights_test_mat, 1, median, na.rm = TRUE)
  
  qlo_pool <- apply(qlo_mat, 1, median, na.rm = TRUE)
  qhi_pool <- apply(qhi_mat, 1, median, na.rm = TRUE)
  cate_pool <- apply(cate_mat, 1, median, na.rm = TRUE)
  
  qhat_pool <- vapply(weights_test_pool, function(wt) {
    weighted_conformal_quantile(scores_calib = scores_pool, alpha = alpha,
                                w_calib = weights_calib_pool, w_test = wt)},
    numeric(1))
  
  tibble(sim = sim, proportion = proportion, conf_level = 1 - alpha,
         method = paste0(method_name, "_median"), person_id = test_idx,
         cate = cate_pool, lower = qlo_pool - qhat_pool, 
         upper = qhi_pool + qhat_pool, true_cate = true_cate)
}
# Pooling stacking + conservative ----
stack_cqr_scores_all <- function(imp_runs, method_name, alpha, sim, proportion,
                                 test_idx, true_cate, 
                                 base_pool = c("stacked", "conservative")) {
  
  base_pool <- match.arg(base_pool)
  nm <- get_cqr_names(alpha)
  
  # run imputation
  outs <- lapply(imp_runs, function(z) z$raw[[method_name]])
  
  # stacking
  scores_stack <- unlist(lapply(outs, function(o) o$cqr_raw[[nm$score_name]]))
  
  weights_calib_stack <- unlist(lapply(outs, function(o) o$cqr_raw$weights_calib))
  weights_test_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw$weights_test))
  
  qlo_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$pred_name]]$q_lower))
  qhi_mat <- do.call(cbind, lapply(outs, function(o) o$cqr_raw[[nm$pred_name]]$q_upper))
  
  cate_mat <- do.call(cbind, lapply(outs, function(o) o$cate_point$cate_hat))
  
  # stacking method
  if (base_pool == "stacked") {
    qlo_pool <- rowMeans(qlo_mat, na.rm = TRUE)
    qhi_pool <- rowMeans(qhi_mat, na.rm = TRUE)
  } else if (base_pool == "conservative") {
    qlo_pool <- apply(qlo_mat, 1, min, na.rm = TRUE)
    qhi_pool <- apply(qhi_mat, 1, max, na.rm = TRUE)
  }
  
  weights_test_pool <- rowMeans(weights_test_mat, na.rm = TRUE)
  cate_pool_vec <- rowMeans(cate_mat, na.rm = TRUE)
  
  qhat_stack <- vapply(weights_test_pool, function(wt) {
    weighted_conformal_quantile(scores_calib = scores_stack, w_test = wt,
                                w_calib = weights_calib_stack, alpha = alpha)},
    numeric(1))
  
  tibble(sim = sim, proportion = proportion, conf_level = 1 - alpha,
    method = paste0(method_name, "_", base_pool), person_id = test_idx,
    cate = cate_pool_vec, lower = qlo_pool - qhat_stack, 
    upper = qhi_pool + qhat_stack, true_cate = true_cate)
}
# iteration ----
M_imp <- 3
imp_runs <- vector("list", M_imp)

for (proportion in c("10", "25", "50")) {
  
  set.seed(30092026 + iter + as.numeric(proportion))
  
  # load data
  load(paste0(path, "/IHDP_", iter, "_", proportion, ".RData"))
    
  # define covariate names
  covariateNames <- setdiff(names(data), c("Trt", "Y", "true_cate"))
  n <- nrow(data)
    
  idx <- sample.int(n)
    
  # subset data
  n_nuis <- floor(0.25 * n)                # size nuisance data frame
  n_quant <- floor(0.5 * n)                # size quantile data frame
  n_calib <- floor(0.15 * n)               # size calibration data frame
  n_test <- n - n_nuis - n_quant - n_calib # size test data frame
    
  nuis_idx  <- idx[1:n_nuis]                                            # nuisance index
  quant_idx <- idx[(n_nuis + 1):(n_nuis + n_quant)]                     # quantile index
  calib_idx <- idx[(n_nuis + n_quant + 1):(n_nuis + n_quant + n_calib)] # calibration index
  test_idx <- idx[(n_nuis + n_quant + n_calib + 1):n]                   # test index
  train_idx <- c(nuis_idx, quant_idx)                                   # train index
    
  # MICE
  ignore_vec <- rep(FALSE, n)
  ignore_vec[c(calib_idx, test_idx)] <- TRUE
    
  # Multiple Imputation
  imp <- mice(data, m = M_imp, ignore = ignore_vec, method = "rf", printFlag = FALSE, 
              seed = 4200 + 10 * iter + as.integer(proportion))
  completed_list <- mice::complete(imp, action = "all")
    
    
  # cqr_dr ----
  # compute pipeline for all imputed datasets
  for (m_imp in seq_len(M_imp)) {

    set.seed(42 + 40 * iter + 2 * as.numeric(proportion))

    imp_runs[[m_imp]] <- run_one(data_imp = completed_list[[m_imp]], iter = iter,
                                 proportion = proportion, imp_id = m_imp,
                                 nuis_idx = nuis_idx, quant_idx = quant_idx,
                                 calib_idx = calib_idx, test_idx = test_idx,
                                 covariateNames = covariateNames, n_cores = n_cores)
  }

  # pooling / stacking
  results_cqr_score_methods <- bind_rows(

    # cqr_dr_rf 95%
    pool_cqr_scores_mean(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                         alpha = 0.05, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    pool_cqr_scores_median(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                           alpha = 0.05, sim = iter, proportion = proportion,
                           test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                         alpha = 0.05, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "stacked"),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                         alpha = 0.05, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "conservative"),

    # cqr_dr_rf 80%
    pool_cqr_scores_mean(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                         alpha = 0.2, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    pool_cqr_scores_median(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                           alpha = 0.2, sim = iter, proportion = proportion,
                           test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                         alpha = 0.2, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "stacked"),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_rf",
                         alpha = 0.2, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "conservative"),

     # cqr_dr_xgb 95%
    pool_cqr_scores_mean(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                         alpha = 0.05, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    pool_cqr_scores_median(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                           alpha = 0.05, sim = iter, proportion = proportion,
                           test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                         alpha = 0.05, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "stacked"),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                         alpha = 0.05, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "conservative"),

    # cqr_dr_xgb 80%
    pool_cqr_scores_mean(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                         alpha = 0.2, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    pool_cqr_scores_median(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                           alpha = 0.2, sim = iter, proportion = proportion,
                           test_idx = test_idx, true_cate = data$true_cate[test_idx]),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                         alpha = 0.2, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "stacked"),

    stack_cqr_scores_all(imp_runs = imp_runs, method_name = "cqr_dr_xgb",
                         alpha = 0.2, sim = iter, proportion = proportion,
                         test_idx = test_idx, true_cate = data$true_cate[test_idx],
                         base_pool = "conservative")
  )
  all_results[[k_idx]] <- results_cqr_score_methods
  k_idx <- k_idx + 1
    
    
  # CAUSAL FOREST ----
  set.seed(42 + 40 * iter + 2 * as.numeric(proportion))
  cf_auto <- causal_forest(data[train_idx, covariateNames], data$Y[train_idx],
                           data$Trt[train_idx], num.threads = n_cores)
    
  cate_cf <- predict(cf_auto, data[test_idx, covariateNames], estimate.variance = TRUE)
  tau_hat <- cate_cf$predictions
  se_hat  <- sqrt(cate_cf$variance.estimates)
    
  # 95% Interval
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.95, 
                                 method = "causal forest", person_id = test_idx, 
                                 cate = tau_hat, lower = tau_hat - 1.96 * se_hat,
                                 upper = tau_hat + 1.96 * se_hat, 
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
    
  # 80% Interval
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.80, 
                                 method = "causal forest", person_id = test_idx, 
                                 cate = tau_hat, lower = tau_hat - 1.282 * se_hat,
                                 upper = tau_hat + 1.282 * se_hat, 
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
    
    
  # BART ----
  set.seed(42 + 40 * iter + 2 * as.numeric(proportion))
  
  prepared <- prepare_missing_indicators(data[train_idx, covariateNames], 
                                         data[test_idx, covariateNames])
  fit <- bartc(confounders = prepared$train, response = data$Y[train_idx], 
               treatment = data$Trt[train_idx], keepTrees = TRUE)
  ite_draws <- predict(fit, newdata = prepared$test, type = "ite")
  bart_cate <- apply(ite_draws, 2, mean)
    
  # 95%
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.95, 
                                 method = "BART", person_id = test_idx, cate = bart_cate, 
                                 lower = apply(ite_draws, 2, quantile, probs = 0.025), 
                                 upper = apply(ite_draws, 2, quantile, probs = 0.975), 
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
    
  # 80%
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.80, 
                                 method = "BART", person_id = test_idx, cate = bart_cate,
                                 lower = apply(ite_draws, 2, quantile, probs = 0.10),
                                 upper = apply(ite_draws, 2, quantile, probs = 0.90),
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
    
  # ZAFFRAN
  set.seed(42 + 40 * iter + 2 * as.numeric(proportion))
  zf <- zaffran(method = "rf", train_df = data[train_idx, ], xvars = covariateNames,
                tvar = "Trt", yvar = "Y", alpha = 0.05,
                nuisance_df = data[nuis_idx, ], quantile_df = data[quant_idx, ],
                calib_df = data[calib_idx, ], test_df = data[test_idx, ])
    
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.95,
                                 method = "zaffran_rf", person_id = test_idx, cate = NA,
                                 lower = zf$res_exact$y_inf,
                                 upper = zf$res_exact$y_sup,
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
    
  zf <- zaffran(method = "rf", train_df = data[train_idx, ], xvars = covariateNames,
                tvar = "Trt", yvar = "Y", alpha = 0.2,
                nuisance_df = data[nuis_idx, ], quantile_df = data[quant_idx, ],
                calib_df = data[calib_idx, ], test_df = data[test_idx, ])
    
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.80,
                                 method = "zaffran_rf", person_id = test_idx, cate = NA,
                                 lower = zf$res_exact$y_inf,
                                 upper = zf$res_exact$y_sup,
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
    
  zf <- zaffran(method = "xgb", train_df = data[train_idx, ], xvars = covariateNames,
                tvar = "Trt", yvar = "Y", alpha = 0.05,
                nuisance_df = data[nuis_idx, ], quantile_df = data[quant_idx, ],
                calib_df = data[calib_idx, ], test_df = data[test_idx, ])
    
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.95,
                                 method = "zaffran_xgb", person_id = test_idx, cate = NA,
                                 lower = zf$res_exact$y_inf,
                                 upper = zf$res_exact$y_sup,
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
    
  zf <- zaffran(method = "xgb", train_df = data[train_idx, ], xvars = covariateNames,
                tvar = "Trt", yvar = "Y", alpha = 0.2,
                nuisance_df = data[nuis_idx, ], quantile_df = data[quant_idx, ],
                calib_df = data[calib_idx, ], test_df = data[test_idx, ])
    
  all_results[[k_idx]] <- tibble(sim = iter, proportion = proportion, conf_level = 0.80,
                                 method = "zaffran_xgb", person_id = test_idx, cate = NA,
                                 lower = zf$res_exact$y_inf,
                                 upper = zf$res_exact$y_sup,
                                 true_cate = data[test_idx, "true_cate"])
  k_idx <- k_idx + 1
}

results <- bind_rows(all_results)
NameOut <- paste0(path, "/results_MAR_iter_", iter, ".RData")
save(results, file = NameOut)
