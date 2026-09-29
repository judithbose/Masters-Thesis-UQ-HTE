# This script is used to generate the tables in Chapters 4.2 (incomplete data), 
# C (incomplete data) and E.

# libraries ----
library(dplyr)
library(readr)
library(purrr)
library(tidyverse)

# paths ----
path <- getwd()
path <- paste0(path, "/Results/Incomplete")

files <- list.files(path = path, pattern = "^data.*\\.RData$", full.names = TRUE)

methods <- c("BART", "Causal Forest", "zaffran_rf", "zaffran_xgb",
  "cqr_dr_rf_mean", "cqr_dr_rf_median", "cqr_dr_rf_stacked", "cqr_dr_rf_conservative",
  "cqr_dr_xgb_mean", "cqr_dr_xgb_median", "cqr_dr_xgb_stacked", "cqr_dr_xgb_conservative")
# function to create tables ----
evaluate_results <- function(file, conf_lvl) {
  
  load(file)
  
  results <- results %>% filter(conf_level == conf_lvl)
  
  # 1. table for CATE estimation
  estimation <- results %>%
    group_by(method, sim, proportion) %>%
    summarise(pehe = sqrt(mean((cate - true_cate)^2, na.rm = TRUE)),
              bias  = mean(cate - true_cate, na.rm = TRUE), .groups = "drop")
  
  # create label
  unique_props <- sort(unique(estimation$proportion))
  prop_levels  <- paste0(unique_props, "%")
  
  estimation <- estimation %>%
    group_by(method, proportion) %>%
    summarise(pehe = mean(pehe, na.rm = TRUE), 
              bias = mean(bias, na.rm = TRUE),.groups = "drop") %>%
    mutate(method = factor(method, levels = methods),
           proportion = factor(paste0(proportion, "%"), levels = prop_levels)) %>%
    arrange(method) %>%
    mutate(across(c(pehe, bias), ~ ifelse(is.nan(.x), "-", sprintf("%.3f", .x))))
  
  estimation_table <- estimation %>%
    select(method, proportion, pehe, bias) %>%
    pivot_wider(names_from  = proportion, values_from = c(pehe, bias),
                names_glue  = "{.value}_{proportion}") %>%
    mutate(method = factor(method, levels = methods)) %>%
    arrange(method)
  
  # 2. table for interval estimation
  interval <- results %>%
    group_by(method, proportion, sim) %>%
    summarise(p_hat_k = mean(true_cate >= lower & true_cate <= upper),
              width_k = mean(upper - lower),
              undershoot_k = mean(upper < true_cate),
              overshoot_k = mean(lower > true_cate), .groups = "drop")
  
  interval <- interval %>% group_by(method, proportion) %>%
    summarise(n_simulations = n(), overall_coverage = mean(p_hat_k),
              mcse_coverage = sd(p_hat_k) / sqrt(n_simulations), mean_width = mean(width_k),
              se_width = sd(width_k) / sqrt(n_simulations), mean_undershoot = mean(undershoot_k),
              mean_overshoot = mean(overshoot_k), .groups = "drop")
  
  # create labels
  unique_props <- sort(unique(interval$proportion))
  prop_levels  <- paste0(unique_props, "%")
  
  interval <- interval %>% 
    mutate(prop_label = factor(paste0(proportion, "%"), levels = prop_levels),
           Coverage = paste0(round(overall_coverage * 100, 1)),
           Width = paste0(round(mean_width, 2)))
  
  # create final table
  interval_table <- interval %>%
    select(method, prop_label, Coverage, Width) %>%
    pivot_wider(names_from  = prop_label, values_from = c(Coverage, Width),
                names_glue  = "{.value}_{prop_label}") %>%
    mutate(method = factor(method, levels = methods)) %>%
    arrange(method)
  
  return(list(estimation_table, interval_table))
}

# loop ----
for (file in files) {
  
  # 95% interval
  evaluation_table <- evaluate_results(file, conf_lvl = 0.95)
  
  # saving table for interval estimation
  output_file <- paste0(path, "/Tables/", tools::file_path_sans_ext(basename(file)), "_95.csv")
  write_csv(evaluation_table[[2]], output_file)
  
  # 80% interval
  evaluation_table <- evaluate_results(file, conf_lvl = 0.8)
  
  # saving table for CATE estimation
  output_file <- paste0(path, "/Tables/CATE_", tools::file_path_sans_ext(basename(file)), ".csv")
  write_csv(evaluation_table[[1]], output_file)
  
  # saving table for interval estimation
  output_file <- paste0(path, "/Tables/", tools::file_path_sans_ext(basename(file)), "_80.csv")
  write_csv(evaluation_table[[2]], output_file)
  
}

# table for 4.2 ----
data <- lapply(files, function(file) {
  env <- new.env()
  load(file, envir = env)
  env$results
})

results <- do.call(rbind, data)

results$method <- recode(results$method, BART = "BART", `Causal Forest` = "Causal Forest",
                         `cqr_dr_rf_conservative` = "cqr_dr_rf", `cqr_dr_rf_mean` = "cqr_dr_rf",
                         `cqr_dr_rf_median` = "cqr_dr_rf", `cqr_dr_rf_stacked` = "cqr_dr_rf",
                         `cqr_dr_xgb_conservative` = "cqr_dr_xgb", `cqr_dr_xgb_mean` = "cqr_dr_xgb",
                         `cqr_dr_xgb_median` = "cqr_dr_xgb", `cqr_dr_xgb_stacked` = "cqr_dr_xgb",
                         `zaffran_rf` = "zaffran_rf", `zaffran_xgb` = "zaffran_xgb")

estimation <- results %>%
  group_by(method, sim) %>%
  summarise(pehe = sqrt(mean((cate - true_cate)^2, na.rm = TRUE)),
            bias  = mean(cate - true_cate, na.rm = TRUE), .groups = "drop")

estimation_table <- estimation %>%
  group_by(method) %>%
  summarise(pehe = mean(pehe, na.rm = TRUE), 
            bias = mean(bias, na.rm = TRUE),.groups = "drop") %>%
  mutate(method = factor(method, levels = c("BART", "Causal Forest", "zaffran_rf", 
                                            "zaffran_xgb", "cqr_dr_rf", "cqr_dr_xgb"))) %>%
  arrange(method) %>%
  mutate(across(c(pehe, bias), ~ ifelse(is.nan(.x), "-", sprintf("%.3f", .x))))

write_csv(estimation_table, paste0(path, "/Tables/CATE.csv"))
