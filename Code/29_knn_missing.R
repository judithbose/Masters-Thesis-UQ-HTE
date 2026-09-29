# This script is used to generate the plots in Chapter D.

# libraries ----
library(FNN)
library(tidyverse)
library(viridis)
library(ggrastr)
library(ggplot2)
library(hexbin)

# paths ----
path1 <- getwd()

simulation_root <- file.path(path1, "Data", "Simulation")
ihdp_root <- file.path(path1, "Data", "IHDP")
results_root <- file.path(path1, "Results", "Incomplete")
plot_root <- file.path(path1, "Results", "Incomplete", "Plots")

simulation_dirs <- list.dirs(simulation_root, recursive = FALSE, full.names = TRUE)
ihdp_dirs <- list.dirs(ihdp_root, recursive = FALSE, full.names = TRUE)
ihdp_dirs <- ihdp_dirs[c(1,3)]
dirs <- c(simulation_dirs, ihdp_dirs)

# settings ----
col <- viridis(4)

methods <- c("BART", "Causal Forest", "cqr_dr_rf", "cqr_dr_xgb", "zaffran_rf", "zaffran_xgb")

cols <- c("Own Approach" = col[3], "Homogeneous Intervals" = col[2],
          "Heterogeneous Intervals" = "black")

k_nachbarn <- 10

# function to create plots ----
analyse_simulation_folder <- function(path, pattern) {
  
  folder_name <- basename(path)
  
  # load simulation results (MCAR)
  results_file <- file.path(results_root, paste0(folder_name, "_MCAR", ".RData"))
  env_results <- new.env()
  load(results_file, envir = env_results)
  results <- env_results$results
  
  results$method <- recode(results$method, BART = "BART", `Causal Forest` = "Causal Forest",
                           `cqr_dr_rf_conservative` = "cqr_dr_rf", `cqr_dr_rf_mean` = "cqr_dr_rf",
                           `cqr_dr_rf_median` = "cqr_dr_rf", `cqr_dr_rf_stacked` = "cqr_dr_rf",
                           `cqr_dr_xgb_conservative` = "cqr_dr_xgb", `cqr_dr_xgb_mean` = "cqr_dr_xgb",
                           `cqr_dr_xgb_median` = "cqr_dr_xgb", `cqr_dr_xgb_stacked` = "cqr_dr_xgb",
                           `zaffran_rf` = "zaffran_rf", `zaffran_xgb` = "zaffran_xgb")
  results$missingness <- "MCAR"
  
  
  # load simulation results (MCAR)
  results_file <- file.path(results_root, paste0(folder_name, "_MAR", ".RData"))
  env_results <- new.env()
  load(results_file, envir = env_results)
  results2 <- env_results$results
  
  results2$method <- recode(results2$method, BART = "BART", `Causal Forest` = "Causal Forest",
                            `cqr_dr_rf_conservative` = "cqr_dr_rf", `cqr_dr_rf_mean` = "cqr_dr_rf",
                            `cqr_dr_rf_median` = "cqr_dr_rf", `cqr_dr_rf_stacked` = "cqr_dr_rf",
                            `cqr_dr_xgb_conservative` = "cqr_dr_xgb", `cqr_dr_xgb_mean` = "cqr_dr_xgb",
                            `cqr_dr_xgb_median` = "cqr_dr_xgb", `cqr_dr_xgb_stacked` = "cqr_dr_xgb",
                            `zaffran_rf` = "zaffran_rf", `zaffran_xgb` = "zaffran_xgb")
  results2$missingness <- "MAR"
  results <- rbind(results, results2)
  

  # load simulated data
  if (pattern == "data"){
    data_files <- list.files(path, pattern = paste0("^data_[0-9]+\\.RData$"), full.names = TRUE)
  } else if (pattern == "IHDP"){
    data_files <- list.files(path, pattern = paste0("^data_[0-9]+_1\\.RData$"), full.names = TRUE)
  }
  
  
  # create knn indicator
  df_analyse <- map_dfr(data_files, function(file) {
    
    # load data
    env_data <- new.env()
    load(file, envir = env_data)
    data <- env_data$data
    
    # extract simulation number
    if (pattern == "data"){
      sim_id <- basename(file) %>% str_remove("data_") %>%
        str_remove("\\.RData$") %>% as.integer()
    } else if (pattern == "IHDP"){
      sim_id <- basename(file) %>% str_remove("data_") %>%
        str_remove("_1\\.RData$") %>% as.integer()
    }

    # rescale X
    covariateNames <- setdiff(names(data), c("Trt", "Y", "true_cate"))
    X_scaled <- scale(data[, covariateNames])
    
    # compute kNN
    knn_ergebnis <- get.knn(X_scaled, k = k_nachbarn)
    tibble(sim = sim_id, person_id = seq_len(nrow(data)),
           knn_density_dist = rowMeans(knn_ergebnis$nn.dist))
  })
  
  # join data
  # res <- results %>%
  #   inner_join(df_analyse, by = c("sim", "person_id")) %>%
  #   mutate(Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
  #                            method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
  #                            TRUE ~ "Heterogeneous Intervals"),
  #          method = factor(method, levels = methods))
  method_order <- c("BART", "zaffran_rf", "cqr_dr_rf", "Causal Forest", "zaffran_xgb", "cqr_dr_xgb")
  
  res <- results %>%
    inner_join(df_analyse, by = c("sim", "person_id")) %>%
    mutate(Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Heterogeneous Intervals"),
           method = factor(method, levels = method_order))
  
  # compute correlation
  cor_by_method <- res %>% group_by(method) %>%
    summarise(correlation = cor(knn_density_dist, upper - lower, use = "complete.obs",
                                method = "pearson"), .groups = "drop")
  
  # create facet labels
  facet_labels <- cor_by_method %>%
    mutate(label = paste0(method, "\nCorrelation \u03c1 = ", sprintf("%.3f", correlation))) %>%
    select(method, label) %>% deframe()
  
  
  p <- ggplot(res, aes(x = knn_density_dist, y = upper - lower, colour = Group)) +
    
    ggrastr::geom_point_rast(alpha = 0.05, size = 0.3, raster.dpi = 200) +
    
    geom_smooth(aes(fill = Group), method = "glm", se = TRUE,
                linewidth = 1, color = "red") +
    
    facet_wrap(~method, scales = "free_y", nrow = 2,
               labeller = labeller(method = facet_labels)) +
    
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    
    labs(title = "",
         x = "Average k-NN Distance", y = "Interval Width", 
         colour = "Method Type", fill = "Method Type") +
    
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(),
          strip.text = element_text(size = 12), legend.position = "bottom",
          legend.title = element_text(), panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey85", linewidth = 0.3))
  
  # save plot
  ggsave(filename = file.path(plot_root, paste0(folder_name, "_density_missing.png")),
         plot = p, width = 12, height = 9, dpi = 300)
  
  return(p)
}

# create plots -----
dir_patterns <- c(rep("data", length(simulation_dirs)), rep("IHDP", length(ihdp_dirs)))
plots <- purrr::map2(dirs, dir_patterns,analyse_simulation_folder)
