# This script is used to generate the plots in Chapter C.

# libraries ----
library(FNN)
library(tidyverse)
library(viridis)
library(ggrastr)
library(ggplot2)
library(hexbin)
library(dplyr)

# paths ----
path1 <- getwd()

simulation_root <- file.path(path1, "Data", "Simulation")
ihdp_root <- file.path(path1, "Data", "IHDP")
results_root <- file.path(path1, "Results", "Complete")
plot_root <- file.path(path1, "Results", "Complete", "Plots")

simulation_dirs <- list.dirs(simulation_root, recursive = FALSE, full.names = TRUE)
ihdp_dirs <- list.dirs(ihdp_root, recursive = FALSE, full.names = TRUE)
ihdp_dirs <- ihdp_dirs[c(1,3)]
dirs <- c(simulation_dirs, ihdp_dirs)

# settings ----
col <- viridis(4)

methods <- c("BART", "Causal Forest", "L&C", "L&C cqr", "dr_rf", "dr_xgb", 
             "cqr_dr_rf", "cqr_dr_xgb")

cols <- c("Own Approach" = col[3], "Homogeneous Intervals" = col[2],
          "Heterogeneous Intervals" = "black")

k_nachbarn <- 10

# function to create plots ----
analyse_simulation_folder <- function(path, pattern) {
  
  folder_name <- basename(path)
  
  # load simulation results
  results_file <- file.path(results_root, paste0(folder_name, ".RData"))
  env_results <- new.env()
  load(results_file, envir = env_results)
  results <- env_results$results
  
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
  res <- results %>%
    inner_join(df_analyse, by = c("sim", "person_id")) %>%
    mutate(Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Heterogeneous Intervals"),
           method = factor(method, levels = methods))
  
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
  ggsave(filename = file.path(plot_root, paste0(folder_name, "_density.png")),
         plot = p, width = 12, height = 9, dpi = 300)
  
}

# create plots -----
dir_patterns <- c(rep("data", length(simulation_dirs)), rep("IHDP", length(ihdp_dirs)))
plots <- purrr::map2(dirs, dir_patterns,analyse_simulation_folder)
