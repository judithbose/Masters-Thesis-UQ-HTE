## This script is used to generate artificial missing values for the simulated data.
## This data is used in Chapter 4.3 of the thesis.

# path
path <- getwd()
path <- paste0(path, "/Code")

# libraries
library(missMethods)
library(missForest)
library(dplyr)

# settings
nsim <- 100
sim_dirs <- list.dirs(paste0(path, "/../Data/Simulation"), recursive = FALSE, 
                      full.names = FALSE)

# MCAR ----
for (sim_dir in sim_dirs) {
  
  out_dir <- paste0(sim_dir, "_MCAR")
  dir.create(file.path(path, "../Data/Simulation_Missing", out_dir),
             recursive = TRUE, showWarnings = FALSE)
  
  for (iter in 1:nsim) {
    for (p in c(0.1, 0.25, 0.5)) {
      
      set.seed(01042026 + iter)
      
      # using delete_MCAR to create MCAR data
      load(file.path(path, "../Data/Simulation", sim_dir, paste0("data_", iter, ".RData")))
      
      data <- delete_MCAR(data, p, 
                          cols_mis = setdiff(names(data), c("Trt", "Y", "true_cate")))
      
      save(data, file = file.path(path, "../Data/Simulation_Missing", out_dir,
                                  paste0("data_", iter, "_", p * 100, ".RData")))
    }
  }
}

# MAR adopted from Thurow et al. (2021) ----

create_mar_vse <- function(missing, rate_ges, datframe, mar_column = 10) {

  mcar_mis <- datframe

  mcar_mis[, missing] <- prodNA(datframe[, missing], rate_ges)
  rates <- apply(mcar_mis, 2, function(x) { sum(is.na(x)) / length(x) })
  
  mar_column_discrete <- dplyr::ntile(datframe[, mar_column], n = round(nrow(datframe) / 20))
  
  mar_mis <- mcar_mis
  for (j in 1:length(missing)) {
    
    current_column <- missing[j]
    
    if (current_column == mar_column) { next }
    
    anz_mis <- sum(is.na(mcar_mis[, current_column]))
    
    levels_mar_column_discrete <- sort(unique(mar_column_discrete))
    
    props <- sort(runif(length(levels_mar_column_discrete)), decreasing = TRUE)

    mar_column_hfgk <- as.vector(table(mar_column_discrete))
    
    tmp <- props * mar_column_hfgk
    tmp <- tmp / sum(tmp)
    
    anz_mis_in_level <- round(tmp * anz_mis)
    
    mar_current <- datframe[, current_column]
    
    for (k in 1:length(anz_mis_in_level)) {
      if (anz_mis_in_level[k] > mar_column_hfgk[k]) {
        diff <- anz_mis_in_level[k] - mar_column_hfgk[k]
        anz_mis_in_level[k] <- mar_column_hfgk[k]
        anz_mis_in_level[k + 1] <- anz_mis_in_level[k + 1] + diff
      }
    }
    
    for (i in 1:length(levels_mar_column_discrete)) {
      indic_cur_mar_column <- which(mar_column_discrete == levels_mar_column_discrete[i])
      getsNA <- sample(indic_cur_mar_column, anz_mis_in_level[i])
      mar_current[getsNA] <- NA
    }
    
    mar_mis[, current_column] <- mar_current
  }

  mar_mis[, mar_column] <- datframe[, mar_column]
  
  return(mar_mis)
}

for (sim_dir in sim_dirs) {
  
  out_dir <- paste0(sim_dir, "_MAR")
  dir.create(file.path(path, "../Data/Simulation_Missing", out_dir),
             recursive = TRUE, showWarnings = FALSE)
  
  for (iter in 1:nsim) {
    for (p in c(0.1, 0.25, 0.5)) {
      
      set.seed(01042026 + iter)
      
      
      load(file.path(path, "../Data/Simulation", sim_dir, paste0("data_", iter, ".RData")))
      
      # using function to create MAR data
      data <- create_mar_vse(datframe = data, missing = setdiff(names(data), c("Trt", "Y", "true_cate")), mar_column = "X1", rate_ges = p)
      
      save(data, file = file.path(path, "../Data/Simulation_Missing", out_dir,
                                  paste0("data_", iter, "_", p * 100, ".RData")))
    }
  }
}