# This script is used to generate the plots in Chapter E.

# libraries ----
library(dplyr)
library(ggplot2)
library(viridis)

# paths ----
path <- getwd()
path <- paste0(path, "/Results/Incomplete")

files <- list.files(path = path, pattern = "^data.*\\.RData$", full.names = TRUE)

# settings ----
col <- viridis(4)
methods <- c("BART", "Causal Forest", "cqr_dr_rf_conservative", "cqr_dr_rf_mean", 
             "cqr_dr_rf_median", "cqr_dr_rf_stacked",
             "cqr_dr_xgb_conservative", "cqr_dr_xgb_mean", 
             "cqr_dr_xgb_median", "cqr_dr_xgb_stacked", "zaffran_rf", "zaffran_xgb")

cols <- c("Own Approach" = col[3], "Homogeneous Intervals" = col[2],
          "Heterogeneous Intervals" = "black")


# function to plot ----
create_plots <- function(file){
  
  load(file)
  
  # Coverage ----
  
  # compute coverage
  cover_df <- results %>%
    mutate(covered = true_cate >= lower & true_cate <= upper) %>%
    group_by(sim, proportion, method, conf_level) %>%
    summarise(coverage = mean(covered, na.rm = TRUE), .groups = "drop") %>%
    
    mutate(Group = case_when(method %in% c("zaffran_rf", "zaffran_xgb") ~ "Heterogeneous Intervals",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Own Approach"),
           
           Subgroup = case_when(method %in% c("BART", "Causal Forest") ~ "Homogeneous",
                                grepl("^cqr_dr_rf", method) ~ "CQR-RF",
                                grepl("^cqr_dr_xgb", method) ~ "CQR-XGB",
                                grepl("^zaffran", method) ~ "Heterogeneous",
                                TRUE ~ "Other"),
           
           suffix_rank = case_when(grepl("_mean$", method) ~ 1,
                                   grepl("_median$", method) ~ 2,
                                   grepl("_conservative$", method) ~ 3,
                                   grepl("_stacked$", method) ~ 4, 
                                   TRUE ~ 5))
  
  # order data
  subgroup_order <- c("Homogeneous", "Heterogeneous", "CQR-RF", "CQR-XGB")
  
  methods_ordered <- cover_df %>%
    distinct(Subgroup, suffix_rank, method) %>%
    arrange(factor(Subgroup, levels = subgroup_order), suffix_rank, method) %>%
    pull(method)
  
  cover_df <- cover_df %>%
    mutate(method = factor(method, levels = methods_ordered))
  
  # create lines within the plot
  main_boundaries <- cover_df %>%
    distinct(method, Group) %>%
    mutate(pos = as.numeric(method)) %>%
    group_by(Group) %>%
    summarise(max_pos = max(pos), .groups = "drop") %>%
    filter(max_pos < length(methods_ordered)) %>%
    pull(max_pos) + 0.5
  
  all_sub_boundaries <- cover_df %>%
    distinct(method, Subgroup) %>%
    mutate(pos = as.numeric(method)) %>%
    group_by(Subgroup) %>%
    summarise(max_pos = max(pos), .groups = "drop") %>%
    filter(max_pos < length(methods_ordered)) %>%
    pull(max_pos) + 0.5
  
  xgb_rf_boundary <- setdiff(all_sub_boundaries, main_boundaries)
  

  # plot 
  p1 <- ggplot(cover_df, aes(x = method, y = coverage, fill = Group, colour = Group)) +
    
    geom_vline(xintercept = xgb_rf_boundary, color = "gray60", linetype = "dotted", linewidth = 0.6) +
    geom_vline(xintercept = main_boundaries, color = "gray45", linetype = "dashed", linewidth = 0.6) +
    
    geom_hline(aes(yintercept = as.numeric(conf_level)), color = "red", linetype = "dashed", linewidth = 0.6) +
    
    geom_boxplot(width = 0.6, outlier.size = 0.8, outlier.alpha = 0.4) +
    
    facet_wrap(~ conf_level + proportion, 
               labeller = labeller(conf_level = c("0.8" = "80% Interval", "0.95" = "95% Interval"),
                                   proportion = c("10" = "Missingness Rate 10%", "25" = "Missingness Rate 25%", "50" = "Missingness Rate 50%"))) +
    
    scale_fill_manual(values = scales::alpha(cols, 0.25)) +
    scale_color_manual(values = cols) +
    
    labs(title = "", x = "Method", y = "Coverage", 
         fill = "Type of Method", colour = "Type of Method") +
    
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  # Width Interval ----
  
  # compute width
  width_df <- results %>%
    mutate(width = upper - lower) %>%
    group_by(sim, proportion, method, conf_level) %>%
    summarise(width = mean(width, na.rm = TRUE), .groups = "drop") %>%
    
    mutate(Group = case_when(method %in% c("zaffran_rf", "zaffran_xgb") ~ "Heterogeneous Intervals",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Own Approach"),
           
           Subgroup = case_when(method %in% c("BART", "Causal Forest") ~ "Homogeneous",
                                grepl("^cqr_dr_rf", method) ~ "CQR-RF",
                                grepl("^cqr_dr_xgb", method) ~ "CQR-XGB",
                                grepl("^zaffran", method) ~ "Heterogeneous",
                                TRUE ~ "Other"),
           
           suffix_rank = case_when(grepl("_mean$", method) ~ 1,
                                   grepl("_median$", method) ~ 2,
                                   grepl("_conservative$", method) ~ 3,
                                   grepl("_stacked$", method) ~ 4, 
                                   TRUE ~ 5))
  
  # order data 
  subgroup_order <- c("Homogeneous", "Heterogeneous", "CQR-RF", "CQR-XGB")
  
  methods_ordered <- width_df %>%
    distinct(Subgroup, suffix_rank, method) %>%
    arrange(factor(Subgroup, levels = subgroup_order), suffix_rank, method) %>%
    pull(method)
  
  width_df <- width_df %>%
    mutate(method = factor(method, levels = methods_ordered))
  
  # create lines within the plot
  main_boundaries <- width_df %>%
    distinct(method, Group) %>%
    mutate(pos = as.numeric(method)) %>%
    group_by(Group) %>%
    summarise(max_pos = max(pos), .groups = "drop") %>%
    filter(max_pos < length(methods_ordered)) %>%
    pull(max_pos) + 0.5
  
  all_sub_boundaries <- width_df %>%
    distinct(method, Subgroup) %>%
    mutate(pos = as.numeric(method)) %>%
    group_by(Subgroup) %>%
    summarise(max_pos = max(pos), .groups = "drop") %>%
    filter(max_pos < length(methods_ordered)) %>%
    pull(max_pos) + 0.5
  
  xgb_rf_boundary <- setdiff(all_sub_boundaries, main_boundaries)
  
  
  # plot 
  p2 <- ggplot(width_df, aes(x = method, y = width, fill = Group, colour = Group)) +
    
    geom_vline(xintercept = xgb_rf_boundary, color = "gray60", linetype = "dotted", linewidth = 0.6) +
    geom_vline(xintercept = main_boundaries, color = "gray45", linetype = "dashed", linewidth = 0.6) +
    
    geom_hline(aes(yintercept = as.numeric(conf_level)), color = "red", linetype = "dashed", linewidth = 0.6) +
    
    geom_boxplot(width = 0.6, outlier.size = 0.8, outlier.alpha = 0.4) +
    
    facet_wrap(~ conf_level + proportion, 
               labeller = labeller(conf_level = c("0.8" = "80% Interval", "0.95" = "95% Interval"),
                                   proportion = c("10" = "Missingness Rate 10%", "25" = "Missingness Rate 25%", "50" = "Missingness Rate 50%"))) +
    
    scale_fill_manual(values = scales::alpha(cols, 0.25)) +
    scale_color_manual(values = cols) +
    
    labs(title = "", x = "Method", y = "Interval Width", 
         fill = "Type of Method", colour = "Type of Method") +
    
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  return(list(p1, p2))
}

# loop ----
for (file in files) {
  
  plt <- create_plots(file)
  
  # save data
  pdf(paste0(path, "/Plots/", tools::file_path_sans_ext(basename(file)), ".pdf"),
      width = 9, height = 6)
  
  print(plt[[1]])
  print(plt[[2]])
  
  dev.off()
}
