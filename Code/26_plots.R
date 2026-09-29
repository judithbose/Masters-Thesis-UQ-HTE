# This script is used to generate the plots in Chapter D.

# libraries ----
library(dplyr)
library(ggplot2)
library(viridis)

# paths ----
path <- getwd()
path <- paste0(path, "/Results/Complete")

files <- list.files(path = path, pattern = "^data.*\\.RData$", full.names = TRUE)
files <- files[c(1:7, 9)]

# settings ----
col <- viridis(4)
methods <- c("BART", "Causal Forest", "L&C", "L&C cqr", "dr_rf", "dr_xgb", "cqr_dr_rf", "cqr_dr_xgb")

cols <- c("Own Approach" = col[3], "Homogeneous Intervals" = col[2],
          "Heterogeneous Intervals" = "black")

# function to plot ----
create_plots <- function(file){
  
  load(file)
  
  # Coverage ----
  cover_df <- results %>%
    mutate(covered = results$true_cate >= lower & results$true_cate <= upper) %>%
    group_by(sim, rep, method, conf_level) %>%
    summarise(coverage = mean(covered, na.rm = TRUE), .groups = "drop")
  
  cover_df <- cover_df %>%
    mutate(Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Heterogeneous Intervals"),
           method = factor(method, levels = methods))
  
  p1 <- ggplot(cover_df, aes(x = method, y = coverage, fill = Group, colour = Group)) +
    
    geom_hline(aes(yintercept = conf_level), color = "red", linetype = "dashed") +
    
    geom_boxplot(width = 0.7) +
    
    facet_wrap(~ conf_level, 
               labeller = labeller(conf_level = c("0.8" = "80% Interval", "0.95" = "95% Interval"))) +
    
    scale_fill_manual(values = scales::alpha(cols, 0.2)) +
    scale_color_manual(values = cols) +
    
    labs(title = "", x = "Method", y = "Coverage", 
         fill = "Type of Method", colour = "Type of Method") +
    
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Width Interval ----
  width_df <- results %>%
    mutate(width = upper - lower) %>%
    group_by(sim, method, conf_level) %>%
    summarise(mean_width = mean(width, na.rm = TRUE), .groups = "drop")
  
  width_df <- width_df %>%
    mutate(Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Heterogeneous Intervals"),
           method = factor(method, levels = methods))
  
  p2 <- ggplot(width_df, aes(x = method, y = mean_width, fill = Group, colour = Group)) +
    
    geom_boxplot(width = 0.7) +
    
    facet_wrap(~ conf_level, 
               labeller = labeller(conf_level = c("0.8" = "80% Interval", "0.95" = "95% Interval"))) +
    
    scale_fill_manual(values = scales::alpha(cols, 0.2)) +
    scale_color_manual(values = cols) +
    
    labs(title = "", x = "Method", y = "Interval Width", 
         fill = "Type of Method", colour = "Type of Method") +
    
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major.x = element_blank())
  
  # Sample Intervals ----
  results <- results %>%
    mutate(Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Heterogeneous Intervals"),
           method = factor(method, levels = methods))
  
  set.seed(15062026)
  results <- results[results$conf_level == "0.8", ]
  
  selected_ids <- results %>%
    filter(sim == 1) %>%
    arrange(true_cate) %>%
    slice(round(seq(1, n(), length.out = 30))) %>%
    pull(person_id)
  
  plot_df <- results %>%
    filter(sim == 1, person_id %in% selected_ids)

  plot_df <- plot_df %>% arrange(method, true_cate) %>%
    group_by(method) %>%
    mutate(person_order = row_number()) %>%
    ungroup()
  
  plot_df$contains_true <- with(plot_df, true_cate >= lower & true_cate <= upper)
  
  df_true  <- subset(plot_df, contains_true)
  df_false <- subset(plot_df, !contains_true)
  
  p4 <- ggplot() +
    
    geom_point(data = plot_df, aes(x = person_order, y = true_cate, colour = Group)) +
    
    geom_errorbar(data = df_true, width = 0.2,
                  aes(x = person_order, ymin = lower, ymax = upper, colour = Group)) +
    
    geom_errorbar(data = df_false, aes(x = person_order, ymin = lower, ymax = upper),
                  colour = "red", width = 0.2, show.legend = FALSE) +
    
    facet_wrap(~ method, nrow = 2, scales = "fixed") +
    coord_flip() +
    
    scale_color_manual(values = cols) +
    labs(x = "Individuals", y = "Interval", title = "", colour = "Type of Method")
  
  
  return(list(p1, p2, p4))
}

# save plots ----
for (file in files) {
  
  plt <- create_plots(file)
  
  pdf(paste0(path, "/Plots/", tools::file_path_sans_ext(basename(file)), ".pdf"),
      width = 9, height = 6)
  
  print(plt[[1]])
  print(plt[[2]])
  print(plt[[3]])
  
  dev.off()
}

