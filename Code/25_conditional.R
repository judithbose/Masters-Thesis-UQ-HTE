# This script is used to generate the plots for evaluating conditional coverage
# in Chapter 4.3.

# libraries ----
library(dplyr)
library(ggplot2)
library(viridis)

# paths ----
path <- getwd()
path <- paste0(path, "/Results/Complete")

# settings ----
col <- viridis(4)
methods <- c("BART", "Causal Forest", "L&C", "L&C cqr", "dr_rf", "dr_xgb", "cqr_dr_rf", "cqr_dr_xgb")

cols <- c("Own Approach" = col[3], "Homogeneous Intervals" = col[2],
          "Heterogeneous Intervals" = "black")

# function to plot ----
plot <- function(results){
  cover_df <- results %>%
    mutate(covered = as.numeric(results$true_cate >= lower & results$true_cate <= upper)) %>%
    mutate(width = results$upper - results$lower) %>% 
    select(c(method, sim, rep, conf_level, covered, width))
  
  cover_df <- cover_df %>%
    group_by(method, sim, conf_level) %>%
    summarise(coverage = mean(covered), .groups = "drop")
  
  cover_df <- cover_df %>%
    mutate(Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
                             method %in% c("BART", "Causal Forest") ~ "Homogeneous Intervals",
                             TRUE ~ "Heterogeneous Intervals"),
           method = factor(method, levels = methods))
  
  p <- ggplot(cover_df, aes(x = method, y = coverage, color = Group)) +
    geom_hline(aes(yintercept = conf_level), color = "red", linetype = "dashed") +
    geom_point(size = 3) +
    scale_color_manual(values = cols) +
    facet_wrap(~ sim + conf_level, nrow = 2, dir = "v",
               labeller = labeller(sim = function(x) paste("Subsetting", x),
                                   conf_level = function(x) paste0(as.numeric(x) * 100, "% Interval"))) +
    labs(title = "Conditional Coverage by Method", x = "Method", 
         y = "Conditional Coverage", colour = "Type of Method") +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# IHDP A ----

# load data
load(paste0(path, "/dataA_con.RData"))

# save plot
pdf(paste0(path, "/Plots/plot_A_con.pdf"), width = 9, height = 6)
plot(results)
dev.off()

# IHDP B ----
load(paste0(path, "/dataB_con.RData"))

pdf(paste0(path, "/Plots/plot_B_con.pdf"), width = 9, height = 6)
plot(results)
dev.off()