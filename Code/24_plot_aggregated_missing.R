# This script is used to generate the plots in Chapters 4.4.

# libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

# path ----
path <- getwd()
path <- paste0(path, "/Results/Incomplete")

file_paths <- list.files(path = path, pattern = "^data.*\\.RData$", full.names = TRUE)

# combine into one data frame ----
df_all <- map_dfr(file_paths, function(file) {
  temp_env <- new.env()
  load(file, envir = temp_env)
  
  obj_name <- ls(temp_env)[1]
  data <- temp_env[[obj_name]]
  
  data <- data %>% mutate(file_name = basename(file))
  return(data)
})

# recode simulation setting name
df_all$file_name <- recode(df_all$file_name, dataB_MAR.RData = "IHDP B & MAR", dataB_MCAR.RData = "IHDP B & MCAR", 
                           dataA_MAR.RData = "IHDP A & MAR", dataA_MCAR.RData = "IHDP A & MCAR", 
                           data_I1200_MAR.RData = "Interaction (n = 1200) & MAR", data_I1200_MCAR.RData = "Interaction (n = 1200) & MCAR",
                           data_I600_MAR.RData = "Interaction (n = 600) & MAR", data_I600_MCAR.RData = "Interaction (n = 600) & MCAR",
                           data_L1200_MAR.RData = "Linear (n = 1200) & MAR", data_L1200_MCAR.RData = "Linear (n = 1200) & MCAR",
                           data_L600_MAR.RData = "Linear (n = 600) & MAR", data_L600_MCAR.RData = "Linear (n = 600) & MCAR",
                           data_N1200_MAR.RData = "Non-Linear (n = 1200) & MAR", data_N1200_MCAR.RData = "Non-Linear (n = 1200) & MCAR",
                           data_N600_MAR.RData = "Non-Linear (n = 600) & MAR", data_N600_MCAR.RData = "Non-Linear (n = 600) & MCAR")

df_all$missingness <- ifelse(df_all$file_name %in% c("IHDP A & MAR", "IHDP B & MAR", "Interaction (n = 1200) & MAR", "Interaction (n = 600) & MAR",
                                                     "Linear (n = 1200) & MAR", "Linear (n = 600) & MAR", "Non-Linear (n = 1200) & MAR", "Non-Linear (n = 600) & MAR"), "MAR", "MCAR")

# compute width and coverage
df_all <- df_all %>%
  mutate(width = upper - lower) %>%
  mutate(covered = true_cate >= lower & true_cate <= upper)

df_all <- df_all %>% group_by(file_name, method, conf_level, missingness) %>%
  summarise(mean_coverage = mean(covered, na.rm = TRUE), 
            mean_width = mean(width, na.rm = TRUE)) %>%
  mutate(scenario = paste(file_name))

df_all$scenario <- recode(df_all$scenario, `IHDP B & MAR` = "IHDP B", `IHDP B & MCAR` = "IHDP B", 
                              `IHDP A & MAR` = "IHDP A", `IHDP A & MCAR` = "IHDP A", 
                              `Interaction (n = 1200) & MAR` = "Interaction (n = 1200)", `Interaction (n = 1200) & MCAR` = "Interaction (n = 1200)",
                              `Interaction (n = 600) & MAR` = "Interaction (n = 600)", `Interaction (n = 600) & MCAR` = "Interaction (n = 600)",
                              
                              `Linear (n = 1200) & MAR` = "Linear (n = 1200)", `Linear (n = 1200) & MCAR` = "Linear (n = 1200)",
                              `Linear (n = 600) & MAR` = "Linear (n = 600)", `Linear (n = 600) & MCAR` = "Linear (n = 600)",
                              
                              `Non-Linear (n = 1200) & MAR` = "Non-Linear (n = 1200)", `Non-Linear (n = 1200) & MCAR` = "Non-Linear (n = 1200)",
                              `Non-Linear (n = 600) & MAR` = "Non-Linear (n = 600)", `Non-Linear (n = 600) & MCAR` = "Non-Linear (n = 600)")


method_order <- c("BART", "Causal Forest", "zaffran_rf", "zaffran_xgb",
                  "cqr_dr_rf_mean", "cqr_dr_rf_median", "cqr_dr_rf_stacked", "cqr_dr_rf_conservative",
                  "cqr_dr_xgb_mean", "cqr_dr_xgb_median", "cqr_dr_xgb_stacked", "cqr_dr_xgb_conservative")

df_all <- df_all %>% mutate(
  
  # group variable
  Group = case_when(method %in% c("zaffran_rf", "zaffran_xgb") ~ "Heterogeneous Approach",
                    method %in% c("BART", "Causal Forest") ~ "Homogeneous Approach",
                    TRUE ~ "Own Approach"),
  
  # order
  Group = factor(Group, levels = c("Homogeneous Approach", "Heterogeneous Approach",
                                   "Own Approach")),
  method = factor(method, levels = rev(method_order)))


plt <- list()
# Conf Level Plot ----
plt[[1]] <- ggplot(df_all %>% filter(conf_level == 0.95, missingness == "MAR"), 
                   aes(x = scenario, y = method, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#8E0152FF",  high = "#276419FF", 
                       midpoint = 0.95, name = "Coverage") +
  facet_wrap( ~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "95% Interval (MAR)",  x = "Setting", y = "Method")

plt[[2]] <- ggplot(df_all %>% filter(conf_level == 0.95, missingness == "MCAR"), 
                   aes(x = scenario, y = method, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#8E0152FF",  high = "#276419FF", 
                       midpoint = 0.95, name = "Coverage") +
  facet_wrap( ~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "95% Interval (MCAR)",  x = "Setting", y = "Method")


plt[[3]] <- ggplot(df_all %>% filter(conf_level == 0.8, missingness == "MAR"), 
                   aes(x = scenario, y = method, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#8E0152FF",  high = "#276419FF", 
                       midpoint = 0.8, name = "Coverage") +
  facet_wrap( ~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "80% Interval (MAR)",  x = "Setting", y = "Method")

plt[[4]] <- ggplot(df_all %>% filter(conf_level == 0.8, missingness == "MCAR"), 
                   aes(x = scenario, y = method, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#8E0152FF",  high = "#276419FF", 
                       midpoint = 0.8, name = "Coverage") +
  facet_wrap( ~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "80% Interval (MCAR)",  x = "Setting", y = "Method")

# Width Plot ----
plt[[5]] <- ggplot(df_all %>% filter(conf_level == 0.95, missingness == "MAR"), 
                   aes(x = scenario, y = method, fill = mean_width)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#276419FF",  high = "#8E0152FF", 
                       midpoint = mean(df_all[df_all$conf_level == 0.95, ]$mean_width),
                       name = "Interval
Width") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "95% Interval (MAR)", x = "Setting", y = "Method")

plt[[6]] <- ggplot(df_all %>% filter(conf_level == 0.95, missingness == "MCAR"), 
                   aes(x = scenario, y = method, fill = mean_width)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#276419FF",  high = "#8E0152FF", 
                       midpoint = mean(df_all[df_all$conf_level == 0.95, ]$mean_width),
                       name = "Interval
Width") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "95% Interval (MCAR)", x = "Setting", y = "Method")

plt[[7]] <- ggplot(df_all %>% filter(conf_level == 0.8, missingness == "MAR"), 
                   aes(x = scenario, y = method, fill = mean_width)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#276419FF",  high = "#8E0152FF", 
                       midpoint = mean(df_all[df_all$conf_level == 0.8, ]$mean_width),
                       name = "Interval
Width") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "80% Interval (MAR)", x = "Setting", y = "Method")

plt[[8]] <- ggplot(df_all %>% filter(conf_level == 0.8, missingness == "MCAR"), 
                   aes(x = scenario, y = method, fill = mean_width)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#276419FF",  high = "#8E0152FF", 
                       midpoint = mean(df_all[df_all$conf_level == 0.8, ]$mean_width),
                       name = "Interval
Width") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "80% Interval (MCAR)", x = "Setting", y = "Method")

pdf(paste0(path, "/Plots/Plot_Aggregated_Missing.pdf"),
    width = 9, height = 6)

plt[[1]] + plt[[2]]
plt[[3]] + plt[[4]]
plt[[5]] + plt[[6]]
plt[[7]] + plt[[8]]

dev.off()

