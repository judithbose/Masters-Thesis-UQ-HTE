# This script is used to generate the plots in Chapters 4.3.

# libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

# path ----
path <- getwd()
path <- paste0(path, "/Results/Complete")

file_paths <- list.files(path = path, pattern = "^data.*\\.RData$", full.names = TRUE)
file_paths <- file_paths[c(1:7, 9)]

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
df_all$file_name <- recode(df_all$file_name, dataB.RData = "IHDP B", data_I1200.RData = "Interaction (n = 1200)",
                           data_I600.RData = "Interaction (n = 600)", data_L1200.RData = "Linear (n = 1200)",
                           data_L600.RData = "Linear (n = 600)", data_N1200.RData = "Non-Linear (n = 1200)",
                           data_N600.RData = "Non-Linear (n = 600)", dataA.RData = "IHDP A")

# compute width and coverage
df_all <- df_all %>%
  mutate(width = upper - lower) %>%
  mutate(covered = true_cate >= lower & true_cate <= upper)

df_all <- df_all %>% group_by(file_name, method, conf_level) %>%
  summarise(mean_coverage = mean(covered, na.rm = TRUE), 
            mean_width = mean(width, na.rm = TRUE)) %>%
  mutate(scenario = paste(file_name))

method_order <- c("BART", "Causal Forest", "L&C", "L&C cqr", "dr_rf", "dr_xgb",
                  "cqr_dr_rf", "cqr_dr_xgb")

df_all <- df_all %>% mutate(
  
  # group variable
  Group = case_when(method %in% c("cqr_dr_rf", "cqr_dr_xgb") ~ "Own Approach",
                    method %in% c("BART", "Causal Forest") ~ "Homogeneous Approach",
                    TRUE ~ "Heterogeneous Approach"),
    
  # order
  Group = factor(Group, levels = c("Homogeneous Approach", "Heterogeneous Approach",
                                   "Own Approach")),
  method = factor(method, levels = rev(method_order)))

plt <- list()
## Conf Level Plot ----
plt[[1]] <- ggplot(df_all %>% filter(conf_level == 0.95), 
       aes(x = scenario, y = method, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#8E0152FF",  high = "#276419FF", 
                       midpoint = 0.95, name = "Coverage") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "95% Interval",  x = "Setting", y = "Method")


plt[[2]] <- ggplot(df_all %>% filter(conf_level == 0.8), 
       aes(x = scenario, y = method, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#8E0152FF",  high = "#276419FF", 
                       midpoint = 0.80, name = "Coverage") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "80% Interval",  x = "Setting", y = "Method")

## Width Plot ----
plt[[3]] <- ggplot(df_all %>% filter(conf_level == 0.95), 
       aes(x = scenario, y = method, fill = log(mean_width))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#276419FF",  high = "#8E0152FF", 
                       midpoint = mean(log(df_all[df_all$conf_level == 0.95, ]$mean_width)),
                       name = "Mean \nlog(Interval Width)") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "95% Interval", x = "Setting", y = "Method")

plt[[4]] <- ggplot(df_all %>% filter(conf_level == 0.80), 
                   aes(x = scenario, y = method, fill = log(mean_width))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#276419FF",  high = "#8E0152FF", 
                       midpoint = mean(log(df_all[df_all$conf_level == 0.80, ]$mean_width)),
                       name = "Mean \nlog(Interval Width)") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "80% Interval", x = "Setting", y = "Method")

# save plots ----
pdf(paste0(path, "/Plots/Plot_Aggregated.pdf"),
    width = 9, height = 6)

plt[[1]] + plt[[2]]
plt[[3]] + plt[[4]]

dev.off()


