# This generates tables and plot in Chapter 5.3.

# libraries ----
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(readr)
library(tidyverse)

# paths ----
path <- getwd()
load(paste0(path, "/Results/Policy/results_iter_100.RData"))

policy <- results %>%
  mutate(contrast = recode(contrast, "1_vs_0" = "tau1", "2_vs_0" = "tau2")) %>%
  pivot_wider(id_cols = c(person_id, sim, conf_level, method),
              names_from = contrast, values_from = c(cate, lower, upper)) %>%
  rename(id = person_id, tau1_hat = cate_tau1, tau2_hat = cate_tau2, 
         tau1_lower = lower_tau1, tau2_lower = lower_tau2,
         tau1_upper = upper_tau1, tau2_upper = upper_tau2)

trt_all <- map_dfr(1:100, function(i) {
  
  env <- new.env()
  
  obj_name <- load(paste0(path, "/Data/Policy/TRUTH_", i, ".RData"), envir = env)
  
  truth_i <- env[[obj_name[1]]]
  
  truth_i %>% mutate(sim = i, person_id = row_number())
})

# create one data set ----
policy_data <- merge(policy, trt_all, by.x = c("sim", "id"), by.y = c("sim", "person_id"))
policy_data <- policy_data %>% rename(tau1_true = tau_1, tau2_true = tau_2,
                                      Y0 = Y_0, Y1 = Y_1, Y2 = Y_2)

policy_data <- policy_data %>%
  filter(method == "cqr_dr_xgb") %>%
  select(-method)

# helper functions ----

# to assign treatment
choose_treatment <- function(score_0, score_1, score_2) {
  
  score_matrix <- cbind(control = score_0, treatment_1 = score_1, 
                        treatment_2 = score_2)
  
  max.col(score_matrix, ties.method = "first") - 1
}

choose_treatment_safe <- function(tau1_hat, tau2_hat, tau1_lower, tau2_lower,
                                  max_allowed_harm = 0) {
  
  eligible_t1 <- tau1_lower > -max_allowed_harm
  eligible_t2 <- tau2_lower > -max_allowed_harm
  
  score_matrix <- cbind(control = 0, treatment_1 = ifelse(eligible_t1, tau1_hat, -Inf),
                        treatment_2 = ifelse(eligible_t2, tau2_hat, -Inf))
  
  max.col(score_matrix, ties.method = "first") - 1L
}

# to define policies 
add_policies <- function(data, lambda_values) {
  
  oracle_treatment <- max.col(as.matrix(data[, c("Y0", "Y1", "Y2")]), 
                              ties.method = "first") - 1
  
  data_out <- data %>% mutate(
    # point, lcb and maximum
    policy_point = choose_treatment(score_0 = 0, score_1 = tau1_hat, score_2 = tau2_hat),
    policy_lcb = choose_treatment(score_0 = 0, score_1 = tau1_lower, score_2 = tau2_lower),
    policy_maximum = choose_treatment(score_0 = 0, score_1 = tau1_upper, score_2 = tau2_upper),
    
    policy_safe = choose_treatment_safe(tau1_hat, tau2_hat, tau1_lower, tau2_lower,
                                        max_allowed_harm = 10),
    
    # baseline policies
    policy_none = 0L, policy_all_t1 = 1L, policy_all_t2 = 2L,
    
    # oracle policy
    policy_oracle = oracle_treatment)
  
  # penalized for different lambda
  for (lambda in lambda_values) {
    
    lambda_name <- str_replace(as.character(lambda), "\\.", "_")
    policy_name <- paste0("policy_penalized_lambda_", lambda_name)
    
    data_out[[policy_name]] <- choose_treatment(score_0 = 0,
      score_1 = data_out$tau1_hat - lambda * (data_out$tau1_upper - data_out$tau1_lower),
      score_2 = data_out$tau2_hat - lambda * (data_out$tau2_upper - data_out$tau2_lower))
  }
  
  data_out
}

# compute policy assignments ----
lambda_values <- c(0, 0.10, 0.25, 0.50, 0.75, 1.00)

policy_data <- add_policies(data = policy_data, lambda_values = lambda_values)

policy_columns <- names(policy_data) %>%
  keep(~ str_detect(.x, "^policy_"))

policy_columns

# evaluate policy ----
evaluate_policy <- function(data, policy_column) {
  
  assigned_treatment <- data[[policy_column]]
  
  y_matrix <- as.matrix(data[, c("Y0", "Y1", "Y2")])
  
  row_id <- seq_len(nrow(data))
  column_id <- assigned_treatment + 1
  
  y_selected <- y_matrix[cbind(row_id, column_id)]
  y_oracle <- apply(y_matrix, 1, max)
  
  regret <- y_oracle - y_selected
  
  harm_indicator <- y_selected < data$Y0
  harm_amount <- pmax(data$Y0 - y_selected, 0)
  
  tibble(policy = policy_column,
         policy_value = mean(y_selected),
         incremental_value_vs_none = mean(y_selected - data$Y0),
         mean_regret = mean(regret),
         harm_rate = mean(harm_indicator))
}

results_by_simulation <- policy_data %>%
  group_by(sim, conf_level) %>%
  group_split() %>%
  map_dfr(function(sim_data) {
    
    current_sim <- unique(sim_data$sim)
    current_conf_level <- unique(sim_data$conf_level)
    
    map_dfr(policy_columns, ~ evaluate_policy(data = sim_data, policy_column = .x)) %>%
      mutate(sim = current_sim, conf_level = current_conf_level)
  })


# table ----
policy_summary <- results_by_simulation %>%
  group_by(conf_level, policy) %>%
  summarise(n_simulations = n(),
    
    policy_value_mean = mean(policy_value, na.rm = TRUE),
    policy_value_sd = sd(policy_value, na.rm = TRUE),
    
    incremental_value_mean = mean(incremental_value_vs_none, na.rm = TRUE),
    incremental_value_sd = sd(incremental_value_vs_none, na.rm = TRUE),
    
    mean_regret_mean = mean(mean_regret, na.rm = TRUE),
    mean_regret_sd = sd(mean_regret, na.rm = TRUE),
    
    harm_rate_mean = mean(harm_rate, na.rm = TRUE),
    harm_rate_sd = sd(harm_rate, na.rm = TRUE),
    
    .groups = "drop") %>%
  arrange(conf_level, desc(policy_value_mean))

main_policies <- c("policy_point", "policy_lcb", "policy_oracle", "policy_maximum",
                   "policy_penalized_lambda_0_25", "policy_safe")

main_results_table <- policy_summary %>%
  filter(conf_level == 0.95, policy %in% main_policies) %>%
  mutate(policy = recode(policy, "policy_point" = "Point",
                         "policy_lcb" = "LCB",
                         "policy_penalized_lambda_0_25" = "Penalized",
                         "policy_oracle" = "Oracle",
                         "policy_maximum" = "Maximum",
                         "policy_safe" = "Safe"),
         policy = factor(policy, levels = c("Oracle", "Point", "LCB", "Maximum",
                                            "Penalized", "Safe")),
         `Policy value` = paste0(round(policy_value_mean * 100, 1),
                                 " (± ", round(policy_value_sd * 100, 2), ")"),
         Regret = paste0(round(mean_regret_mean, 2), " (± ", round(mean_regret_sd, 4), ")"),
         Harm = paste0(round(harm_rate_mean, 2), " (± ", round(harm_rate_sd, 4), ")")) %>%
  arrange(policy) %>%
  select(Policy = policy, `Policy value`, Regret, Harm)

write_csv(main_results_table, paste0(path, "/Results/Policy/Tables/Policy_95.csv"))

# treatment shares ----
treatment_assignments <- policy_data %>%
  select(sim, id, conf_level, all_of(policy_columns)) %>%
  pivot_longer(cols = all_of(policy_columns), names_to = "policy", values_to = "assigned_treatment") %>%
  mutate(assigned_treatment = recode(as.character(assigned_treatment), 
                                     "0" = "Control", "1" = "Treatment 1", "2" = "Treatment 2"))

treatment_shares_by_simulation <- treatment_assignments %>%
  group_by(sim, conf_level, policy, assigned_treatment) %>%
  summarise(treatment_share = n() / sum(n()), .groups = "drop")

treatment_share_summary <- treatment_shares_by_simulation %>%
  group_by(policy, assigned_treatment) %>%
  summarise(treatment_share = sum(treatment_share), .groups = "drop") %>%
  complete(policy, assigned_treatment = c("Control", "Treatment 1", "Treatment 2"),
           fill = list(n = 0)) %>%
  group_by(policy) %>%
  mutate(anteil = treatment_share / sum(treatment_share, na.rm = TRUE)) %>%
  ungroup()

plot_share <- treatment_share_summary %>%
  filter(policy %in% c("policy_point", "policy_lcb", "policy_penalized_lambda_0_25",
                       "policy_oracle", "policy_maximum", "policy_safe")) %>%
  mutate(policy = recode(policy, "policy_point" = "Point", "policy_lcb" = "LCB",
                         "policy_penalized_lambda_0_25" = "Penalized",
                         "policy_oracle" = "Oracle",
                         "policy_maximum" = "Maximum",
                         "policy_safe" = "Safe"),
         policy = factor(policy, levels = c("Oracle", "Point", "LCB", "Maximum", "Penalized", "Safe")),) %>%
  ggplot(aes(x = policy, y = anteil, fill = assigned_treatment)) +
  geom_col(position = "stack") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Share of assigned individuals", fill = "Assigned treatment") +
  theme_minimal()

pdf(paste0(path, "/Results/Policy/Plots/plot_share.pdf"), width = 9, height = 6)
plot_share
dev.off()

# lambda sensitivity analysis ----
lambda_sensitivity <- policy_summary %>%
  filter(conf_level == 0.95, str_detect(policy, "^policy_penalized_lambda_")) %>%
  mutate(lambda = policy %>% str_remove("^policy_penalized_lambda_") %>%
           str_replace("_", ".") %>% as.numeric()) %>%
  arrange(lambda)

lambda_sensitivity_long <- lambda_sensitivity %>%
  select(lambda, policy_value_mean, mean_regret_mean, harm_rate_mean) %>%
  pivot_longer(cols = -lambda, names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, "policy_value_mean" = "Policy value",
                         "mean_regret_mean" = "Mean regret",
                         "harm_rate_mean" = "Harm rate"))

plot_lambda <- ggplot(lambda_sensitivity_long, aes(x = lambda, y = value)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = expression(lambda), y = NULL) +
  theme_minimal()

pdf(paste0(path, "/Results/Policy/Plots/plot_lambda.pdf"), width = 9, height = 6)
plot_lambda
dev.off()

# regret distribution ----
individual_policy_results <- treatment_assignments %>%
  left_join(policy_data %>% select(sim, id, conf_level, Y0, Y1, Y2), 
            by = c("sim", "id", "conf_level")) %>%
  mutate(y_selected = case_when(assigned_treatment == "Control" ~ Y0,
                                assigned_treatment == "Treatment 1" ~ Y1,
                                assigned_treatment == "Treatment 2" ~ Y2),
         y_oracle = pmax(Y0, Y1, Y2), individual_regret = y_oracle - y_selected,
         harmed = y_selected < Y0)

plot_regret_distribution <- individual_policy_results %>%
  filter(conf_level == 0.95, policy %in% c("policy_point", "policy_lcb", "policy_maximum",
                                           "policy_penalized_lambda_0_25", "policy_safe")) %>%
  mutate(policy = recode(policy, "policy_point" = "Point", "policy_lcb" = "LCB",
      "policy_penalized_lambda_0_25" = "Penalized", "policy_maximum" = "Maximum",
      "policy_safe" = "Safe")) %>%
  ggplot(aes(x = policy, y = individual_regret)) +
  geom_boxplot(outlier.alpha = 0.15) +
  labs(x = NULL, y = "Individual regret relative to oracle") +
  guides(fill = "none") +
  theme_minimal()


pdf(paste0(path, "/Results/Policy/Plots/plot_regret.pdf"), width = 9, height = 6)
plot_regret_distribution
dev.off()