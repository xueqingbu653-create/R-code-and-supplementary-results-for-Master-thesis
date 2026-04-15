## Differences across treatments (ICBT vs SSRI)

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(broom.mixed)
library(performance)
library(ggplot2)
library(see)
library(stringr)

## 1. Data preparation
df <- combined_treat

df <- df %>%
  mutate(
    subject = factor(subject),
    sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
    treatment = factor(treatment, levels = c("ICBT", "SSRI"))
  )

# Missing values check
print(colSums(is.na(df)))

## 2. Centering / standardization
cont_vars <- c(
  "age",
  "SP",
  "NRS",
  "apnea",
  "sleepiness",
  "KSQ_total",
  "KSS_pre",
  "LSAS_pre"
)

for (v in cont_vars) {
  df[[paste0(v, "_c")]] <- as.numeric(scale(df[[v]], center = TRUE, scale = FALSE))
  df[[paste0(v, "_z")]] <- as.numeric(scale(df[[v]], center = TRUE, scale = TRUE))
}

## 3. Long format for LMM
df_long <- df %>%
  pivot_longer(
    cols = c(LSAS_pre, LSAS_post),
    names_to = "time",
    values_to = "LSAS"
  ) %>%
  mutate(
    time = factor(
      time,
      levels = c("LSAS_pre", "LSAS_post"),
      labels = c("pre", "post")
    )
  )

print(head(df_long))
print(table(df_long$time, df_long$treatment))

## 4. Descriptive statistics
desc_table <- df_long %>%
  group_by(treatment, time) %>%
  summarise(
    n = sum(!is.na(LSAS)),
    mean_LSAS = mean(LSAS, na.rm = TRUE),
    sd_LSAS = sd(LSAS, na.rm = TRUE),
    .groups = "drop"
  )

print(desc_table)

## 5. Main LMM
## LSAS ~ time * treatment * sleep_predictor + age + sex + (1 | subject)
run_lmm <- function(data_long, sleep_var_z) {
  
  formula_txt <- paste0(
    "LSAS ~ time * treatment * ", sleep_var_z,
    " + age_z + sex + (1 | subject)"
  )
  
  fit <- lmer(
    as.formula(formula_txt),
    data = data_long,
    REML = FALSE,
    na.action = na.omit
  )
  
  return(fit)
}

sleep_predictors <- c(
  "SP_z",
  "NRS_z",
  "apnea_z",
  "sleepiness_z",
  "KSQ_total_z",
  "KSS_pre_z"
)

model_list <- lapply(sleep_predictors, function(x) run_lmm(df_long, x))
names(model_list) <- sleep_predictors

## Inspect one model
print(summary(model_list$SP_z))
print(anova(model_list$SP_z))
check_model(model_list$SP_z)

## 6. Extract fixed effects
all_fixed_results <- bind_rows(
  lapply(names(model_list), function(pred) {
    broom.mixed::tidy(model_list[[pred]], effects = "fixed") %>%
      mutate(predictor = pred)
  })
)

print(all_fixed_results)

## Three-way interactions only
three_way_results <- all_fixed_results %>%
  filter(str_detect(term, "^timepost:treatmentSSRI:")) %>%
  select(predictor, term, estimate, std.error, statistic, p.value) %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "fdr"),
    p_Holm = p.adjust(p.value, method = "holm"),
    p_Bonferroni = p.adjust(p.value, method = "bonferroni")
  ) %>%
  arrange(p.value)

print(three_way_results)

## Export
write.csv(all_fixed_results, "all_fixed_effects_results.csv", row.names = FALSE)
write.csv(three_way_results, "three_way_interaction_results.csv", row.names = FALSE)

## 7. Sensitivity analysis: change score
df_change <- df %>%
  mutate(
    LSAS_improve = LSAS_pre - LSAS_post
  )

run_change_model <- function(data, sleep_var_z) {
  formula_txt <- paste0(
    "LSAS_improve ~ treatment * ", sleep_var_z,
    " + LSAS_pre_z + age_z + sex"
  )
  
  lm(as.formula(formula_txt), data = data, na.action = na.omit)
}

change_model_list <- lapply(sleep_predictors, function(x) run_change_model(df_change, x))
names(change_model_list) <- sleep_predictors

## Extract results
change_results <- bind_rows(
  lapply(names(change_model_list), function(pred) {
    broom::tidy(change_model_list[[pred]]) %>%
      mutate(predictor = pred)
  })
)

print(change_results)

## treatment × sleep interaction
change_interaction_results <- change_results %>%
  filter(str_detect(term, "^treatmentSSRI:")) %>%
  select(predictor, term, estimate, std.error, statistic, p.value) %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "fdr"),
    p_Holm = p.adjust(p.value, method = "holm"),
    p_Bonferroni = p.adjust(p.value, method = "bonferroni")
  ) %>%
  arrange(p.value)

print(change_interaction_results)

## Export
write.csv(change_results, "change_model_results.csv", row.names = FALSE)
write.csv(change_interaction_results, "change_model_interaction_results.csv", row.names = FALSE)

## Inspect one change-score model
print(summary(change_model_list$SP_z))

## 8. Overall treatment effect model
fit_overall <- lmer(
  LSAS ~ time * treatment + age_z + sex + (1 | subject),
  data = df_long,
  REML = FALSE,
  na.action = na.omit
)

print(summary(fit_overall))
print(anova(fit_overall))
check_model(fit_overall)

## Estimated marginal means
emm_overall <- emmeans(fit_overall, ~ time * treatment)
print(emm_overall)

## Pre-post contrasts within each treatment
emm_overall_pairs <- pairs(emm_overall, by = "treatment")
print(emm_overall_pairs)

## Compare treatments within each time point
emm_treatment_by_time <- pairs(emm_overall, by = "time")
print(emm_treatment_by_time)

write.csv(as.data.frame(emm_overall), "EMMeans_time_by_treatment.csv", row.names = FALSE)
write.csv(as.data.frame(emm_overall_pairs), "EMMeans_pre_post_contrasts_by_treatment.csv", row.names = FALSE)
write.csv(as.data.frame(emm_treatment_by_time), "EMMeans_treatment_contrasts_within_time.csv", row.names = FALSE)

## 9. Plot data
plot_data <- df_long %>%
  group_by(treatment, time) %>%
  summarise(
    n = sum(!is.na(LSAS)),
    mean_LSAS = mean(LSAS, na.rm = TRUE),
    sd_LSAS = sd(LSAS, na.rm = TRUE),
    se = sd_LSAS / sqrt(n),
    ci_low = mean_LSAS - qt(0.975, df = n - 1) * se,
    ci_high = mean_LSAS + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

summary_data <- plot_data

## Export descriptive tables
write.csv(desc_table, "descriptive_table_LSAS.csv", row.names = FALSE)
write.csv(plot_data, "plot_data_basic_figure.csv", row.names = FALSE)
write.csv(summary_data, "plot_data_upgraded_figure.csv", row.names = FALSE)

#### visualization ####
## 10. Figure 1: basic mean plot
p_basic <- ggplot(
  plot_data,
  aes(x = time, y = mean_LSAS, group = treatment, color = treatment, linetype = treatment)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.08,
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c(
      "ICBT" = "#1F78B4",
      "SSRI" = "#D55E00"
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Change in LSAS from Pre- to Post-treatment",
    x = "Time",
    y = "LSAS score",
    color = "Treatment",
    linetype = "Treatment"
  )

print(p_basic)

ggsave(
  "Figure_basic_LSAS_pre_post_by_treatment.png",
  p_basic,
  width = 6.5,
  height = 4.5,
  dpi = 300
)


## 11. Figure 2: trajectories + group means + 95% CI
p_upgraded <- ggplot(df_long, aes(x = time, y = LSAS, group = subject)) +
  geom_line(color = "grey70", alpha = 0.35, linewidth = 0.5) +
  geom_point(color = "grey70", alpha = 0.35, size = 1.2) +
  geom_ribbon(
    data = summary_data,
    aes(x = time, ymin = ci_low, ymax = ci_high, fill = treatment, group = treatment),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  geom_line(
    data = summary_data,
    aes(x = time, y = mean_LSAS, color = treatment, group = treatment),
    inherit.aes = FALSE,
    linewidth = 1.4
  ) +
  geom_point(
    data = summary_data,
    aes(x = time, y = mean_LSAS, color = treatment, group = treatment),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "ICBT" = "#1F78B4",
      "SSRI" = "#D55E00"
    )
  ) +
  scale_fill_manual(
    values = c(
      "ICBT" = "#1F78B4",
      "SSRI" = "#D55E00"
    )
  ) +
  
  facet_wrap(~ treatment) +
  theme_bw(base_size = 14) +
  labs(
    title = "LSAS Trajectories from Pre- to Post-treatment",
    subtitle = "Thin lines represent individual participants; thick lines represent group means with 95% CIs",
    x = "Time",
    y = "LSAS score",
    color = "Treatment",
    fill = "Treatment"
  )

print(p_upgraded)

ggsave(
  "Figure_upgraded_LSAS_trajectories.png",
  p_upgraded,
  width = 9,
  height = 5,
  dpi = 300
)

## 12. Figure 3: combined treatment figure
p_combined <- ggplot() +
  geom_line(
    data = df_long,
    aes(x = time, y = LSAS, group = subject),
    color = "grey70",
    alpha = 0.35,
    linewidth = 0.45
  ) +
  geom_ribbon(
    data = summary_data,
    aes(x = time, ymin = ci_low, ymax = ci_high, fill = treatment, group = treatment),
    alpha = 0.18
  ) +
  geom_line(
    data = summary_data,
    aes(x = time, y = mean_LSAS, group = treatment, color = treatment, linetype = treatment),
    linewidth = 1.4
  ) +
  geom_point(
    data = summary_data,
    aes(x = time, y = mean_LSAS, group = treatment, color = treatment),
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "ICBT" = "#1F78B4",
      "SSRI" = "#D55E00"
    )
  ) +
  scale_fill_manual(
    values = c(
      "ICBT" = "#1F78B4",
      "SSRI" = "#D55E00"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "ICBT" = "solid",
      "SSRI" = "dashed"
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Pre- to Post-treatment Change in LSAS by Treatment",
    subtitle = "Individual trajectories with group means and 95% CIs",
    x = "Time",
    y = "LSAS score",
    color = "Treatment",
    fill = "Treatment",
    linetype = "Treatment"
  )

print(p_combined)

ggsave(
  "Figure_combined_LSAS_pre_post.png",
  p_combined,
  width = 7,
  height = 5,
  dpi = 300
)

## 13. Model fit summary table
model_fit_info <- bind_rows(
  lapply(names(model_list), function(pred) {
    tibble(
      predictor = pred,
      AIC = AIC(model_list[[pred]]),
      BIC = BIC(model_list[[pred]]),
      logLik = as.numeric(logLik(model_list[[pred]]))
    )
  })
)

print(model_fit_info)
write.csv(model_fit_info, "model_fit_information.csv", row.names = FALSE)



