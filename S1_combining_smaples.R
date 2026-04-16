#### S1: Difference Test for Sample 1 and Sample 2 ####
#### Cohort effect for SAD ####
df1 <- cohort_SAD
sleep_vars <- c("SP", "NRS", "apnea", "sleepiness", "KSQ_total", "KSS_pre")

# Factor
df1$sex <- factor(df1$sex)
df1$sample <- factor(df1$sample)

library(dplyr)
library(purrr)
library(car)
library(ggplot2)
library(tidyr)
library(ggdist)
library(broom)
library(effectsize)
library(writexl)

# Normality test
shapiro_results_sad <- map_dfr(sleep_vars, function(var) {
  df1 %>%
    group_by(sample) %>%
    summarise(
      variable = var,
      W = shapiro.test(na.omit(.data[[var]]))$statistic,
      p_value = shapiro.test(na.omit(.data[[var]]))$p.value,
      .groups = "drop"
    )
})

shapiro_results_sad

# Example QQ plot for KSS_pre
qqnorm(df1$KSQ_total[df1$sample == "1"], main = "QQ Plot: KSS_pre, Sample 1")
qqline(df1$KSQ_total[df1$sample == "1"])

qqnorm(df1$KSS_pre[df1$sample == "2"], main = "QQ Plot: KSS_pre, Sample 2")
qqline(df1$KSS_pre[df1$sample == "2"])

# Homogeneity of variance
levene_results_sad <- lapply(sleep_vars, function(v) {
  formula <- as.formula(paste(v, "~ sample"))
  out <- leveneTest(formula, data = df1)
  data.frame(
    variable = v,
    F = out$`F value`[1],
    p_value = out$`Pr(>F)`[1]
  )
})

levene_results_sad <- do.call(rbind, levene_results_sad)
levene_results_sad

# Group comparisons
w_apnea_sad <- wilcox.test(apnea ~ sample, data = df1, exact = FALSE)

tt_SP_sad <- t.test(SP ~ sample, data = df1, var.equal = TRUE)
tt_NRS_sad <- t.test(NRS ~ sample, data = df1, var.equal = TRUE)
tt_sleep_sad <- t.test(sleepiness ~ sample, data = df1, var.equal = TRUE)
tt_KSQ_sad <- t.test(KSQ_total ~ sample, data = df1, var.equal = TRUE)
tt_KSS_sad <- t.test(KSS_pre ~ sample, data = df1, var.equal = TRUE)

# Effect sizes
g_SP_sad <- hedges_g(SP ~ sample, data = df1)
g_NRS_sad <- hedges_g(NRS ~ sample, data = df1)
g_sleep_sad <- hedges_g(sleepiness ~ sample, data = df1)
g_KSQ_sad <- hedges_g(KSQ_total ~ sample, data = df1)
g_KSS_sad <- hedges_g(KSS_pre ~ sample, data = df1)

rb_apnea_sad <- rank_biserial(apnea ~ sample, data = df1, ci = 0.95)

# Multiple comparisons correction
p_vals_sad <- c(
  apnea = w_apnea_sad$p.value,
  SP = tt_SP_sad$p.value,
  NRS = tt_NRS_sad$p.value,
  sleepiness = tt_sleep_sad$p.value,
  KSQ_total = tt_KSQ_sad$p.value,
  KSS_pre = tt_KSS_sad$p.value
)

results_sad_p <- data.frame(
  variable = names(p_vals_sad),
  p_raw = p_vals_sad,
  p_FDR = p.adjust(p_vals_sad, method = "fdr"),
  p_Holm = p.adjust(p_vals_sad, method = "holm"),
  p_Bonferroni = p.adjust(p_vals_sad, method = "bonferroni")
)

results_sad_p

#### Visualization ####
df_long1 <- df1 %>%
  select(sample, SP, NRS, apnea, sleepiness, KSQ_total, KSS_pre) %>%
  pivot_longer(-sample, names_to = "variable", values_to = "value")

p1 <- ggplot(df_long1, aes(x = sample, y = value, fill = sample)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.5) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Comparison of SAD groups across two samples", y = "Score") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("cohort_comparison_SAD_violin_box.pdf", plot = p1, width = 10, height = 8)

p1_rain <- ggplot(df_long1, aes(x = sample, y = value, fill = sample)) +
  stat_halfeye(
    adjust = 0.6,
    width = 0.6,
    .width = 0,
    justification = -0.25,
    point_colour = NA,
    alpha = 0.6
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.08,
    size = 1,
    alpha = 0.45
  ) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Sample", y = "Score") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("cohort_comparison_SAD_raincloud.pdf", plot = p1_rain, width = 10, height = 8)

# Organize SAD results into table
results_full_sad <- data.frame(
  variable = c("apnea", "SP", "NRS", "sleepiness", "KSQ_total", "KSS_pre"),
  test = c("Wilcoxon", "t-test", "t-test", "t-test", "t-test", "t-test"),
  statistic = c(
    unname(w_apnea_sad$statistic),
    unname(tt_SP_sad$statistic),
    unname(tt_NRS_sad$statistic),
    unname(tt_sleep_sad$statistic),
    unname(tt_KSQ_sad$statistic),
    unname(tt_KSS_sad$statistic)
  ),
  df = c(
    NA,
    unname(tt_SP_sad$parameter),
    unname(tt_NRS_sad$parameter),
    unname(tt_sleep_sad$parameter),
    unname(tt_KSQ_sad$parameter),
    unname(tt_KSS_sad$parameter)
  ),
  p_raw = unname(p_vals_sad),
  p_FDR = unname(p.adjust(p_vals_sad, method = "fdr")),
  p_Holm = unname(p.adjust(p_vals_sad, method = "holm")),
  p_Bonferroni = unname(p.adjust(p_vals_sad, method = "bonferroni")),
  effect_size = c(
    rb_apnea_sad$r_rank_biserial,
    g_SP_sad$Hedges_g,
    g_NRS_sad$Hedges_g,
    g_sleep_sad$Hedges_g,
    g_KSQ_sad$Hedges_g,
    g_KSS_sad$Hedges_g
  ),
  CI_low = c(
    rb_apnea_sad$CI_low,
    g_SP_sad$CI_low,
    g_NRS_sad$CI_low,
    g_sleep_sad$CI_low,
    g_KSQ_sad$CI_low,
    g_KSS_sad$CI_low
  ),
  CI_high = c(
    rb_apnea_sad$CI_high,
    g_SP_sad$CI_high,
    g_NRS_sad$CI_high,
    g_sleep_sad$CI_high,
    g_KSQ_sad$CI_high,
    g_KSS_sad$CI_high
  )
)

results_full_sad <- results_full_sad %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

results_full_sad
write_xlsx(results_full_sad, "cohort_results_SAD.xlsx")

#### Cohort effect for HC ####
df2 <- cohort_HC
sleep_vars <- c("SP", "NRS", "apnea", "sleepiness", "KSQ_total", "KSS_pre")

# Factor
df2$sex <- factor(df2$sex)
df2$sample <- factor(df2$sample)

# Normality test
shapiro_results_hc <- map_dfr(sleep_vars, function(var) {
  df2 %>%
    group_by(sample) %>%
    summarise(
      variable = var,
      W = shapiro.test(na.omit(.data[[var]]))$statistic,
      p_value = shapiro.test(na.omit(.data[[var]]))$p.value,
      .groups = "drop"
    )
})

shapiro_results_hc

# Homogeneity of variance
levene_results_hc <- lapply(sleep_vars, function(v) {
  formula <- as.formula(paste(v, "~ sample"))
  out <- leveneTest(formula, data = df2)
  data.frame(
    variable = v,
    F = out$`F value`[1],
    p_value = out$`Pr(>F)`[1]
  )
})

levene_results_hc <- do.call(rbind, levene_results_hc)
levene_results_hc

qqnorm(df1$KSQ_total[df1$sample == "1"],
       main = "QQ Plot: KSQ_total (Sample 1)")
qqline(df1$KSQ_total[df1$sample == "1"])
qqnorm(df1$KSQ_total[df1$sample == "2"],
       main = "QQ Plot: KSQ_total (Sample 2)")
qqline(df1$KSQ_total[df1$sample == "2"])

# Group comparisons
tt_SP_hc <- t.test(SP ~ sample, data = df2, var.equal = TRUE)
tt_KSQ_hc <- t.test(KSQ_total ~ sample, data = df2, var.equal = TRUE)

w_NRS_hc <- wilcox.test(NRS ~ sample, data = df2, exact = FALSE)
w_apnea_hc <- wilcox.test(apnea ~ sample, data = df2, exact = FALSE)
w_sleep_hc <- wilcox.test(sleepiness ~ sample, data = df2, exact = FALSE)
w_KSS_hc <- wilcox.test(KSS_pre ~ sample, data = df2, exact = FALSE)

# Effect sizes
# t-test → Hedges g
g_SP_hc <- hedges_g(SP ~ sample, data = df2)
g_KSQ_hc <- hedges_g(KSQ_total ~ sample, data = df2)

# Wilcoxon → rank-biserial
rb_NRS_hc <- rank_biserial(NRS ~ sample, data = df2, ci = 0.95)
rb_apnea_hc <- rank_biserial(apnea ~ sample, data = df2, ci = 0.95)
rb_sleep_hc <- rank_biserial(sleepiness ~ sample, data = df2, ci = 0.95)
rb_KSS_hc <- rank_biserial(KSS_pre ~ sample, data = df2, ci = 0.95)

# Multiple comparisons correction
p_vals_hc <- c(
  SP = tt_SP_hc$p.value,
  KSQ_total = tt_KSQ_hc$p.value,
  NRS = w_NRS_hc$p.value,
  apnea = w_apnea_hc$p.value,
  sleepiness = w_sleep_hc$p.value,
  KSS_pre = w_KSS_hc$p.value
)

results_hc_p <- data.frame(
  variable = names(p_vals_hc),
  p_raw = unname(p_vals_hc),
  p_FDR = unname(p.adjust(p_vals_hc, method = "fdr")),
  p_Holm = unname(p.adjust(p_vals_hc, method = "holm")),
  p_Bonferroni = unname(p.adjust(p_vals_hc, method = "bonferroni"))
)

results_hc_p

#### Visualization ####
df_long2 <- df2 %>%
  select(sample, SP, NRS, apnea, sleepiness, KSQ_total, KSS_pre) %>%
  pivot_longer(-sample, names_to = "variable", values_to = "value")

p2 <- ggplot(df_long2, aes(x = sample, y = value, fill = sample)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.5) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Comparison of HC groups across two samples", y = "Score") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("cohort_comparison_HC_violin_box.pdf", plot = p2, width = 10, height = 8)

p2_rain <- ggplot(df_long2, aes(x = sample, y = value, fill = sample)) +
  stat_halfeye(
    adjust = 0.6,
    width = 0.6,
    .width = 0,
    justification = -0.25,
    point_colour = NA,
    alpha = 0.6
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.08,
    size = 1,
    alpha = 0.45
  ) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Sample", y = "Score") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("cohort_comparison_HC_raincloud.pdf", plot = p2_rain, width = 10, height = 8)

# Organize HC results into table
results_full_hc <- data.frame(
  variable = c("SP", "KSQ_total", "NRS", "apnea", "sleepiness", "KSS_pre"),
  
  test = c("t-test", "t-test", "Wilcoxon", "Wilcoxon", "Wilcoxon", "Wilcoxon"),
  
  statistic = c(
    unname(tt_SP_hc$statistic),
    unname(tt_KSQ_hc$statistic),
    unname(w_NRS_hc$statistic),
    unname(w_apnea_hc$statistic),
    unname(w_sleep_hc$statistic),
    unname(w_KSS_hc$statistic)
  ),
  
  df = c(
    unname(tt_SP_hc$parameter),
    unname(tt_KSQ_hc$parameter),
    NA,
    NA,
    NA,
    NA
  ),
  
  p_raw = c(
    unname(tt_SP_hc$p.value),
    unname(tt_KSQ_hc$p.value),
    unname(w_NRS_hc$p.value),
    unname(w_apnea_hc$p.value),
    unname(w_sleep_hc$p.value),
    unname(w_KSS_hc$p.value)
  ),
  
  p_FDR = c(
    unname(results_hc_p$p_FDR[results_hc_p$variable == "SP"]),
    unname(results_hc_p$p_FDR[results_hc_p$variable == "KSQ_total"]),
    unname(results_hc_p$p_FDR[results_hc_p$variable == "NRS"]),
    unname(results_hc_p$p_FDR[results_hc_p$variable == "apnea"]),
    unname(results_hc_p$p_FDR[results_hc_p$variable == "sleepiness"]),
    unname(results_hc_p$p_FDR[results_hc_p$variable == "KSS_pre"])
  ),
  
  p_Holm = c(
    unname(results_hc_p$p_Holm[results_hc_p$variable == "SP"]),
    unname(results_hc_p$p_Holm[results_hc_p$variable == "KSQ_total"]),
    unname(results_hc_p$p_Holm[results_hc_p$variable == "NRS"]),
    unname(results_hc_p$p_Holm[results_hc_p$variable == "apnea"]),
    unname(results_hc_p$p_Holm[results_hc_p$variable == "sleepiness"]),
    unname(results_hc_p$p_Holm[results_hc_p$variable == "KSS_pre"])
  ),
  
  p_Bonferroni = c(
    unname(results_hc_p$p_Bonferroni[results_hc_p$variable == "SP"]),
    unname(results_hc_p$p_Bonferroni[results_hc_p$variable == "KSQ_total"]),
    unname(results_hc_p$p_Bonferroni[results_hc_p$variable == "NRS"]),
    unname(results_hc_p$p_Bonferroni[results_hc_p$variable == "apnea"]),
    unname(results_hc_p$p_Bonferroni[results_hc_p$variable == "sleepiness"]),
    unname(results_hc_p$p_Bonferroni[results_hc_p$variable == "KSS_pre"])
  ),
  
  effect_size = c(
    g_SP_hc$Hedges_g,
    g_KSQ_hc$Hedges_g,
    rb_NRS_hc$r_rank_biserial,
    rb_apnea_hc$r_rank_biserial,
    rb_sleep_hc$r_rank_biserial,
    rb_KSS_hc$r_rank_biserial
  ),
  
  CI_low = c(
    g_SP_hc$CI_low,
    g_KSQ_hc$CI_low,
    rb_NRS_hc$CI_low,
    rb_apnea_hc$CI_low,
    rb_sleep_hc$CI_low,
    rb_KSS_hc$CI_low
  ),
  
  CI_high = c(
    g_SP_hc$CI_high,
    g_KSQ_hc$CI_high,
    rb_NRS_hc$CI_high,
    rb_apnea_hc$CI_high,
    rb_sleep_hc$CI_high,
    rb_KSS_hc$CI_high
  )
)

results_full_hc <- results_full_hc %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

results_full_hc

write_xlsx(results_full_hc, "cohort_results_HC.xlsx")

# Combined visualization
df_all <- combined_data

df_long_all <- df_all %>%
  select(group, sample, SP, NRS, apnea, sleepiness, KSQ_total, KSS_pre) %>%
  pivot_longer(
    -c(group, sample),
    names_to = "variable",
    values_to = "value"
  )

df_long_all$sample <- factor(df_long_all$sample, levels = c("1", "2"))

df_long_all$group <- factor(
  df_long_all$group,
  levels = c(1, 0),
  labels = c("SAD", "HC")
)

p_all <- ggplot(df_long_all, aes(x = sample, y = value, fill = sample)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.4) +
  facet_grid(group ~ variable, scales = "free_y") +
  scale_fill_manual(values = c("1" = "#E64B35", "2" = "#4DBBD5")) +
  labs(
    x = "Sample",
    y = "Score",
    title = "Comparison of sleep variables across samples and groups"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

ggsave(
  "cohort_combined_violin1.pdf",
  plot = p_all,
  width = 12,
  height = 8
)

#### Descriptive analysis for combined sample ####
df <- combined_data
library(dplyr)
desc_table <- df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    SP_mean = mean(SP, na.rm = TRUE),
    SP_sd = sd(SP, na.rm = TRUE),
    NRS_mean = mean(NRS, na.rm = TRUE),
    NRS_sd = sd(NRS, na.rm = TRUE),
    apnea_mean = mean(apnea, na.rm = TRUE),
    apnea_sd = sd(apnea, na.rm = TRUE),
    sleepiness_mean = mean(sleepiness, na.rm = TRUE),
    sleepiness_sd = sd(sleepiness, na.rm = TRUE),
    KSQ_total_mean = mean(KSQ_total, na.rm = TRUE),
    KSQ_total_sd = sd(KSQ_total, na.rm = TRUE),
    KSS_pre_mean = mean(KSS_pre, na.rm = TRUE),
    KSS_pre_sd = sd(KSS_pre, na.rm = TRUE)
  )
write.csv(desc_table, "descriptive_statistics.csv", row.names = FALSE)

## sex
sex_table <- df %>%
  group_by(group, sex) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(group) %>%
  mutate(percent = round(n/sum(n)*100,1))
write.csv(sex_table, "sex_distribution.csv", row.names = FALSE)


## difference test
by(df$age, df$group, shapiro.test)
by(df$SP, df$group, shapiro.test)
by(df$NRS, df$group, shapiro.test)
by(df$apnea, df$group, shapiro.test)
by(df$sleepiness, df$group, shapiro.test)
by(df$KSQ_total, df$group, shapiro.test)
by(df$KSS_pre, df$group, shapiro.test)

# Mann–Whitney U tests for continuous variables
p_values <- data.frame(
  Variable = c("Age", "SP", "NRS", "Apnea", "Sleepiness", "KSQ_total", "KSS_pre"),
  p_value = c(
    wilcox.test(age ~ group, data = df)$p.value,
    wilcox.test(SP ~ group, data = df)$p.value,
    wilcox.test(NRS ~ group, data = df)$p.value,
    wilcox.test(apnea ~ group, data = df)$p.value,
    wilcox.test(sleepiness ~ group, data = df)$p.value,
    wilcox.test(KSQ_total ~ group, data = df)$p.value,
    wilcox.test(KSS_pre ~ group, data = df)$p.value
  )
)

# Format p values
p_values$p_value <- ifelse(
  p_values$p_value < 0.001,
  "<0.001",
  format(round(p_values$p_value, 3), nsmall = 3)
)

# Sex chi-square test
sex_test <- chisq.test(table(df$group, df$sex))

sex_result <- data.frame(
  Variable = "Sex",
  p_value = ifelse(
    sex_test$p.value < 0.001,
    "<0.001",
    format(round(sex_test$p.value, 3), nsmall = 3)
  )
)


results <- rbind(p_values, sex_result)

library(writexl)
write_xlsx(results, "group_difference_tests.xlsx")



