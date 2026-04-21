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
# Normality check
by(df$age, df$group, shapiro.test)
by(df$SP, df$group, shapiro.test)
by(df$NRS, df$group, shapiro.test)
by(df$apnea, df$group, shapiro.test)
by(df$sleepiness, df$group, shapiro.test)
by(df$KSQ_total, df$group, shapiro.test)
by(df$KSS_pre, df$group, shapiro.test)

# Heretogeneity check
library(car)
df$group <- as.factor(df$group)
leveneTest(age ~ group, data = df)
leveneTest(SP ~ group, data = df)
leveneTest(NRS ~ group, data = df)
leveneTest(apnea ~ group, data = df)
leveneTest(sleepiness ~ group, data = df)
leveneTest(KSQ_total ~ group, data = df)
leveneTest(KSS_pre ~ group, data = df)


#t-test
vars <- c("age", "SP", "NRS", "apnea", "sleepiness", "KSQ_total", "KSS_pre")

p_values_t <- lapply(vars, function(v) {
  
  if (v == "SP") {
    res <- t.test(as.formula(paste(v, "~ group")), data = df)
  } else {
    res <- t.test(as.formula(paste(v, "~ group")), data = df, var.equal = TRUE)
  }
  
  data.frame(
    Variable = v,
    t_value = res$statistic,
    df = res$parameter,
    p_value = res$p.value
  )
})

p_values_t <- do.call(rbind, p_values_t)
print(p_values_t)

# Sex chi-square test
sex_test <- chisq.test(table(df$group, df$sex))
sex_test
sex_result <- data.frame(
  Variable = "Sex",
  p_value = ifelse(
    sex_test$p.value < 0.001,
    "<0.001",
    format(round(sex_test$p.value, 3), nsmall = 3)
  )
)

sex_result

#### Binary logistic regression analysis for combined sample ####
pkgs <- c(
  "openxlsx", "ggplot2", "pROC", "dplyr", "tidyr",
  "caret", "e1071", "kernlab", "gridExtra"
)

df <- combined_data

sleep_vars <- c("SP", "NRS", "apnea", "sleepiness", "KSQ_total", "KSS_pre")
covars <- c("age", "sex")
all_needed <- c("group", sleep_vars, covars)

df <- df[, all_needed]
df <- na.omit(df)

# standardization and centering
to_scale <- c(sleep_vars, "age")
df[, to_scale] <- scale(df[, to_scale], center = TRUE, scale = TRUE)

# group code
df$group_num <- ifelse(df$group %in% c(1, "1", "SAD"), 1, 0)
df$group_fac01 <- factor(df$group_num, levels = c(0, 1))
df$group_fac <- factor(df$group_num, levels = c(0, 1), labels = c("HC", "SAD"))
df$sex <- factor(df$sex)

str(df$group_fac01)
table(df$group_fac)

# Binary logistic regression
run_logistic <- function(var_list, data) {
  
  results_list <- lapply(var_list, function(v) {
    
    formula <- as.formula(paste("group_fac01 ~", v, "+ age + sex"))
    model <- glm(formula, data = data, family = binomial)
    
    coef_table <- summary(model)$coefficients
    cis <- suppressMessages(confint(model))
    
    data.frame(
      variable = v,
      term = rownames(coef_table),
      beta = coef_table[, "Estimate"],
      SE = coef_table[, "Std. Error"],
      z = coef_table[, "z value"],
      p = coef_table[, "Pr(>|z|)"],
      lowerCI_logit = cis[, 1],
      upperCI_logit = cis[, 2],
      row.names = NULL
    )
  })
  
  results_all <- do.call(rbind, results_list)
  
  main_results <- do.call(rbind, lapply(seq_along(var_list), function(i) {
    v <- var_list[i]
    tmp <- results_list[[i]]
    tmp[tmp$term %in% c(v, "age") | grepl("^sex", tmp$term), ]
  }))
  
  main_results$OR <- exp(main_results$beta)
  main_results$lowerCI <- exp(main_results$lowerCI_logit)
  main_results$upperCI <- exp(main_results$upperCI_logit)
  
  main_results$p_FDR <- p.adjust(main_results$p, method = "BH")
  main_results$p_Holm <- p.adjust(main_results$p, method = "holm")
  main_results$p_Bonferroni <- p.adjust(main_results$p, method = "bonferroni")
  
  main_results$OR_95CI <- paste0(
    sprintf("%.2f", main_results$OR),
    " (",
    sprintf("%.2f", main_results$lowerCI),
    "-",
    sprintf("%.2f", main_results$upperCI),
    ")"
  )
  
  format_p <- function(p) {
    ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))
  }
  
  final_table <- data.frame(
    Variable = main_results$variable,
    OR_95CI = main_results$OR_95CI,
    P_raw = main_results$p,
    FDR_raw = main_results$p_FDR,
    Holm_raw = main_results$p_Holm,
    Bonf_raw = main_results$p_Bonferroni,
    P_value = format_p(main_results$p),
    FDR = format_p(main_results$p_FDR),
    Holm = format_p(main_results$p_Holm),
    Bonferroni = format_p(main_results$p_Bonferroni),
    OR_num = main_results$OR,
    lower = main_results$lowerCI,
    upper = main_results$upperCI,
    beta = main_results$beta,
    stringsAsFactors = FALSE
  )
  
  list(main = final_table, full = results_all)
}

sleep_results <- run_logistic(sleep_vars, df)

openxlsx::write.xlsx(
  list(
    Sleep_main_results = sleep_results$main,
    Sleep_full_results = sleep_results$full
  ),
  file = "Binary_logistic_results_combined_sample.xlsx",
  overwrite = TRUE
)

cor(
  df[, intersect(c("SP", "NRS", "sleepiness", "KSQ_total", "KSS_pre"), colnames(df))],
  use = "complete.obs"
)

# Bonferroni-significant variables
sig_sleep_vars <- c("SP", "NRS", "sleepiness", "KSS_pre")
sig_sleep_vars

if (length(sig_sleep_vars) == 0) {
  stop("No Bonferroni-significant sleep variables were identified.")
}

formula_multi <- as.formula(
  paste("group_fac01 ~", paste(c(sig_sleep_vars, "age", "sex"), collapse = " + "))
)

model_multi <- glm(formula_multi, data = df, family = binomial)
summary(model_multi)

## final table
coef_table <- summary(model_multi)$coefficients
cis <- suppressMessages(confint(model_multi))
final_table <- data.frame(
  Variable = rownames(coef_table),
  Beta = round(coef_table[, "Estimate"], 3),
  OR = round(exp(coef_table[, "Estimate"]), 3),
  CI_lower = round(exp(cis[, 1]), 3),
  CI_upper = round(exp(cis[, 2]), 3),
  p = round(coef_table[, "Pr(>|z|)"], 3),
  row.names = NULL
)
final_table

final_table$CI_95 <- paste0(
  final_table$CI_lower,
  "–",
  final_table$CI_upper
)
final_table <- final_table[, c("Variable","Beta","OR","CI_95","p")]
final_table

#### visualization ####
# Forest plot
sleep_results$main$Sig3 <- with(
  sleep_results$main,
  ifelse(Bonf_raw < 0.05, "Bonferroni",
         ifelse(P_raw < 0.05, "P < 0.05", "NS"))
)

sleep_results$main$Sig3 <- factor(
  sleep_results$main$Sig3,
  levels = c("Bonferroni", "P < 0.05", "NS")
)

plot_forest_pub <- function(dat, title_text, base_size = 12, show_bonf = TRUE) {
  
  dat <- dat[order(dat$OR_num, decreasing = FALSE), ]
  dat$Variable <- factor(dat$Variable, levels = dat$Variable)
  
  p <- ggplot(dat, aes(x = Variable, y = OR_num, color = Sig3)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, color = "grey45") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18, linewidth = 0.7) +
    geom_point(size = 3) +
    coord_flip() +
    scale_y_log10() +
    scale_color_manual(values = c(
      "Bonferroni" = "#D62728",
      "P < 0.05" = "#FF7F0E",
      "NS" = "#4C72B0"
    )) +
    labs(
      title = title_text,
      x = NULL,
      y = "Odds ratio (95% CI)",
      color = NULL
    ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(face = "bold"),
      legend.position = "top"
    )
  
  if (show_bonf) {
    p <- p +
      geom_text(
        data = subset(dat, Bonf_raw < 0.05),
        aes(label = paste0("Bonf=", Bonferroni)),
        hjust = -0.1,
        size = 3.4,
        color = "black"
      ) +
      expand_limits(y = max(dat$upper) * 1.5)
  }
  
  return(p)
}

p_sleep <- plot_forest_pub(
  sleep_results$main,
  "Sleep-related variables",
  base_size = 12,
  show_bonf = TRUE
)

print(p_sleep)

ggsave(
  "Forestplot_sleep_combined_sample.pdf",
  p_sleep,
  width = 6,
  height = 4.5,
  device = cairo_pdf
)

## ROC for each significant variables (age+sex)
library(pROC)
roc_list_raw <- list()
roc_list_smooth <- list()
cutoff_list <- list()

auc_values <- numeric(length(sig_sleep_vars))
names(auc_values) <- sig_sleep_vars


for (i in seq_along(sig_sleep_vars)) {
  
  var <- sig_sleep_vars[i]
  
  formula_i <- as.formula(
    paste("group_fac01 ~", paste(c(var, "age", "sex"), collapse = " + "))
  )
  
  model_i <- glm(formula_i, data = df, family = binomial)
  prob_i <- predict(model_i, type = "response")
  
  roc_i_raw <- roc(df$group_num, prob_i, ci = TRUE, quiet = TRUE)
  roc_i_smooth <- smooth(roc_i_raw)
  
  best_i <- coords(
    roc_i_raw,
    "best",
    ret = c("threshold","sensitivity","specificity"),
    transpose = FALSE
  )
  
  cutoff_list[[i]] <- data.frame(
    variable_raw = var,
    threshold = best_i$threshold,
    sensitivity = best_i$sensitivity,
    specificity = best_i$specificity
  )
  
  auc_i <- as.numeric(auc(roc_i_raw))
  auc_values[var] <- auc_i
  
  roc_df_raw <- data.frame(
    FPR = 1 - roc_i_raw$specificities,
    TPR = roc_i_raw$sensitivities,
    variable_raw = var
  )
  
  roc_df_smooth <- data.frame(
    FPR = 1 - roc_i_smooth$specificities,
    TPR = roc_i_smooth$sensitivities,
    variable_raw = var
  )
  
  roc_list_raw[[i]] <- roc_df_raw
  roc_list_smooth[[i]] <- roc_df_smooth
}

roc_all_raw <- bind_rows(roc_list_raw)
roc_all_smooth <- bind_rows(roc_list_smooth)

auc_df <- data.frame(
  variable_raw = names(auc_values),
  auc = as.numeric(auc_values)
) %>%
  arrange(desc(auc))

cutoff_df <- bind_rows(cutoff_list)
roc_summary <- auc_df %>%
  left_join(cutoff_df, by = "variable_raw")

roc_summary


label_map <- setNames(
  paste0(auc_df$variable_raw, " (AUC = ", sprintf("%.3f", auc_df$auc), ")"),
  auc_df$variable_raw
)

roc_all_raw$variable <- factor(
  label_map[roc_all_raw$variable_raw],
  levels = label_map[auc_df$variable_raw]
)

roc_all_smooth$variable <- factor(
  label_map[roc_all_smooth$variable_raw],
  levels = label_map[auc_df$variable_raw]
)

palette_top <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0E442"
)
palette_use <- palette_top[seq_len(length(levels(roc_all_smooth$variable)))]


# raw ROC
p_multi_roc_sleep_raw <- ggplot(roc_all_raw, aes(x = FPR, y = TPR, color = variable)) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8, color = "grey70") +
  scale_color_manual(values = palette_use) +
  coord_equal() +
  labs(
    title = "ROC curves for individual sleep variables",
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    legend.position = "right",
    legend.text = element_text(size = 10.5),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_multi_roc_sleep_raw)
ggsave(
  "ROC_sleep_individual_combined_sample.pdf",
  plot = p_multi_roc_sleep_raw,
  width = 6.5,
  height = 5.5,
  device = cairo_pdf
)

## smoothed ROC
p_multi_roc_sleep <- ggplot(roc_all_smooth, aes(x = FPR, y = TPR, color = variable)) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8, color = "grey70") +
  scale_color_manual(values = palette_use) +
  coord_equal() +
  labs(
    title = "Smoothed ROC curves for individual sleep variables",
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    legend.position = "right",
    legend.text = element_text(size = 10.5),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_multi_roc_sleep)

ggsave(
  "ROC_sleep_individual_smoothed_combined_sample.pdf",
  plot = p_multi_roc_sleep,
  width = 6.5,
  height = 5.5,
  device = cairo_pdf
)

## ROC for all significant variables
formula_sig_sleep <- as.formula(
  paste("group_fac01 ~", paste(c(sig_sleep_vars, "age", "sex"), collapse = " + "))
)

model_sig_sleep <- glm(formula_sig_sleep, data = df, family = binomial)
summary(model_sig_sleep)

prob_sig_sleep <- predict(model_sig_sleep, type = "response")

roc_sleep <- roc(df$group_num, prob_sig_sleep, ci = TRUE, quiet = TRUE)

best_sleep <- coords(
  roc_sleep,
  "best",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

best_sleep

auc_sleep <- as.numeric(auc(roc_sleep))
auc_ci_sleep <- ci.auc(roc_sleep)

roc_df <- data.frame(
  FPR = 1 - roc_sleep$specificities,
  TPR = roc_sleep$sensitivities
)

p_roc_sleep <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(linewidth = 1.0, color = "#1F78B4") +
  geom_abline(intercept = 0, slope = 1,
              linetype = 2, linewidth = 0.8, color = "grey60") +
  coord_equal() +
  labs(
    title = "ROC curve for significant sleep variables",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  annotate(
    "label",
    x = 0.65, y = 0.15,
    label = paste0(
      "AUC = ", sprintf("%.3f", auc_sleep),
      "\n95% CI: ", sprintf("%.3f", auc_ci_sleep[1]),
      "–", sprintf("%.3f", auc_ci_sleep[3]),
      "\nN variables = ", length(sig_sleep_vars)
    ),
    size = 4.2,
    label.size = 0.3
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

print(p_roc_sleep)

ggsave(
  "ROC_sleep_significant_combined_sample.pdf",
  plot = p_roc_sleep,
  width = 5.5,
  height = 5,
  device = cairo_pdf
)

#### Machine learning: SVM ####
set.seed(123)

library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(openxlsx)
library(kernlab)
library(e1071)
library(grid)

## 1. Data preparation
analysis_vars <- c("SP", "NRS", "sleepiness", "KSS_pre", "age", "sex")
df_ml <- df[, c("group_fac", "group_num", analysis_vars)]

## Ensure outcome factor order is correct:
## HC = negative class, SAD = positive class
df_ml$group_fac <- factor(df_ml$group_fac, levels = c("HC", "SAD"))

## Optional: make sure sex is coded appropriately
## If sex is binary numeric (0/1), you may keep it numeric.
## If you want sex treated as categorical, uncomment the next line:
# df_ml$sex <- factor(df_ml$sex)

## Remove rows with missing values
df_ml <- na.omit(df_ml)


## 2. Fixed repeated CV splits
set.seed(123)
folds <- createMultiFolds(df_ml$group_fac, k = 10, times = 10)
index_out <- lapply(folds, function(ind) setdiff(seq_len(nrow(df_ml)), ind))

## 3. Fixed tuning grid for SVM
## Use only predictor matrix for sigma estimation
x_mat <- model.matrix(
  as.formula(paste("~", paste(analysis_vars, collapse = " + "))),
  data = df_ml
)[, -1, drop = FALSE]

sigma_est <- as.numeric(kernlab::sigest(x_mat, frac = 1)[2])

svm_grid <- expand.grid(
  sigma = sigma_est,
  C = 2^(-2:2)
)

## 4. Fixed seeds for caret
## repeatedcv: 5 folds x 10 repeats = 50 resamples
n_resamples <- 10 * 10
n_models_svm <- nrow(svm_grid)

set.seed(123)
seeds_svm <- vector("list", length = n_resamples + 1)
for (i in seq_len(n_resamples)) {
  seeds_svm[[i]] <- sample.int(100000, n_models_svm)
}
seeds_svm[[n_resamples + 1]] <- sample.int(100000, 1)

## Logistic regression has 1 model per resample
set.seed(456)
seeds_glm <- vector("list", length = n_resamples + 1)
for (i in seq_len(n_resamples)) {
  seeds_glm[[i]] <- sample.int(100000, 1)
}
seeds_glm[[n_resamples + 1]] <- sample.int(100000, 1)

## 5. TrainControl objects
ctrl_svm <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  index = folds,
  indexOut = index_out,
  seeds = seeds_svm,
  allowParallel = FALSE
)

ctrl_glm <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  index = folds,
  indexOut = index_out,
  seeds = seeds_glm,
  allowParallel = FALSE
)

## 6. Train SVM model
set.seed(123)
svm_fit <- train(
  as.formula(paste("group_fac ~", paste(analysis_vars, collapse = " + "))),
  data = df_ml,
  method = "svmRadial",
  metric = "ROC",
  trControl = ctrl_svm,
  preProcess = c("center", "scale"),
  tuneGrid = svm_grid
)

print(svm_fit)
svm_fit$bestTune
## Extract out-of-fold predictions for the best SVM tuning parameters
best_pred <- svm_fit$pred[
  svm_fit$pred$sigma == svm_fit$bestTune$sigma &
    svm_fit$pred$C == svm_fit$bestTune$C, ]

## Confusion matrix and performance metrics
cm_svm <- confusionMatrix(best_pred$pred, best_pred$obs, positive = "SAD")
print(cm_svm)

acc_svm  <- unname(cm_svm$overall["Accuracy"])
sens_svm <- unname(cm_svm$byClass["Sensitivity"])
spec_svm <- unname(cm_svm$byClass["Specificity"])

## ROC for SVM
roc_svm <- roc(
  response = best_pred$obs,
  predictor = best_pred$SAD,
  levels = c("HC", "SAD"),
  quiet = TRUE
)

auc_svm <- as.numeric(auc(roc_svm))
ci_svm <- ci.auc(roc_svm)

auc_svm
ci_svm
acc_svm
sens_svm
spec_svm

## 7. Plot SVM ROC
p_svm_roc <- ggplot(
  data.frame(
    FPR = 1 - roc_svm$specificities,
    TPR = roc_svm$sensitivities
  ),
  aes(x = FPR, y = TPR)
) +
  geom_line(linewidth = 1.0, color = "#D55E00") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey60") +
  coord_equal() +
  annotate(
    "label",
    x = 0.65, y = 0.15,
    label = paste0(
      "SVM (RBF)\n",
      "AUC = ", sprintf("%.3f", auc_svm),
      "\n95% CI: ", sprintf("%.3f", ci_svm[1]),
      "–", sprintf("%.3f", ci_svm[3])
    ),
    size = 4.2,
    label.size = 0.3
  ) +
  labs(
    title = "Cross-validated ROC curve for SVM",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  )

print(p_svm_roc)

ggsave(
  "ROC_SVM_cv.pdf",
  p_svm_roc,
  width = 5.5,
  height = 5,
  device = cairo_pdf
)

## 8. Train logistic regression using the same CV splits
set.seed(123)
logit_fit <- train(
  as.formula(paste("group_fac ~", paste(analysis_vars, collapse = " + "))),
  data = df_ml,
  method = "glm",
  family = binomial,
  metric = "ROC",
  trControl = ctrl_glm
)

print(logit_fit)

logit_pred <- logit_fit$pred

roc_logit <- roc(
  response = logit_pred$obs,
  predictor = logit_pred$SAD,
  levels = c("HC", "SAD"),
  quiet = TRUE
)

auc_logit <- as.numeric(auc(roc_logit))
ci_logit <- ci.auc(roc_logit)

auc_logit
ci_logit

## DeLong test
delong_test <- roc.test(roc_logit, roc_svm, method = "delong")

## 9. Logistic vs SVM ROC comparison plot
label_logit <- paste0(
  "Logistic regression\n",
  "AUC = ", sprintf("%.3f", auc_logit), "\n",
  "95% CI ", sprintf("%.3f", ci_logit[1]), "–", sprintf("%.3f", ci_logit[3])
)

label_svm <- paste0(
  "SVM (RBF)\n",
  "AUC = ", sprintf("%.3f", auc_svm), "\n",
  "95% CI ", sprintf("%.3f", ci_svm[1]), "–", sprintf("%.3f", ci_svm[3])
)

roc_df_logit <- data.frame(
  FPR = 1 - roc_logit$specificities,
  TPR = roc_logit$sensitivities,
  Model = label_logit
)

roc_df_svm <- data.frame(
  FPR = 1 - roc_svm$specificities,
  TPR = roc_svm$sensitivities,
  Model = label_svm
)

roc_df_all <- rbind(roc_df_logit, roc_df_svm)

p_cv_roc <- ggplot(roc_df_all, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey60") +
  coord_equal() +
  scale_color_manual(
    values = setNames(
      c("#0072B2", "#D55E00"),
      c(label_logit, label_svm)
    )
  ) +
  labs(
    title = paste0(
      "Cross-validated ROC curves\n",
      "DeLong test: p = ", sprintf("%.3f", delong_test$p.value)
    ),
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )

print(p_cv_roc)

ggsave(
  "ROC_logistic_vs_SVM_cv.pdf",
  plot = p_cv_roc,
  width = 6,
  height = 5,
  device = cairo_pdf
)

## 10. Calibration plot
make_calibration_df <- function(obs, pred, n_bins = 10) {
  d <- data.frame(obs = obs, pred = pred)
  
  ## protect against duplicated quantile breaks
  q_breaks <- quantile(
    d$pred,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  )
  q_breaks <- unique(q_breaks)
  
  if (length(q_breaks) < 3) {
    stop("Not enough unique predicted probabilities to create calibration bins.")
  }
  
  d$bin <- cut(
    d$pred,
    breaks = q_breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  d %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(pred, na.rm = TRUE),
      obs_rate = mean(obs == "SAD", na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}

cal_logit <- make_calibration_df(logit_pred$obs, logit_pred$SAD, n_bins = 10)
cal_svm <- make_calibration_df(best_pred$obs, best_pred$SAD, n_bins = 10)

brier_logit <- mean((ifelse(logit_pred$obs == "SAD", 1, 0) - logit_pred$SAD)^2)
brier_svm <- mean((ifelse(best_pred$obs == "SAD", 1, 0) - best_pred$SAD)^2)

label_logit_cal <- paste0(
  "Logistic regression\n",
  "Brier = ", sprintf("%.3f", brier_logit)
)

label_svm_cal <- paste0(
  "SVM (RBF)\n",
  "Brier = ", sprintf("%.3f", brier_svm)
)

cal_logit$Model <- label_logit_cal
cal_svm$Model <- label_svm_cal

cal_df <- bind_rows(cal_logit, cal_svm)

p_cal <- ggplot(cal_df, aes(x = mean_pred, y = obs_rate, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey60") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(
    values = setNames(
      c("#0072B2", "#D55E00"),
      c(label_logit_cal, label_svm_cal)
    )
  ) +
  labs(
    title = "Calibration plot",
    x = "Mean predicted probability",
    y = "Observed event rate",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10, lineheight = 1.1)
  )

print(p_cal)

ggsave(
  "Calibration_logistic_vs_SVM_cv.pdf",
  plot = p_cal,
  width = 6,
  height = 5,
  device = cairo_pdf
)

## 11. Combined figure: ROC + calibration
p_combo <- gridExtra::grid.arrange(p_cv_roc, p_cal, ncol = 2)

ggsave(
  "Figure_logistic_vs_SVM_ROC_calibration.pdf",
  p_combo,
  width = 11,
  height = 5,
  device = cairo_pdf
)

## 12. SVM decision boundary in PCA space
## IMPORTANT:
## This is for visualization only, based on the full dataset.
## It is NOT a cross-validated boundary.

set.seed(123)

svm_model <- svm(
  as.formula(paste("group_fac ~", paste(analysis_vars, collapse = " + "))),
  data = df_ml,
  kernel = "radial",
  probability = TRUE,
  scale = TRUE
)

pred_class <- predict(svm_model, newdata = df_ml)

vars_pca <- df_ml[, analysis_vars[analysis_vars != "sex"], drop = FALSE]
vars_scaled <- scale(vars_pca)
pca <- prcomp(vars_scaled, center = FALSE, scale. = FALSE)

plot_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = df_ml$group_fac,
  pred = pred_class
)

plot_df$y <- ifelse(plot_df$group == "SAD", 1, 0)
plot_df$sv <- FALSE
plot_df$sv[svm_model$index] <- TRUE

## Approximate decision boundary in PCA space for visualization
glm2d <- glm(
  y ~ PC1 + PC2 + I(PC1^2) + I(PC2^2) + PC1:PC2,
  data = plot_df,
  family = binomial
)

x_seq <- seq(min(plot_df$PC1), max(plot_df$PC1), length.out = 300)
y_seq <- seq(min(plot_df$PC2), max(plot_df$PC2), length.out = 300)

grid <- expand.grid(
  PC1 = x_seq,
  PC2 = y_seq
)

grid$prob <- predict(glm2d, newdata = grid, type = "response")

auc_lab <- paste0(
  "SVM (RBF)\n",
  "CV AUC = ", sprintf("%.3f", auc_svm), "\n",
  "95% CI: ", sprintf("%.3f", ci_svm[1]), "–", sprintf("%.3f", ci_svm[3])
)

p_svm_pca <- ggplot(plot_df, aes(PC1, PC2)) +
  stat_ellipse(
    aes(color = group),
    level = 0.95,
    linewidth = 1.0,
    alpha = 0.3
  ) +
  geom_contour(
    data = grid,
    aes(x = PC1, y = PC2, z = prob),
    breaks = 0.5,
    color = "black",
    linewidth = 1
  ) +
  geom_point(
    aes(color = group, shape = pred),
    size = 2.6,
    alpha = 0.75
  ) +
  geom_point(
    data = subset(plot_df, sv),
    aes(PC1, PC2),
    shape = 21,
    size = 4.2,
    stroke = 0.7,
    fill = NA,
    color = "black"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.02,
    vjust = 1.1,
    label = auc_lab,
    size = 5
  ) +
  scale_color_manual(values = c("HC" = "#D55E00", "SAD" = "#0072B2")) +
  scale_shape_manual(values = c("HC" = 16, "SAD" = 17)) +
  labs(
    title = "SVM decision boundary in PCA space",
    subtitle = "Visualization based on the model fitted to the full dataset",
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
    color = "True group",
    shape = "Predicted"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_svm_pca)

ggsave(
  "svm_pca_boundary_plot.pdf",
  p_svm_pca,
  width = 10,
  height = 8,
  device = cairo_pdf
)

## 13. Summary output tables
summary_table <- data.frame(
  Model = c("Logistic regression", "SVM (RBF)"),
  AUC = c(round(auc_logit, 3), round(auc_svm, 3)),
  CI_lower = c(round(ci_logit[1], 3), round(ci_svm[1], 3)),
  CI_upper = c(round(ci_logit[3], 3), round(ci_svm[3], 3)),
  Brier = c(round(brier_logit, 3), round(brier_svm, 3))
)

summary_table$AUC_95CI <- paste0(
  sprintf("%.3f", summary_table$AUC),
  " (",
  sprintf("%.3f", summary_table$CI_lower),
  "–",
  sprintf("%.3f", summary_table$CI_upper),
  ")"
)

print(summary_table)

final_table <- data.frame(
  Model = c("Logistic regression", "SVM (RBF)"),
  AUC_95CI = c(
    paste0(sprintf("%.3f", auc_logit), " (", sprintf("%.3f", ci_logit[1]), "–", sprintf("%.3f", ci_logit[3]), ")"),
    paste0(sprintf("%.3f", auc_svm), " (", sprintf("%.3f", ci_svm[1]), "–", sprintf("%.3f", ci_svm[3]), ")")
  ),
  Brier = c(round(brier_logit, 3), round(brier_svm, 3)),
  DeLong_p = c(round(delong_test$p.value, 3), NA)
)

print(final_table)

svm_perf_table <- data.frame(
  Performance_Metric = c("AUC (95% CI)", "Accuracy", "Sensitivity", "Specificity"),
  Value = c(
    paste0(sprintf("%.3f", auc_svm), " (", sprintf("%.3f", ci_svm[1]), "–", sprintf("%.3f", ci_svm[3]), ")"),
    sprintf("%.3f", acc_svm),
    sprintf("%.3f", sens_svm),
    sprintf("%.3f", spec_svm)
  )
)

print(svm_perf_table)

openxlsx::write.xlsx(
  list(
    "Model_comparison" = summary_table,
    "Final_summary" = final_table,
    "SVM_performance" = svm_perf_table
  ),
  file = "Model_comparison_summary.xlsx",
  overwrite = TRUE
)

## 14. RBF kernel visualization
## 14. Pairwise RBF landscape visualization
## Two feature pairs:
## 1) SP + NRS
## 2) sleepiness + KSS_pre

make_rbf_plots <- function(data, vars, group_var = "group_fac", gamma_val, file_prefix) {
  
  ## keep required columns
  plot_df <- data[, c(group_var, vars)]
  plot_df <- na.omit(plot_df)
  
  ## standardize the two variables
  x_scaled <- scale(plot_df[, vars, drop = FALSE])
  
  points_df <- data.frame(
    x = x_scaled[, 1],
    y = x_scaled[, 2],
    group = plot_df[[group_var]]
  )
  
  ## RBF function
  rbf_kernel <- function(x, y, cx, cy, gamma) {
    exp(-gamma * ((x - cx)^2 + (y - cy)^2))
  }
  
  ## grid
  x_seq <- seq(min(points_df$x) - 1, max(points_df$x) + 1, length.out = 180)
  y_seq <- seq(min(points_df$y) - 1, max(points_df$y) + 1, length.out = 180)
  grid_df <- expand.grid(x = x_seq, y = y_seq)
  
  ## summed RBF bumps
  z_sum <- rep(0, nrow(grid_df))
  for (i in seq_len(nrow(points_df))) {
    z_sum <- z_sum + rbf_kernel(
      x = grid_df$x,
      y = grid_df$y,
      cx = points_df$x[i],
      cy = points_df$y[i],
      gamma = gamma_val
    )
  }
  
  grid_df$z <- z_sum / nrow(points_df)
  
  ## 2D plot
  p_2d <- ggplot(grid_df, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = z), color = "white", alpha = 0.7, linewidth = 0.35) +
    geom_point(
      data = points_df,
      aes(x = x, y = y, color = group),
      inherit.aes = FALSE,
      size = 2.3,
      alpha = 0.9
    ) +
    scale_color_manual(values = c("HC" = "#D55E00", "SAD" = "#0072B2")) +
    labs(
      title = paste0("RBF kernel landscape: ", vars[1], " + ", vars[2]),
      subtitle = "Multiple Gaussian bumps centered at observed data points",
      x = paste0(vars[1], " (scaled)"),
      y = paste0(vars[2], " (scaled)"),
      fill = "Kernel\nintensity",
      color = "Group"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(color = "black")
    )
  
  print(p_2d)
  
  ggsave(
    paste0(file_prefix, "_2D.pdf"),
    plot = p_2d,
    width = 7,
    height = 5.5,
    device = cairo_pdf
  )
  
  ## 3D plot
  z_mat <- matrix(grid_df$z, nrow = length(x_seq), ncol = length(y_seq))
  
  pdf(paste0(file_prefix, "_3D.pdf"), width = 7, height = 6)
  persp(
    x_seq, y_seq, z_mat,
    theta = 35,
    phi = 28,
    expand = 0.7,
    col = "steelblue",
    shade = 0.35,
    border = NA,
    ticktype = "detailed",
    xlab = paste0(vars[1], " (scaled)"),
    ylab = paste0(vars[2], " (scaled)"),
    zlab = "Kernel intensity",
    main = paste0("RBF kernel landscape (3D): ", vars[1], " + ", vars[2])
  )
  dev.off()
}

## use the same sigma estimated earlier
gamma_val <- sigma_est

## Plot 1: SP + NRS
make_rbf_plots(
  data = df_ml,
  vars = c("SP", "NRS"),
  gamma_val = gamma_val,
  file_prefix = "RBF_landscape_SP_NRS"
)

## Plot 2: sleepiness + KSS_pre
make_rbf_plots(
  data = df_ml,
  vars = c("sleepiness", "KSS_pre"),
  gamma_val = gamma_val,
  file_prefix = "RBF_landscape_sleepiness_KSSpre"
)


## 15. Predicted probability distribution plot for SVM

p_svm_prob <- ggplot(best_pred, aes(x = obs, y = SAD, fill = obs)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = obs), width = 0.08, height = 0, size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = c("HC" = "#D55E00", "SAD" = "#0072B2")) +
  scale_color_manual(values = c("HC" = "#D55E00", "SAD" = "#0072B2")) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") +
  labs(
    title = "Distribution of SVM-predicted probabilities",
    subtitle = "Predicted probability of SAD by true group",
    x = "True group",
    y = "Predicted probability of SAD"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

print(p_svm_prob)

ggsave(
  "SVM_predicted_probability_distribution.pdf",
  plot = p_svm_prob,
  width = 5.5,
  height = 5,
  device = cairo_pdf
)
