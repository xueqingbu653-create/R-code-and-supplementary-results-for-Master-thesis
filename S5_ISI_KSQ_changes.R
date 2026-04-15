## KSQ_pre_follow up ####

df <- KSQ_followup

df[df == "Missing"] <- NA
df[,3:ncol(df)] <- lapply(df[,3:ncol(df)], function(x) as.numeric(as.character(x)))
df$i[is.na(df$i)] <- median(df$i, na.rm = TRUE)
colSums(is.na(df))
write.csv(df, "KSQ_followup_completed.csv", row.names = FALSE)

#### KSQ pre-post comparison ####
df <- KSQ_pre_post

library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(emmeans)
library(ggplot2)
library(effectsize)
library(performance)

long_data <- df %>%
  pivot_longer(
    cols = c(SP_pre, SP_post,
             NRS_pre, NRS_post,
             apnea_pre, apnea_post,
             sleepiness_pre, sleepiness_post,
             KSQ_pre, KSQ_post),
    names_to = c("variable", "time"),
    names_sep = "_",
    values_to = "score"
  )

long_data$time <- factor(long_data$time, levels = c("pre", "post"))
long_data$variable <- factor(
  long_data$variable,
  levels = c("SP", "NRS", "apnea", "sleepiness", "KSQ")
)

## LMM
model_all <- lmer(score ~ time * variable + (1 | subject), data = long_data)

summary(model_all)
anova(model_all)

# Overall model R2
r2_all <- performance::r2_nakagawa(model_all)
print(r2_all)

r2_all_df <- data.frame(
  R2_marginal = round(r2_all$R2_marginal, 3),
  R2_conditional = round(r2_all$R2_conditional, 3)
)

write.csv(r2_all_df, "LMM_prepost_all_variables_R2.csv", row.names = FALSE)

# Fixed effects table for the overall model
results_all <- broom.mixed::tidy(model_all, effects = "fixed", conf.int = TRUE) %>%
  mutate(across(c(estimate, std.error, conf.low, conf.high, p.value), ~ round(., 3)))

print(results_all)
write.csv(results_all, "LMM_prepost_all_variables_results.csv", row.names = FALSE)


## Separate LMM for each variable

run_lmm <- function(data, var_name) {
  
  # Reshape data from wide to long format
  d <- data %>%
    select(subject, all_of(paste0(var_name, "_pre")), all_of(paste0(var_name, "_post"))) %>%
    pivot_longer(
      cols = -subject,
      names_to = "time",
      values_to = "score"
    ) %>%
    drop_na(score)
  
  # Recode time as pre/post
  d$time <- factor(
    d$time,
    levels = c(paste0(var_name, "_pre"), paste0(var_name, "_post")),
    labels = c("pre", "post")
  )
  
  # Fit separate LMM for each variable
  m <- lmer(score ~ time + (1 | subject), data = d)
  
  # Fixed effects table
  res <- broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>%
    mutate(variable = var_name) %>%
    select(variable, term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(across(c(estimate, std.error, conf.low, conf.high, p.value), ~ round(., 3)))
  
  # Estimated marginal means
  emm <- emmeans(m, ~ time)
  emm_df <- as.data.frame(emm) %>%
    mutate(variable = var_name) %>%
    relocate(variable) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  # Pre-post contrast
  contrast_df <- as.data.frame(pairs(emm, reverse = TRUE)) %>%
    mutate(variable = var_name) %>%
    relocate(variable) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  # LMM-based Cohen's d
  eff_df <- as.data.frame(
    eff_size(
      emm,
      sigma = sigma(m),
      edf = df.residual(m),
      method = "revpairwise"
    )
  ) %>%
    mutate(variable = var_name) %>%
    relocate(variable) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  # Paired-samples Cohen's dz
  wide_d <- data %>%
    select(
      subject,
      pre = all_of(paste0(var_name, "_pre")),
      post = all_of(paste0(var_name, "_post"))
    ) %>%
    drop_na(pre, post)
  
  paired_d_df <- as.data.frame(
    effectsize::cohens_d(
      wide_d$post,
      wide_d$pre,
      paired = TRUE,
      ci = 0.95
    )
  ) %>%
    mutate(variable = var_name) %>%
    relocate(variable) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  # Model R2
  r2_obj <- performance::r2_nakagawa(m)
  r2_df <- data.frame(
    variable = var_name,
    R2_marginal = round(r2_obj$R2_marginal, 3),
    R2_conditional = round(r2_obj$R2_conditional, 3)
  )
  
  list(
    model = m,
    results = res,
    emmeans = emm_df,
    contrast = contrast_df,
    effect_size = eff_df,
    paired_d = paired_d_df,
    r2 = r2_df,
    data = d
  )
}

# Variables to analyze
vars <- c("SP", "NRS", "apnea", "sleepiness", "KSQ")

# Run separate LMMs
lmm_list <- lapply(vars, function(v) run_lmm(df, v))
names(lmm_list) <- vars


## Keep the original tables

# Fixed effects table for each variable
results_each <- bind_rows(lapply(lmm_list, function(x) x$results))
print(results_each)
write.csv(results_each, "LMM_prepost_each_variable_results.csv", row.names = FALSE)

# Estimated marginal means for each variable
emm_all <- bind_rows(lapply(lmm_list, function(x) x$emmeans))
print(emm_all)
write.csv(emm_all, "LMM_prepost_emmeans.csv", row.names = FALSE)


## Additional output tables
# Collect all results (no separate output)
contrast_all <- bind_rows(lapply(lmm_list, function(x) x$contrast))
eff_all <- bind_rows(lapply(lmm_list, function(x) x$effect_size))
paired_d_all <- bind_rows(lapply(lmm_list, function(x) x$paired_d))
r2_each <- bind_rows(lapply(lmm_list, function(x) x$r2))

## Final summary table
summary_table <- contrast_all %>%
  left_join(
    eff_all %>%
      select(variable, effect.size, lower.CL, upper.CL),
    by = "variable"
  ) %>%
  left_join(
    paired_d_all %>%
      select(variable, Cohens_d, CI_low, CI_high),
    by = "variable"
  ) %>%
  left_join(
    r2_each,
    by = "variable"
  ) %>%
  rename(
    mean_diff = estimate,
    t_value = t.ratio,
    p_value = p.value,
    cohens_d_lmm = effect.size,
    cohens_d_lmm_low = lower.CL,
    cohens_d_lmm_high = upper.CL,
    cohens_dz = Cohens_d,
    cohens_dz_low = CI_low,
    cohens_dz_high = CI_high
  ) %>%
  select(
    variable,
    contrast,
    mean_diff,
    SE,
    df,
    t_value,
    p_value,
    cohens_d_lmm,
    cohens_d_lmm_low,
    cohens_d_lmm_high,
    cohens_dz,
    cohens_dz_low,
    cohens_dz_high,
    R2_marginal,
    R2_conditional
  )

print(summary_table)

write.csv(
  summary_table,
  "LMM_prepost_summary_table.csv",
  row.names = FALSE
)
## Plot individual trajectories + mean change
p <- ggplot(long_data, aes(x = time, y = score, group = subject)) +
  geom_line(color = "grey70", linewidth = 0.4, alpha = 0.5) +
  geom_point(color = "grey70", size = 1.2, alpha = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "#2C7FB8", linewidth = 1.3) +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               color = "#2C7FB8", size = 3) +
  stat_summary(aes(group = 1), fun.data = mean_cl_boot, geom = "errorbar",
               color = "#2C7FB8", width = 0.08, linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(
    title = "Pre-post changes in sleep variables",
    x = "Time",
    y = "Score"
  ) +
  scale_x_discrete(labels = c("Pre", "Post")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p)

ggsave(
  "LMM_prepost_sleep_variables_plot.pdf",
  plot = p,
  width = 10,
  height = 6
)

ggsave(
  "LMM_prepost_sleep_variables_plot.tiff",
  plot = p,
  width = 10,
  height = 6,
  dpi = 600,
  compression = "lzw"
)

#### ISI changes ####
df <- sample1_SAD

library(tidyverse)
library(effectsize)

isi_long <- df %>%
  pivot_longer(cols = c(isi1, isi18, isi112),
               names_to = "time",
               values_to = "isi")

isi_long$time <- factor(isi_long$time,
                        levels = c("isi1","isi18","isi112"),
                        labels = c("pre","post","followup"))

library(lme4)
library(lmerTest)

model <- lmer(isi ~ time + age + sex + (1 | subject), data = isi_long)

summary(model)
anova(model)

library(broom.mixed)
library(dplyr)
table_lmm <- tidy(model, effects = "fixed", conf.int = TRUE) %>%
  select(term, estimate, std.error, conf.low, conf.high, p.value)

table_lmm
write.csv(table_lmm, "lmm_isi_results.csv", row.names = FALSE)

## Post hoc comparisons and effect sizes
library(emmeans)

# Estimated marginal means
emm_isi <- emmeans(model, ~ time)
emm_isi_df <- as.data.frame(emm_isi) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

print(emm_isi_df)
write.csv(emm_isi_df, "isi_emmeans.csv", row.names = FALSE)

# Pairwise comparisons between time points
pairs_isi <- pairs(emm_isi, adjust = "bonferroni")
pairs_isi_df <- as.data.frame(pairs_isi) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

print(pairs_isi_df)
write.csv(pairs_isi_df, "isi_pairwise_comparisons.csv", row.names = FALSE)

# Cohen's d for pairwise comparisons
eff_isi <- as.data.frame(
  eff_size(
    emm_isi,
    sigma = sigma(model),
    edf = df.residual(model),
    method = "pairwise"
  )
) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

print(eff_isi)
write.csv(eff_isi, "isi_effect_sizes_cohens_d.csv", row.names = FALSE)

# Final summary table
isi_summary_table <- pairs_isi_df %>%
  left_join(
    eff_isi %>%
      select(contrast, effect.size, lower.CL, upper.CL),
    by = "contrast"
  ) %>%
  rename(
    mean_diff = estimate,
    t_value = t.ratio,
    p_value = p.value,
    cohens_d = effect.size,
    cohens_d_low = lower.CL,
    cohens_d_high = upper.CL
  ) %>%
  select(
    contrast,
    mean_diff,
    SE,
    df,
    t_value,
    p_value,
    cohens_d,
    cohens_d_low,
    cohens_d_high
  )

print(isi_summary_table)
write.csv(isi_summary_table, "isi_summary_table.csv", row.names = FALSE)

## visualization
library(ggplot2)

p_isi <- ggplot(isi_long, aes(x = time, y = isi, group = subject)) +
  
  geom_line(color = "grey70", linewidth = 0.4, alpha = 0.5) +
  geom_point(color = "grey70", size = 1.2, alpha = 0.5) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    color = "#2C7FB8",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    color = "#2C7FB8",
    size = 3
  ) +
  
  # 95% CI
  stat_summary(
    aes(group = 1),
    fun.data = mean_cl_boot,
    geom = "errorbar",
    color = "#2C7FB8",
    width = 0.1,
    linewidth = 0.9
  ) +
  
  labs(
    x = "Time",
    y = "ISI score",
    title = "Individual and mean changes in ISI over time"
  ) +
  
  theme_classic(base_size = 14)

ggsave(
  "isi_change_plot.png",
  plot = p_isi,
  width = 6,
  height = 4,
  dpi = 300
)


