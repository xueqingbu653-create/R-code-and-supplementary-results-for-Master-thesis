
# FULL ANALYSIS PIPELINE
# Study 1: Sleep subtype identification and validation
# Study 2: Clinical utility of sleep subtypes
# Study 3: Neural correlates of sleep subtypes

df <- KSQ_combined

library(dplyr)
library(tidyr)
library(ggplot2)
library(mclust)
library(lme4)
library(emmeans)
library(purrr)


# Sleep variables used for subtype identification
sleep_vars <- c("SP", "NRS", "apnea", "sleepiness")

# ROI variables for Study 3
roi_vars <- c(
  "Amygdala_02","Caudata_02","Cingulate_Ant_02","Cingulate_Post_02",
  "Frontal_Mid_02","Frontal_Sup_02",
  "Frontal_Sup_Medial_02","Hypothalamus_02","Insula_02",
  "LC_02","N_Acc_02","OFCmed_02","Pallidum_02",
  "Precuneus_02","Putamen_02","Raphe_D_02","Raphe_M_02",
  "Thalamus_02","VTA_02"
)


#### STUDY 1: Sleep subtype identification (full sample) ####

# Create clustering dataset
# Only retain variables needed for clustering so that HC cases
# are not dropped because of missing LSAS_post
cluster_dat <- df %>%
  dplyr::mutate(
    row_id = dplyr::row_number(),
    group_label = factor(group, levels = c(0, 1), labels = c("HC", "SAD"))
  ) %>%
  dplyr::select(
    row_id,
    subject,
    sample,
    group,
    group_label,
    dplyr::all_of(sleep_vars)
  ) %>%
  tidyr::drop_na(dplyr::all_of(sleep_vars))

cat("N included in Study 1 clustering =", nrow(cluster_dat), "\n")
print(table(cluster_dat$group_label))

# Standardize clustering variables
X_full <- cluster_dat %>%
  dplyr::select(dplyr::all_of(sleep_vars))

Xz_full <- scale(X_full)

print(summary(X_full))
print(cor(X_full, use = "pairwise.complete.obs"))


# Gaussian mixture model clustering

set.seed(123)

gmm_full <- mclust::Mclust(
  data = Xz_full,
  G = 1:5
)

summary(gmm_full)

cat("Best model =", gmm_full$modelName, "\n")
cat("Best number of classes =", gmm_full$G, "\n")

# Save BIC plot
bic_mat <- gmm_full$BIC

bic_long <- as.data.frame(as.table(bic_mat)) %>%
  rename(
    G = Var1,
    model = Var2,
    BIC = Freq
  ) %>%
  mutate(
    G = as.numeric(as.character(G))
  ) %>%
  filter(!is.na(BIC))

p_bic <- ggplot(bic_long, aes(x = G, y = BIC, color = model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Number of components",
    y = "BIC",
    color = "Model"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "study_gmm_BIC_plot.png",
  p_bic,
  width = 8,
  height = 6,
  dpi = 300
)


# Extract subtype assignment

cluster_dat$subtype <- factor(gmm_full$classification)

# Posterior probabilities
post_probs <- as.data.frame(gmm_full$z)
colnames(post_probs) <- paste0("prob_class_", seq_len(ncol(post_probs)))

cluster_dat <- dplyr::bind_cols(cluster_dat, post_probs)

# Maximum posterior probability
cluster_dat$max_post_prob <- apply(gmm_full$z, 1, max)

# Entropy-like individual uncertainty
cluster_dat$entropy <- -rowSums(gmm_full$z * log(gmm_full$z + 1e-10))

# Save subtype assignment table
write.csv(
  cluster_dat,
  "study1_fullsample_subtype_assignments.csv",
  row.names = FALSE
)

# Subtype size and classification quality

subtype_size <- cluster_dat %>%
  dplyr::count(subtype) %>%
  dplyr::mutate(percent = round(100 * n / sum(n), 2))

print(subtype_size)

write.csv(
  subtype_size,
  "study1_subtype_sizes.csv",
  row.names = FALSE
)

classification_summary <- cluster_dat %>%
  dplyr::summarise(
    mean_max_post_prob = mean(max_post_prob),
    sd_max_post_prob = sd(max_post_prob),
    mean_entropy = mean(entropy),
    sd_entropy = sd(entropy)
  )

print(classification_summary)

write.csv(
  classification_summary,
  "study1_classification_quality.csv",
  row.names = FALSE
)


# Subtype profiles
# Raw-value subtype profile
profile_raw <- cluster_dat %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sleep_vars), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sleep_vars),
    names_to = "dimension",
    values_to = "mean_raw"
  )

write.csv(
  profile_raw,
  "study1_subtype_profiles_raw.csv",
  row.names = FALSE
)

# Z-score subtype profile for plotting
df_z <- as.data.frame(Xz_full)
df_z$subtype <- cluster_dat$subtype

profile_z <- df_z %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sleep_vars), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sleep_vars),
    names_to = "dimension",
    values_to = "z_mean"
  )

profile_z$dimension <- factor(
  profile_z$dimension,
  levels = sleep_vars,
  labels = c("SP", "NRS", "Apnea", "Sleepiness")
)

p_profile <- ggplot(
  profile_z,
  aes(x = dimension, y = z_mean, group = subtype, color = subtype)
) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    title = "Study 1: Sleep subtype profiles (z-scores)",
    x = NULL,
    y = "Standardized mean (z-score)",
    color = "Subtype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_profile)

ggsave(
  "study1_subtype_profiles.png",
  p_profile,
  width = 7,
  height = 5.5
)

# Rename subtype labels after profile inspection
# Adjust these labels if your profile interpretation changes
cluster_dat$subtype_label <- dplyr::recode(
  cluster_dat$subtype,
  `1` = "Apnea-dominant",
  `2` = "Insomnia-sleepiness",
  `3` = "Low disturbance"
)

# Bootstrap stability
set.seed(123)

best_G <- gmm_full$G
best_model <- gmm_full$modelName
orig_class <- gmm_full$classification

boot_gmm_stability <- function(X, B, G, modelName, orig_class) {
  n <- nrow(X)
  ari_values <- numeric(B)
  
  for (b in 1:B) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]
    
    fit_b <- tryCatch(
      mclust::Mclust(Xb, G = G, modelNames = modelName),
      error = function(e) NULL
    )
    
    if (!is.null(fit_b)) {
      ari_values[b] <- mclust::adjustedRandIndex(
        orig_class[idx],
        fit_b$classification
      )
    } else {
      ari_values[b] <- NA
    }
  }
  
  return(ari_values)
}

boot_ari <- boot_gmm_stability(
  X = Xz_full,
  B = 500,
  G = best_G,
  modelName = best_model,
  orig_class = orig_class
)

stability_summary <- data.frame(
  mean_ARI = mean(boot_ari, na.rm = TRUE),
  sd_ARI = sd(boot_ari, na.rm = TRUE),
  median_ARI = median(boot_ari, na.rm = TRUE),
  min_ARI = min(boot_ari, na.rm = TRUE),
  max_ARI = max(boot_ari, na.rm = TRUE)
)

print(stability_summary)

write.csv(
  stability_summary,
  "study1_bootstrap_stability_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(ARI = boot_ari),
  "study1_bootstrap_ARI_values.csv",
  row.names = FALSE
)

p_boot <- ggplot(
  data.frame(ARI = boot_ari),
  aes(x = ARI)
) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(
    title = "Bootstrap stability of GMM subtypes",
    x = "Adjusted Rand Index (ARI)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_boot)

ggsave(
  "study1_bootstrap_ARI_histogram.png",
  p_boot,
  width = 7,
  height = 5.5
)


# Subtype distribution by group and by sample
# SAD vs HC
subtype_by_group <- table(cluster_dat$group_label, cluster_dat$subtype_label)
print(subtype_by_group)

chisq_res <- chisq.test(subtype_by_group)
print(chisq_res)
# Post-hoc: standardized residuals
std_res <- chisq_res$stdres
print(round(std_res, 3))

chi_sq <- unname(chisq_res$statistic)
df_chi <- unname(chisq_res$parameter)
p_chi <- chisq_res$p.value

cat(sprintf("Chi-square = %.3f, df = %d, p = %.3f\n", chi_sq, df_chi, p_chi))

pct_by_group <- prop.table(subtype_by_group, margin = 1) * 100

subtype_table_pct <- matrix(
  paste0(
    subtype_by_group,
    " (",
    sprintf("%.1f", pct_by_group),
    "%)"
  ),
  nrow = nrow(subtype_by_group),
  dimnames = dimnames(subtype_by_group)
)

subtype_table_pct <- as.data.frame.matrix(subtype_table_pct)
subtype_table_pct <- tibble::rownames_to_column(subtype_table_pct, var = "Group")

print(subtype_table_pct)

write.csv(
  subtype_table_pct,
  "study1_subtype_by_group_with_percent.csv",
  row.names = FALSE
)

chisq_summary <- data.frame(
  chisq = round(chi_sq, 3),
  df = df_chi,
  p = ifelse(p_chi < .001, "< .001", sub("^0", "", sprintf("%.3f", p_chi)))
)

write.csv(
  chisq_summary,
  "study1_subtype_by_group_chisq.csv",
  row.names = FALSE
)

# Distribution plot
plot_dat_group <- cluster_dat %>%
  dplyr::count(group_label, subtype) %>%
  dplyr::group_by(group_label) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

p_group_dist <- ggplot(
  plot_dat_group,
  aes(x = group_label, y = prop, fill = subtype)
) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = c(
      "1" = "#E69F00",
      "2" = "#009E73",
      "3" = "#0072B2"
    )
  ) +
  labs(
    title = "Distribution of sleep subtypes in SAD and HC",
    x = "Group",
    y = "Proportion",
    fill = "Subtype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_group_dist)

ggsave(
  "study1_subtype_distribution_SAD_HC.png",
  p_group_dist,
  width = 7,
  height = 5.5
)

# sample1 vs sample2
subtype_by_sample <- table(cluster_dat$sample, cluster_dat$subtype)
print(subtype_by_sample)
print(chisq.test(subtype_by_sample))

write.csv(
  as.data.frame.matrix(subtype_by_sample),
  "study1_subtype_by_sample_crosstab.csv"
)


# Merge subtype back to main dataset
analysis_dat <- df %>%
  dplyr::mutate(
    row_id = dplyr::row_number(),
    group_label = factor(group, levels = c(0, 1), labels = c("HC", "SAD"))
  ) %>%
  dplyr::left_join(
    cluster_dat %>%
      dplyr::select(row_id, subtype, subtype_label),
    by = "row_id"
  ) %>%
  dplyr::mutate(
    treatment = dplyr::case_when(
      sample == 1 & group == 1 ~ "ICBT",
      sample == 2 & group == 1 ~ "SSRI",
      TRUE ~ NA_character_
    ),
    treatment = factor(treatment, levels = c("ICBT", "SSRI"))
  )

#### sample characteristics

table1_dat <- df %>%
  mutate(
    group_label = factor(group, levels = c(0,1), labels = c("HC","SAD"))
  )

# mean (SD) functions
fmt_mean_sd1 <- function(x){
  sprintf("%.1f (%.1f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
}

fmt_mean_sd3 <- function(x){
  sprintf("%.3f (%.3f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
}

# n (%)
fmt_n_pct <- function(x){
  n <- sum(x, na.rm=TRUE)
  pct <- 100 * mean(x, na.rm=TRUE)
  sprintf("%d (%.1f%%)", n, pct)
}

# APA p value format
fmt_p <- function(p){
  if (is.na(p)) return("")
  if (p < .001) return("< .001")
  sub("^0","", sprintf("%.3f", p))
}

table1 <- bind_rows(
  
  data.frame(
    Variable = "N",
    Full = as.character(nrow(table1_dat)),
    HC = as.character(sum(table1_dat$group_label=="HC")),
    SAD = as.character(sum(table1_dat$group_label=="SAD")),
    p = ""
  ),
  
  data.frame(
    Variable = "Age, mean (SD)",
    Full = fmt_mean_sd1(table1_dat$age),
    HC = fmt_mean_sd1(table1_dat$age[table1_dat$group_label=="HC"]),
    SAD = fmt_mean_sd1(table1_dat$age[table1_dat$group_label=="SAD"]),
    p = fmt_p(t.test(age ~ group_label, data=table1_dat)$p.value)
  ),
  
  data.frame(
    Variable = "Male, n (%)",
    Full = fmt_n_pct(table1_dat$sex==1),
    HC = fmt_n_pct(table1_dat$sex[table1_dat$group_label=="HC"]==1),
    SAD = fmt_n_pct(table1_dat$sex[table1_dat$group_label=="SAD"]==1),
    p = fmt_p(chisq.test(table(table1_dat$group_label, table1_dat$sex))$p.value)
  ),
  
  data.frame(
    Variable = "SP, mean (SD)",
    Full = fmt_mean_sd3(table1_dat$SP),
    HC = fmt_mean_sd3(table1_dat$SP[table1_dat$group_label=="HC"]),
    SAD = fmt_mean_sd3(table1_dat$SP[table1_dat$group_label=="SAD"]),
    p = fmt_p(t.test(SP ~ group_label, data=table1_dat)$p.value)
  ),
  
  data.frame(
    Variable = "NRS, mean (SD)",
    Full = fmt_mean_sd3(table1_dat$NRS),
    HC = fmt_mean_sd3(table1_dat$NRS[table1_dat$group_label=="HC"]),
    SAD = fmt_mean_sd3(table1_dat$NRS[table1_dat$group_label=="SAD"]),
    p = fmt_p(t.test(NRS ~ group_label, data=table1_dat)$p.value)
  ),
  
  data.frame(
    Variable = "Apnea, mean (SD)",
    Full = fmt_mean_sd3(table1_dat$apnea),
    HC = fmt_mean_sd3(table1_dat$apnea[table1_dat$group_label=="HC"]),
    SAD = fmt_mean_sd3(table1_dat$apnea[table1_dat$group_label=="SAD"]),
    p = fmt_p(t.test(apnea ~ group_label, data=table1_dat)$p.value)
  ),
  
  data.frame(
    Variable = "Sleepiness, mean (SD)",
    Full = fmt_mean_sd3(table1_dat$sleepiness),
    HC = fmt_mean_sd3(table1_dat$sleepiness[table1_dat$group_label=="HC"]),
    SAD = fmt_mean_sd3(table1_dat$sleepiness[table1_dat$group_label=="SAD"]),
    p = fmt_p(t.test(sleepiness ~ group_label, data=table1_dat)$p.value)
  )
)

print(table1)

write.csv(
  table1,
  "study1_table1_sample_characteristics.csv",
  row.names = FALSE
)


#### STUDY 2: Clinical utility (SAD only) ####
sad_dat <- analysis_dat %>%
  dplyr::filter(group == 1)

cat("N in SAD clinical dataset =", nrow(sad_dat), "\n")
print(table(sad_dat$treatment, useNA = "ifany"))
print(table(sad_dat$subtype, useNA = "ifany"))


# Baseline severity
# This tests whether sleep subtypes differ in baseline social anxiety severity
lm_baseline <- lm(
  LSAS_pre ~ subtype + age + sex,
  data = sad_dat
)

print(summary(lm_baseline))
print(anova(lm_baseline))

baseline_emm <- emmeans::emmeans(lm_baseline, pairwise ~ subtype)
print(baseline_emm)

# Treatment improvement (change score)
# This tests whether sleep subtypes differ in overall symptom improvement
sad_dat <- sad_dat %>%
  dplyr::mutate(delta_LSAS = LSAS_pre - LSAS_post)

lm_delta <- lm(
  delta_LSAS ~ subtype + age + sex + treatment,
  data = sad_dat
)

print(summary(lm_delta))
print(anova(lm_delta))

# test subtype-by-treatment interaction
lm_delta_interaction <- lm(
  delta_LSAS ~ subtype * treatment + age + sex,
  data = sad_dat
)

print(summary(lm_delta_interaction))
print(anova(lm_delta_interaction))

## Longitudinal symptom trajectory
# This tests whether symptom change over time differs by subtype
library(lmerTest)
sad_long <- sad_dat %>%
  dplyr::select(subject, subtype, subtype_label, treatment, age, sex, LSAS_pre, LSAS_post) %>%
  tidyr::pivot_longer(
    cols = c(LSAS_pre, LSAS_post),
    names_to = "time",
    values_to = "LSAS"
  ) %>%
  dplyr::mutate(
    time = factor(time,
                  levels = c("LSAS_pre", "LSAS_post"),
                  labels = c("pre", "post"))
  )

lmm_subtype <- lmerTest::lmer(
  LSAS ~ time * subtype + age + sex + (1 | subject),
  data = sad_long
)

print(summary(lmm_subtype))
print(anova(lmm_subtype))

## test whether subtype effects differ by treatment modality
lmm_subtype_treatment <- lmerTest::lmer(
  LSAS ~ time * subtype * treatment + age + sex + (1 | subject),
  data = sad_long
)

print(summary(lmm_subtype_treatment))
print(anova(lmm_subtype_treatment))


## Plot treatment trajectories by subtype
library(emmeans)

emm_traj <- emmeans(
  lmm_subtype,
  ~ time * subtype
)

plot_dat <- as.data.frame(emm_traj)

plot_dat$subtype_label <- factor(
  plot_dat$subtype,
  levels = c(1,2,3),
  labels = c("Apnea-dominant", "Insomnia-sleepiness", "Low disturbance")
)

p_traj <- ggplot(
  plot_dat,
  aes(x = time, y = emmean, group = subtype_label, color = subtype_label)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.08
  ) +
  scale_color_manual(
  values = c(
    "Apnea-dominant" = "#E69F00",
    "Insomnia-sleepiness" = "#009E73",
    "Low disturbance" = "#0072B2"
  ),
  breaks = c(
    "Apnea-dominant",
    "Insomnia-sleepiness",
    "Low disturbance"
  )
) +
  labs(
    title = "LSAS change over time by sleep subtype",
    x = "Time",
    y = "Estimated LSAS",
    color = "Sleep subtype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_traj)

ggsave(
  "study2_LSAS_trajectory_by_subtype.png",
  p_traj,
  width = 7,
  height = 5.5
)


#### STUDY 3: Neural correlates of sleep subtypes ####
# Imaging available only in sample1

img_dat <- analysis_dat %>%
  dplyr::filter(sample == 1)

img_sad <- img_dat %>%
  dplyr::filter(group == 1)

cat("N in sample1 imaging dataset =", nrow(img_dat), "\n")
cat("N in sample1 SAD imaging dataset =", nrow(img_sad), "\n")

# Multivariate overall test (MANOVA)
# Restrict to complete ROI cases for MANOVA
img_sad_complete <- img_sad %>%
  tidyr::drop_na(dplyr::all_of(roi_vars), subtype, age, sex)

roi_formula_text <- paste(roi_vars, collapse = ", ")
manova_formula <- as.formula(
  paste0("cbind(", roi_formula_text, ") ~ subtype + age + sex")
)

fit_manova <- manova(manova_formula, data = img_sad_complete)
print(summary(fit_manova, test = "Pillai"))

# ROI-wise subtype comparisons
roi_results <- purrr::map_dfr(roi_vars, function(roi) {
  
  roi_dat <- img_sad %>%
    dplyr::select(subtype, age, sex, dplyr::all_of(roi)) %>%
    tidyr::drop_na()
  
  formula <- stats::as.formula(
    paste0(roi, " ~ subtype + age + sex")
  )
  
  fit <- lm(formula, data = roi_dat)
  aov_tab <- anova(fit)
  
  data.frame(
    ROI = roi,
    N = nrow(roi_dat),
    df1 = aov_tab["subtype", "Df"],
    df2 = aov_tab["Residuals", "Df"],
    F_value = aov_tab["subtype", "F value"],
    p_value = aov_tab["subtype", "Pr(>F)"]
  )
})

roi_results$p_fdr <- p.adjust(roi_results$p_value, method = "fdr")
roi_results$p_bonf <- p.adjust(roi_results$p_value, method = "bonferroni")
roi_results$p_fwe <- p.adjust(roi_results$p_value, method = "holm")

print(roi_results)

write.csv(
  roi_results,
  "study3_roi_subtype_results.csv",
  row.names = FALSE
)


# Optional post hoc tests for significant ROIs
# Example: run post hoc tests for ROIs with nominal p < .05
sig_roi_nominal <- roi_results %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::pull(ROI)

posthoc_list <- lapply(sig_roi_nominal, function(roi) {
  roi_dat <- img_sad %>%
    dplyr::select(subtype_label, age, sex, dplyr::all_of(roi)) %>%
    tidyr::drop_na()
  
  formula <- stats::as.formula(
    paste0(roi, " ~ subtype_label + age + sex")
  )
  
  fit <- lm(formula, data = roi_dat)
  emm <- emmeans::emmeans(fit, pairwise ~ subtype_label)
  
  list(
    ROI = roi,
    model = summary(fit),
    posthoc = emm
  )
})

# Print post hoc results if any nominally significant ROI exists
if (length(posthoc_list) > 0) {
  print(posthoc_list)
}


# PCA-based multivariate brain pattern analysis

brain_dat <- img_sad %>%
  dplyr::select(subject, subtype, subtype_label, age, sex, dplyr::all_of(roi_vars)) %>%
  tidyr::drop_na()

brain_mat <- brain_dat %>%
  dplyr::select(dplyr::all_of(roi_vars))

pca_brain <- prcomp(brain_mat, scale. = TRUE)

print(summary(pca_brain))

pca_scores <- as.data.frame(pca_brain$x)
pca_scores$subtype <- brain_dat$subtype
pca_scores$subtype_label <- brain_dat$subtype_label
pca_scores$age <- brain_dat$age
pca_scores$sex <- brain_dat$sex

# Test subtype effect on first principal component
lm_pc1 <- lm(
  PC1 ~ subtype + age + sex,
  data = pca_scores
)

print(summary(lm_pc1))
print(anova(lm_pc1))

# Plot PCA brain pattern
p_pca <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = subtype_label)
) +
  geom_point(size = 2.8, alpha = 0.85) +
  stat_ellipse(level = 0.68, linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "Insomnia-sleepiness" = "#E69F00",
      "Apnea-dominant" = "#009E73",
      "Low disturbance" = "#0072B2"
    )
  ) +
  labs(
    title = "PCA of ROI SDBOLD by sleep subtype",
    x = "PC1",
    y = "PC2",
    color = "Sleep subtype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_pca)

ggsave(
  "study3_PCA_brain_pattern_by_subtype.png",
  p_pca,
  width = 7,
  height = 5.5
)


# Optional case-control exploratory analysis in sample1
# This exploratory analysis tests whether subtype and case-control group
# are jointly associated with ROI variability in sample1
img_dat_complete <- img_dat %>%
  tidyr::drop_na(subtype, age, sex)

roi_results_group <- purrr::map_dfr(roi_vars, function(roi) {
  
  roi_dat <- img_dat_complete %>%
    dplyr::select(group_label, subtype, age, sex, dplyr::all_of(roi)) %>%
    tidyr::drop_na()
  
  formula <- stats::as.formula(
    paste0(roi, " ~ group_label + subtype + age + sex")
  )
  
  fit <- lm(formula, data = roi_dat)
  aov_tab <- anova(fit)
  
  data.frame(
    ROI = roi,
    F_group = ifelse("group_label" %in% rownames(aov_tab), aov_tab["group_label", "F value"], NA),
    p_group = ifelse("group_label" %in% rownames(aov_tab), aov_tab["group_label", "Pr(>F)"], NA),
    F_subtype = ifelse("subtype" %in% rownames(aov_tab), aov_tab["subtype", "F value"], NA),
    p_subtype = ifelse("subtype" %in% rownames(aov_tab), aov_tab["subtype", "Pr(>F)"], NA)
  )
})

roi_results_group$p_group_fdr <- p.adjust(roi_results_group$p_group, method = "fdr")
roi_results_group$p_subtype_fdr <- p.adjust(roi_results_group$p_subtype, method = "fdr")

write.csv(
  roi_results_group,
  "study3_roi_group_subtype_exploratory_results.csv",
  row.names = FALSE
)

print(roi_results_group)

# Save key analysis datasets

write.csv(cluster_dat, "cluster_dat_with_subtypes.csv", row.names = FALSE)
write.csv(analysis_dat, "analysis_dat_with_subtypes.csv", row.names = FALSE)
write.csv(sad_dat, "study2_sad_dat.csv", row.names = FALSE)
write.csv(img_sad, "study3_img_sad_dat.csv", row.names = FALSE)





#### SAD group ####
df_sad <- KSQ_SAD_combined

dat_sad <- df_sad %>%
  dplyr::select(SP, NRS, apnea, sleepiness) %>%
  mutate(.row_id = row_number()) %>%
  na.omit()

dat_vars_sad <- dat_sad %>%
  dplyr::select(SP, NRS, apnea, sleepiness)

# Check sample size and correlations
cat("SAD sample N =", nrow(dat_vars_sad), "\n")
print(summary(dat_vars_sad))
print(cor(dat_vars_sad, use = "pairwise.complete.obs"))

# Standardize
Xz_sad <- scale(dat_vars_sad)

# Hierarchical clustering
d_sad <- dist(Xz_sad, method = "euclidean")
hc_sad <- hclust(d_sad, method = "ward.D2")

# Evaluate optimal number of clusters
sil_width_sad <- sapply(2:6, function(k) {
  cl_tmp <- cutree(hc_sad, k = k)
  mean(silhouette(cl_tmp, d_sad)[, 3])
})

sil_summary_sad <- data.frame(k = 2:6, silhouette = sil_width_sad)
print(sil_summary_sad)

best_k_sad <- sil_summary_sad$k[which.max(sil_summary_sad$silhouette)]
cat("Best k for SAD sample =", best_k_sad, "\n")

write.csv(sil_summary_sad, "SAD_silhouette_summary.csv", row.names = FALSE)

# Silhouette plot
sil_sad <- fviz_nbclust(
  as.data.frame(Xz_sad),
  FUNcluster = hcut,
  method = "silhouette"
)

p_sil_sad <- ggplot(sil_sad$data, aes(x = sil_sad$data[[1]], y = sil_sad$data[[2]])) +
  geom_line(aes(group = 1), linewidth = 1.8, color = "#2C7FB8") +   # 粗线
  geom_point(size = 3.5, color = "#2C7FB8") +                       # 大点
  geom_vline(
    xintercept = sil_sad$nbclust,
    linetype = "dashed",
    linewidth = 1.1,
    color = "#2C7FB8"
  ) +
  theme_minimal() +
  labs(
    title = "Optimal number of clusters (SAD group)",
    x = "Number of clusters k",
    y = "Average silhouette width"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

print(p_sil_sad)

ggsave(
  "SAD_silhouette_plot.pdf",
  plot = p_sil_sad,
  width = 8,
  height = 6
)
# Final cluster solution
cl_sad <- cutree(hc_sad, k = best_k_sad)

dat_cluster_sad <- dat_sad %>%
  mutate(cluster_hc = factor(cl_sad))

print(table(dat_cluster_sad$cluster_hc))
write.csv(as.data.frame(table(dat_cluster_sad$cluster_hc)),
          "SAD_cluster_size.csv",
          row.names = FALSE)

# K-means validation
set.seed(123)
km_sad <- kmeans(Xz_sad, centers = best_k_sad, nstart = 50)

dat_cluster_sad$cluster_km <- factor(km_sad$cluster)

hc_km_table_sad <- table(dat_cluster_sad$cluster_hc, dat_cluster_sad$cluster_km)
print(hc_km_table_sad)
write.csv(as.data.frame.matrix(hc_km_table_sad),
          "SAD_hc_kmeans_crosstab.csv")

# Cluster profile plot
df_z_sad <- as.data.frame(Xz_sad)
df_z_sad$cluster <- factor(cl_sad)

cluster_n_sad <- table(df_z_sad$cluster)
cluster_labels_sad <- paste0("Cluster ", names(cluster_n_sad), " (n=", cluster_n_sad, ")")

df_z_sad$cluster <- factor(
  df_z_sad$cluster,
  levels = names(cluster_n_sad),
  labels = cluster_labels_sad
)

profile_sad <- df_z_sad %>%
  group_by(cluster) %>%
  summarise(
    across(c(apnea, NRS, sleepiness, SP), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -cluster,
    names_to = "dimension",
    values_to = "z_mean"
  )

profile_sad$dimension <- factor(
  profile_sad$dimension,
  levels = c("apnea", "NRS", "sleepiness", "SP"),
  labels = c("Apnea", "NRS", "Sleepiness", "SP")
)

write.csv(profile_sad, "SAD_cluster_profile_means.csv", row.names = FALSE)

p_profile_sad <- ggplot(
  profile_sad,
  aes(x = dimension, y = z_mean, group = cluster, color = cluster)
) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    title = "SAD cluster profiles",
    x = NULL,
    y = "Standardized mean (z-score)",
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.text = element_text(size = 9),
    legend.key.size = grid::unit(0.6, "cm"),
    legend.position = "right"
  )

print(p_profile_sad)

ggsave(
  "SAD_cluster_profile.pdf",
  plot = p_profile_sad,
  width = 7,
  height = 5.5
)

ggsave(
  "SAD_cluster_profile.tiff",
  plot = p_profile_sad,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# PCA plot
p_pca_sad <- fviz_cluster(
  list(data = Xz_sad, cluster = cl_sad),
  geom = "point",
  ellipse.type = "convex",
  pointsize = 2.5,
  palette = "jco",
  ggtheme = theme_classic(base_size = 14)
) +
  labs(
    title = "SAD PCA cluster plot",
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_pca_sad)

ggsave(
  "SAD_PCA_cluster_plot.pdf",
  plot = p_pca_sad,
  width = 8,
  height = 6
)

# Difference tests
df_cluster_sad <- df_sad[dat_sad$.row_id, ]
df_cluster_sad$cluster <- factor(cl_sad)

library(rstatix)
mw_results_sad <- df_cluster_sad %>%
  pivot_longer(
    cols = all_of(sleep_vars),
    names_to = "sleep_var",
    values_to = "value"
  ) %>%
  group_by(sleep_var) %>%
  wilcox_test(value ~ cluster)

print(mw_results_sad)

write.csv(
  mw_results_sad,
  "SAD_mann_whitney_results.csv",
  row.names = FALSE
)
getwd()

# Violin plot
long_df_sad <- dat_cluster_sad %>%
  pivot_longer(
    cols = c(SP, NRS, apnea, sleepiness),
    names_to = "dimension",
    values_to = "score"
  )

long_df_sad$dimension <- factor(
  long_df_sad$dimension,
  levels = c("SP", "NRS", "apnea", "sleepiness"),
  labels = c("SP", "NRS", "Apnea", "Sleepiness")
)

p_violin_sad <- ggplot(long_df_sad, aes(x = cluster_hc, y = score, fill = cluster_hc)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  facet_wrap(~ dimension, scales = "free", nrow = 2) +
  labs(
    title = "SAD: distribution of clustering variables by cluster",
    x = "Cluster",
    y = "Score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_violin_sad)

ggsave(
  "SAD_violin_boxplot_clusters.pdf",
  plot = p_violin_sad,
  width = 8,
  height = 6
)
getwd()




#### HC group ####
df_hc <- KSQ_HC_combined

dat_hc <- df_hc %>%
  dplyr::select(SP, NRS, apnea, sleepiness) %>%
  mutate(.row_id = row_number()) %>%
  na.omit()

dat_vars_hc <- dat_hc %>%
  dplyr::select(SP, NRS, apnea, sleepiness)

# Check sample size and correlations
cat("HC sample N =", nrow(dat_vars_hc), "\n")
print(summary(dat_vars_hc))
print(cor(dat_vars_hc, use = "pairwise.complete.obs"))

# Standardize
Xz_hc <- scale(dat_vars_hc)

# Hierarchical clustering
d_hc <- dist(Xz_hc, method = "euclidean")
hc_hc <- hclust(d_hc, method = "ward.D2")

# Evaluate optimal number of clusters
sil_width_hc <- sapply(2:6, function(k) {
  cl_tmp <- cutree(hc_hc, k = k)
  mean(silhouette(cl_tmp, d_hc)[, 3])
})

sil_summary_hc <- data.frame(k = 2:6, silhouette = sil_width_hc)
print(sil_summary_hc)

best_k_hc <- sil_summary_hc$k[which.max(sil_summary_hc$silhouette)]
cat("Best k for HC sample =", best_k_hc, "\n")

write.csv(sil_summary_hc, "HC_silhouette_summary.csv", row.names = FALSE)

# Silhouette plot
sil_hc <- fviz_nbclust(
  as.data.frame(Xz_hc),
  FUNcluster = hcut,
  method = "silhouette"
)

p_sil_hc <- ggplot(sil_hc$data, aes(x = sil_hc$data[[1]], y = sil_hc$data[[2]])) +
  geom_line(aes(group = 1), linewidth = 1.8, color = "#2C7FB8") +
  geom_point(size = 3.5, color = "#2C7FB8") +
  geom_vline(
    xintercept = sil_hc$nbclust,
    linetype = "dashed",
    linewidth = 1.1,
    color = "#2C7FB8"
  ) +
  theme_minimal() +
  labs(
    title = "Optimal number of clusters",
    x = "Number of clusters k",
    y = "Average silhouette width"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

print(p_sil_hc)

ggsave(
  "HC_silhouette_plot.pdf",
  plot = p_sil_hc,
  width = 8,
  height = 6
)

# Final cluster solution
cl_hc <- cutree(hc_hc, k = best_k_hc)

dat_cluster_hc <- dat_hc %>%
  mutate(cluster_hc = factor(cl_hc))

print(table(dat_cluster_hc$cluster_hc))

write.csv(
  as.data.frame(table(dat_cluster_hc$cluster_hc)),
  "HC_cluster_size.csv",
  row.names = FALSE
)

# K-means validation
set.seed(123)
km_hc <- kmeans(Xz_hc, centers = best_k_hc, nstart = 50)

dat_cluster_hc$cluster_km <- factor(km_hc$cluster)

hc_km_table_hc <- table(dat_cluster_hc$cluster_hc, dat_cluster_hc$cluster_km)
print(hc_km_table_hc)

write.csv(
  as.data.frame.matrix(hc_km_table_hc),
  "HC_hc_kmeans_crosstab.csv"
)

# Cluster profile plot
df_z_hc <- as.data.frame(Xz_hc)
df_z_hc$cluster <- factor(cl_hc)

cluster_n_hc <- table(df_z_hc$cluster)
cluster_labels_hc <- paste0("Cluster ", names(cluster_n_hc), " (n=", cluster_n_hc, ")")

df_z_hc$cluster <- factor(
  df_z_hc$cluster,
  levels = names(cluster_n_hc),
  labels = cluster_labels_hc
)

profile_hc <- df_z_hc %>%
  group_by(cluster) %>%
  summarise(
    across(c(apnea, NRS, sleepiness, SP), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -cluster,
    names_to = "dimension",
    values_to = "z_mean"
  )

profile_hc$dimension <- factor(
  profile_hc$dimension,
  levels = c("apnea", "NRS", "sleepiness", "SP"),
  labels = c("Apnea", "NRS", "Sleepiness", "SP")
)

write.csv(profile_hc, "HC_cluster_profile_means.csv", row.names = FALSE)

p_profile_hc <- ggplot(
  profile_hc,
  aes(x = dimension, y = z_mean, group = cluster, color = cluster)
) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    title = "HC group cluster profiles",
    x = NULL,
    y = "Standardized mean (z-score)",
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.text = element_text(size = 9),
    legend.key.size = grid::unit(0.6, "cm"),
    legend.position = "right"
  )

print(p_profile_hc)

ggsave(
  "HC_cluster_profile.pdf",
  plot = p_profile_hc,
  width = 7,
  height = 5.5
)

ggsave(
  "HC_cluster_profile.tiff",
  plot = p_profile_hc,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# PCA plot
p_pca_hc <- fviz_cluster(
  list(data = Xz_hc, cluster = cl_hc),
  geom = "point",
  ellipse.type = "convex",
  pointsize = 2.5,
  palette = "jco",
  ggtheme = theme_classic(base_size = 14)
) +
  labs(
    title = "HC group PCA cluster plot",
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_pca_hc)

ggsave(
  "HC_PCA_cluster_plot.pdf",
  plot = p_pca_hc,
  width = 8,
  height = 6
)

# Difference tests
df_cluster_hc <- df_hc[dat_hc$.row_id, ]
df_cluster_hc$cluster <- factor(cl_hc)

sleep_vars <- c("SP", "NRS", "apnea", "sleepiness")

library(rstatix)
kw_results_hc <- df_cluster_hc %>%
  pivot_longer(
    cols = all_of(sleep_vars),
    names_to = "sleep_var",
    values_to = "value"
  ) %>%
  group_by(sleep_var) %>%
  kruskal_test(value ~ cluster)

print(kw_results_hc)

write.csv(
  kw_results_hc,
  "HC_kruskal_results.csv",
  row.names = FALSE
)

posthoc_results_hc <- df_cluster_hc %>%
  pivot_longer(
    cols = all_of(sleep_vars),
    names_to = "sleep_var",
    values_to = "value"
  ) %>%
  group_by(sleep_var) %>%
  dunn_test(value ~ cluster, p.adjust.method = "bonferroni")

print(posthoc_results_hc)

write.csv(
  posthoc_results_hc,
  "HC_dunn_posthoc_results.csv",
  row.names = FALSE
)

# Violin plot
long_df_hc <- dat_cluster_hc %>%
  pivot_longer(
    cols = c(SP, NRS, apnea, sleepiness),
    names_to = "dimension",
    values_to = "score"
  )

long_df_hc$dimension <- factor(
  long_df_hc$dimension,
  levels = c("SP", "NRS", "apnea", "sleepiness"),
  labels = c("SP", "NRS", "Apnea", "Sleepiness")
)

p_violin_hc <- ggplot(long_df_hc, aes(x = cluster_hc, y = score, fill = cluster_hc)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  facet_wrap(~ dimension, scales = "free", nrow = 2) +
  labs(
    title = "HC group: distribution of clustering variables by cluster",
    x = "Cluster",
    y = "Score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_violin_hc)

ggsave(
  "HC_violin_boxplot_clusters.pdf",
  plot = p_violin_hc,
  width = 8,
  height = 6
)

