#### Q1a: Binary logistic regression analysis (sleep variables predict SAD diagnosis, without age and sex) ####

library(dplyr)
library(openxlsx)
library(pROC)
library(ggplot2)
library(pscl)

df1 <- combined_data

sleep_vars1 <- c(
  "SP", "NRS", "apnea",
  "sleepiness", "KSQ_total", "KSS_pre"
)

covars1 <- c("sample")

all_needed1 <- c("group", sleep_vars1, covars1)

df1 <- df1[, all_needed1]
df1 <- na.omit(df1)

#### Standardization (ONLY sleep variables, NO age anymore)
df1[, sleep_vars1] <- scale(df1[, sleep_vars1])

#### Group coding
df1$group_num1 <- ifelse(df1$group %in% c(1, "1", "SAD"), 1, 0)
df1$group_fac01_1 <- factor(df1$group_num1, levels = c(0, 1))
df1$group_fac1 <- factor(df1$group_num1, levels = c(0, 1),
                         labels = c("HC", "SAD"))

df1$sample <- factor(df1$sample)


#### Single-variable logistic regression (without age and sex)
run_logistic1 <- function(var_list, data) {
  
  results_list <- lapply(var_list, function(v) {
    
    form <- as.formula(
      paste("group_fac01_1 ~", v, "+ sample")
    )
    
    model <- glm(form, data = data, family = binomial)
    
    coef_table <- summary(model)$coefficients
    cis <- suppressMessages(confint(model))
    
    data.frame(
      variable = v,
      term = rownames(coef_table),
      beta = coef_table[, "Estimate"],
      p = coef_table[, "Pr(>|z|)"],
      lower = cis[, 1],
      upper = cis[, 2]
    )
  })
  
  out <- do.call(rbind, results_list)
  out <- out[out$term %in% var_list, ]
  
  out$OR <- exp(out$beta)
  out$lower <- exp(out$lower)
  out$upper <- exp(out$upper)
  
  out$P_raw <- out$p
  out$Bonf_raw <- p.adjust(out$p, method = "bonferroni")
  
  return(list(main = out))
}

## Run single variable analysis
sleep_results1 <- run_logistic1(sleep_vars1, df1)
sleep_results1
## Forest plot
sleep_results1$main$Sig3 <- with(
  sleep_results1$main,
  ifelse(Bonf_raw < 0.05, "Bonferroni",
         ifelse(P_raw < 0.05, "P < 0.05", "NS"))
)

sleep_results1$main$Sig3 <- factor(
  sleep_results1$main$Sig3,
  levels = c("Bonferroni", "P < 0.05", "NS")
)


plot_forest_pub1 <- function(dat, title_text, base_size = 12, show_bonf = TRUE) {
  
  dat <- dat[order(dat$OR, decreasing = FALSE), ]
  dat$Variable <- factor(dat$variable, levels = dat$variable)
  
  p <- ggplot(dat, aes(x = Variable, y = OR, color = Sig3)) +
    
    geom_hline(yintercept = 1,
               linetype = "dashed",
               linewidth = 0.6,
               color = "grey45") +
    
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.18,
                  linewidth = 0.7) +
    
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
  
  ## Bonferroni
  if (show_bonf) {
    p <- p +
      geom_text(
        data = subset(dat, Bonf_raw < 0.05),
        aes(label = paste0("Bonf=", sprintf("%.3f", Bonf_raw))),
        hjust = -0.1,
        size = 3.4,
        color = "black"
      ) +
      expand_limits(y = max(dat$upper) * 1.5)
  }
  
  return(p)
}


p_sleep1 <- plot_forest_pub1(
  sleep_results1$main,
  "Sleep-related variables (sample only)",
  base_size = 12,
  show_bonf = TRUE
)

print(p_sleep1)

ggsave(
  "Forestplot_sleep_sample_only.pdf",
  p_sleep1,
  width = 6,
  height = 4.5,
  device = cairo_pdf
)

#### ROC for each variable (without age and sex)

# Define significant variables based on p < 0.05 from single variable results
sig_sleep_vars1 <- sleep_results1$main %>%
  filter(P_raw < 0.05, variable != "(Intercept)") %>%
  pull(variable) %>%
  unique()

# If no variables are significant at p < 0.05, use all or top ones
if(length(sig_sleep_vars1) == 0) {
  sig_sleep_vars1 <- sleep_vars1
  message("No variables with p < 0.05, using all sleep variables")
}

roc_list_raw1 <- list()
roc_list_smooth1 <- list()
cutoff_list1 <- list()

auc_values1 <- numeric(length(sig_sleep_vars1))
names(auc_values1) <- sig_sleep_vars1


for (i in seq_along(sig_sleep_vars1)) {
  
  var <- sig_sleep_vars1[i]
  
  formula_i <- as.formula(
    paste("group_fac01_1 ~", paste(c(var, "sample"), collapse = " + "))
  )
  
  model_i <- glm(formula_i, data = df1, family = binomial)
  prob_i <- predict(model_i, type = "response")
  
  roc_i_raw <- roc(df1$group_num1, prob_i, ci = TRUE, quiet = TRUE)
  roc_i_smooth <- smooth(roc_i_raw)
  
  best_i <- coords(
    roc_i_raw,
    "best",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  cutoff_list1[[i]] <- data.frame(
    variable_raw = var,
    threshold = best_i$threshold,
    sensitivity = best_i$sensitivity,
    specificity = best_i$specificity
  )
  
  auc_i <- as.numeric(auc(roc_i_raw))
  auc_values1[var] <- auc_i
  
  roc_list_raw1[[i]] <- data.frame(
    FPR = 1 - roc_i_raw$specificities,
    TPR = roc_i_raw$sensitivities,
    variable_raw = var
  )
  
  roc_list_smooth1[[i]] <- data.frame(
    FPR = 1 - roc_i_smooth$specificities,
    TPR = roc_i_smooth$sensitivities,
    variable_raw = var
  )
}

roc_all_raw1 <- bind_rows(roc_list_raw1)
roc_all_smooth1 <- bind_rows(roc_list_smooth1)

auc_df1 <- data.frame(
  variable_raw = names(auc_values1),
  auc = as.numeric(auc_values1)
) %>%
  arrange(desc(auc))

cutoff_df1 <- bind_rows(cutoff_list1)

roc_summary1 <- auc_df1 %>%
  left_join(cutoff_df1, by = "variable_raw")

roc_summary1

label_map1 <- setNames(
  paste0(auc_df1$variable_raw,
         " (AUC = ", sprintf("%.3f", auc_df1$auc), ")"),
  auc_df1$variable_raw
)

roc_all_raw1$variable <- factor(
  label_map1[roc_all_raw1$variable_raw],
  levels = label_map1[auc_df1$variable_raw]
)

roc_all_smooth1$variable <- factor(
  label_map1[roc_all_smooth1$variable_raw],
  levels = label_map1[auc_df1$variable_raw]
)

palette_top <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0E442"
)

palette_use1 <- palette_top[seq_len(length(levels(roc_all_raw1$variable)))]


p_multi_roc_sleep_raw1 <- ggplot(
  roc_all_raw1,
  aes(x = FPR, y = TPR, color = variable)
) +
  geom_line(linewidth = 0.8) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey70"
  ) +
  scale_color_manual(values = palette_use1) +
  coord_equal() +
  labs(
    title = "ROC curves for individual sleep variables (sample-adjusted only)",
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

print(p_multi_roc_sleep_raw1)

ggsave(
  "ROC_sleep_individual_sample_only.pdf",
  plot = p_multi_roc_sleep_raw1,
  width = 6.5,
  height = 5.5,
  device = cairo_pdf
)

#### Multivariable model (without adjusting age and sex) ####

# Use variables that were significant in single variable analysis
# Or specify manually if needed
sig_sleep_vars_unadj <- c("SP", "NRS", "sleepiness", "KSS_pre")

formula_multi_unadj <- as.formula(
  paste(
    "group_fac01_1 ~",
    paste(c(sig_sleep_vars_unadj, "sample"),
          collapse = " + ")
  )
)

model_multi_unadj <- glm(
  formula_multi_unadj,
  data = df1,
  family = binomial
)

summary(model_multi_unadj)

#### Model significance

chi_squared_unadj <- model_multi_unadj$null.deviance -
  model_multi_unadj$deviance

p_value_unadj <- pchisq(
  chi_squared_unadj,
  df = model_multi_unadj$df.null -
    model_multi_unadj$df.residual,
  lower.tail = FALSE
)

p_value_unadj

#### McFadden pseudo-R²

pseudo_r2_unadj <- pR2(model_multi_unadj)

mcfadden_r2_unadj <- pseudo_r2_unadj["McFadden"]

mcfadden_r2_unadj

#### Events Per Variable (EPV)

events_unadj <- sum(df1$group_num1 == 1)

u_unadj <- length(coef(model_multi_unadj)) - 1

EPV_unadj <- events_unadj / u_unadj

EPV_unadj

#### Final regression table

coef_table_unadj <- summary(model_multi_unadj)$coefficients

cis_unadj <- suppressMessages(
  confint(model_multi_unadj)
)

final_table_unadj <- data.frame(
  Variable = rownames(coef_table_unadj),
  
  Beta = round(
    coef_table_unadj[, "Estimate"],
    3
  ),
  
  OR = round(
    exp(coef_table_unadj[, "Estimate"]),
    3
  ),
  
  CI_lower = round(
    exp(cis_unadj[, 1]),
    3
  ),
  
  CI_upper = round(
    exp(cis_unadj[, 2]),
    3
  ),
  
  p = round(
    coef_table_unadj[, "Pr(>|z|)"],
    3
  )
)

final_table_unadj$CI_95 <- paste0(
  final_table_unadj$CI_lower,
  "–",
  final_table_unadj$CI_upper
)

final_table_unadj <- final_table_unadj[, c(
  "Variable",
  "Beta",
  "OR",
  "CI_95",
  "p"
)]

final_table_unadj

#### ROC for multivariable model

model_sig_sleep_unadj <- glm(
  formula_multi_unadj,
  data = df1,
  family = binomial
)

prob_sig_sleep_unadj <- predict(
  model_sig_sleep_unadj,
  type = "response"
)

roc_sleep_unadj <- roc(
  df1$group_num1,  # Note: group_num1
  prob_sig_sleep_unadj,
  ci = TRUE
)

auc_sleep_unadj <- as.numeric(
  auc(roc_sleep_unadj)
)

auc_ci_sleep_unadj <- ci.auc(
  roc_sleep_unadj
)

roc_sleep_unadj

#### Best cutoff

best_sleep_unadj <- coords(
  roc_sleep_unadj,
  "best",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  transpose = FALSE
)

best_sleep_unadj

#### Brier Score

brier_sleep_unadj <- mean(
  (prob_sig_sleep_unadj - df1$group_num1)^2
)

brier_sleep_unadj

#### ROC visualization

roc_df_unadj <- data.frame(
  FPR = 1 - roc_sleep_unadj$specificities,
  TPR = roc_sleep_unadj$sensitivities
)

p_roc_sleep_unadj <- ggplot(
  roc_df_unadj,
  aes(x = FPR, y = TPR)
) +
  
  geom_line(
    linewidth = 1.0,
    color = "#D55E00"
  ) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    linewidth = 0.8,
    color = "grey60"
  ) +
  
  coord_equal() +
  
  labs(
    title = "ROC curve (without age/sex adjustment)",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  
  annotate(
    "label",
    x = 0.65,
    y = 0.15,
    
    label = paste0(
      "AUC = ",
      sprintf("%.3f", auc_sleep_unadj),
      
      "\n95% CI: ",
      sprintf("%.3f", auc_ci_sleep_unadj[1]),
      "–",
      sprintf("%.3f", auc_ci_sleep_unadj[3]),
      
      "\nEPV = ",
      round(EPV_unadj, 2),
      
      "\nMcFadden R² = ",
      round(as.numeric(mcfadden_r2_unadj), 3)
    ),
    
    size = 4.2,
    label.size = 0.3
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      face = "bold"
    )
  )

print(p_roc_sleep_unadj)

ggsave(
  "ROC_multivariable_logistic_unadjusted.png",
  p_roc_sleep_unadj,
  width = 5.5,
  height = 5,
  device = "png"
)



#### Q1b: Brain variables predict SAD diagnosis (without age and sex) ####

## Q1b_brain variables predict group membership

## Univariate Firth + Bayesian Multivariable Prediction
## Objective: Evaluate whether 36 ROI BOLD variability measures
## can discriminate SAD from HC, and assess their joint predictive value

df2 <- sample1

df2$group <- as.numeric(as.character(df2$group))

# 19 ROI variables
roi_vars <- c(
  "Amygdala_02","Caudata_02","Cingulate_Ant_02","Cingulate_Post_02",
  "Frontal_Mid_02","Frontal_Sup_02",
  "Frontal_Sup_Medial_02","Hypothalamus_02","Insula_02",
  "LC_02","N_Acc_02","OFCmed_02","Pallidum_02",
  "Precuneus_02","Putamen_02","Raphe_D_02","Raphe_M_02",
  "Thalamus_02","VTA_02"
)

# Keep only variables needed for analysis (NO age, NO sex)
df_analysis <- df2[, c("group", roi_vars)]

# Remove rows with missing values
df_analysis <- na.omit(df_analysis)

# Standardize all ROI variables (NO age anymore)
to_scale <- roi_vars
df_analysis[, to_scale] <- scale(df_analysis[, to_scale])

# Check final sample size
cat("Final N =", nrow(df_analysis), "\n")
cat("Number of ROI =", length(roi_vars), "\n")


#### Univariate Firth logistic regression (NO covariates) ####
## Model: group ~ ROI (without age and sex)
library(logistf)
library(dplyr)
library(openxlsx)

run_firth_univariate <- function(var_list, data, outcome = "group") {
  
  results_list <- lapply(var_list, function(v) {
    
    form <- as.formula(paste(outcome, "~", v))
    
    fit <- logistf(formula = form, data = data)
    
    beta <- coef(fit)[v]
    ci <- confint(fit)[v, ]
    pval <- fit$prob[v]
    
    data.frame(
      ROI = v,
      beta = unname(beta),
      OR = exp(unname(beta)),
      lowerCI = exp(ci[1]),
      upperCI = exp(ci[2]),
      p = unname(pval),
      stringsAsFactors = FALSE
    )
  })
  
  res <- bind_rows(results_list)
  res$FDR <- p.adjust(res$p, method = "BH")
  res$Holm <- p.adjust(res$p, method = "holm")
  res$Bonferroni <- p.adjust(res$p, method = "bonferroni")
  
  res$OR_95CI <- paste0(
    sprintf("%.2f", res$OR), " (",
    sprintf("%.2f", res$lowerCI), "-",
    sprintf("%.2f", res$upperCI), ")"
  )
  
  res <- res %>%
    arrange(p)
  
  return(res)
}

firth_results <- run_firth_univariate(
  var_list = roi_vars,
  data = df_analysis,
  outcome = "group"
)

print(firth_results)

write.xlsx(
  firth_results,
  file = "Firth_univariate_results_noCovars.xlsx",
  overwrite = TRUE
)


#### Power analysis for univariate Firth logistic regression (without age and sex) ####

# Base probability of group
base_prob <- mean(df_analysis$group)
n_total <- nrow(df_analysis)

# For univariate model without covariates, r2.other.x = 0
r2_covariates <- 0

power_results_by_roi <- firth_results %>%
  rowwise() %>%
  mutate(
    PostHoc_Power = pwrss.z.logistic(
      p0 = base_prob,
      odds.ratio = OR, 
      n = n_total,
      r2.other.x = r2_covariates,
      alpha = 0.05,
      alternative = "not equal"
    )$power,
    
    Required_N = pwrss.z.logistic(
      p0 = base_prob,
      odds.ratio = OR,
      power = 0.80,
      r2.other.x = r2_covariates,
      alpha = 0.05,
      alternative = "not equal"
    )$n
  ) %>%
  select(ROI, OR, p, FDR, PostHoc_Power, Required_N) %>%
  mutate(
    PostHoc_Power_percent = round(PostHoc_Power * 100, 1),
    Sufficient_Power = PostHoc_Power >= 0.80
  )

print(power_results_by_roi)

write.xlsx(power_results_by_roi, "PostHoc_Power_by_ROI_noCovars.xlsx")

# Sensitivity power analysis
or_range <- seq(1.1, 5.0, by = 0.1)
power_vals <- sapply(or_range, function(or) {
  pwrss.z.logistic(
    p0 = base_prob,
    odds.ratio = or,
    n = n_total,
    r2.other.x = r2_covariates,
    alpha = 0.05,
    alternative = "not equal",
    verbose = FALSE
  )$power
})

# Find the smallest OR achieving 80% power
target_idx <- which(power_vals >= 0.80)[1]

if(!is.na(target_idx)) {
  min_OR <- or_range[target_idx]
  cat("Minimum detectable OR (80% power):", round(min_OR, 2), "\n")
} else {
  cat("Cannot detect OR ≤ 5.0 with 80% power\n")
  cat("Maximum power achieved:", round(max(power_vals) * 100, 1), "% at OR =", 
      or_range[which.max(power_vals)], "\n")
}
getwd()

#### Volcano plot for univariate Firth results
library(ggplot2)
library(ggrepel)

volcano_data <- firth_results
volcano_data$logOR <- log(volcano_data$OR)
volcano_data$neglogP <- -log10(volcano_data$p)
volcano_data$Sig <- ifelse(volcano_data$p < 0.05, "P < 0.05", "NS")

# Label the top 5 ROIs ranked by raw p-value
top_vars <- volcano_data %>%
  arrange(p) %>%
  slice(1:5)

p_volcano <- ggplot(volcano_data, aes(x = logOR, y = neglogP, color = Sig)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = top_vars,
    aes(label = ROI),
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("P < 0.05" = "#D55E00", "NS" = "#4C72B0")) +
  labs(
    title = "Volcano plot of univariate Firth logistic regression (no covariates)",
    x = "log(Odds Ratio)",
    y = "-log10(p-value)",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = c(0.9, 0.88)
  )

print(p_volcano)

ggsave("Volcano_ROI_Firth_noCovars.pdf", p_volcano, width = 6.5, height = 5.5)


#### FULL BAYESIAN MULTIVARIABLE LOGISTIC REGRESSION (NO COVARIATES) ####
## Model: group ~ all ROIs (without age and sex)
## Prior: Regularized Horseshoe (for shrinkage & variable selection)
library(brms)
library(bayesplot)
library(tidybayes)
library(pROC)
library(ggplot2)
library(openxlsx)

#### 1. Prepare Data and Formula ####

# Ensure data is ready
df_analysis <- na.omit(df_analysis)  # Already done, but safe check

# Formula with all predictors (NO age, NO sex)
formula_bayes <- as.formula(
  paste("group ~", paste(roi_vars, collapse = " + "))
)

#### 2. Set Priors ####

# Regularized Horseshoe prior (recommended)
priors_horseshoe <- c(
  prior(horseshoe(scale_global = 0.05, scale_slab = 2), class = "b"),
  prior(normal(0, 2.5), class = "Intercept")
)

cat("========== Full Bayesian Model with Regularized Horseshoe Prior ==========\n")
cat("Prior settings:\n")
cat("  - scale_global = 0.05 (controls overall sparsity)\n")
cat("  - scale_slab = 2 (allows for large effects when supported by data)\n")
cat("  - Intercept: Normal(0, 2.5)\n")
cat("  - NO covariates (age and sex excluded)\n\n")

#### 3. Fit Full Bayesian Model ####

cat("Fitting model with 4 chains, 8000 iterations...\n")
cat("This may take several minutes depending on your data size...\n\n")

fit_bayes_full <- brm(
  formula = formula_bayes,
  data = df_analysis,
  family = bernoulli(link = "logit"),
  prior = priors_horseshoe,
  chains = 4,
  iter = 8000,
  warmup = 4000,
  thin = 2,
  cores = 4,
  seed = 1234,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  save_pars = save_pars(all = TRUE)
)

cat("Model fitting complete!\n\n")

#### 4. Model Diagnostics ####

cat("========== Model Diagnostics ==========\n")

# R-hat values
rhat_values <- summary(fit_bayes_full)$fixed[, "Rhat"]
max_rhat <- max(rhat_values)
cat("Maximum R-hat:", round(max_rhat, 3), "\n")
if(max_rhat < 1.01) {
  cat("✓ Good convergence (all R-hat < 1.01)\n")
} else {
  cat("⚠️ Some parameters may not have converged\n")
}

# Effective sample size
neff_ratio <- summary(fit_bayes_full)$fixed[, "Bulk_ESS"] / 8000
min_neff <- min(neff_ratio)
cat("Minimum ESS ratio:", round(min_neff, 2), "\n")
if(min_neff > 0.1) {
  cat("✓ Adequate effective sample size\n")
}

# Trace plots for key parameters
roi_params <- paste0("b_", roi_vars)

# ROI 1-6
trace_roi_1 <- mcmc_trace(fit_bayes_full, pars = roi_params[1:6]) +
  ggtitle("Trace Plots: ROIs 1-6 (no covariates)")

ggsave("Trace_Plots_ROI_1-6_noCovars.png", trace_roi_1, width = 12, height = 8)

# ROI 7-12
trace_roi_2 <- mcmc_trace(fit_bayes_full, pars = roi_params[7:12]) +
  ggtitle("Trace Plots: ROIs 7-12 (no covariates)")

ggsave("Trace_Plots_ROI_7-12_noCovars.png", trace_roi_2, width = 12, height = 8)

# ROI 13-19
trace_roi_3 <- mcmc_trace(fit_bayes_full, pars = roi_params[13:19]) +
  ggtitle("Trace Plots: ROIs 13-19 (no covariates)")

ggsave("Trace_Plots_ROI_13-19_noCovars.png", trace_roi_3, width = 12, height = 8)
cat("✓ Trace plots saved\n")


#### 5. Posterior Estimates ####

cat("\n========== Posterior Estimates ==========\n")

library(tibble)
fixed_effects <- summary(fit_bayes_full)$fixed

posterior_summary <- as.data.frame(fixed_effects) %>%
  rownames_to_column("Parameter") %>%
  mutate(
    OR = exp(Estimate),
    OR_lower = exp(`l-95% CI`),   
    OR_upper = exp(`u-95% CI`),  
    # Probability of being positive (for effects > 0)
    Prob_positive = ifelse(Parameter == "b_Intercept", NA, 
                           1 - pnorm(0, mean = Estimate, sd = Est.Error)),
    # Probability of meaningful effect (|logOR| > log(1.2))
    Prob_meaningful = ifelse(Parameter == "b_Intercept",
                             NA,
                             pmax(
                               pnorm(log(1.2), mean = Estimate, sd = Est.Error, 
                                     lower.tail = FALSE),
                               pnorm(-log(1.2), mean = Estimate, sd = Est.Error, 
                                     lower.tail = TRUE)
                             ))
  )

# Display results for ROIs only
cat("\nTop 10 ROIs by absolute effect size:\n")
posterior_summary %>%
  filter(Parameter %in% paste0("b_", roi_vars)) %>%
  arrange(desc(abs(Estimate))) %>%
  head(10) %>%
  select(Parameter, OR, OR_lower, OR_upper, Prob_meaningful) %>%
  print(digits = 3)

# Save full summary
write.xlsx(posterior_summary, "Bayesian_Posterior_Summary_Horseshoe_noCovars.xlsx")

#### 6. Variable Importance (Posterior Inclusion Probabilities) ####

cat("\n========== Variable Importance ==========\n")

# Extract posterior draws for all coefficients (simple method)
# Get all b_ parameters except Intercept
all_params <- variables(fit_bayes_full)
b_params <- all_params[grep("^b_", all_params)]
b_params_no_intercept <- b_params[b_params != "b_Intercept"]

# Extract posterior draws
posterior_draws <- as_draws_matrix(fit_bayes_full, variable = b_params_no_intercept)

cat("Extracted", ncol(posterior_draws), "parameters\n")

# Calculate inclusion probability (|beta| > 0.1 on log-odds scale)
inclusion_threshold <- 0.1
inclusion_prob <- apply(abs(posterior_draws) > inclusion_threshold, 2, mean)

# Probability of direction (consistent sign)
prob_direction <- apply(posterior_draws, 2, function(x) max(mean(x > 0), mean(x < 0)))

# Create importance summary
variable_importance <- data.frame(
  Parameter = colnames(posterior_draws),
  Inclusion_Prob = inclusion_prob,
  Prob_Direction = prob_direction,
  Median_LogOR = apply(posterior_draws, 2, median),
  OR = exp(apply(posterior_draws, 2, median)),
  OR_lower = exp(apply(posterior_draws, 2, function(x) quantile(x, 0.025))),
  OR_upper = exp(apply(posterior_draws, 2, function(x) quantile(x, 0.975)))
) %>%
  arrange(desc(Inclusion_Prob))

cat("Variables with inclusion probability > 0.5:\n")
important_vars <- variable_importance %>% filter(Inclusion_Prob > 0.5)
if(nrow(important_vars) > 0) {
  print(important_vars[, c("Parameter", "Inclusion_Prob", "OR", "OR_lower", "OR_upper")])
} else {
  cat("  No variables with inclusion probability > 0.5\n")
  cat("  Top 5 variables:\n")
  print(head(variable_importance[, c("Parameter", "Inclusion_Prob", "OR")], 5))
}

# Save variable importance
write.xlsx(variable_importance, "Bayesian_Variable_Importance_Horseshoe_noCovars.xlsx")

#### 7. Posterior Predictive Checks ####

cat("\n========== Posterior Predictive Checks ==========\n")

# Generate posterior predictions
pred_bayes <- posterior_predict(fit_bayes_full, ndraws = 100)
observed <- df_analysis$group

# PPC plot
ppc_plot <- ppc_dens_overlay(observed, pred_bayes[1:50, ]) +
  ggtitle("Posterior Predictive Check - Density (no covariates)")

ggsave("Bayesian_PPC_Density_noCovars.png", ppc_plot, width = 8, height = 6)

# Bayesian p-value
bayes_p_value <- mean(apply(pred_bayes, 1, function(x) mean(x)) > mean(observed))
cat("Bayesian p-value:", round(bayes_p_value, 3), "\n")
if(bayes_p_value > 0.025 & bayes_p_value < 0.975) {
  cat("✓ Model fits the data adequately\n")
}

#### 8. Model Performance Metrics ####

cat("\n========== Model Performance ==========\n")

# Predicted probabilities
pred_prob_full <- fitted(fit_bayes_full, scale = "response")[, "Estimate"]
pred_ci_lower <- fitted(fit_bayes_full, scale = "response")[, "Q2.5"]
pred_ci_upper <- fitted(fit_bayes_full, scale = "response")[, "Q97.5"]

# AUC
roc_obj_full <- roc(observed, pred_prob_full, quiet = TRUE)
auc_full <- as.numeric(auc(roc_obj_full))
auc_ci_full <- ci.auc(roc_obj_full)

cat("AUC:", round(auc_full, 3), "\n")
cat("AUC 95% CI:", round(auc_ci_full[1], 3), "-", round(auc_ci_full[3], 3), "\n")

# Brier score
brier_full <- mean((observed - pred_prob_full)^2)
cat("Brier score:", round(brier_full, 4), "\n")

# Bayesian R-squared
bayes_r2 <- bayes_R2(fit_bayes_full)
cat("Bayesian R-squared:", round(mean(bayes_r2), 3), 
    "| 95% CI:", round(quantile(bayes_r2, 0.025), 3), "-", round(quantile(bayes_r2, 0.975), 3), "\n")

#### 9. ROC Curve ####

roc_df <- data.frame(
  FPR = 1 - roc_obj_full$specificities,
  TPR = roc_obj_full$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(linewidth = 1.2, color = "#1F78B4") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  coord_equal() +
  annotate(
    "label",
    x = 0.65, y = 0.15,
    label = paste0(
      "AUC = ", sprintf("%.3f", auc_full),
      "\n95% CI: ", sprintf("%.3f", auc_ci_full[1]), "–", sprintf("%.3f", auc_ci_full[3])
    ),
    size = 4
  ) +
  labs(
    title = "ROC Curve - Full Bayesian Model (Horseshoe Prior, no covariates)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_classic(base_size = 13)

ggsave("Bayesian_ROC_Curve_Horseshoe_noCovars.pdf", p_roc, width = 6, height = 6)
ggsave("Bayesian_ROC_Curve_Horseshoe_noCovars.png", p_roc, width = 6, height = 6)


#### 10. LOO Cross-Validation ####

cat("\n========== Model Comparison ==========\n")

# LOO for full model
loo_full <- loo(fit_bayes_full, moment_match = TRUE)
print(loo_full)

# Null model for comparison
fit_null <- brm(
  group ~ 1,
  data = df_analysis,
  family = bernoulli(link = "logit"),
  prior = prior(normal(0, 2.5), class = "Intercept"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234
)

loo_null <- loo(fit_null)
loo_compare <- loo_compare(loo_full, loo_null)
print(loo_compare)

#### 11. Save All Results ####

# Compile all results
final_results <- list(
  Model_Summary = summary(fit_bayes_full),
  Posterior_Estimates = posterior_summary,
  Variable_Importance = variable_importance,
  Model_Performance = data.frame(
    Metric = c("AUC", "AUC_CI_lower", "AUC_CI_upper", 
               "Brier_Score", "Bayesian_R2_mean", "Bayesian_R2_lower", "Bayesian_R2_upper"),
    Value = c(auc_full, auc_ci_full[1], auc_ci_full[3],
              brier_full, mean(bayes_r2), quantile(bayes_r2, 0.025), quantile(bayes_r2, 0.975))
  ),
  LOO_Comparison = as.data.frame(loo_compare),
  Predictions = data.frame(
    observed = observed,
    predicted = pred_prob_full,
    lower_ci = pred_ci_lower,
    upper_ci = pred_ci_upper
  )
)

# Save to Excel
write.xlsx(
  list(
    Readme = data.frame(
      Note = "Full Bayesian Model with Regularized Horseshoe Prior (NO covariates)",
      Date = Sys.Date(),
      Prior_scale_global = 0.05,
      Prior_scale_slab = 2,
      Chains = 4,
      Iterations = 8000,
      Covariates_included = "None"
    ),
    Posterior_Summary = posterior_summary,
    Variable_Importance = variable_importance,
    Model_Performance = final_results$Model_Performance,
    LOO_Comparison = final_results$LOO_Comparison,
    Predictions = final_results$Predictions
  ),
  file = "Full_Bayesian_Results_Horseshoe_noCovars.xlsx",
  overwrite = TRUE
)

# Save R objects
saveRDS(fit_bayes_full, "Full_Bayesian_Model_Horseshoe_noCovars.rds")
saveRDS(final_results, "Full_Bayesian_Results_Horseshoe_noCovars.rds")

cat("\n========== Complete! ==========\n")
cat("All results saved:\n")
cat("  - Full_Bayesian_Results_Horseshoe_noCovars.xlsx (main results)\n")
cat("  - Full_Bayesian_Model_Horseshoe_noCovars.rds (model object)\n")
cat("  - Bayesian_ROC_Curve_Horseshoe_noCovars.pdf\n")
cat("  - Bayesian_PPC_Density_noCovars.png\n")
cat("  - Trace_Plots_ROI_*_noCovars.png\n")

#### Q2a, Q2b: sleep and brain variables predict ICBT treatment outcome ####

library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(writexl)
library(performance)
library(simr)

#### Sample 1 (NO age, NO sex) ####
df3 <- sample1_SAD

## reshape to long format
df_long <- df3 %>% 
  pivot_longer(
    cols = c("LSAS_pre","LSAS_post"),
    names_to = "time",
    names_prefix = "LSAS_",
    values_to = "LSAS"
  )

## factor
df_long$time <- factor(df_long$time, levels = c("pre","post"))

sleep_vars <- c("SP","NRS","apnea","sleepiness","KSQ_total","KSS_pre")

brain_vars <- c(
  "Amygdala_02","Caudata_02","Cingulate_Ant_02","Cingulate_Post_02",
  "Frontal_Mid_02","Frontal_Sup_02","Frontal_Sup_Medial_02",
  "Hypothalamus_02","Insula_02",
  "LC_02","N_Acc_02","OFCmed_02","Pallidum_02",
  "Precuneus_02","Putamen_02","Raphe_D_02","Raphe_M_02",
  "Thalamus_02","VTA_02"
)

all_predictors <- c(sleep_vars, brain_vars)

## standardize variables (NO age)
to_center <- all_predictors

for (v in to_center){
  new_name <- paste0(v,"_c")
  df_long[[new_name]] <- scale(df_long[[v]], center=TRUE, scale=TRUE)
}

sleep_vars_c <- paste0(sleep_vars,"_c")
brain_vars_c <- paste0(brain_vars,"_c")


## Function to run LMM and extract main + interaction (NO covariates)

run_lmm <- function(var_c, data){
  
  formula_txt <- paste0(
    "LSAS ~ ", var_c," * time + (1|subject)"
  )
  
  model_formula <- as.formula(formula_txt)
  
  out <- tryCatch({
    
    model <- lmer(model_formula, data=data, REML=FALSE)
    
    coef_tab <- summary(model)$coefficients
    rn <- rownames(coef_tab)
    
    r2 <- performance::r2(model)
    
    int1 <- paste0(var_c,":timepost")
    int2 <- paste0("timepost:",var_c)
    
    interaction_name <- if(int1 %in% rn){
      int1
    } else if(int2 %in% rn){
      int2
    } else{
      NA_character_
    }
    
    main_name <- var_c
    time_name <- "timepost"
    
    terms_keep <- c(main_name, time_name, interaction_name)
    terms_keep <- terms_keep[terms_keep %in% rn]
    terms_keep <- na.omit(terms_keep)
    
    ci <- suppressMessages(
      confint(model, parm=terms_keep, method="Wald")
    )
    
    res <- lapply(terms_keep, function(tt){
      
      effect_type <- if(tt==main_name){
        "main_predictor"
      } else if(tt==time_name){
        "main_time"
      } else{
        "interaction"
      }
      
      data.frame(
        variable = gsub("_c$","",var_c),
        effect = effect_type,
        term = tt,
        beta = coef_tab[tt,"Estimate"],
        se   = coef_tab[tt,"Std. Error"],
        CI_lower = ci[tt,1],
        CI_upper = ci[tt,2],
        p = coef_tab[tt,"Pr(>|t|)"],
        
        R2_marginal   = r2$R2_marginal,
        R2_conditional= r2$R2_conditional,
        
        n = nrow(model.frame(model)),
        stringsAsFactors = FALSE
      )
      
    }) %>% bind_rows()
    
    res
    
  }, error=function(e){
    
    data.frame(
      variable = gsub("_c$","",var_c),
      effect = NA,
      term = NA,
      beta = NA,
      se = NA,
      CI_lower = NA,
      CI_upper = NA,
      p = NA,
      R2_marginal = NA,
      R2_conditional = NA,
      n = NA
    )
    
  })
  
  return(out)
}


## Sleep models
sleep_results <- lapply(
  sleep_vars_c,
  run_lmm,
  data=df_long
) %>%
  bind_rows() %>%
  mutate(group="sleep")

## Brain models
brain_results <- lapply(
  brain_vars_c,
  run_lmm,
  data=df_long
) %>%
  bind_rows() %>%
  mutate(group="brain")


## Multiple comparison correction (interaction only)
sleep_results <- sleep_results %>%
  mutate(
    p_fdr = p.adjust(p, method = "BH"), 
    p_bonf = p.adjust(p, method = "bonferroni"), 
    p_fwe = p.adjust(p, method = "holm") 
  )

brain_results <- brain_results %>%
  mutate(
    p_fdr = p.adjust(p, method = "BH"),
    p_bonf = p.adjust(p, method = "bonferroni"),
    p_fwe = p.adjust(p, method = "holm") 
  )

all_results <- bind_rows(sleep_results, brain_results)

write_xlsx(
  list(
    sleep_results = sleep_results,
    brain_results = brain_results,
    all_results   = all_results
  ),
  path = "LMM_sample1_noCovars_main_and_interaction.xlsx"
)

#### Statistical power (NO covariates) ####
## sleep variables

## LMM_SP
model_SP <- lmer(
  LSAS ~ SP_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

## check fixed effects
names(fixef(model_SP))

## Power analysis: interaction effect
## beta = 0.2

fixef(model_SP)["SP_c:timepost"] <- 0.2

power_beta_SP <- powerSim(
  model_SP,
  test = fixed("SP_c:timepost", "t"),
  nsim = 100
)

print(power_beta_SP)

## LMM_NRS
model_NRS <- lmer(
  LSAS ~ NRS_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

## check fixed effects
names(fixef(model_NRS))

## Power analysis: interaction effect
## beta = 0.2

fixef(model_NRS)["NRS_c:timepost"] <- 0.2

power_beta_NRS <- powerSim(
  model_NRS,
  test = fixed("NRS_c:timepost", "t"),
  nsim = 100
)

print(power_beta_NRS)

## LMM_apnea
model_apnea <- lmer(
  LSAS ~ apnea_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_apnea)["apnea_c:timepost"] <- 0.2

power_apnea <- powerSim(
  model_apnea,
  test = fixed("apnea_c:timepost","t"),
  nsim = 100
)

print(power_apnea)

## LMM_sleepiness
model_sleepiness <- lmer(
  LSAS ~ sleepiness_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_sleepiness)["sleepiness_c:timepost"] <- 0.2

power_sleepiness <- powerSim(
  model_sleepiness,
  test = fixed("sleepiness_c:timepost","t"),
  nsim = 100
)

print(power_sleepiness)

## LMM_KSQ_total
model_KSQ <- lmer(
  LSAS ~ KSQ_total_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_KSQ)["KSQ_total_c:timepost"] <- 0.2

power_KSQ <- powerSim(
  model_KSQ,
  test = fixed("KSQ_total_c:timepost","t"),
  nsim = 100
)

print(power_KSQ)

## LMM_KSS
model_KSS <- lmer(
  LSAS ~ KSS_pre_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_KSS)["KSS_pre_c:timepost"] <- 0.2

power_KSS <- powerSim(
  model_KSS,
  test = fixed("KSS_pre_c:timepost","t"),
  nsim = 100
)

print(power_KSS)

#### Brain variables (NO covariates) ####

## Amygdala
model_Amygdala <- lmer(
  LSAS ~ Amygdala_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Amygdala)["Amygdala_02_c:timepost"] <- 0.2

power_Amygdala <- powerSim(model_Amygdala, test=fixed("Amygdala_02_c:timepost","t"), nsim=100)
print(power_Amygdala)

## Caudata
model_Caudata <- lmer(
  LSAS ~ Caudata_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Caudata)["Caudata_02_c:timepost"] <- 0.2
power_Caudata <- powerSim(model_Caudata, test=fixed("Caudata_02_c:timepost","t"), nsim=100)
print(power_Caudata)

## Cingulate_Ant
model_CingA <- lmer(
  LSAS ~ Cingulate_Ant_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_CingA)["Cingulate_Ant_02_c:timepost"] <- 0.2
power_CingA <- powerSim(model_CingA, test=fixed("Cingulate_Ant_02_c:timepost","t"), nsim=100)
print(power_CingA)

## Cingulate_Post
model_CingP <- lmer(
  LSAS ~ Cingulate_Post_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_CingP)["Cingulate_Post_02_c:timepost"] <- 0.2
power_CingP <- powerSim(model_CingP, test=fixed("Cingulate_Post_02_c:timepost","t"), nsim=100)
print(power_CingP)

## Frontal_Mid
model_FM <- lmer(
  LSAS ~ Frontal_Mid_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_FM)["Frontal_Mid_02_c:timepost"] <- 0.2
power_FM <- powerSim(model_FM, test=fixed("Frontal_Mid_02_c:timepost","t"), nsim=100)
print(power_FM)

## Frontal_Sup
model_FS <- lmer(
  LSAS ~ Frontal_Sup_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_FS)["Frontal_Sup_02_c:timepost"] <- 0.2
power_FS <- powerSim(model_FS, test=fixed("Frontal_Sup_02_c:timepost","t"), nsim=100)
print(power_FS)

## Frontal_Sup_Medial
model_FSM <- lmer(
  LSAS ~ Frontal_Sup_Medial_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_FSM)["Frontal_Sup_Medial_02_c:timepost"] <- 0.2
power_FSM <- powerSim(model_FSM, test=fixed("Frontal_Sup_Medial_02_c:timepost","t"), nsim=100)
print(power_FSM)

## Hypothalamus
model_Hyp <- lmer(
  LSAS ~ Hypothalamus_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Hyp)["Hypothalamus_02_c:timepost"] <- 0.2
power_Hyp <- powerSim(model_Hyp, test=fixed("Hypothalamus_02_c:timepost","t"), nsim=100)
print(power_Hyp)

## Insula
model_Insula <- lmer(
  LSAS ~ Insula_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Insula)["Insula_02_c:timepost"] <- 0.2
power_Insula <- powerSim(model_Insula, test=fixed("Insula_02_c:timepost","t"), nsim=100)
print(power_Insula)

## LC
model_LC <- lmer(
  LSAS ~ LC_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_LC)["LC_02_c:timepost"] <- 0.2
power_LC <- powerSim(model_LC, test=fixed("LC_02_c:timepost","t"), nsim=100)
print(power_LC)

## N_Acc
model_NAcc <- lmer(
  LSAS ~ N_Acc_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_NAcc)["N_Acc_02_c:timepost"] <- 0.2

power_NAcc <- powerSim(
  model_NAcc,
  test = fixed("N_Acc_02_c:timepost","t"),
  nsim = 100
)

print(power_NAcc)

## OFCmed
model_OFCmed <- lmer(
  LSAS ~ OFCmed_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_OFCmed)["OFCmed_02_c:timepost"] <- 0.2

power_OFCmed <- powerSim(
  model_OFCmed,
  test = fixed("OFCmed_02_c:timepost","t"),
  nsim = 100
)

print(power_OFCmed)

## Pallidum
model_Pallidum <- lmer(
  LSAS ~ Pallidum_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Pallidum)["Pallidum_02_c:timepost"] <- 0.2

power_Pallidum <- powerSim(
  model_Pallidum,
  test = fixed("Pallidum_02_c:timepost","t"),
  nsim = 100
)

print(power_Pallidum)

## Precuneus
model_Precuneus <- lmer(
  LSAS ~ Precuneus_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Precuneus)["Precuneus_02_c:timepost"] <- 0.2

power_Precuneus <- powerSim(
  model_Precuneus,
  test = fixed("Precuneus_02_c:timepost","t"),
  nsim = 100
)

print(power_Precuneus)

## Putamen
model_Putamen <- lmer(
  LSAS ~ Putamen_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Putamen)["Putamen_02_c:timepost"] <- 0.2

power_Putamen <- powerSim(
  model_Putamen,
  test = fixed("Putamen_02_c:timepost","t"),
  nsim = 100
)

print(power_Putamen)

## Raphe_D
model_Raphe_D <- lmer(
  LSAS ~ Raphe_D_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Raphe_D)["Raphe_D_02_c:timepost"] <- 0.2

power_Raphe_D <- powerSim(
  model_Raphe_D,
  test = fixed("Raphe_D_02_c:timepost","t"),
  nsim = 100
)

print(power_Raphe_D)

## Raphe_M
model_Raphe_M <- lmer(
  LSAS ~ Raphe_M_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Raphe_M)["Raphe_M_02_c:timepost"] <- 0.2

power_Raphe_M <- powerSim(
  model_Raphe_M,
  test = fixed("Raphe_M_02_c:timepost","t"),
  nsim = 100
)

print(power_Raphe_M)

## Thalamus
model_Thalamus <- lmer(
  LSAS ~ Thalamus_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Thalamus)["Thalamus_02_c:timepost"] <- 0.2

power_Thalamus <- powerSim(
  model_Thalamus,
  test = fixed("Thalamus_02_c:timepost","t"),
  nsim = 100
)

print(power_Thalamus)

## VTA
model_VTA <- lmer(
  LSAS ~ VTA_02_c * time + (1|subject),
  data = df_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_VTA)["VTA_02_c:timepost"] <- 0.2

power_VTA <- powerSim(
  model_VTA,
  test = fixed("VTA_02_c:timepost","t"),
  nsim = 100
)

print(power_VTA)



#### Q2c: sleep variables predict SSRI treatment outcome ####
df4 <- sample2_SAD

## predictors for sample2
sleep_vars_s2 <- c("SP", "NRS", "apnea", "sleepiness", "KSQ_total", "KSS_pre")

brain_vars_s2 <- c(
  "Midbrain", "Rapha", "ACC", "PCC", "Thalamus",
  "Insula", "Striatum", "Caudate", "Putamen",
  "Hippocampus", "Amygdala", "Hypothalamus",
  "OFC", "MedInfFrontal", "SupFrontal"
)

all_predictors_s2 <- c(sleep_vars_s2, brain_vars_s2)

## reshape to long format
df4_long <- df4 %>%
  pivot_longer(
    cols = c("LSASSR_pre", "LSASSR_post"),
    names_to = "time",
    names_prefix = "LSASSR_",
    values_to = "LSASSR"
  )

## factor
df4_long$time <- factor(df4_long$time, levels = c("pre", "post"))

## standardize variables (NO age)
to_center_s2 <- all_predictors_s2

for (v in to_center_s2) {
  new_name <- paste0(v, "_c")
  df4_long[[new_name]] <- as.numeric(scale(df4_long[[v]], center = TRUE, scale = TRUE))
}

sleep_vars_c_s2 <- paste0(sleep_vars_s2, "_c")
brain_vars_c_s2 <- paste0(brain_vars_s2, "_c")

## function to run LMM and extract main predictor + main time + interaction (NO covariates)
run_lmm_s2 <- function(var_c, data) {
  
  formula_txt <- paste0(
    "LSASSR ~ ", var_c, " * time + (1|subject)"
  )
  
  model_formula <- as.formula(formula_txt)
  
  out <- tryCatch({
    
    model <- lmer(model_formula, data = data, REML = FALSE)
    
    coef_tab <- summary(model)$coefficients
    rn <- rownames(coef_tab)
    
    ## R2
    r2 <- performance::r2(model)
    
    main_name <- var_c
    time_name <- "timepost"
    
    int1 <- paste0(var_c, ":timepost")
    int2 <- paste0("timepost:", var_c)
    
    interaction_name <- if (int1 %in% rn) {
      int1
    } else if (int2 %in% rn) {
      int2
    } else {
      NA_character_
    }
    
    terms_keep <- c(main_name, time_name, interaction_name)
    terms_keep <- na.omit(terms_keep)
    terms_keep <- terms_keep[terms_keep %in% rn]
    
    ci <- suppressMessages(
      confint(model, parm = terms_keep, method = "Wald")
    )
    
    res <- lapply(terms_keep, function(tt) {
      
      effect_type <- if (tt == main_name) {
        "main_predictor"
      } else if (tt == time_name) {
        "main_time"
      } else {
        "interaction"
      }
      
      data.frame(
        group = NA_character_,
        variable = gsub("_c$", "", var_c),
        effect = effect_type,
        term = tt,
        beta = coef_tab[tt, "Estimate"],
        se = coef_tab[tt, "Std. Error"],
        CI_lower = ci[tt, 1],
        CI_upper = ci[tt, 2],
        p = coef_tab[tt, "Pr(>|t|)"],
        R2_marginal = r2$R2_marginal,
        R2_conditional = r2$R2_conditional,
        n = nrow(model.frame(model)),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
    
    res
    
  }, error = function(e) {
    
    message("Error in variable: ", var_c)
    message(e$message)
    
    data.frame(
      group = NA_character_,
      variable = gsub("_c$", "", var_c),
      effect = NA_character_,
      term = NA_character_,
      beta = NA_real_,
      se = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      p = NA_real_,
      R2_marginal = NA_real_,
      R2_conditional = NA_real_,
      n = NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  
  return(out)
}

## sleep models
sleep_results_s2 <- lapply(
  sleep_vars_c_s2,
  run_lmm_s2,
  data = df4_long
) %>%
  bind_rows() %>%
  mutate(group = "sleep") %>%
  group_by(effect) %>%
  mutate(
    p_fdr  = p.adjust(p, method = "BH"),
    p_bonf = p.adjust(p, method = "bonferroni"),
    p_fwe  = p.adjust(p, method = "holm")
  ) %>%
  ungroup() %>%
  select(group, variable, effect, term, beta, se, CI_lower, CI_upper,
         p, p_fdr, p_bonf, p_fwe, R2_marginal, R2_conditional, n)

## brain models
brain_results_s2 <- lapply(
  brain_vars_c_s2,
  run_lmm_s2,
  data = df4_long
) %>%
  bind_rows() %>%
  mutate(group = "brain") %>%
  group_by(effect) %>%
  mutate(
    p_fdr  = p.adjust(p, method = "BH"),
    p_bonf = p.adjust(p, method = "bonferroni"),
    p_fwe  = p.adjust(p, method = "holm")
  ) %>%
  ungroup() %>%
  select(group, variable, effect, term, beta, se, CI_lower, CI_upper,
         p, p_fdr, p_bonf, p_fwe, R2_marginal, R2_conditional, n)

## all results
all_results_s2 <- bind_rows(sleep_results_s2, brain_results_s2)

## export
write_xlsx(
  list(
    sleep_results_sample2 = sleep_results_s2,
    brain_results_sample2 = brain_results_s2,
    all_results_sample2   = all_results_s2
  ),
  path = "LMM_sample2_noCovars_main_time_interaction_with_R2.xlsx"
)

#### Statistical power for sample 2 (NO covariates) ####

## SP
model_SP_s2 <- lmer(
  LSASSR ~ SP_c * time + (1|subject),
  data = df4_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_SP_s2)["SP_c:timepost"] <- 0.2

power_SP_s2 <- powerSim(
  model_SP_s2,
  test = fixed("SP_c:timepost","t"),
  nsim = 100
)

print(power_SP_s2)

## NRS
model_NRS_s2 <- lmer(
  LSASSR ~ NRS_c * time + (1|subject),
  data = df4_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_NRS_s2)["NRS_c:timepost"] <- 0.2

power_NRS_s2 <- powerSim(
  model_NRS_s2,
  test = fixed("NRS_c:timepost","t"),
  nsim = 100
)

print(power_NRS_s2)

## apnea
model_apnea_s2 <- lmer(
  LSASSR ~ apnea_c * time + (1|subject),
  data = df4_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_apnea_s2)["apnea_c:timepost"] <- 0.2

power_apnea_s2 <- powerSim(
  model_apnea_s2,
  test = fixed("apnea_c:timepost","t"),
  nsim = 100
)

print(power_apnea_s2)

## sleepiness
model_sleepiness_s2 <- lmer(
  LSASSR ~ sleepiness_c * time + (1|subject),
  data = df4_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_sleepiness_s2)["sleepiness_c:timepost"] <- 0.2

power_sleepiness_s2 <- powerSim(
  model_sleepiness_s2,
  test = fixed("sleepiness_c:timepost","t"),
  nsim = 100
)

print(power_sleepiness_s2)

## KSQ_total
model_KSQ_s2 <- lmer(
  LSASSR ~ KSQ_total_c * time + (1|subject),
  data = df4_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_KSQ_s2)["KSQ_total_c:timepost"] <- 0.2

power_KSQ_s2 <- powerSim(
  model_KSQ_s2,
  test = fixed("KSQ_total_c:timepost","t"),
  nsim = 100
)

print(power_KSQ_s2)

## KSS
model_KSS_s2 <- lmer(
  LSASSR ~ KSS_pre_c * time + (1|subject),
  data = df4_long,
  REML = FALSE,
  control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_KSS_s2)["KSS_pre_c:timepost"] <- 0.2

power_KSS_s2 <- powerSim(
  model_KSS_s2,
  test = fixed("KSS_pre_c:timepost","t"),
  nsim = 100
)

print(power_KSS_s2)

#### Q3: Differences across treatments (ICBT vs SSRI) ####
## Differences across treatments (ICBT vs SSRI) - NO covariates

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
library(simr)

#### Data preparation ####
df5 <- combined_treat

df5 <- df5 %>%
  mutate(
    subject = factor(subject),
    treatment = factor(treatment, levels = c("ICBT", "SSRI"))
  )

# Missing values check
print(colSums(is.na(df5)))

## Centering / standardization (NO age, NO sex)
cont_vars <- c(
  "SP",
  "NRS",
  "apnea",
  "sleepiness",
  "KSQ_total",
  "KSS_pre",
  "LSAS_pre"
)

for (v in cont_vars) {
  df5[[paste0(v, "_c")]] <- as.numeric(scale(df5[[v]], center = TRUE, scale = FALSE))
  df5[[paste0(v, "_z")]] <- as.numeric(scale(df5[[v]], center = TRUE, scale = TRUE))
}

## Long format for LMM
df5_long <- df5 %>%
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

print(head(df5_long))
print(table(df5_long$time, df5_long$treatment))

## Descriptive statistics
desc_table <- df5_long %>%
  group_by(treatment, time) %>%
  summarise(
    n = sum(!is.na(LSAS)),
    mean_LSAS = mean(LSAS, na.rm = TRUE),
    sd_LSAS = sd(LSAS, na.rm = TRUE),
    .groups = "drop"
  )

print(desc_table)

#### Main LMM (NO age, NO sex) ####
## LSAS ~ time * treatment * sleep_predictor + (1 | subject)
run_lmm <- function(data_long, sleep_var_z) {
  
  formula_txt <- paste0(
    "LSAS ~ time * treatment * ", sleep_var_z,
    " + (1 | subject)"
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

model_list <- lapply(sleep_predictors, function(x) run_lmm(df5_long, x))
names(model_list) <- sleep_predictors

## Inspect one model
print(summary(model_list$SP_z))
print(anova(model_list$SP_z))
check_model(model_list$SP_z)

model_fit_info <- bind_rows(
  lapply(names(model_list), function(pred) {
    
    model <- model_list[[pred]]
    
    r2 <- as.data.frame(performance::r2_nakagawa(model))
    
    tibble(
      predictor = pred,
      AIC = AIC(model),
      BIC = BIC(model),
      logLik = as.numeric(logLik(model)),
      R2_marginal = r2$R2_marginal[1],
      R2_conditional = r2$R2_conditional[1]
    )
    
  })
)

print(model_fit_info)

## Extract fixed effects
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
write.csv(all_fixed_results, "all_fixed_effects_results_noCovars.csv", row.names = FALSE)
write.csv(three_way_results, "three_way_interaction_results_noCovars.csv", row.names = FALSE)
write.csv(model_fit_info, "model_fit_information_with_R2_noCovars.csv", row.names = FALSE)

#### Three-way interaction statistical power (NO covariates) ####

## SP
m <- model_list$SP_z
set.seed(123)
power_SP <- powerSim(
  m,
  fixed("timepost:treatmentSSRI:SP_z"),
  nsim = 100
)
print(power_SP)

## NRS
m <- model_list$NRS_z
set.seed(123)
power_NRS <- powerSim(
  m,
  fixed("timepost:treatmentSSRI:NRS_z"),
  nsim = 100
)
print(power_NRS)

## apnea
m <- model_list$apnea_z
set.seed(123)
power_apnea <- powerSim(
  m,
  fixed("timepost:treatmentSSRI:apnea_z"),
  nsim = 100
)
print(power_apnea)

## sleepiness
m <- model_list$sleepiness_z
set.seed(123)
power_sleepiness <- powerSim(
  m,
  fixed("timepost:treatmentSSRI:sleepiness_z"),
  nsim = 100
)
print(power_sleepiness)

## KSQ_total
m <- model_list$KSQ_total_z
set.seed(123)
power_KSQ_total <- powerSim(
  m,
  fixed("timepost:treatmentSSRI:KSQ_total_z"),
  nsim = 100
)
print(power_KSQ_total)

## KSS_pre
m <- model_list$KSS_pre_z
set.seed(123)
power_KSS_pre <- powerSim(
  m,
  fixed("timepost:treatmentSSRI:KSS_pre_z"),
  nsim = 100
)
print(power_KSS_pre)

#### Minimum Detectable Effect (MDE) - NO covariates ####

## Collect all MDE results
mde_all <- bind_rows(lapply(names(model_list), function(pred) {
  m <- model_list[[pred]]
  coef_table <- summary(m)$coefficients
  row_name <- paste0("timepost:treatmentSSRI:", pred)
  
  # Check if the interaction term exists
  if(row_name %in% rownames(coef_table)) {
    est <- coef_table[row_name, "Estimate"]
    se  <- coef_table[row_name, "Std. Error"]
    
    data.frame(
      predictor = pred,
      interaction = row_name,
      estimate = est,
      SE = se,
      MDE_95 = 1.96 * se,
      MDE_80 = 0.84 * se
    )
  } else {
    data.frame(
      predictor = pred,
      interaction = row_name,
      estimate = NA,
      SE = NA,
      MDE_95 = NA,
      MDE_80 = NA
    )
  }
}))

print(mde_all)

# Save MDE results
write.csv(mde_all, "MDE_three_way_interaction_noCovars.csv", row.names = FALSE)







