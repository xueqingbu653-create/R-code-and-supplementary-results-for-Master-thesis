## Q1b_brain variables predict group membership

## Univariate Firth + Bayesian Multivariable Prediction
## Objective: Evaluate whether 36 ROI BOLD variability measures
## can discriminate SAD from HC, and assess their joint predictive value

df <- sample1

df$group <- as.numeric(as.character(df$group))

# Covariates
df$sex <- factor(df$sex)

# 19 ROI variables
roi_vars <- c(
  "Amygdala_02","Caudata_02","Cingulate_Ant_02","Cingulate_Post_02",
  "Frontal_Mid_02","Frontal_Sup_02",
  "Frontal_Sup_Medial_02","Hypothalamus_02","Insula_02",
  "LC_02","N_Acc_02","OFCmed_02","Pallidum_02",
  "Precuneus_02","Putamen_02","Raphe_D_02","Raphe_M_02",
  "Thalamus_02","VTA_02"
)

# Keep only variables needed for analysis
df_analysis <- df[, c("group", "age", "sex", roi_vars)]

# Remove rows with missing values
df_analysis <- na.omit(df_analysis)

# Standardize all ROI variables and age
to_scale <- c("age", roi_vars)
df_analysis[, to_scale] <- scale(df_analysis[, to_scale])

# Check final sample size
cat("Final N =", nrow(df_analysis), "\n")
cat("Number of ROI =", length(roi_vars), "\n")


#### Univariate Firth logistic regression ####
## Model: group ~ ROI + age + sex
library(logistf)
library(dplyr)
library(openxlsx)

run_firth_univariate <- function(var_list, data, outcome = "group", covariates = c("age", "sex")) {
  
  results_list <- lapply(var_list, function(v) {
    
    rhs <- paste(c(v, covariates), collapse = " + ")
    form <- as.formula(paste(outcome, "~", rhs))
    
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
  outcome = "group",
  covariates = c("age", "sex")
)

print(firth_results)

write.xlsx(
  firth_results,
  file = "Firth_univariate_results.xlsx",
  overwrite = TRUE
)


#### Power analysis for univariate Firth logistic regression ####
## Using standard glm (with same predictors) for power calculation
## Small effect assumed: OR = 1.5 (log = 0.405)

## post hoc statistical power
#### Power analysis for univariate Firth logistic regression ####
library(pwrss)
null_model <- glm(group ~ age + sex, data = df_analysis, family = binomial)
r2_covariates <- 1 - (null_model$deviance / null_model$null.deviance)

n_total <- nrow(df_analysis)
base_prob <- mean(df_analysis$group)

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

write.xlsx(power_results_by_roi, "PostHoc_Power_by_ROI.xlsx")

## Sevnsitive power analysis
# Search for minimum detectable OR by testing a range of values
or_range <- seq(1.1, 5.0, by = 0.1)
power_vals <- sapply(or_range, function(or) {
  pwrss.z.logistic(
    p0 = base_prob,              # Baseline probability
    odds.ratio = or,             # OR to test
    n = n_total,                 # Total sample size
    r2.other.x = r2_covariates,  # Variance explained by covariates
    alpha = 0.05,                # Significance level
    alternative = "not equal",   # Two-tailed test
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
    title = "Volcano plot of univariate Firth logistic regression",
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

ggsave("Volcano_ROI_Firth_adjusted.pdf", p_volcano, width = 6.5, height = 5.5)


#### FULL BAYESIAN MULTIVARIABLE LOGISTIC REGRESSION ####
## Model: group ~ all ROIs + age + sex
## Prior: Regularized Horseshoe (for shrinkage & variable selection)
install.packages("tidybayes")
library(brms)
library(bayesplot)
library(tidybayes)
library(pROC)
library(ggplot2)
library(openxlsx)

#### 1. Prepare Data and Formula ####

# Ensure data is ready
df_analysis <- na.omit(df_analysis)  # Already done, but safe check

# Formula with all predictors
formula_bayes <- as.formula(
  paste("group ~ age + sex +", paste(roi_vars, collapse = " + "))
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
cat("  - Intercept: Normal(0, 2.5)\n\n")

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
key_params_corrected <- c(
  "b_Intercept",                                
  "b_age",                                          
  "b_sexMale",                                    
  paste0("b_", roi_vars[1:5])                   
)

roi_params <- paste0("b_", roi_vars)
# ROI 1-6
trace_roi_1 <- mcmc_trace(fit_bayes_full, pars = roi_params[1:6]) +
  ggtitle("Trace Plots: ROIs 1-6")

ggsave("Trace_Plots_ROI_1-6.png", trace_roi_1, width = 12, height = 8)

# ROI 7-12
trace_roi_2 <- mcmc_trace(fit_bayes_full, pars = roi_params[7:12]) +
  ggtitle("Trace Plots: ROIs 7-12")

ggsave("Trace_Plots_ROI_7-12.png", trace_roi_2, width = 12, height = 8)

# ROI 13-19
trace_roi_3 <- mcmc_trace(fit_bayes_full, pars = roi_params[13:19]) +
  ggtitle("Trace Plots: ROIs 13-19")

ggsave("Trace_Plots_ROI_13-19.png", trace_roi_3, width = 12, height = 8)
cat("✓ Trace plots saved: Bayesian_Trace_Plots.png\n")


#### 5. Posterior Estimates ####

cat("\n========== Posterior Estimates ==========\n")

library(tibble)
fixed_effects <- summary(fit_bayes_full)$fixed
colnames(fixed_effects)


# Extract and transform estimates
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
    Prob_meaningful = ifelse(Parameter %in% c("b_Intercept", "b_age", "b_sexMale"),
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
write.xlsx(posterior_summary, "Bayesian_Posterior_Summary_Horseshoe.xlsx")

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
write.xlsx(variable_importance, "Bayesian_Variable_Importance_Horseshoe.xlsx")

#### 7. Posterior Predictive Checks ####

cat("\n========== Posterior Predictive Checks ==========\n")

# Generate posterior predictions
pred_bayes <- posterior_predict(fit_bayes_full, ndraws = 100)
observed <- df_analysis$group

# PPC plot
ppc_plot <- ppc_dens_overlay(observed, pred_bayes[1:50, ]) +
  ggtitle("Posterior Predictive Check - Density")

ggsave("Bayesian_PPC_Density.png", ppc_plot, width = 8, height = 6)

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
    title = "ROC Curve - Full Bayesian Model (Horseshoe Prior)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_classic(base_size = 13)

ggsave("Bayesian_ROC_Curve_Horseshoe.pdf", p_roc, width = 6, height = 6)
ggsave("Bayesian_ROC_Curve_Horseshoe.png", p_roc, width = 6, height = 6)


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
      Note = "Full Bayesian Model with Regularized Horseshoe Prior",
      Date = Sys.Date(),
      Prior_scale_global = 0.05,
      Prior_scale_slab = 2,
      Chains = 4,
      Iterations = 8000
    ),
    Posterior_Summary = posterior_summary,
    Variable_Importance = variable_importance,
    Model_Performance = final_results$Model_Performance,
    LOO_Comparison = final_results$LOO_Comparison,
    Predictions = final_results$Predictions
  ),
  file = "Full_Bayesian_Results_Horseshoe.xlsx",
  overwrite = TRUE
)

# Save R objects
saveRDS(fit_bayes_full, "Full_Bayesian_Model_Horseshoe.rds")
saveRDS(final_results, "Full_Bayesian_Results_Horseshoe.rds")

cat("\n========== Complete! ==========\n")
cat("All results saved:\n")
cat("  - Full_Bayesian_Results_Horseshoe.xlsx (main results)\n")
cat("  - Full_Bayesian_Model_Horseshoe.rds (model object)\n")
cat("  - Bayesian_ROC_Curve_Horseshoe.pdf\n")
cat("  - Bayesian_Forest_Plot_Horseshoe.png\n")
cat("  - Bayesian_Trace_Plots.png\n")
cat("  - Bayesian_PPC_Density.png\n")


