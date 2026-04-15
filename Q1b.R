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

## Volcano plot for univariate Firth results
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


#### Bayesian multivariable logistic regression ####
## Model: group ~ all 36 ROIs + age + sex
## Weakly informative priors are used for shrinkage
#options(download.file.method = "libcurl")
#options(repos = c(CRAN = "https://cran.r-project.org"))
#install.packages("posterior", dependencies = TRUE)
#install.packages("brms", dependencies = TRUE)
#install.packages("rstantools", dependencies = TRUE)

library(brms)

formula_bayes <- as.formula(
  paste("group ~ age + sex +", paste(roi_vars, collapse = " + "))
)

priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 2.5), class = "Intercept")
)

fit_bayes <- brm(
  formula = formula_bayes,
  data = df_analysis,
  family = bernoulli(link = "logit"),
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

summary(fit_bayes)


## Extract predicted probabilities from Bayesian model
pred_prob <- fitted(fit_bayes, scale = "response")[, "Estimate"]

pred_df <- data.frame(
  observed = df_analysis$group,
  predicted = pred_prob
)

head(pred_df)


## ROC curve and AUC
library(pROC)

roc_obj <- roc(response = pred_df$observed, predictor = pred_df$predicted)
auc_val <- as.numeric(auc(roc_obj))
auc_ci <- ci.auc(roc_obj)

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(linewidth = 1.2, color = "#1F78B4") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  coord_equal() +
  annotate(
    "label",
    x = 0.65, y = 0.15,
    label = paste0(
      "AUC = ", sprintf("%.3f", auc_val),
      "\n95% CI: ", sprintf("%.3f", auc_ci[1]), "–", sprintf("%.3f", auc_ci[3])
    ),
    size = 4
  ) +
  labs(
    title = "ROC curve for Bayesian multivariable model",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

print(p_roc)

ggsave("ROC_Bayesian_multivariable_adjusted.pdf", p_roc, width = 5.5, height = 5)


## Calibration plot
## Compare predicted probabilities with observed event rates
calibration_df <- pred_df %>%
  mutate(decile = ntile(predicted, 10)) %>%
  group_by(decile) %>%
  summarise(
    mean_pred = mean(predicted),
    obs_rate = mean(observed),
    n = n(),
    .groups = "drop"
  )

p_cal <- ggplot(calibration_df, aes(x = mean_pred, y = obs_rate)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  xlim(0, 1) + ylim(0, 1) +
  labs(
    title = "Calibration plot",
    x = "Mean predicted probability",
    y = "Observed event rate"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

print(p_cal)

ggsave("Calibration_Bayesian_multivariable_adjusted.pdf", p_cal, width = 5.5, height = 5)


## Brier score
## Lower values indicate better probabilistic prediction
brier_score <- mean((pred_df$observed - pred_df$predicted)^2)
cat("Brier score =", round(brier_score, 4), "\n")

## 1×10 repeated stratified cross-validation
## Internal validation of model performance

library(caret)

set.seed(1234)

## Create 10 stratified folds
## returnTrain = TRUE means each element contains training indices
folds <- createFolds(df_analysis$group, k = 10, returnTrain = TRUE)

cv_auc <- numeric(length(folds))
cv_brier <- numeric(length(folds))

for (i in seq_along(folds)) {
  
  ## createFolds(..., returnTrain = TRUE) returns training indices
  train_idx <- folds[[i]]
  
  ## define test indices as remaining observations
  test_idx <- setdiff(seq_len(nrow(df_analysis)), train_idx)
  
  train_data <- df_analysis[train_idx, ]
  test_data  <- df_analysis[test_idx, ]
  
  ## Fit Bayesian logistic regression on training data
  fit_cv <- brm(
    formula = formula_bayes,
    data = train_data,
    family = bernoulli(link = "logit"),
    prior = priors,
    chains = 2,
    iter = 2000,
    warmup = 1000,
    cores = 2,
    seed = 100 + i,
    refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  
  ## Predict probabilities for test set
  test_prob <- fitted(fit_cv, newdata = test_data, scale = "response")[, "Estimate"]
  
  ## Compute AUC
  roc_cv <- roc(test_data$group, test_prob, quiet = TRUE)
  cv_auc[i] <- as.numeric(auc(roc_cv))
  
  ## Compute Brier score
  cv_brier[i] <- mean((test_data$group - test_prob)^2)
}

## Report 10-fold cross-validated performance
cat(
  "10-fold CV AUC: ",
  round(mean(cv_auc), 3),
  " ± ",
  round(sd(cv_auc), 3),
  "\n",
  sep = ""
)

cat(
  "10-fold CV Brier: ",
  round(mean(cv_brier), 4),
  " ± ",
  round(sd(cv_brier), 4),
  "\n",
  sep = ""
)

## Export results
write.xlsx(
  list(
    Firth_univariate_results = firth_results,
    Bayesian_predictions = pred_df,
    Calibration_data = calibration_df,
    CV_AUC = data.frame(Fold = 1:length(cv_auc), AUC = cv_auc),
    CV_Brier = data.frame(Fold = 1:length(cv_brier), Brier = cv_brier)
  ),
  file = "Analysis_PlanA_Firth_Bayesian_adjusted.xlsx",
  overwrite = TRUE
)
getwd()
