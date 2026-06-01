### Partial Least Squares Discriminant Analysis (PLS-DA) Results for Brain Variables

df <- sample1
library(pls)
library(caret)
library(pROC)
roi_vars <- c(
  "Amygdala_02","Caudata_02","Cingulate_Ant_02","Cingulate_Post_02",
  "Frontal_Mid_02","Frontal_Sup_02",
  "Frontal_Sup_Medial_02","Hypothalamus_02","Insula_02",
  "LC_02","N_Acc_02","OFCmed_02","Pallidum_02",
  "Precuneus_02","Putamen_02","Raphe_D_02","Raphe_M_02",
  "Thalamus_02","VTA_02"
)

# 1. PREPARE DATA
# Assuming your dataframe is named 'df'

# Convert sex to numeric: Male = 1, Female = 0
df$sex_num <- ifelse(df$sex == "Male", 1, 0)

# Create predictor matrix (19 ROIs + age + sex_num)
X <- df[, c(roi_vars, "age", "sex_num")]

# Create outcome variable (group: 0 or 1)
Y <- df$group

# 2. SPLIT DATA INTO TRAINING AND TEST SETS
set.seed(123)  # For reproducibility

# 70% training, 30% testing
train_index <- createDataPartition(Y, p = 0.7, list = FALSE)

X_train_raw <- X[train_index, ]
X_test_raw <- X[-train_index, ]
Y_train <- Y[train_index]
Y_test <- Y[-train_index]

# 3. STANDARDIZE PREDICTORS (scale using training set parameters)
# Important: Use mean and sd from training set to avoid data leakage
train_mean <- colMeans(X_train_raw)
train_sd <- apply(X_train_raw, 2, sd)

X_train <- scale(X_train_raw, center = train_mean, scale = train_sd)
X_test <- scale(X_test_raw, center = train_mean, scale = train_sd)

# 4. PERFORM PLS-DA WITH CROSS-VALIDATION 
# Fit PLS model with up to 10 components and cross-validation
set.seed(123)
pls_model <- plsr(Y_train ~ X_train, 
                  ncomp = 10, 
                  validation = "CV",  # Cross-validation
                  segments = 10)      # 10-fold CV

# 5. SELECT OPTIMAL NUMBER OF COMPONENTS 
# View cross-validation results
summary(pls_model)

# Plot RMSEP (Root Mean Square Error of Prediction)
plot(RMSEP(pls_model), main = "RMSEP vs Number of Components")

# Extract RMSEP values
rmsep_values <- RMSEP(pls_model)$val[1, 1, ]
ncomp_opt <- which.min(rmsep_values)
cat("Optimal number of components:", ncomp_opt, "\n")

# Alternative: Use one-standard-error rule
rmsep_min <- min(rmsep_values)
rmsep_se <- RMSEP(pls_model)$se[1, 1, ]
ncomp_1se <- which(rmsep_values <= rmsep_min + rmsep_se[ncomp_opt])[1]
cat("Optimal components (1-SE rule):", ncomp_1se, "\n")

# Use the chosen optimal components
final_ncomp <- ncomp_opt  # or ncomp_1se

# 6. REFIT MODEL WITH OPTIMAL COMPONENTS
pls_final <- plsr(Y_train ~ X_train, ncomp = final_ncomp)

# 7. MAKE PREDICTIONS ON TEST SET
# Predict continuous values
pred_test_continuous <- predict(pls_final, newdata = X_test, ncomp = final_ncomp)
pred_test_continuous <- as.vector(pred_test_continuous)

# Convert to binary class predictions (threshold = 0.5)
pred_test_class <- ifelse(pred_test_continuous > 0.5, 1, 0)

# 8. EVALUATE MODEL PERFORMANCE
# Confusion matrix
conf_matrix <- confusionMatrix(as.factor(pred_test_class), as.factor(Y_test))
print(conf_matrix)

# Calculate metrics manually
accuracy <- mean(pred_test_class == Y_test)
sensitivity <- conf_matrix$byClass["Sensitivity"]  # True positive rate
specificity <- conf_matrix$byClass["Specificity"]  # True negative rate
precision <- conf_matrix$byClass["Precision"]
f1_score <- conf_matrix$byClass["F1"]

cat("\n=== Model Performance ===\n")
cat("Accuracy:", round(accuracy, 3), "\n")
cat("Sensitivity (Recall):", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")
cat("Precision:", round(precision, 3), "\n")
cat("F1 Score:", round(f1_score, 3), "\n")

# 9. ROC CURVE AND AUC 
roc_curve <- roc(Y_test, pred_test_continuous)
auc_value <- auc(roc_curve)

cat("AUC:", round(auc_value, 3), "\n")

# Plot ROC curve
library(ggplot2)
library(pROC)
optimal_fpr <- 1 - optimal_specificity
# Calculate ROC and confidence intervals
roc_obj_full <- roc(Y_test, pred_test_continuous)
auc_full <- auc(roc_obj_full)
auc_ci_full <- ci.auc(roc_obj_full, conf.level = 0.95)

# Create data frame for ggplot
roc_df <- data.frame(
  FPR = 1 - roc_obj_full$specificities,
  TPR = roc_obj_full$sensitivities
)

# Find optimal cutoff (Youden's index)
optimal_idx <- which.max(roc_obj_full$sensitivities + roc_obj_full$specificities - 1)
optimal_cutoff <- roc_obj_full$thresholds[optimal_idx]
optimal_sensitivity <- roc_obj_full$sensitivities[optimal_idx]
optimal_specificity <- roc_obj_full$specificities[optimal_idx]

# Professional ROC curve
p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  # Main ROC line
  geom_line(size = 1.3, color = "#1F78B4") +
  
  # Confidence interval ribbon (optional, if you have bootstrap)
  # geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper), alpha = 0.2, fill = "#1F78B4") +
  
  # Diagonal reference line (chance level)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "grey50", size = 0.8) +
  
  # Add optimal cutoff point
  geom_point(x = optimal_fpr, y = optimal_sensitivity,
             size = 5, color = "#E31A23", shape = 18, stroke = 1.5) +
  
  # Add annotation for optimal cutoff
  annotate("text", 
           x = optimal_fpr - 0.08,         
           y = optimal_sensitivity + 0.08,  
           label = paste0("Optimal cutoff = ", round(optimal_cutoff, 3)),
           size = 3.5, color = "#E31A23", fontface = "italic")  +
  
  # Add AUC annotation
  annotate("label",
           x = 0.65, y = 0.2,
           label = paste0(
             "AUC = ", sprintf("%.3f", auc_full),
             "\n95% CI: ", sprintf("%.3f", auc_ci_full[1]), 
             "–", sprintf("%.3f", auc_ci_full[3])
           ),
           size = 4.5, 
           fill = "white", 
           color = "black",
           label.size = 0.3,
           label.padding = unit(0.5, "lines")) +
  
  # Labels and title
  labs(
    title = "PLS-DA ROC Curve",
    subtitle = paste0("Predicting Group Membership from 19 ROI Volumes"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  
  # Set axis limits
  xlim(0, 1) +
  ylim(0, 1) +
  
  # Professional theme
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
  ) +
  
  # Add minor gridlines
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))

# Display the plot
print(p_roc)

# Save high-resolution versions
ggsave("PLS_Figures/Figure2_ROC_Curve_Professional.png", 
       p_roc, width = 7, height = 6.5, dpi = 300, bg = "white")

ggsave("PLS_Figures/Figure2_ROC_Curve_Professional.pdf", 
       p_roc, width = 7, height = 6.5, device = "pdf")

# Also save as TIFF for journal submission
tiff("PLS_Figures/Figure2_ROC_Curve_Professional.tiff", 
     width = 7, height = 6.5, units = "in", res = 300, compression = "lzw")
print(p_roc)
dev.off()

cat("  Saved: Figure2_ROC_Curve_Professional.png/pdf/tiff\n")


# 10. VARIABLE IMPORTANCE (VIP SCORES) - CORRECTED VERSION

# Function to calculate VIP scores (only need ONE definition)
calculate_vip <- function(pls_object, ncomp) {
  # Ensure ncomp is valid
  if (ncomp < 1) ncomp <- 1
  if (ncomp > ncol(pls_object$loading.weights)) {
    ncomp <- ncol(pls_object$loading.weights)
  }
  
  # Extract loading weights and Y loadings
  w <- pls_object$loading.weights[, 1:ncomp, drop = FALSE]
  q <- pls_object$Yloadings[, 1:ncomp, drop = FALSE]
  
  # Calculate VIP
  SS <- q^2 * colSums(w^2)
  VIP <- sqrt(nrow(w) * apply(SS, 1, sum) / sum(SS))
  
  return(VIP)
}

# Check if final_ncomp is valid
if (!exists("final_ncomp") || final_ncomp < 1) {
  final_ncomp <- pls_final$ncomp
  cat("Using model's ncomp:", final_ncomp, "\n")
}

# Calculate VIP scores
vip_scores <- calculate_vip(pls_final, final_ncomp)

# Check length before assigning names
cat("VIP scores length:", length(vip_scores), "\n")
cat("X_train columns:", ncol(X_train), "\n")

# Assign names (only if lengths match)
if (length(vip_scores) == ncol(X_train)) {
  names(vip_scores) <- colnames(X_train)
} else {
  # If lengths don't match, use the available variables
  cat("Warning: Length mismatch. Using available variables.\n")
  names(vip_scores) <- colnames(X_train)[1:length(vip_scores)]
}

# Sort by importance
vip_sorted <- sort(vip_scores, decreasing = TRUE)
cat("\n=== Top 10 Most Important Variables (VIP scores) ===\n")
print(round(vip_sorted[1:min(10, length(vip_sorted))], 3))

# Plot VIP scores
par(mar = c(10, 4, 4, 2))
barplot(vip_sorted, 
        las = 2, 
        col = ifelse(vip_sorted > 1, "steelblue", "lightgray"),
        main = paste("Variable Importance in Projection (VIP) -", final_ncomp, "components"),
        ylab = "VIP Score",
        cex.names = 0.7)
abline(h = 1, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("VIP > 1 (Important)", "VIP <= 1"), 
       fill = c("steelblue", "lightgray"), cex = 0.8)

# 11. PLOT MODEL COEFFICIENTS
# Extract regression coefficients
coefficients <- coef(pls_final, ncomp = final_ncomp, intercept = TRUE)
coef_values <- as.vector(coefficients[,,1])
names(coef_values) <- c("Intercept", colnames(X_train))

# Plot top positive and negative coefficients
coef_sorted <- sort(coef_values[-1], decreasing = TRUE)  # Exclude intercept
par(mar = c(8, 4, 4, 2))
barplot(head(coef_sorted, 10), 
        las = 2, 
        col = ifelse(head(coef_sorted, 10) > 0, "darkgreen", "darkred"),
        main = "Top 10 Coefficients (Excluding Intercept)",
        ylab = "Coefficient Value",
        cex.names = 0.8)
abline(h = 0, lwd = 2)

# 12. CROSS-VALIDATION PREDICTIONS ON TRAINING SET
# Get cross-validated predictions on training set
cv_predictions <- pls_model$validation$pred[, , final_ncomp]
cv_predictions <- as.vector(cv_predictions)
cv_class <- ifelse(cv_predictions > 0.5, 1, 0)

# Cross-validated performance
cv_accuracy <- mean(cv_class == Y_train)
cat("\n=== Cross-Validation Performance (Training Set) ===\n")
cat("CV Accuracy:", round(cv_accuracy, 3), "\n")

# 13. ADDITIONAL DIAGNOSTIC PLOTS
# Plot predicted vs actual
# Close any existing graphics devices
graphics.off()

# Open a new graphics device with larger size
dev.new(width = 10, height = 10)  # For R GUI
# Or use: windows(width = 10, height = 10)  # For Windows
# Or use: x11(width = 10, height = 10)      # For Linux

# Set up 2x2 plot layout
par(mfrow = c(2, 2))
par(mar = c(4, 4, 3, 2))  # Adjust margins: bottom, left, top, right

# Plot 1: Training set: Actual vs Predicted
plot(Y_train, as.vector(predict(pls_final, ncomp = final_ncomp)),
     main = "Training Set: Actual vs Predicted",
     xlab = "Actual Group", ylab = "Predicted Value",
     col = ifelse(Y_train == 1, "blue", "red"), pch = 19,
     xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
abline(h = 0.5, lty = 2, col = "gray")
abline(0, 1, lty = 1, col = "black")

# Plot 2: Test set: Actual vs Predicted
plot(Y_test, pred_test_continuous,
     main = "Test Set: Actual vs Predicted",
     xlab = "Actual Group", ylab = "Predicted Value",
     col = ifelse(Y_test == 1, "blue", "red"), pch = 19,
     xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
abline(h = 0.5, lty = 2, col = "gray")
abline(0, 1, lty = 1, col = "black")

# Plot 3: Histogram of predicted probabilities by group
hist(pred_test_continuous[Y_test == 0], 
     col = rgb(1, 0, 0, 0.5), 
     main = "Test Set: Predicted Probabilities",
     xlab = "Predicted Probability", 
     xlim = c(0, 1),
     ylim = c(0, max(table(cut(pred_test_continuous, breaks = 20)))),
     breaks = 20)
hist(pred_test_continuous[Y_test == 1], 
     col = rgb(0, 0, 1, 0.5), 
     add = TRUE,
     breaks = 20)
legend("topright", legend = c("Group 0", "Group 1"), 
       fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), cex = 0.8)

# Plot 4: RMSEP plot
plot(RMSEP(pls_model), main = "RMSEP vs Number of Components")
abline(v = final_ncomp, col = "red", lty = 2)

# 14. SAVE RESULTS (OPTIONAL)
results <- list(
  model = pls_final,
  optimal_components = final_ncomp,
  vip_scores = vip_scores,
  test_accuracy = accuracy,
  test_auc = auc_value,
  confusion_matrix = conf_matrix,
  roc_curve = roc_curve
)

# Save to RDS file
saveRDS(results, "pls_da_results.rds")

# Save VIP scores to CSV
vip_df <- data.frame(Variable = names(vip_scores), VIP = vip_scores)
vip_df <- vip_df[order(vip_df$VIP, decreasing = TRUE), ]
write.csv(vip_df, "vip_scores.csv", row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat("Results saved to 'pls_da_results.rds' and 'vip_scores.csv'\n")


