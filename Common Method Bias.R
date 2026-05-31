#### Common Method Bias Assessment (Marker Variable Technique) ####
# Using sample1 as df

library(dplyr)

df <- sample1

# 1. Prepare data (baseline variables only)
df_cmb <- df %>%
  select(
    age,                  
    SP, NRS, apnea, sleepiness, KSQ_total, KSS_pre, 
    isi1,                 
    LSAS_pre              
  ) %>%
  na.omit() 

# If isi1 is not available in sample1, use this instead:
# df_cmb <- df %>% select(age, SP, NRS, apnea, sleepiness, KSQ_total, KSS_pre, LSAS_pre) %>% na.omit()

# Check data overview
cat("Sample size:", nrow(df_cmb), "\n")
cat("Variables included:", paste(names(df_cmb), collapse = ", "), "\n\n")

# 2. Compute correlation matrix
cor_matrix <- cor(df_cmb, use = "pairwise.complete.obs")

# 3. Extract correlations between marker variable (age) and other variables
marker_cors <- cor_matrix["age", setdiff(names(df_cmb), "age")]

# 4. Compute mean absolute correlation
mean_abs_cor <- mean(abs(marker_cors), na.rm = TRUE)

# 5. Output results
cat("========== Marker Variable Technique Results ==========\n")
cat("Marker variable: age\n")
cat("Core variables:", paste(setdiff(names(df_cmb), "age"), collapse = ", "), "\n\n")
cat("Correlations between age and each variable:\n")
print(round(marker_cors, 3))
cat("\nMean absolute correlation (|r̄|) =", round(mean_abs_cor, 3), "\n\n")

# 6. Interpretation criteria
if (mean_abs_cor < 0.2) {
  cat("✅ Conclusion: Common method bias is unlikely to be a major concern (mean |r| < 0.20)\n")
} else if (mean_abs_cor < 0.3) {
  cat("⚠️ Conclusion: Some common method bias may be present, but the magnitude is modest\n")
} else {
  cat("❌ Conclusion: Common method bias could be a concern; results should be interpreted with caution\n")
}

#### Sample 2 ####
df <- sample2

# Prepare data 
df_cmb <- df %>%
  select(
    age,
    SP, NRS, apnea, sleepiness, KSQ_total, KSS_pre,
    `LSASSR_tot_pre`,    
    LSAS_pre
  )  

cat("Sample size:", nrow(df_cmb), "\n")  
cat("Variables:", paste(names(df_cmb), collapse = ", "), "\n\n")

# Correlation matrix (pairwise complete 
cor_matrix <- cor(df_cmb, use = "pairwise.complete.obs")

# Marker variable correlations
marker_cors <- cor_matrix["age", c("SP", "NRS", "apnea", "sleepiness", "KSQ_total", "KSS_pre", "LSASSR_tot_pre", "LSAS_pre")]

# Mean absolute correlation
mean_abs_cor <- mean(abs(marker_cors), na.rm = TRUE)

# Results
cat("========== Marker Variable Technique Results ==========\n")
cat("Marker variable: age\n\n")
cat("Correlations between age and each variable:\n")
print(round(marker_cors, 3))
cat("\nMean absolute correlation (|r̄|) =", round(mean_abs_cor, 3), "\n\n")

# Interpretation
if (mean_abs_cor < 0.2) {
  cat("✅ Conclusion: Common method bias is unlikely to be a major concern (mean |r| < 0.20)\n")
} else if (mean_abs_cor < 0.3) {
  cat("⚠️ Conclusion: Some common method bias may be present, but the magnitude is modest\n")
} else {
  cat("❌ Conclusion: Common method bias could be a concern; results should be interpreted with caution\n")
}


