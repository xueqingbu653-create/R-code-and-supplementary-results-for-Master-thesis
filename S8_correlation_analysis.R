# Correlation analyses

#### SAD group ####

df <- sample1_SAD
sleep_vars <- c("SP","NRS","apnea","sleepiness", "KSQ_total", "KSS_pre")
brain_vars <- c(
  "Amygdala_02","Caudata_02","Cingulate_Ant_02","Cingulate_Post_02",
  "Frontal_Mid_02","Frontal_Sup_02", "Frontal_Sup_Medial_02",
  "Hypothalamus_02","Insula_02",
  "LC_02","N_Acc_02","OFCmed_02","Pallidum_02",
  "Precuneus_02","Putamen_02","Raphe_D_02","Raphe_M_02",
  "Thalamus_02","VTA_02"
)

#Normality test
library(purrr)
library(dplyr)

all_vars <- c(sleep_vars, brain_vars)

shapiro_results <- map_dfr(all_vars, function(v){
  test <- shapiro.test(df[[v]])
  
  data.frame(
    variable = v,
    W = test$statistic,
    p_value = test$p.value
  )
})

shapiro_results

#Spearman correlation 
cor_results <- expand.grid(
  sleep = sleep_vars,
  brain = brain_vars
) %>%
  rowwise() %>%
  mutate(
    cor_test = list(cor.test(
      df[[sleep]],
      df[[brain]],
      method = "spearman"
    )),
    rho = unname(cor_test$estimate),
    p = cor_test$p.value
  ) %>%
  ungroup() %>%
  select(sleep, brain, rho, p)

#Multiple comparisons correction
cor_results <- cor_results %>%
  mutate(
    p_bonf = p.adjust(p, method = "bonferroni"),
    p_fdr  = p.adjust(p, method = "fdr"),
    p_fwe  = p.adjust(p, method = "holm")
  )

library(writexl)
write_xlsx(cor_results, "sleep_brain_SAD_correlations.xlsx")

#### full sample ####
df1 <- sample1
all_vars <- c(sleep_vars, brain_vars)

shapiro_results1 <- map_dfr(all_vars, function(v){
  test <- shapiro.test(df1[[v]])
  
  data.frame(
    variable = v,
    W = test$statistic,
    p_value = test$p.value
  )
})

shapiro_results1

#Spearman correlation 
cor_results1 <- expand.grid(
  sleep = sleep_vars,
  brain = brain_vars
) %>%
  rowwise() %>%
  mutate(
    cor_test = list(cor.test(
      df1[[sleep]],
      df1[[brain]],
      method = "spearman"
    )),
    rho = unname(cor_test$estimate),
    p = cor_test$p.value
  ) %>%
  ungroup() %>%
  select(sleep, brain, rho, p)

#Multiple comparisons correction
cor_results1 <- cor_results1 %>%
  mutate(
    p_bonf = p.adjust(p, method = "bonferroni"),
    p_fdr  = p.adjust(p, method = "fdr"),
    p_fwe  = p.adjust(p, method = "holm")
  )

library(writexl)
write_xlsx(cor_results1, "sleep_brain_full_correlations.xlsx")

### Brain-Brain correlations

brain_pairs <- combn(brain_vars, 2, simplify = FALSE)

brain_cor_results <- map_dfr(brain_pairs, function(pair){
  
  var1 <- pair[1]
  var2 <- pair[2]
  
  # Pearson
  pearson_test <- cor.test(
    df1[[var1]],
    df1[[var2]],
    method = "pearson"
  )
  
  # Spearman
  spearman_test <- cor.test(
    df1[[var1]],
    df1[[var2]],
    method = "spearman"
  )
  
  data.frame(
    brain1 = var1,
    brain2 = var2,
    
    # Pearson
    pearson_r = unname(pearson_test$estimate),
    pearson_p = pearson_test$p.value,
    
    # Spearman
    spearman_rho = unname(spearman_test$estimate),
    spearman_p = spearman_test$p.value
  )
})

# Multiple comparison correction

brain_cor_results <- brain_cor_results %>%
  mutate(
    
    # Pearson corrections
    pearson_bonf = p.adjust(pearson_p, method = "bonferroni"),
    pearson_fdr  = p.adjust(pearson_p, method = "fdr"),
    pearson_holm = p.adjust(pearson_p, method = "holm"),
    
    # Spearman corrections
    spearman_bonf = p.adjust(spearman_p, method = "bonferroni"),
    spearman_fdr  = p.adjust(spearman_p, method = "fdr"),
    spearman_holm = p.adjust(spearman_p, method = "holm")
  )

brain_cor_results
write.csv(
  brain_cor_results,
  "brain_brain_correlations.csv",
  row.names = FALSE
)

#### HC group ####
df2 <- sample1_HC

all_vars <- c(sleep_vars, brain_vars)

shapiro_results2 <- map_dfr(all_vars, function(v){
  test <- shapiro.test(df2[[v]])
  
  data.frame(
    variable = v,
    W = test$statistic,
    p_value = test$p.value
  )
})

shapiro_results2

#Spearman correlation 
cor_results2 <- expand.grid(
  sleep = sleep_vars,
  brain = brain_vars
) %>%
  rowwise() %>%
  mutate(
    cor_test = list(cor.test(
      df2[[sleep]],
      df2[[brain]],
      method = "spearman"
    )),
    rho = unname(cor_test$estimate),
    p = cor_test$p.value
  ) %>%
  ungroup() %>%
  select(sleep, brain, rho, p)

#Multiple comparisons correction
cor_results2 <- cor_results2 %>%
  mutate(
    p_bonf = p.adjust(p, method = "bonferroni"),
    p_fdr  = p.adjust(p, method = "fdr"),
    p_fwe  = p.adjust(p, method = "holm")
  )

write_xlsx(cor_results2, "sleep_brain_HC_correlations.xlsx")





