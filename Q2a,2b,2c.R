## LMM (Q2a, 2b, 2c Sample 1, Sample 2)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(writexl)

#### Sample 1 ####
df <- sample1_SAD

## reshape to long format
df_long <- df %>% 
  pivot_longer(
    cols = c("LSAS_pre","LSAS_post"),
    names_to = "time",
    names_prefix = "LSAS_",
    values_to = "LSAS"
  )

## factor
df_long$time <- factor(df_long$time, levels = c("pre","post"))
df_long$sex  <- factor(df_long$sex, levels = c("Female","Male"))

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

## standardize variables
to_center <- c("age", all_predictors)

for (v in to_center){
  new_name <- paste0(v,"_c")
  df_long[[new_name]] <- scale(df_long[[v]], center=TRUE, scale=TRUE)
}

sleep_vars_c <- paste0(sleep_vars,"_c")
brain_vars_c <- paste0(brain_vars,"_c")


## Function to run LMM and extract main + interaction

run_lmm <- function(var_c, data){
  
  formula_txt <- paste0(
    "LSAS ~ ", var_c," * time + age_c + sex + (1|subject)"
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
  path = "LMM_sample1_main_and_interaction.xlsx"
)
getwd()

#### Sample 2 ####
df1 <- sample2_SAD

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
df1_long <- df1 %>%
  pivot_longer(
    cols = c("LSASSR_pre", "LSASSR_post"),
    names_to = "time",
    names_prefix = "LSASSR_",
    values_to = "LSASSR"
  )

## factor
df1_long$time <- factor(df1_long$time, levels = c("pre", "post"))
df1_long$sex  <- factor(df1_long$sex)

## standardize variables
to_center_s2 <- c("age", all_predictors_s2)

for (v in to_center_s2) {
  new_name <- paste0(v, "_c")
  df1_long[[new_name]] <- as.numeric(scale(df1_long[[v]], center = TRUE, scale = TRUE))
}

sleep_vars_c_s2 <- paste0(sleep_vars_s2, "_c")
brain_vars_c_s2 <- paste0(brain_vars_s2, "_c")

## function to run LMM and extract main predictor + main time + interaction
run_lmm_s2 <- function(var_c, data) {
  
  formula_txt <- paste0(
    "LSASSR ~ ", var_c, " * time + age_c + sex + (1|subject)"
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
  data = df1_long
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
  data = df1_long
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
  path = "LMM_sample2_main_time_interaction_with_R2.xlsx"
)


