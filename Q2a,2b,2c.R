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

#### Statistical power ####
#### sleep variables

library(lme4)
library(simr)

## LMM_SP
model_SP <- lmer(
  LSAS ~ SP_c * time + age_c + sex + (1|subject),
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
  LSAS ~ NRS_c * time + age_c + sex + (1|subject),
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
  LSAS ~ apnea_c * time + age_c + sex + (1|subject),
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
  LSAS ~ sleepiness_c * time + age_c + sex + (1|subject),
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
  LSAS ~ KSQ_total_c * time + age_c + sex + (1|subject),
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
  LSAS ~ KSS_pre_c * time + age_c + sex + (1|subject),
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

#### Brain variables

## Amygdala
model_Amygdala <- lmer(
  LSAS ~ Amygdala_02_c * time + age_c + sex + (1|subject),
  data=df_long, REML=FALSE,
  control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

fixef(model_Amygdala)["Amygdala_02_c:timepost"] <- 0.2

power_Amygdala <- powerSim(model_Amygdala, test=fixed("Amygdala_02_c:timepost","t"), nsim=100)
print(power_Amygdala)

## Caudata
model_Caudata <- lmer(LSAS ~ Caudata_02_c*time + age_c + sex + (1|subject),
                      data=df_long, REML=FALSE,
                      control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_Caudata)["Caudata_02_c:timepost"] <- 0.2
power_Caudata <- powerSim(model_Caudata, test=fixed("Caudata_02_c:timepost","t"), nsim=100)
print(power_Caudata)

## Cingulate_Ant
model_CingA <- lmer(LSAS ~ Cingulate_Ant_02_c*time + age_c + sex + (1|subject),
                    df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_CingA)["Cingulate_Ant_02_c:timepost"] <- 0.2
power_CingA <- powerSim(model_CingA, test=fixed("Cingulate_Ant_02_c:timepost","t"), nsim=100)
print(power_CingA)

## Cingulate_Post
model_CingP <- lmer(LSAS ~ Cingulate_Post_02_c*time + age_c + sex + (1|subject),
                    df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_CingP)["Cingulate_Post_02_c:timepost"] <- 0.2
power_CingP <- powerSim(model_CingP, test=fixed("Cingulate_Post_02_c:timepost","t"), nsim=100)
print(power_CingP)

## Frontal_Mid
model_FM <- lmer(LSAS ~ Frontal_Mid_02_c*time + age_c + sex + (1|subject),
                 df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_FM)["Frontal_Mid_02_c:timepost"] <- 0.2
power_FM <- powerSim(model_FM, test=fixed("Frontal_Mid_02_c:timepost","t"), nsim=100)
print(power_FM)

## Frontal_Sup
model_FS <- lmer(LSAS ~ Frontal_Sup_02_c*time + age_c + sex + (1|subject),
                 df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_FS)["Frontal_Sup_02_c:timepost"] <- 0.2
power_FS <- powerSim(model_FS, test=fixed("Frontal_Sup_02_c:timepost","t"), nsim=100)
print(power_FS)

## Frontal_Sup_Medial
model_FSM <- lmer(LSAS ~ Frontal_Sup_Medial_02_c*time + age_c + sex + (1|subject),
df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_FSM)["Frontal_Sup_Medial_02_c:timepost"] <- 0.2
power_FSM <- powerSim(model_FSM, test=fixed("Frontal_Sup_Medial_02_c:timepost","t"), nsim=100)
print(power_FSM)

## Hypothalamus
model_Hyp <- lmer(LSAS ~ Hypothalamus_02_c*time + age_c + sex + (1|subject),
                  df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_Hyp)["Hypothalamus_02_c:timepost"] <- 0.2
power_Hyp <- powerSim(model_Hyp, test=fixed("Hypothalamus_02_c:timepost","t"), nsim=100)
print(power_Hyp)

## Insula
model_Insula <- lmer(LSAS ~ Insula_02_c*time + age_c + sex + (1|subject),
                     df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_Insula)["Insula_02_c:timepost"] <- 0.2
power_Insula <- powerSim(model_Insula, test=fixed("Insula_02_c:timepost","t"), nsim=100)
print(power_Insula)

## LC
model_LC <- lmer(LSAS ~ LC_02_c*time + age_c + sex + (1|subject),
                 df_long, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fixef(model_LC)["LC_02_c:timepost"] <- 0.2
power_LC <- powerSim(model_LC, test=fixed("LC_02_c:timepost","t"), nsim=100)
print(power_LC)

## N_Acc
model_NAcc <- lmer(
  LSAS ~ N_Acc_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ OFCmed_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ Pallidum_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ Precuneus_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ Putamen_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ Raphe_D_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ Raphe_M_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ Thalamus_02_c * time + age_c + sex + (1|subject),
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
  LSAS ~ VTA_02_c * time + age_c + sex + (1|subject),
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

#### Statistical power for sample 2 #####
## SP
model_SP_s2 <- lmer(
  LSASSR ~ SP_c * time + age_c + sex + (1|subject),
  data = df1_long,
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
  LSASSR ~ NRS_c * time + age_c + sex + (1|subject),
  data = df1_long,
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
  LSASSR ~ apnea_c * time + age_c + sex + (1|subject),
  data = df1_long,
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
  LSASSR ~ sleepiness_c * time + age_c + sex + (1|subject),
  data = df1_long,
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
  LSASSR ~ KSQ_total_c * time + age_c + sex + (1|subject),
  data = df1_long,
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
  LSASSR ~ KSS_pre_c * time + age_c + sex + (1|subject),
  data = df1_long,
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

