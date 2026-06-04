#backward selection procedure
#Name in old directory: 03_robustness_appendixD_best_secbest.R

rm(list = ls())

library('data.table')
library('dplyr')
library('clubSandwich')
library('stargazer')
library('Hmisc')
library('readstata13')
library('car')
library('miceadds')
library('multiwayvcov')
library('estimatr')
library('hdm')
library('lme4')
library('stringr')
library('pracma')
library('ggplot2')
library('hash')
library('Rfast')
library('ggsci')
library('scales')
library('this.path')

setwd(this.path::here())
source("ecma_directory.R")

datadir <- paste0(ecmadir,"Data/")
path_input_data <- paste0(datadir,"Input Data/")
path_prepared_data <- paste0(datadir,"Prepared Data/")
path_output_data <- paste0(datadir, "Output Data/")

path_figures <- paste0(ecmadir,"Figures/")
path_tables <- paste0(ecmadir,"Tables/")
path_functions <- paste0(ecmadir,"Code/Helper_functions/")

setwd(path_functions)


source("pooling_functions.R")

treatment_profiles <- list(c("random", "noIncentive", "noReminder"), #seeds only
                           c("trusted", "noIncentive", "noReminder"),
                           c("gossip", "noIncentive", "noReminder"),
                           
                           c("noSeed", "slope", "noReminder"), #incentives only
                           c("noSeed", "flat", "noReminder"),
                           
                           c("noSeed", "noIncentive", "SMS"), #SMS only
                           
                           c("random", "slope", "noReminder"), #seeds and incentives (slopes) only
                           c("trusted", "slope", "noReminder"),
                           c("gossip", "slope", "noReminder"),
                           
                           c("random", "flat", "noReminder"), #seeds and incentives (flats) only
                           c("trusted", "flat", "noReminder"),
                           c("gossip", "flat", "noReminder"),
                           
                           c("random", "noIncentive", "SMS"), #seeds and SMS only
                           c("trusted", "noIncentive", "SMS"),
                           c("gossip", "noIncentive", "SMS"),
                           
                           c("noSeed", "slope", "SMS"), #incentives (slope) and SMS only
                           c("noSeed", "flat", "SMS"), #incentives (flats) and SMS only
                           
                           c("random", "slope", "SMS"), #seeds and incentives (slopes) and SMS
                           c("trusted", "slope", "SMS"),
                           c("gossip", "slope", "SMS"),
                           
                           c("random", "flat", "SMS"), #seeds and incentives (flats) and SMS
                           c("trusted", "flat", "SMS"),
                           c("gossip", "flat", "SMS"))



#outcome = "shot_Measles1"
outcome = "shots_per_dollar"

#--READ VILLAGEXMONTH dataset

villagexmonth_level <- fread(paste0(path_prepared_data,"/Tablet_VillageXMonth_Costs.csv"),header = TRUE, sep = ",", data.table = FALSE)

#--SELECT SEEDS RISK EXPT
villagexmonth_level <- villagexmonth_level %>% filter(seedsrisk == 1)

#--SELECT ONLY FIRST IMPLEMENTATION
villagexmonth_level <- villagexmonth_level %>% filter(first_implementation == 1)


#--CREATE ALL POLICIES -- 
source('create_sp_variables.R')

#District-Time FE
villagexmonth_level$fes <- group_indices(villagexmonth_level, id_district, created_year, created_month)


#create FES dummies
fes_dummies <- data.frame(lme4::dummy(villagexmonth_level$fes))

villagexmonth_level <- cbind(villagexmonth_level, fes_dummies)

#-- LASSO SELECTION
variables <- grep("^SP", names(villagexmonth_level), value = TRUE)
variables <- variables[variables != "SP_noSeedXnoIncentiveXnoReminder"]

variables_expanded <- c(variables, colnames(fes_dummies))

current_variables <- variables_expanded

deselect_list <- c()
deselect_pval <- c()

pval_cutoff_list <- logspace(-13,-2,n = 12)

#one time OLS model to choose p

full_sp_formula <- as.formula(paste0(outcome,"~",paste0(variables_expanded,collapse = "+")))
full_model_ols <- estimatr::lm_robust(formula = full_sp_formula, data = villagexmonth_level, weights = village_population, se_type = "classical")

source("pval_lambda_mapping_functions.R")

pval_cutoff_lb <- 10^(-13)
pval_cutoff_ub <- 10^(-2)

lambda_extremes <- get_lambda_from_p(c(pval_cutoff_lb,pval_cutoff_ub))
lambda_grid <- linspace(lambda_extremes[2],lambda_extremes[1],20)

if (outcome == "shot_Measles1") {
  lambda_grid <- linspace(0.15,0.51,20)
} else if (outcome == "shots_per_dollar") {
  lambda_grid <- linspace(0.00045,0.0015,20)
}

pval_cutoff_list <- get_p_from_lambda(lambda_grid)

best_pol_list <- c()
best_pol_est_list <- c()
best_pol_CI_lower <- c()
best_pol_CI_upper <- c()

#Store best and second best estimates pre-WC adjustment 
best_pol_list_preWC <- c()
best_pol_est_list_preWC <- c()
best_pol_se_preWC <- c()
best_pol_p_preWC <- c()

sec_best_pol_list_preWC <- c()
sec_best_pol_est_list_preWC <- c()
sec_best_pol_se_preWC <- c()
sec_best_pol_p_preWC <- c()


source("map_key_policy_names.R")


for (pval_cutoff in pval_cutoff_list) {
  current_variables <- variables_expanded
  
  deselect_list <- c()
  deselect_pval <- c()
  
  #Initialize the BE loop
  current_sp_formula <- as.formula(paste0(outcome,"~",paste0(current_variables,collapse = "+")))
  current_model_ols <- estimatr::lm_robust(formula = current_sp_formula, data = villagexmonth_level, weights = village_population, se_type = "classical")
  current_max_pval <- max(current_model_ols$p.value)
  
  while (current_max_pval > pval_cutoff) {
    
    deselect_name <- names(current_model_ols$p.value[which.max(current_model_ols$p.value)])
    deselect_list <- c(deselect_list, deselect_name)
    deselect_pval <- c(deselect_pval, current_model_ols$p.value[which.max(current_model_ols$p.value)])
    
    current_variables <- current_variables[current_variables != deselect_name]
    current_sp_formula <- as.formula(paste0(outcome,"~",paste0(current_variables,collapse = "+")))
    current_model_ols <- estimatr::lm_robust(formula = current_sp_formula, data = villagexmonth_level, weights = village_population, se_type = "classical")
    current_max_pval <- max(current_model_ols$p.value)

  }
  
  support_SP <- variables_expanded[!(variables_expanded %in% deselect_list)]
  
  toRemove <- grep(pattern = "ntercept",x = support_SP) #don't want intercept to be lasso selected!
  if (length(toRemove) > 0) {
    support_SP <- support_SP[-toRemove]
  }
  
  fes_chosen <- grep("^X", support_SP, value = TRUE)
  
  #Final support policies
  support_SP_policies <- grep("^X", support_SP, value = TRUE, invert = TRUE)
  
  # ######################
  #-- AUTOMATED POOLING
  # ######################

  final_data <- villagexmonth_level

  for (tp_name in treatment_profiles) {
    tp_support <- get_relevant_sp_in_tp(support_SP, tp_name)
    final_data <- add_pooled_policies_tp(final_data, tp_support, tp_name)
  }

  pooled_policies <- grep("^POOLED_", colnames(final_data), value = TRUE)


  # ######################
  #-- POST LASSO
  # ######################

  variables_pooled_expanded <- c(pooled_policies,fes_chosen)

  formula_pl <- as.formula(paste0(outcome,"~",paste0(variables_pooled_expanded ,collapse = "+")))
  model_pl <- estimatr::lm_robust(formula = formula_pl, data = final_data, clusters = id_sc, weights = village_population, se_type = "CR0")

  pl_effects <- model_pl$coefficients[pooled_policies]
  pl_pval <-model_pl$p.value[pooled_policies]

  ranked_policy_names <- names(sort(pl_effects, decreasing = TRUE))
  pol_best_name <- ranked_policy_names[1]
  pol_2nd_name <- ranked_policy_names[2]

  #-- Store post LASSO results
  best_pol_list_preWC     <- c(best_pol_list_preWC, policy_name_mapping[[pol_best_name]])
  best_pol_est_list_preWC <- c(best_pol_est_list_preWC,max(pl_effects))
  best_pol_se_preWC       <- c(best_pol_se_preWC,model_pl$std.error[pooled_policies][pol_best_name])
  best_pol_p_preWC        <- c(best_pol_p_preWC,pl_pval[pol_best_name])

  sec_best_pol_list_preWC     <- c(sec_best_pol_list_preWC, policy_name_mapping[[pol_2nd_name]])
  sec_best_pol_est_list_preWC <- c(sec_best_pol_est_list_preWC,nth(pl_effects, 2, descending = T))#second best effect
  sec_best_pol_se_preWC       <- c(sec_best_pol_se_preWC,model_pl$std.error[pooled_policies][pol_2nd_name])
  sec_best_pol_p_preWC       <- c(sec_best_pol_p_preWC,pl_pval[pol_2nd_name])



  # #########################
  #-- ANDREWS WC ADJUSTMENT
  # #########################

  source('inference_on_winners_functions.R')

  alpha = 0.05
  beta = 0.005

  ntreat <- length(pooled_policies) + 1

  best_scaled_effect <- sqrt(nobs(model_pl)) * pl_effects[pol_best_name]
  trunc_scaled_effect <- max(0, sqrt(nobs(model_pl)) * pl_effects[pol_2nd_name]) #in case
  var_around_best <- nobs(model_pl) * (model_pl$std.error[pol_best_name])^2

  hybrid_results_scaled <- get_hybrid_Y_alpha_beta_custom(best_scaled_effect, trunc_scaled_effect, var_around_best, ntreat, alpha, beta)
  unbiased_results_scaled <- get_perfectly_unbiased_custom(best_scaled_effect, trunc_scaled_effect, var_around_best, ntreat, alpha)


  hybrid_results <- (1/sqrt(nobs(model_pl))) * hybrid_results_scaled
  unbiased_results <- (1/sqrt(nobs(model_pl))) * unbiased_results_scaled


  #-- Store WC adjusted results
  best_pol_list <- c(best_pol_list, policy_name_mapping[[pol_best_name]])
  best_pol_est_list <- c(best_pol_est_list, hybrid_results[1])
  best_pol_CI_lower <- c(best_pol_CI_lower, hybrid_results[2])
  best_pol_CI_upper <- c(best_pol_CI_upper, hybrid_results[3])



  print(pval_cutoff)
  print(policy_name_mapping[[pol_best_name]])
  print(hybrid_results)

  
}



# #########################
#-- OUTPUT GRAPHS
# #########################

setwd(path_figures)

#-- Post WC Adjustment
df_vis <- data.frame(best_policy = best_pol_list,
                     #p_cutoff = pval_cutoff_list,
                     lambda_val = lambda_grid,
                     wc_adj_estimate = best_pol_est_list,
                     L =best_pol_CI_lower,
                     U =best_pol_CI_upper)

policy_switch_count <- sum(tail(df_vis$best_policy, -1) != head(df_vis$best_policy, -1))
policy_selection_summary <- df_vis %>% 
  count(best_policy, name = "times_selected") %>% 
  mutate(share_selected = times_selected / nrow(df_vis)) %>% 
  arrange(desc(times_selected), best_policy)

policy_stability_table_dir <- paste0(path_tables, "WC_adjusted_estimates")
dir.create(policy_stability_table_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(policy_selection_summary, paste0(policy_stability_table_dir, "/policy_selection_stability_", outcome, ".csv"), row.names = FALSE)

print(policy_selection_summary)
print(paste0("Policy switches across lambda grid: ", policy_switch_count))


p_WC <- ggplot(df_vis, aes(x =  lambda_val , y = wc_adj_estimate, color= best_policy)) +
  geom_errorbar(aes(ymax = U, ymin = L),lwd = 1) + 
  geom_point(size = 3) +
  ggtitle(ifelse(outcome == "shot_Measles1", "A2: Immunizations", "A1: Immunizations/$")) + 
  xlab("Lambda") + 
  ylab("WC adjusted Estimate") + 
  scale_colour_discrete(name = "Best policy") + 
  #scale_x_continuous(breaks = lambda_grid) + 
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(face = "bold", size = 17),
        legend.position = "bottom")


ggsave(paste0("Regularization_Path/lambda_robustness_",outcome,".pdf"), plot=p_WC, width = 12, height = 8)


plot_policy_stability <- ggplot(df_vis, aes(x = lambda_val, y = factor(best_policy, levels = policy_selection_summary$best_policy), color = best_policy, group = 1)) +
  geom_step(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 3) +
  ggtitle(ifelse(outcome == "shot_Measles1", "Policy stability: Immunizations", "Policy stability: Immunizations/$")) +
  xlab("Lambda") +
  ylab("Selected best policy") +
  labs(caption = paste0("Policy switches across lambda grid: ", policy_switch_count)) +
  scale_colour_discrete(name = "Best policy") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(face = "bold", size = 17),
        legend.position = "bottom")

ggsave(paste0("Regularization_Path/lambda_policy_stability_",outcome,".pdf"), plot = plot_policy_stability, width = 12, height = 8)


#-- Pre WC Adjustment
df_best <- data.frame(Policy = best_pol_list_preWC,
                      lambda_val = lambda_grid,
                      Estimate = best_pol_est_list_preWC,
                      SE =  best_pol_se_preWC ,
                      p_val = best_pol_p_preWC,
                      Rank = "Best policy") 

df_sec_best <- data.frame(Policy = sec_best_pol_list_preWC,
                          lambda_val = lambda_grid,
                          Estimate = sec_best_pol_est_list_preWC,
                          SE = sec_best_pol_se_preWC,
                          p_val = sec_best_pol_p_preWC,
                          Rank = "Second best")


df_pre_WC <- rbind(df_best,df_sec_best)

#Measles shots
dodge_width = switch(outcome, "shot_Measles1" = 0.015, "shots_per_dollar" = 0.00005)
colvalues = switch(outcome, "shot_Measles1" = c("#B79F00","#F8766D","#00BFC4", "#00BA38", "#619CFF", "#F564E3"), "shots_per_dollar" = hue_pal()(7))
legend_rows = switch(outcome, "shot_Measles1" = c(2,3), "shots_per_dollar" = c(2,4))

plot_pre_WC <- ggplot(df_pre_WC,aes(x = lambda_val, y = Estimate, group = Estimate, color= Policy)) +
  geom_point(aes(shape = Rank), size = 3, position=position_dodge(width=dodge_width))  +
  geom_errorbar(aes(ymax =Estimate+1.96*SE , ymin =Estimate-1.96*SE),position=position_dodge(width=dodge_width),lwd = 1)  + 
  ggtitle(ifelse(outcome == "shot_Measles1", "B2: Immunizations", "B1: Immunizations/$")) + 
  xlab("Lambda") + 
  ylab("pre-WC adjusted Estimate") + 
  scale_color_manual(values = colvalues,guide = guide_legend(nrow = legend_rows[2])) + 
  scale_shape_manual(values = c("Best policy" = 17, "Second best" = 16), guide = guide_legend(nrow = legend_rows[1]))+
  #scale_x_continuous(breaks = lambda_grid) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(face = "bold", size = 17),
        legend.position = "bottom")

ggsave(paste0("Regularization_Path/lambda_robustness_preWC_",outcome,".pdf"), plot=plot_pre_WC, width = 12, height = 8)
