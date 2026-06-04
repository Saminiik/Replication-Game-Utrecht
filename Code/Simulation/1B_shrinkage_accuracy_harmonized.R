#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE:  Compare performance of the OLS estimator and the Puffer one
# AUTHOR:		Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	02/11/2021
# LAST MODIFIED: 
# STATUS: 	Draft

#Old name: 1B_shrinkage_accuracy_harmonized.R

#######################################################################################


#################################################################################
#
# 0.ENVIRONMENT SETUP
#
#################################################################################

rm(list = ls())

library('tidyverse')
library('pracma')
library('gtools')
library('numbers')
library('glmnet')
library('purrr')
library('hdm')
library('pracma')
library("scales")
library('data.table')
library('comprehenr')

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

set.seed(25)

source('crossprod_functions.R')
source('simulation_functions.R')


#################################################################################
#
# I.GENERATE ORACLE INFORMATION FOR SIMULATIONS
#
#################################################################################
R = c(5,5,3)
M = 3

stdev <- 1
epsilon = 2.3*stdev

#Generate policy vectors
policy_vectors <- create_policy_vectors(R,M)
G <- create_sp_to_unique_transformation(R,M)

#Omit the pure control
policy_vectors <- policy_vectors[-1] #omit_control_policy
G <- G[-1,-1]

pos_treat_effects <- linspace(1,5,M) #policy effects

# pos_indices_universe <- c()
# 
# # for (i in 1:length(policy_vectors)) {
# #  candidate <- policy_vectors[[i]]
# #  if ((sum(candidate == 1) == 0) & (sum(candidate) != R*M)) { #all treatment arms are on, but don't pick the very highest intensity (intuition: to fail irrepresentability)
# #    pos_indices_universe <- c(pos_indices_universe,i)
# #  }
# # }
# 
# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

pos_indices_universe <- seq(1, length(policy_vectors))

#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################
nobs_list <- round(linspace(1000,10000,10)) #nobs_list <- round(logspace(3,4,10))
ncomb = 5
nsim <- 20 #20 to replicate the paper


pol_sets <- list()
for (i in 1:ncomb) {
  pol_sets[[i]] <- sample(pos_indices_universe,M, replace = FALSE)
}


#Initialize measures over n

#MSE
mse_ols_over_n <- c()
mse_hybrid_ols_over_n <- c()
mse_puffer_over_n <- c()
mse_hybrid_puffer_over_n <- c()

#Best policy inclusion
ols_some_best_acc_over_n <- c()
ols_min_best_acc_over_n <- c()
puffer_some_best_acc_over_n <- c()
puffer_min_best_acc_over_n <- c()

#List with information on coefficients and amount of shrinkage
coef_shrinkage_data_list <- list(list()) #nested list
coef_shrinkage_input_sims <- list()     #to feed the nested list

best_policy_random_rate_over_n <- c()
for (i in 1:length(nobs_list)) {
  nobs = nobs_list[i]
  cat("\nN = ", nobs)
  
  #Initialize accuracies over configurations
  #MSE
  mse_ols_over_conf <- c()
  mse_hybrid_ols_over_conf <- c()
  mse_puffer_over_conf <- c()
  mse_hybrid_puffer_over_conf <- c()
  
  #Best policy inclusion
  ols_some_best_acc_over_conf <- c()
  ols_min_best_acc_over_conf <- c()
  puffer_some_best_acc_over_conf <- c()
  puffer_min_best_acc_over_conf <- c()
  
  best_policy_random_rate_over_conf <- c()
  for (comb in 1:ncomb) {
    cat("\n\tConfiguration = ", comb)
    
    pos_indices <- pol_sets[[comb]] #just taking a random combination
    #pos_indices <- c(pos_indices, Position(function(x) identical(x, c(1,2,2)), policy_vectors)) 
    
    
    #derived quantities
    sp_treat_effects_pop <- rep.int(0, length(policy_vectors)) #vector of smart pooling effects
    sp_treat_effects_pop[pos_indices] <- pos_treat_effects * stdev #zero except at pos_indices - effects scaled by stdev
    unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies
    
    support_pop <- policy_vectors[pos_indices]
    #true_alpha_support <- policy_vectors[pos_indices]
    support_pop_binned <- retrieve_binned_support(support_pop, R, M)
    (collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M))
    
    (collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M)) 
    
    
    ## The true best policies of the model:
    best_policy <- names(unique_treat_effects_pop[,1])[which(unique_treat_effects_pop == max(unique_treat_effects_pop))]
    efficient_best_policy <- best_policy[1]
    efficient_best_policy_coef <- unique_treat_effects_pop[which(names(unique_treat_effects_pop[,1]) == efficient_best_policy)]
    
    best_policy_random_rate_over_conf <- c(best_policy_random_rate_over_conf,1/length(best_policy)) #theoretical random selection rate
  
    # #True best pooled policy
    # best_pol_ind_pop <- which(collapsed_pop_trueeffects == max(collapsed_pop_trueeffects)) 
    # best_pol_pop <- collapsed_pop_info[best_pol_ind_pop] #this may be more than one policy, so keep it in a list
    # 
    # #True best minimum dosage policy
    # best_pol_min <- min(do.call(paste0, as.data.frame(best_pol_pop)))
    
  
    #Initialize accuracies over simulations
   
    #MSE
    mse_ols_over_sim <- c()
    mse_hybrid_ols_over_sim <- c()
    mse_puffer_over_sim <- c()
    mse_hybrid_puffer_over_sim <- c()
    
    #Best policy inclusion
    ols_some_best_acc_over_sim <- c()
    ols_min_best_acc_over_sim <- c()
    puffer_some_best_acc_over_sim <- c()
    puffer_min_best_acc_over_sim <- c()
    
    #Coefficients
    ols_coef_list <- c()
    ols_hybrid_coef_list <- c()
    puffer_coef_list <- c()
    puffer_hybrid_coef_list <- c()
    
    #Shrinkage
    ols_prct_diff_list <- c()
    ols_abs_diff_list <- c()
    puffer_prct_diff_list <- c()
    puffer_abs_diff_list <- c()
    
    
    for (sim in 1:nsim) {
      cat("\n\t\tSimulation = ", sim)
    
      #DATA _ ORIGINAL ASSIGNMENTS
      original_treatments <- create_original_treatment_assignments(R,M,nobs)
      
      #DATA - SP MATRIX
      sp_matrix <- create_sp_matrix(R,M, original_treatments, policy_vectors) #no need to omit control since we feed policy_vectors
      
      beta_matrix <- sp_matrix %*% solve(G) # Y= X alpha + e becomes Y = X G^-1 G alpha + e - THIS IS THE UNIQUE SLIM POLICIES
      
      #outcome generation
      y <- rep.int(0,nobs)
      
      for (j in 1:length(policy_vectors)) {
        y <- y + sp_treat_effects_pop[j]* sp_matrix[,j]
      }
      
      y  = y + rnorm(nobs,0,epsilon)
      
      
      #------------------------------------------------#
      # OLS ON BETA MATRIX (unique policy regression)
      #------------------------------------------------#
      
      beta_df_pop <- as.data.frame(cbind(y, beta_matrix))
      #beta_df_pop <- beta_df_pop[,colSums(beta_df_pop) > 0] # Just to check if a policy is not asigned in all the data set 
      
      #OLS model
      model_beta_ols <- estimatr::lm_robust(formula = as.formula("y~ ."), data = beta_df_pop)
      model_term <- model_beta_ols$term[-1] %>% str_remove("`") %>% str_remove("`")
      
      #Get best coefficients / best policies
      beta_coefficients <- model_beta_ols$coefficients[-1]
      best_coef_ols <- max(beta_coefficients)
      hybrid_coef_ols <- get_Andrews_estimates_custom(model_beta_ols, type = "hybrid", alpha = 0.05 ,beta = 0.005) #adjust for Winner's Curse
     
      beta_best_policies <- model_term[which(beta_coefficients == max(beta_coefficients))]
    
      # Check best policy inclusion:
      includes_some_best_ols <- 0 
      includes_min_best_ols  <- 0
    
      if(beta_best_policies %in% best_policy){
        includes_some_best_ols = 1
        }
      if(beta_best_policies == efficient_best_policy){
       includes_min_best_ols = 1
      } 
      
      
      #Store Results
      
      #Coefficients
      ols_coef_list <- c(ols_coef_list,best_coef_ols)
      ols_hybrid_coef_list <- c(ols_hybrid_coef_list,hybrid_coef_ols[1])
      
      #Best policy inclusion
      ols_some_best_acc_over_sim <- c(ols_some_best_acc_over_sim,includes_some_best_ols)
      ols_min_best_acc_over_sim <- c(ols_min_best_acc_over_sim,includes_min_best_ols)
      
      #Shrinkage
      ols_prct_diff_list <- c(ols_prct_diff_list, (best_coef_ols-hybrid_coef_ols[1])/best_coef_ols)
      ols_abs_diff_list <- c(ols_abs_diff_list, best_coef_ols-hybrid_coef_ols[1])
      
      #MSE
      deviation_ols <- best_coef_ols - max(collapsed_pop_trueeffects)
      deviation_hybrid_ols <- hybrid_coef_ols - max(collapsed_pop_trueeffects)
      mse_ols_over_sim <- c(mse_ols_over_sim,deviation_ols)
      mse_hybrid_ols_over_sim <- c(mse_hybrid_ols_over_sim,deviation_hybrid_ols)
      
      
      #--------------------------#
      # BACKWARD ELIMINATION
      #--------------------------#
     
      support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors)
      support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)
      
      collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
      collapsed_data_df_backward<- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)
      
      #Post LASSO on the unique pooled policies
      model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward)
      
      backward_coefficients <- model_data_pl_backward$coefficients[-1]
      backward_model_term <-  model_data_pl_backward$term[-1] %>% str_remove("`") %>% str_remove("`")
      
      
      best_coef_backward <- max(backward_coefficients)
      hybrid_coef_backward <- get_Andrews_estimates_custom(model_data_pl_backward, type = "hybrid", alpha = 0.05 ,beta = 0.005)
      
      backward_best_policies <- backward_model_term[which(backward_coefficients == max(backward_coefficients))]
      
      # Check best policy inclusion:
      includes_some_best_puffer <- 0 
      includes_min_best_puffer  <- 0
      
      if (any(unlist(str_split(backward_best_policies,"X")) %in% best_policy)){
        includes_some_best_puffer = 1
      }
      if (efficient_best_policy %in%  unlist(str_split(backward_best_policies,"X"))){
        includes_min_best_puffer = 1
      }
      
      
      #Store results
      
      #Coefficients
      puffer_coef_list <- c(puffer_coef_list,best_coef_backward)
      puffer_hybrid_coef_list <- c(puffer_hybrid_coef_list,hybrid_coef_backward[1])
      
      #Best policy inclusion
      puffer_some_best_acc_over_sim <- c(puffer_some_best_acc_over_sim, includes_some_best_puffer)
      puffer_min_best_acc_over_sim <-c(puffer_min_best_acc_over_sim, includes_min_best_puffer)
      
      #Shrinkage
      puffer_prct_diff_list <- c(puffer_prct_diff_list, (best_coef_backward - hybrid_coef_backward[1])/best_coef_backward) #% shrinkage 
      puffer_abs_diff_list<- c(puffer_abs_diff_list, best_coef_backward - hybrid_coef_backward[1]) #shrinkage in absolute terms
     
      #MSE
      deviation_puffer <- best_coef_backward - max(collapsed_pop_trueeffects)
      deviation_hybrid_ols <- hybrid_coef_backward - max(collapsed_pop_trueeffects)
      mse_puffer_over_sim <- c(mse_puffer_over_sim,deviation_puffer)
      mse_hybrid_puffer_over_sim <- c(mse_hybrid_puffer_over_sim, mse_puffer_over_sim)
    
      
      cat("\n\t\t\tOLS Some best = ", includes_some_best_ols)
      cat("\n\t\t\tOLS Min best = ", includes_min_best_ols)
      cat("\n\t\t\tPuffer Some best = ", includes_some_best_puffer)
      cat("\n\t\t\tPuffer Min best = ", includes_min_best_puffer)
      
    }
    
    
    coef_shrinkage_input_sims[[comb]] = list("OLS_best" = ols_coef_list, "OLS_hybrid" = ols_hybrid_coef_list,
                                           "SP_best" = puffer_coef_list, "SP_hybrid" = puffer_hybrid_coef_list,
                                           "OLS_prct_diff" = ols_prct_diff_list, "SP_prct_diff" = puffer_prct_diff_list,
                                           "OLS_abs_diff" = ols_abs_diff_list , "SP_abs_diff" = puffer_abs_diff_list,
                                           "True_best" = max(collapsed_pop_trueeffects))
    

    
    
    #Store results
    #MSE
    mse_ols_over_conf <- c(mse_ols_over_conf,mean(mse_ols_over_sim^2))
    mse_hybrid_ols_over_conf <- c(mse_hybrid_ols_over_conf,mean(mse_hybrid_ols_over_sim^2))
    mse_puffer_over_conf <- c(mse_puffer_over_conf,mean(mse_puffer_over_sim^2))
    mse_hybrid_puffer_over_conf <- c(mse_hybrid_puffer_over_conf,mean(mse_hybrid_puffer_over_sim^2))
    
    #Best policy inclusion
    ols_some_best_acc_over_conf <- c(ols_some_best_acc_over_conf,mean(ols_some_best_acc_over_sim))
    ols_min_best_acc_over_conf <- c(ols_min_best_acc_over_conf,mean(ols_min_best_acc_over_sim))
    puffer_some_best_acc_over_conf <- c(puffer_some_best_acc_over_conf,mean(puffer_some_best_acc_over_sim))
    puffer_min_best_acc_over_conf <- c(puffer_min_best_acc_over_conf,mean(puffer_min_best_acc_over_sim))
  }
  
  #Store results
  #MSE
  mse_ols_over_n <- c(mse_ols_over_n,mean(mse_ols_over_conf))
  mse_hybrid_ols_over_n <- c(mse_hybrid_ols_over_n,mean(mse_hybrid_ols_over_conf))
  mse_puffer_over_n <- c(mse_puffer_over_n,mean(mse_puffer_over_conf))
  mse_hybrid_puffer_over_n <- c(mse_hybrid_puffer_over_n,mean(mse_hybrid_puffer_over_conf))
  
  #Best policy inclusion
  ols_some_best_acc_over_n <- c(ols_some_best_acc_over_n,mean(ols_some_best_acc_over_conf))
  ols_min_best_acc_over_n <- c(ols_min_best_acc_over_n,mean(ols_min_best_acc_over_conf))
  puffer_some_best_acc_over_n <- c(puffer_some_best_acc_over_n,mean(puffer_some_best_acc_over_conf))
  puffer_min_best_acc_over_n <- c(puffer_min_best_acc_over_n,mean(puffer_min_best_acc_over_conf))
  
  coef_shrinkage_data_list[[i]] <- coef_shrinkage_input_sims
  
  best_policy_random_rate_over_n = c(best_policy_random_rate_over_n, mean(best_policy_random_rate_over_conf))
  
}








#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################


output_wide = data.frame("mse_ols_over_n" = mse_ols_over_n, "mse_hybrid_ols_over_n" = mse_hybrid_ols_over_n,
                         "mse_puffer_over_n"  = mse_puffer_over_n, "mse_hybrid_puffer_over_n" = mse_hybrid_puffer_over_n,
                         "ols_some_best_acc_over_n" = ols_some_best_acc_over_n, "ols_min_best_acc_over_n"  = ols_min_best_acc_over_n,
                         "puffer_some_best_acc_over_n" = puffer_some_best_acc_over_n,"puffer_min_best_acc_over_n" = puffer_min_best_acc_over_n,
                         "theoretical_random_rate" = best_policy_random_rate_over_n, "nobs" = nobs_list, "ncomb" = ncomb, "nsim" = nsim)


output_long = data.frame()

for (i in 1:length(coef_shrinkage_data_list)){
  for (c in 1:ncomb){
    new_df <- as.data.frame(coef_shrinkage_data_list[[i]][[c]]) %>% mutate("n" = rep(nobs_list[i],nsim),"comb" = rep(c,nsim),
                                                                           "sim" = seq(1:nsim),"true_best" = max(collapsed_pop_trueeffects))
    output_long = rbind(output_long, new_df)
  }
}


setwd(path_output_data)

write.csv(output_wide, "Simulation_Data/1B_OLS_wide_harmonized_any.csv")
write.csv(output_long, "Simulation_Data/1B_OLS_long_harmonized_any.csv")

#################################################################################
#
# END
#
#################################################################################
