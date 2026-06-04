#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Compare MSE Naive LASSO and SP estimator on equation 3.2 (alpha space)
# AUTHOR:		Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	02/11/2021
# LAST MODIFIED: 
# STATUS: 	Draft

#######################################################################################



#################################################################################
#
# 0.SET UP
#
#################################################################################


rm(list = ls())

library('selectiveInference')
library('tidyverse')
library('pracma')
library('gtools')
library('numbers')
library('glmnet')
library('purrr')
library('hdm')
library('pracma')
library('ggpubr')
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

source('crossprod_functions.R')
source('simulation_functions.R')

set.seed(1000)


#################################################################################
#
# I.GENERATE ORACLE INFORMATION FOR SIMULATIONS
#
#################################################################################

R = c(5,5,3)
M = 3

pos_treat_effects <- linspace(1,5,M) #policy treatment effects
stdev <- 1
epsilon <- 2.3*stdev

#Create policy vectors and omit pure control
policy_vectors <- create_policy_vectors(R,M)[-1]
G <- create_sp_to_unique_transformation(R,M)[-1,-1]


#Setup which policies may be effective and sample indices of effective policies
pos_indices_universe <- c()

# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

#pos_indices <- sample(pos_indices_universe,M, replace = FALSE)
pos_indices <- sample(1:length(policy_vectors),M, replace = FALSE) #just taking a random combination
#pos_indices <- c(pos_indices, Position(function(x) identical(x, c(2,2,2)), policy_vectors))


#Treatment Effects
sp_treat_effects_pop <- rep.int(0, length(policy_vectors))  #vector of smart pooling effects
sp_treat_effects_pop[pos_indices] <- pos_treat_effects*stdev #zero except at pos_indices scaled by stdev
unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies

support_pop <- policy_vectors[pos_indices] 
support_pop_binned <- retrieve_binned_support(support_pop, R, M)
collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M)

collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M) 

#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################

# Setting up simulation variables
#nobs_list <- round(logspace(3,4,10))
nobs_list <- linspace(1000,10000,10)
nsim <- 20 #for paper

#Initialize accuracy measures over n
deviations_SP_hybrid_over_n <- c()
deviations_ols_hybrid_over_n <- c()
deviations_debiased_lasso_over_n <- c()


coef_shrinkage_data_list <- list() #List with information on coefficients and amount of shrinkage

ind = 1 #this is used to fill data in the coef_shrinkage list

for (nobs in nobs_list) {
  cat("\n N = ", nobs)
  
  #Initialize accuracy measures over sim
  deviations_SP_hybrid_over_sims <- c()
  deviations_ols_hybrid_over_sims <- c()
  deviations_debiased_lasso_over_sims <- c()

  
  #Store hybrid coefficients and amount of shrinkage measures
  ols_coef_list <- c()
  backward_coef_list <- c()
  deblasso_coef_list <- c()
  
  hybrid_ols_coef_list <- c()
  hybrid_backward_coef_list <- c()
  hybrid_deblasso_coef_list <- c()
  
  ols_prct_diff_list <- c()
  backward_prct_diff_list <- c()
  debiased_lasso_prct_diff_list <- c()
  
  ols_abs_diff_list <- c()
  backward_abs_diff_list <- c()
  debiased_lasso_abs_diff_list <- c()
  
  #j = 1
  for (sim in 1:nsim) {
    cat("\n\t Simulation = ", sim)
    
    #DATA _ ORIGINAL ASSIGNMENTS
    original_treatments <- create_original_treatment_assignments(R,M,nobs)
    
    #DATA - SP MATRIX
    sp_matrix <- create_sp_matrix(R,M, original_treatments)
    
    # OMIT THE pure control
    sp_matrix <- sp_matrix[,-1]
    beta_matrix <- sp_matrix %*% solve(G) # Y= X alpha + e becomes Y = X G^-1 G alpha + e - THIS IS THE UNIQUE SLIM POLICIES
    
    #outcome generation
    y <- rep.int(0,nobs)
    
    for (i in 1:length(policy_vectors)) {
      y <- y + sp_treat_effects_pop[i] * sp_matrix[,i] #keep it releative to standard deviation - best policy is 2*stdev
    }
    
    y  = y + rnorm(nobs,0,epsilon)
    
    #--------------------------#
    # I.BACKWARD ELIMINATION
    #--------------------------#

    support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors)
    support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)
   
    collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
    collapsed_data_df_backward<- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)

    #Post LASSO on the unique pooled policies
    model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward)

    #Get the best policy
    best_pol_ind_data_backward <- which(model_data_pl_backward$coefficients[-1] == max(model_data_pl_backward$coefficients[-1])) #this will always be a single policy
    best_pol_data_backward <- collapsed_data_info_backward[[best_pol_ind_data_backward]]


    best_coef_backward <- max(model_data_pl_backward$coefficients[-1])
    hybrid_coef_backward <- get_Andrews_estimates_custom(model_data_pl_backward, type = "hybrid", alpha = 0.05 ,beta = 0.005)

    backward_coef_list <- c(backward_coef_list,best_coef_backward)
    hybrid_backward_coef_list <- c(hybrid_backward_coef_list,hybrid_coef_backward[1])
    
    
    #------------------------------------------------#
    # II. OLS ON BETA MATRIX (unique policy regression)
    #------------------------------------------------#
    
    beta_df_pop <- as.data.frame(cbind(y, beta_matrix))
    
    #OLS model
    model_beta_ols <- estimatr::lm_robust(formula = as.formula("y~ ."), data = beta_df_pop)
    
    #get best policy
    best_coef_ols <- max(model_beta_ols$coefficients[-1])
    hybrid_coef_ols <- get_Andrews_estimates_custom(model_beta_ols, type = "hybrid", alpha = 0.05 ,beta = 0.005) #adjust for Winner's Curse
    
    ols_coef_list <- c(ols_coef_list, best_coef_ols)
    hybrid_ols_coef_list <- c(hybrid_ols_coef_list, hybrid_coef_ols[1])
    
    #-------------------------------------------------------------#
    # III. DEBIASED LASSO (fully desparsified - on beta space)
    #-------------------------------------------------------------#
    
    #get coefficients:
    debiased_lasso <- get_debiased_lasso(y, beta_matrix)
    coefficients_debiased_lasso <- debiased_lasso[[1]] # Coefficient of the debiased support
    model_lasso <- debiased_lasso[[2]] # the naive model lasso
    M_JM <- debiased_lasso[[3]] # The matrix M based on Javanmard and Montanari

    # Transformation of the new coefficients :
    coefficients_debiased_lasso <- as.vector(coefficients_debiased_lasso)
    names(coefficients_debiased_lasso) = names(model_lasso$coefficients[-1]) %>% str_remove("`") %>% str_remove("`")

    # Andrews estimate (in the debiased case):
    hybrid_results <- get_Andrews_estimates_custom(coefficients_debiased_lasso, type = "hybrid", alpha = 0.05 ,beta = 0.005, debiased = TRUE, x = beta_matrix, y = y, M_JM = M_JM, nobs_andrews = nobs)

    pl_effects <- coefficients_debiased_lasso
    pol_best_name <- names(which(pl_effects == max(pl_effects)))
    
    

    best_debiased_lasso = pl_effects[pol_best_name]

    deblasso_coef_list <- c(deblasso_coef_list,best_debiased_lasso)
    hybrid_deblasso_coef_list <- c(hybrid_deblasso_coef_list,hybrid_results[1])



    #Store Results
    deviation_SP_hybrid <- hybrid_coef_backward[1]  - max(collapsed_pop_trueeffects)
    deviation_ols_hybrid <- hybrid_coef_ols[1] - max(collapsed_pop_trueeffects)
    deviation_debiased_lasso <- hybrid_results[1] - max(collapsed_pop_trueeffects)

    deviations_SP_hybrid_over_sims <- c(deviations_SP_hybrid_over_sims,deviation_SP_hybrid )
    deviations_ols_hybrid_over_sims <- c(deviations_ols_hybrid_over_sims, deviation_ols_hybrid)
    deviations_debiased_lasso_over_sims <- c(deviations_debiased_lasso_over_sims, deviation_debiased_lasso)

    #Amount of shrinkage in % and absolute terms
    backward_prct_diff_list <- c(backward_prct_diff_list, (best_coef_backward-hybrid_coef_backward[1])/best_coef_backward) #% shrinkage
    ols_prct_diff_list <- c(ols_prct_diff_list, (best_coef_ols - hybrid_coef_ols[1])/best_coef_ols)
    debiased_lasso_prct_diff_list <- c(debiased_lasso_prct_diff_list, (best_debiased_lasso - hybrid_results[1])/best_debiased_lasso)

    backward_abs_diff_list <- c(backward_abs_diff_list, best_coef_backward-hybrid_coef_backward[1]) #shrinkage in absolute terms
    ols_abs_diff_list <- c(ols_abs_diff_list, best_coef_ols - hybrid_coef_ols[1])
    debiased_lasso_abs_diff_list <- c(debiased_lasso_abs_diff_list, best_debiased_lasso - hybrid_results[1])
    
    cat("\n\t\t TVA % diff = ", (best_coef_backward-hybrid_coef_backward[1])/best_coef_backward)
    cat("\n\t\t Debiased % diff = ", (best_debiased_lasso - hybrid_results[1])/best_debiased_lasso)
    
  } 
  
  deviations_SP_hybrid_over_n <- c(deviations_SP_hybrid_over_n, mean(deviations_SP_hybrid_over_sims^2))
  deviations_ols_hybrid_over_n <- c(deviations_ols_hybrid_over_n, mean(deviations_ols_hybrid_over_sims ^2))
  deviations_debiased_lasso_over_n <- c(deviations_debiased_lasso_over_n, mean(deviations_debiased_lasso_over_sims^2))
  
  coef_shrinkage_data_list[[ind]] = list("OLS_best" = ols_coef_list, "OLS_hybrid" = hybrid_ols_coef_list, 
                                         "SP_best" = backward_coef_list, "SP_hybrid" = hybrid_backward_coef_list, 
                                         "DLasso_best" = deblasso_coef_list, "DLasso_hybrid" = hybrid_deblasso_coef_list,
                                         "OLS_prct_diff" = ols_prct_diff_list, "SP_prct_diff" = backward_prct_diff_list, "DL_prct_diff" = debiased_lasso_prct_diff_list,
                                         "OLS_abs_diff" = ols_abs_diff_list , "SP_abs_diff" = backward_abs_diff_list, "DL_abs_diff" = debiased_lasso_abs_diff_list)
  
                                         
                                         
  ind = ind+1
  
}


coef_shrinkage_df = as.data.frame(coef_shrinkage_data_list[[1]]) %>% mutate("n" = rep(nobs_list[1],nsim))

for (i in 2:length(coef_shrinkage_data_list)){
  new_df <- as.data.frame(coef_shrinkage_data_list[[i]]) %>% mutate("n" = rep(nobs_list[i],nsim))
  coef_shrinkage_df = rbind(coef_shrinkage_df, new_df)
}

#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################

setwd(path_output_data)
output2 = data.frame("deviations_SP_hybrid_over_n" = deviations_SP_hybrid_over_n,
                    "deviations_ols_hybrid_over_n" = deviations_ols_hybrid_over_n,
                    "deviations_debiased_lasso_over_n" = deviations_debiased_lasso_over_n,
                    "nobs_list" = nobs_list,
                    "nsim" = nsim)


write.csv(coef_shrinkage_df, "Simulation_Data/4A_Debiased_LASSO_any.csv")
write.csv(output2, "Simulation_Data/4A_Debiased_LASSO_mse_any.csv")



#################################################################################
#
# END
#
#################################################################################
