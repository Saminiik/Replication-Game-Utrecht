#######################################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Compare MSE Naive LASSO on 3.1 (beta space) to SP estimator on 3.2 (alpha space)
#           conditional on selecting the correct support
# AUTHOR:	Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	02/11/2021
# LAST MODIFIED: 
# STATUS: 	Draft

#######################################################################################################



#################################################################################
#
# 0.SET UP
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
library('foreign')
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

#Experiment structure
R = c(5,5,3)
M = 3



#Given R and M create all possible policy combinations and omit pure control
policy_vectors <- create_policy_vectors(R,M)[-1]
G <- create_sp_to_unique_transformation(R,M)[-1,-1]

#Setup which policies may be effective and sample indices of effective policies
# pos_indices_universe <- c()
# 
# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

pos_indices_universe <- seq(1:length(policy_vectors))
pos_indices <- sample(pos_indices_universe,M, replace = FALSE)

#pos_indices <- c(pos_indices, Position(function(x) identical(x, c(2,2,2)), policy_vectors))#just taking a randomc combination


#Noise
stdev <- 1
epsilon <- 2.3*stdev

#Treatment Effects
pos_treat_effects <- linspace(1,5,M) #policy effects

sp_treat_effects_pop <- rep.int(0, length(policy_vectors))     #vector of smart pooling effects
sp_treat_effects_pop[pos_indices] <- pos_treat_effects*stdev   #zero except at pos_indices, scaled by standard deviation
unique_treat_effects_pop <- G %*% sp_treat_effects_pop         #note this is for the slim policies

support_pop <- policy_vectors[pos_indices] 
support_pop_binned <- retrieve_binned_support(support_pop, R, M)
collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M)

collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M) 


#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################

#Set up simulation variables
nobs_list <- round(logspace(3,4,10))
nsim<- 20

#Initialize accuracy measures over n 
deviations_SP_hybrid_over_n <- c()
deviations_beta_hybrid_over_n <- c()


for (nobs in nobs_list) {
   cat("\nN = ", nobs)
  
  #Initialize accuracy measures over simulations
  deviations_SP_hybrid_over_sims <- c()
  deviations_beta_hybrid_over_sims <- c()
  
 for (sim in 1:nsim) {
    cat("\n\tSimulation = ", sim)
   
   #DATA _ ORIGINAL ASSIGNMENTS
   original_treatments <- create_original_treatment_assignments(R,M,nobs)
   
   #DATA - SP MATRIX
   sp_matrix <- create_sp_matrix(R,M, original_treatments)
   
   # Omit THE pure control
   sp_matrix <- sp_matrix[,-1]
   beta_matrix <- sp_matrix %*% solve(G) # Y= X alpha + e becomes Y = X G^-1 G alpha + e - THIS IS THE UNIQUE SLIM POLICIES
   
   #outcome generation
   y <- rep.int(0,nobs)
   
   for (i in 1:length(policy_vectors)) {
     y <- y + sp_treat_effects_pop[i] * sp_matrix[,i] #keep it relative to standard deviation - best policy is 2*stdev
   }
   
   y  = y + rnorm(nobs,0,epsilon)
   
   
   #-------------------------------------------------------#
   # I.NAIVE LASSO (on Beta matrix) with correct support
   #-------------------------------------------------------#
   
   unique_policy_support_pop <- get_unique_policy_support_pop(collapsed_pop_info)
   unique_policy_support_pop_ind <- get_unique_policy_support_pop(collapsed_pop_info, get_ind = TRUE)
   
   beta_supported_df_pop <- cbind(y,beta_matrix[,unique_policy_support_pop_ind])
   beta_supported_df_pop <- as.data.frame(beta_supported_df_pop)
   
   model_beta_pop <- estimatr::lm_robust(formula = as.formula("y~ ."), data = beta_supported_df_pop)
   
   # Andrews estimate beta model "ideal" - Andrews estimate if we took right support of beta model
   hybrid_beta_pop <- get_Andrews_estimates_custom(model_beta_pop, type = "hybrid", alpha = 0.05 ,beta = 0.005) 
   
   best_policy_name_beta_pop <- colnames(beta_supported_df_pop)[which.max(model_beta_pop$coef)]
   best_policy_beta_pl_effect <- max(model_beta_pop$coef)
   
   
   #-----------------------------------------------#
   # I.BACKWARD ELIMINATION with perfect pooling
   #-----------------------------------------------#
   
   collapsed_pop_df <- get_collapsed_df(collapsed_pop_info, y, beta_matrix, R, M)
    
   model_pop_pl <- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_pop_df) #post lasso with perfect pooling
    
   # Andrews estimate smart pool "ideal" - Andrews estimate if we pooled perfectly
   hybrid_SP_pop <- get_Andrews_estimates_custom(model_pop_pl, type = "hybrid", alpha = 0.05 ,beta = 0.005)
   
   best_policy_name_SP_pop <- colnames(collapsed_pop_df)[which.max(model_pop_pl$coef)]
   best_policy_SP_pl_effect <- max(model_pop_pl$coef)

   #Store Results
   deviation_SP_hybrid <- hybrid_SP_pop[1] - max(collapsed_pop_trueeffects)
   deviation_beta_hybrid <- hybrid_beta_pop[1] - max(collapsed_pop_trueeffects)
   
   deviations_SP_hybrid_over_sims <- c(deviations_SP_hybrid_over_sims,deviation_SP_hybrid)
   deviations_beta_hybrid_over_sims <- c(deviations_beta_hybrid_over_sims, deviation_beta_hybrid)
   
   #Print Some results
   # print("post lasso beta")
   # print(best_policy_beta_pl_effect)
   # print("hybrid beta")
   # print(hybrid_beta_pop[1])
   # 
   # print("post lasso SP")
   # print(best_policy_SP_pl_effect)
   # print("hybrid SP")
   # print(hybrid_SP_pop[1])
   
 } 
  
  deviations_SP_hybrid_over_n <- c(deviations_SP_hybrid_over_n, mean(deviations_SP_hybrid_over_sims^2))
  deviations_beta_hybrid_over_n <- c(deviations_beta_hybrid_over_n, mean(deviations_beta_hybrid_over_sims^2))

  # print("deviation SP hybrid")
  # print(mean(deviations_SP_hybrid_over_sims^2))
  # # print("deviation SP unbiased")
  # # print(mean(deviations_SP_unbiased_over_sims))
  # print("deviation beta hybrid")
  # print(mean(deviations_beta_hybrid_over_sims^2))
  # # print("deviation beta unbiased")
  # # print(mean(deviations_beta_unbiased_over_sims))

}




#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################

setwd(path_output_data)
output = data.frame("deviations_beta_hybrid_over_n"=deviations_beta_hybrid_over_n,
                    "deviations_SP_hybrid_over_n" = deviations_SP_hybrid_over_n,
                    "nobs_list" = nobs_list, "nsim" = nsim)



#write.csv(output, "Simulation_Data/2A_LASSO_beta_MSE_conditional.csv")
write.csv(output, "Simulation_Data/2A_LASSO_beta_MSE_conditional_any.csv")



#################################################################################
#
# END
#
#################################################################################