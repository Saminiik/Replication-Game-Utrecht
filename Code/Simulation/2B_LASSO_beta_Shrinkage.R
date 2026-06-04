#####################################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Show that LASSO on beta space (3.1) imposes stronger shrinkage when using WC adjustment
# AUTHOR:		Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	02/11/2021
# LAST MODIFIED: 
# STATUS: 	Draft

#####################################################################################################



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
library('progress')
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


#Define universe of indices of policies which can be effective
pos_indices_universe <- seq(1:length(policy_vectors))

# pos_indices_universe <- c()
# 
# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

pos_indices <- sample(pos_indices_universe,M, replace = FALSE) #just taking a random combination 
#pos_indices <- c(pos_indices, Position(function(x) identical(x, c(2,2,2)), policy_vectors))


#Treatment Effects
sp_treat_effects_pop <- rep.int(0, length(policy_vectors))  #vector of smart pooling effects
sp_treat_effects_pop[pos_indices] <- pos_treat_effects*stdev #zero except at pos_indices scaled by stdev
unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies

support_pop <- policy_vectors[pos_indices] 
support_pop_binned <- retrieve_binned_support(support_pop, R, M)
collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M)

collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M) 

#collapsed_pop_info
#collapsed_pop_trueeffects
#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################

# Setting up simulation variables
nobs_list <- round(linspace(1000,10000,10))
nsim <- 20 #20 to replicate the paper 


#List with information on coefficients and amount of shrinkage
coef_shrinkage_data_list <- list() 

ind = 1 #this is used to fill data in the coef_shrinkage list

for (nobs in nobs_list) {
  cat("\nN = ", nobs)
  
  #Store hybrid coefficients and amount of shrinkage measures
  Victor_vanilla_coef_list <- c()
  backward_coef_list <- c()
  hybrid_Victor_vanilla_coef_list <- c()
  hybrid_backward_coef_list <- c()
  Victor_vanilla_prct_diff_list <- c()
  backward_prct_diff_list <- c()
  Victor_vanilla_abs_diff_list <- c()
  backward_abs_diff_list <- c()
  
  for (sim in 1:nsim) {
    cat("\n\tSimulation = ", sim)
    
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
    
    #-------------------------------#
    # II.NAIVE LASSO (on beta matrix)
    #-------------------------------#
    
    support_data_Victor_vanilla_beta_i <- retrieve_vanilla_support_Victor(y, beta_matrix, policy_vectors, get_ind = TRUE)
    support_data_Victor_vanilla_beta <- retrieve_vanilla_support_Victor(y, beta_matrix, policy_vectors, get_ind = FALSE)
    
    beta_supported_df <- cbind(y,beta_matrix[,support_data_Victor_vanilla_beta_i])
    beta_supported_df <- as.data.frame(beta_supported_df)
    
    #Post LASSO
    model_beta <- estimatr::lm_robust(formula = as.formula("y~ ."), data = beta_supported_df)
    
    #get best policy
    best_coef_Victor_vanilla <- max(model_beta$coefficients[-1])
    hybrid_coef_Victor_vanilla_beta<- get_Andrews_estimates_custom(model_beta, type = "hybrid", alpha = 0.05 ,beta = 0.005) #adjust for Winner's Curse
    
    Victor_vanilla_coef_list <- c(Victor_vanilla_coef_list,best_coef_Victor_vanilla)
    hybrid_Victor_vanilla_coef_list <- c(hybrid_Victor_vanilla_coef_list,hybrid_coef_Victor_vanilla_beta[1])
    
  
    #Store amount of shrinkage in % and absolute terms 
    backward_prct_diff_list <- c(backward_prct_diff_list, (best_coef_backward-hybrid_coef_backward[1])/best_coef_backward) #% shrinkage 
    Victor_vanilla_prct_diff_list <- c(Victor_vanilla_prct_diff_list, (best_coef_Victor_vanilla-hybrid_coef_Victor_vanilla_beta[1])/best_coef_Victor_vanilla)
    
    backward_abs_diff_list <- c(backward_abs_diff_list, best_coef_backward-hybrid_coef_backward[1]) #shrinkage in absolute terms
    Victor_vanilla_abs_diff_list <- c(Victor_vanilla_abs_diff_list, best_coef_Victor_vanilla-hybrid_coef_Victor_vanilla_beta[1])
    
    
    #Print Some Results
    # print("post ols backward")
    # print(best_coef_backward)
    # print("hybrid backward")
    # print(hybrid_coef_backward[1])
    # print("post ols lasso")
    # print( best_coef_Victor_vanilla)
    # print("hybrid lasso ")
    # print(hybrid_coef_Victor_vanilla_beta[1])
    
  } 
  
  coef_shrinkage_data_list[[ind]] = list("Lasso_best" = Victor_vanilla_coef_list, "Lasso_hybrid" = hybrid_Victor_vanilla_coef_list,
                                         "SP_best" = backward_coef_list, "SP_hybrid" = hybrid_backward_coef_list,
                                         "Lasso_prct_diff" = Victor_vanilla_prct_diff_list, "SP_prct_diff" = backward_prct_diff_list,
                                         "Lasso_abs_diff" = Victor_vanilla_abs_diff_list , "SP_abs_diff" = backward_abs_diff_list)
  
  ind = ind+1
  
}


coef_shrinkage_df = as.data.frame(coef_shrinkage_data_list[[1]]) %>% mutate("n" = rep(nobs_list[1],nsim),
                                                                            "sim" = seq(1:nsim))

for (i in 2:length(coef_shrinkage_data_list)){
  new_df <- as.data.frame(coef_shrinkage_data_list[[i]]) %>% mutate("n" = rep(nobs_list[i],nsim),
                                                                    "sim" = seq(1:nsim))
  coef_shrinkage_df = rbind(coef_shrinkage_df, new_df)
}

#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################

setwd(path_output_data)

#write.csv(coef_shrinkage_df, "/Simulation_Data/2B_LASSO_beta_Shrinkage.csv")


write.csv(coef_shrinkage_df, "Simulation_Data/2B_LASSO_beta_Shrinkage_any.csv")




#################################################################################
#
# END
#
#################################################################################
