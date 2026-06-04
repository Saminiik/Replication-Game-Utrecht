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

library('tidyverse')
library('pracma')
library('gtools')
library('numbers')
library('glmnet')
library('purrr')
library('hdm')
library('pracma')

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

set.seed(25)


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
# pos_indices_universe <- c()
# 
# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

pos_indices_universe <- seq(1:length(policy_vectors))

pos_indices <- sample(pos_indices_universe,M, replace = FALSE) #just taking a random combination 
#pos_indices <- c(pos_indices, Position(function(x) identical(x, c(4,2,2)), policy_vectors))


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
nobs_list <- round(logspace(3,4,10))
nsim <- 20 

#Initialize accuracy measures over n
deviations_SP_hybrid_over_n <- c()
deviations_Victor_vanilla_hybrid_over_n <- c()

support_acc_backward_over_n <- c()
support_acc_Victor_vanilla_over_n <- c()

coef_shrinkage_data_list <- list() #List with information on coefficients and amount of shrinkage

ind = 1 #this is used to fill data in the coef_shrinkage list

for (nobs in nobs_list) {
  cat("\nN = ", nobs)
  
  #Initialize accuracy measures over sim
  deviations_SP_hybrid_over_sims <- c()
  deviations_Victor_vanilla_hybrid_over_sims <- c()
  support_acc_Victor_vanilla_over_sims <- c()
  support_acc_backward_over_sims <- c()
  
  #Store hybrid coefficients and amount of shrinkage measures
  hybrid_Victor_vanilla_coef_list <- c()
  hybrid_backward_coef_list <- c()
  Victor_vanilla_prct_diff_list <- c()
  backward_prct_diff_list <- c()
  Victor_vanilla_abs_diff_list <- c()
  backward_abs_diff_list <- c()
  
  #j = 1
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

    hybrid_backward_coef_list <- c(hybrid_backward_coef_list,hybrid_coef_backward[1])
    
    #-------------------------------#
    # II.NAIVE LASSO (on SP matrix)
    #-------------------------------#

    support_data_Victor_vanilla <- retrieve_vanilla_support_Victor(y, sp_matrix, policy_vectors)
    support_data_binned_Victor_vanilla <- retrieve_binned_support(support_data_Victor_vanilla, R, M)

    collapsed_data_info_Victor_vanilla <- create_collapsed_policies(support_data_binned_Victor_vanilla, R, M)
    collapsed_data_df_Victor_vanilla <- get_collapsed_df(collapsed_data_info_Victor_vanilla, y, beta_matrix, R, M)

    #Post LASSO on the unique pooled policies
    model_data_Victor_vanilla <- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_Victor_vanilla)

    #get best policy
    best_pol_ind_data_Victor_vanilla <- which(model_data_Victor_vanilla$coefficients[-1] == max(model_data_Victor_vanilla$coefficients[-1])) #this will always be a single policy
    best_pol_data_Victor_vanilla <- collapsed_data_info_Victor_vanilla[[best_pol_ind_data_Victor_vanilla]]

    best_coef_Victor_vanilla <- max(model_data_Victor_vanilla$coefficients[-1])
    hybrid_coef_Victor_vanilla<- get_Andrews_estimates_custom(model_data_Victor_vanilla, type = "hybrid", alpha = 0.05 ,beta = 0.005)

    hybrid_Victor_vanilla_coef_list <- c(hybrid_Victor_vanilla_coef_list,hybrid_coef_Victor_vanilla[1])
    
    #-------------------------------#
    # ADD SUPPORT ACCURACY MEASURES
    #-------------------------------#
    
    support_acc_backward <- length(intersect(support_pop, support_data_backward))/length(union(support_pop, support_data_backward))
    
    if (length(support_data_Victor_vanilla) == 0) { #degenerate case
      
      support_acc_Victor_vanilla <- 0
      
    } else { #not degenerate, something collected in support
      
      #Measure of support accuracy for naive OLS
      support_acc_Victor_vanilla <- length(intersect(support_pop, support_data_Victor_vanilla))/length(union(support_pop, support_data_Victor_vanilla))
    }
    
    
    #Store Results
    deviation_SP_hybrid <- hybrid_coef_backward[1]  - max(collapsed_pop_trueeffects)
    deviation_Victor_vanilla_hybrid <- hybrid_coef_Victor_vanilla[1] - max(collapsed_pop_trueeffects)
    
    deviations_SP_hybrid_over_sims <- c(deviations_SP_hybrid_over_sims,deviation_SP_hybrid )
    deviations_Victor_vanilla_hybrid_over_sims <- c(deviations_Victor_vanilla_hybrid_over_sims, deviation_Victor_vanilla_hybrid)
    
    support_acc_backward_over_sims <- c(support_acc_backward_over_sims,  support_acc_backward)
    support_acc_Victor_vanilla_over_sims <- c(support_acc_Victor_vanilla_over_sims, support_acc_Victor_vanilla)
    
    #Amount of shrinkage in % and absolute terms 
    backward_prct_diff_list <- c(backward_prct_diff_list, deviation_SP_hybrid/max(collapsed_pop_trueeffects))
    Victor_vanilla_prct_diff_list <- c(Victor_vanilla_prct_diff_list, deviation_Victor_vanilla_hybrid/max(collapsed_pop_trueeffects))
  
    backward_abs_diff_list <- c(backward_abs_diff_list,abs(deviation_SP_hybrid))
    Victor_vanilla_abs_diff_list <- c(Victor_vanilla_abs_diff_list, abs(deviation_Victor_vanilla_hybrid))
    
    
    # #Print Some Results
    # print("True Support")
    # print(support_pop)
    # print("Backward Support")
    # print(support_data_backward)
    # print("LASSO support")
    # print(support_data_Victor_vanilla)

    cat("\n\t\tSupport acc backward = ",support_acc_backward)
    cat("\n\t\tSupport accuracy LASSO = ",support_acc_Victor_vanilla)

    cat("\n\t\tpost ols backward = ",best_coef_backward)
    cat("\n\t\thybrid backward = ",hybrid_coef_backward[1])
    cat("\n\t\tpost ols lasso = ",best_coef_Victor_vanilla)
    cat("\n\t\thybrid lasso = ",hybrid_coef_Victor_vanilla[1])
    
  } 
  
  deviations_SP_hybrid_over_n <- c(deviations_SP_hybrid_over_n, mean(deviations_SP_hybrid_over_sims^2))
  deviations_Victor_vanilla_hybrid_over_n <- c(deviations_Victor_vanilla_hybrid_over_n, mean(deviations_Victor_vanilla_hybrid_over_sims ^2))
  
  support_acc_backward_over_n <- c( support_acc_backward_over_n, mean(support_acc_backward_over_sims))
  support_acc_Victor_vanilla_over_n <- c(support_acc_Victor_vanilla_over_n, mean(support_acc_Victor_vanilla_over_sims))
  
  coef_shrinkage_data_list[[ind]] = list("Lasso_hybrid" = hybrid_Victor_vanilla_coef_list, "SP_hybrid" = hybrid_backward_coef_list,
                                       "Lasso_prct_diff" = Victor_vanilla_prct_diff_list, "SP_prct_diff" = backward_prct_diff_list,
                                       "Lasso_abs_diff" = Victor_vanilla_abs_diff_list , "SP_abs_diff" = backward_abs_diff_list)
  
  ind = ind+1
  
  #print(deviations_SP_hybrid_over_n)
  #print("deviation SP hybrid")
  #print(mean(deviations_SP_hybrid_over_sims^2))
  # print("deviation SP unbiased")
  # print(mean(deviations_SP_unbiased_over_sims))
  #print("deviation beta hybrid")
  #print(mean(deviations_beta_hybrid_over_sims^2))
  # print("deviation beta unbiased")
  # print(mean(deviations_beta_unbiased_over_sims))
  
}


coef_shrinkage_df = as.data.frame(coef_shrinkage_data_list[[1]]) %>% mutate("n" = rep(nobs_list[1],nsim))

for (i in 2:length(coef_shrinkage_data_list)){
  new_df <- as.data.frame(coef_shrinkage_data_list[[i]]) %>% mutate("n" = rep(nobs_list[i],nsim))
  coef_shrinkage_df = rbind(coef_shrinkage_df, new_df)
}

#################################################################################
#
# III.PLOT RESULTS
#
#################################################################################

setwd(path_output_data)
output <- data.frame("deviations_Victor_vanilla_hybrid_over_n" = deviations_Victor_vanilla_hybrid_over_n,
                     "deviations_SP_hybrid_over_n" = deviations_SP_hybrid_over_n,
                     "nobs_list" = nobs_list, "nsim" = nsim)


write.csv(output, "Simulation_Data/2D_LASSO_alpha_MSE_any.csv")


#################################################################################
#
# END
#
#################################################################################


