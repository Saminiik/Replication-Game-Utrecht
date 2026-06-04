######################################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Show that Puffered LASSO outperforms Naive LASSO on support selection (alpha space)
# AUTHOR:		Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	02/11/2021
# LAST MODIFIED: 03/01/2022 by Apara Venkat
# STATUS: 	Draft

######################################################################################################


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
library("BBSSL")
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
source('sslasso_functions.R')



#################################################################################
#
# I.GENERATE ORACLE INFORMATION FOR SIMULATIONS
#
#################################################################################
R = c(5,5,3)
M = 3

stdev <- 1
epsilon = 2.3*stdev

#Create policy vectors given R and M and omit pure control
policy_vectors  <- create_policy_vectors(R,M)[-1]
G <- create_sp_to_unique_transformation(R,M)[-1,-1]

pos_treat_effects <- linspace(1,5,M) #policy effects from 1 to 5 (in the end its going to be a multiple of standard deviation)
# pos_indices_universe <- c()
# 
# #Two possible conditional random logics
# 
# # for (i in 1:length(policy_vectors)) {
# #   candidate <- policy_vectors[[i]]
# #   if ((sum(candidate == 1) == 0) & (sum(candidate) != R*M)) { #all treatment arms are on, but don't pick the very highest intensity (intuition: to fail irrepresentability). 1 means off (lowest intensity)
# #     pos_indices_universe <- c(pos_indices_universe,i)
# #   }
# # }
# 
# 
# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

pos_indices_universe <- seq(1:length(policy_vectors))

#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################

## SETTING UP SIMULATION VARIABLES ##

nobs_list <- round(logspace(3,4,10))
#nobs = nobs_list[3]
ncomb = 5 #for the paper     # nb of support configurations i.e combinations of effective policies, each simulation having a sample size of n
nsim = 20 #for the paper

#construct the combinations out of pos_indices_universe which is the set of policies where all arms are on but excluding the one with highest intensity
pol_sets <- list()
for (i in 1:ncomb) {
  pol_sets[[i]] <- sample(pos_indices_universe,M, replace = FALSE)   
}


# Parameters for the spike and slab priors (lambda1 is defined at function call)
lambda0 <- 1e5


#Initialise support accuracy measures over n
support_acc_n_backward <- c()
support_acc_n_bbssl <- c()
support_acc_subset_n_backward <- c()
support_acc_subset_n_bbssl <- c()

#Initialize best policy inclusion measures over n
cover_some_best_n_backward <- c()
cover_some_best_n_bbssl <- c()
cover_min_n_backward <- c()
cover_min_n_bbssl <- c()

#Initialize MSE measures over n
deviations_SP_hybrid_over_n <- c()
deviations_sslasso_over_n <- c()

for (nobs in nobs_list) {
  cat("N =", nobs, "\n---\n")
  
  #Initialise accuracy measures over conf
  support_acc_n_list_backward <- c()
  support_acc_subset_n_list_backward <- c()
  support_acc_n_list_bbssl <- c()
  support_acc_subset_n_list_bbssl <- c()
  
  #Initialize best policy inclusion measures over conf
  cover_some_best_config_backward <- c()
  cover_some_best_config_bbssl <- c()
  cover_min_config_backward <- c()
  cover_min_config_bbssl <- c()
  
  #Initialize MSE measures over conf
  deviations_SP_hybrid_over_conf <- c()
  deviations_sslasso_over_conf <- c()
  
  for (comb in 1:ncomb) {
    cat("Combination:", comb, "\n")
    
    pos_indices <- pol_sets[[comb]]#just taking a random combination (comb will range from 1 to ncomb. Combinations in pol_sets were randomly computed)
    #pos_indices <- c(pos_indices, Position(function(x) identical(x, c(2,2,2)), policy_vectors))
    
    sp_treat_effects_pop <- rep.int(0, length(policy_vectors)) #vector of smart pooling effects
    sp_treat_effects_pop[pos_indices] <- pos_treat_effects * stdev #zero except at pos_indices / effects are scaled by stdev
    
    unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies
    
    support_pop <- policy_vectors[pos_indices]  #True support
    support_pop_binned <- retrieve_binned_support(support_pop, R, M)  #Organize support into profiles
    collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M) #Generate intersection of profiles and complements to get final pooled policies
    
    collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M) #compute the coefficients of pooled policies as sum of alphas
    
    #True best pooled policy
    best_pol_ind_pop <- which(collapsed_pop_trueeffects == max(collapsed_pop_trueeffects)) 
    best_pol_pop <- collapsed_pop_info[best_pol_ind_pop] #this may be more than one policy, so keep it in a list
    
    #True best minimum dosage policy
    best_pol_min <- min(do.call(paste0, as.data.frame(best_pol_pop)))
    
    #Initialize performance measures over sim
    support_acc_comb_list_backward <- c()
    support_acc_comb_subset_list_backward <- c()
    support_acc_comb_list_bbssl <- c()
    support_acc_comb_subset_list_bbssl <- c()
    
    cover_some_best_sim_backward  <-c()
    cover_min_sim_backward <-c()
    cover_some_best_sim_bbssl <- c()
    cover_min_sim_bbssl <-c()
    
    deviations_SP_hybrid_over_sims <- c()
    deviations_sslasso_over_sims <- c()
    
    for (sim in 1:nsim) {
      cat("\tSimulation:", sim, "\n")
      
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
        y <- y + sp_treat_effects_pop[i] * sp_matrix[,i] 
      }
      y  = y + rnorm(nobs,0,epsilon)
      
      
      #--------------------------#
      # I.BACKWARD ELIMINATION
      #--------------------------#
      cat("\t\tPuffer via backward elimination\n")

      #Estimated Support from polling/pruning
      support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors)

      #Measure of accuracy of the estimated support (cardinal of intersection / cardinal of union)
      support_acc_backward <- length(intersect(support_pop, support_data_backward))/length(union(support_pop, support_data_backward))

      support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)

      collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
      collapsed_data_df_backward<- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)

      #Post OLS on the unique pooled policies
      model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward)

      #Get the best policy
      best_pol_ind_data_backward <- which(model_data_pl_backward$coefficients[-1] == max(model_data_pl_backward$coefficients[-1])) #this will always be a single policy
      best_pol_data_backward <- collapsed_data_info_backward[[best_pol_ind_data_backward]]

      best_coef_backward <- max(model_data_pl_backward$coefficients[-1])
      hybrid_coef_backward <- get_Andrews_estimates_custom(model_data_pl_backward, type = "hybrid", alpha = 0.05 ,beta = 0.005)
      
      #Check if best policy is a subset , and thickness captures to what extent is its a subset
      is_subset_backward <- 0
      includes_some_best_backward <- 0 
      includes_min_best_backward  <- 0
      
      for (j in 1:length(best_pol_pop)) {
        if (any(do.call(paste0,as.data.frame(best_pol_pop[[j]])) %in% do.call(paste0, as.data.frame(best_pol_data_backward)))){  #Some best policy was selected
          includes_some_best_backward <- 1
        }
        
        if (best_pol_min %in% do.call(paste0, as.data.frame(best_pol_data_backward))){        #The minimum dosage best policy was selected
          includes_min_best_backward <- 1
        }
      }
      
      #Store Results
      support_acc_comb_list_backward <- c(support_acc_comb_list_backward, support_acc_backward)
      
      cover_some_best_sim_backward <- c(cover_some_best_sim_backward, includes_some_best_backward)
      cover_min_sim_backward <-c(cover_min_sim_backward, includes_min_best_backward)
      
      deviation_SP_hybrid <- hybrid_coef_backward[1]  - max(collapsed_pop_trueeffects)
      deviations_SP_hybrid_over_sims <- c(deviations_SP_hybrid_over_sims,deviation_SP_hybrid)
      
      #Print results to keep track
      # cat("\t\t\tSupport accuracy =  ",support_acc_backward, "\n")
      # cat("\t\t\tSome best =  ",includes_some_best_backward, "\n")
      # cat("\t\t\tMin best =  ",includes_min_best_backward, "\n")
      # print(best_pol_pop)
      # print(best_pol_data_backward)
    
      #-------------------------------------------------#
      # II. BAYESIAN BOOTSTRAP SPIKE AND SLAB LASSO
      #-------------------------------------------------#
      cat("\t\tBayesian bootstrap spike and slab lasso\n")
      
      #Run BB_SSL with Puffer pre-conditioning of the design matrix
      #Puffer (based on Rohe 2014)
      X_matrix <- sp_matrix
      X.svd <- svd(X_matrix)
      U = X.svd$u
      D = diag(X.svd$d)

      puffer_F <-  U %*% solve(D) %*% t(U)

      X_PT <- puffer_F %*% X_matrix
      y_PT <- puffer_F %*% y

      #Puffer_Lasso to feed to Puffer_SSL
      reg_lasso_PT <- rlasso(y_PT~X_PT, penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE,lambda.start = NULL))
      reg_sslasso_PT <- SSLASSO_2(X_PT,y_PT, initial.beta = reg_lasso_PT$beta, lambda0 = seq(1, lambda0, length.out = 100), lambda1 = 10^(-nobs/100), max.iter =5000)

      #plot(reg_sslasso_PT)
      selected <- rowSums(reg_sslasso_PT$select)
      support_ind <- which(selected >= 95)    #take the variables that were selected 95% of the time
      support_data_bbssl <- policy_vectors[support_ind]
      
      support_data_binned_bbssl <- retrieve_binned_support(support_data_bbssl, R, M)

      collapsed_data_info_bbssl <- create_collapsed_policies(support_data_binned_bbssl, R, M)
      collapsed_data_df_bbssl <- get_collapsed_df(collapsed_data_info_bbssl, y, beta_matrix, R, M)

      #Post LASSO on the unique pooled policies
      model_data_bbssl <- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_bbssl)
      
      #Get best policy
      best_pol_ind_data_bbssl <- which(model_data_bbssl$coefficients[-1] == max(model_data_bbssl$coefficients[-1]))
      best_pol_data_bbssl <- collapsed_data_info_bbssl[[best_pol_ind_data_bbssl]]

      best_coef_bbssl <- max(model_data_bbssl$coefficients[-1])
      hybrid_coef_bbssl<- get_Andrews_estimates_custom(model_data_bbssl, type = "hybrid", alpha = 0.05 ,beta = 0.005)
      
      #special circumstances: choose nothing in support:
      if (length(support_data_bbssl) == 0) { #degenerate case
          
          support_acc_bbssl <- 0
          
      } else { #not degenerate, something collected in support
          
          #Measure of support accuracy for SSLASSO
          support_acc_bbssl <- length(intersect(support_pop, support_data_bbssl))/length(union(support_pop, support_data_bbssl))
      }
      
    
      includes_some_best_bbssl <- 0 
      includes_min_best_bbssl  <- 0
      
      
      for (j in 1:length(best_pol_pop)) {
        if (any(do.call(paste0,as.data.frame(best_pol_pop[[j]])) %in% do.call(paste0, as.data.frame(best_pol_data_bbssl)))){  #Some best policy was selected
          includes_some_best_bbssl <- 1
        }
        
        if (best_pol_min %in% do.call(paste0, as.data.frame(best_pol_data_bbssl))){        #The minimum dosage best policy was selected
          includes_min_best_bbssl <- 1
        }
      }
      
      #Store Results
      support_acc_comb_list_bbssl <- c(support_acc_comb_list_bbssl, support_acc_bbssl)
      cover_some_best_sim_bbssl <- c(cover_some_best_sim_bbssl, includes_some_best_bbssl)
      cover_min_sim_bbssl <-c(cover_min_sim_bbssl, includes_min_best_bbssl)
      
      deviation_sslasso <- best_coef_bbssl - max(collapsed_pop_trueeffects)
      deviations_sslasso_over_sims <- c(deviations_sslasso_over_sims, deviation_sslasso)
      
      #Print results to keep track
      # cat("\t\t\tSupport contains ", sum(selected >= 95), "policies\n")
      # cat("\t\t\tSupport accuracy =  ",support_acc_bbssl, "\n")
      # cat("\t\t\tSome best =  ",includes_some_best_bbssl, "\n")
      # cat("\t\t\tMin best =  ",includes_min_best_bbssl, "\n")
      # print(best_pol_pop)
      # print(best_pol_data_bbssl)
      
    }
    
    support_acc_n_list_backward <- c(support_acc_n_list_backward, mean(support_acc_comb_list_backward))
    support_acc_n_list_bbssl <- c(support_acc_n_list_bbssl, mean(support_acc_comb_list_bbssl))
    
    cover_some_best_config_backward <- c(cover_some_best_config_backward, mean(cover_some_best_sim_backward))
    cover_some_best_config_bbssl <- c(cover_some_best_config_bbssl, mean(cover_some_best_sim_bbssl))
    cover_min_config_backward <- c(cover_min_config_backward, mean(cover_min_sim_backward ))
    cover_min_config_bbssl <- c(cover_min_config_bbssl, mean(cover_min_sim_bbssl))
    
    deviations_SP_hybrid_over_conf <- c(deviations_SP_hybrid_over_conf, mean(deviations_SP_hybrid_over_sims^2))
    deviations_sslasso_over_conf <- c(deviations_sslasso_over_conf, mean(deviations_sslasso_over_sims^2))
    
  }
  
  support_acc_n_backward <- c(support_acc_n_backward, mean(support_acc_n_list_backward))
  support_acc_n_bbssl <- c(support_acc_n_bbssl, mean(support_acc_n_list_bbssl))
  
  cover_some_best_n_backward <- c(cover_some_best_n_backward, mean(cover_some_best_config_backward))
  cover_some_best_n_bbssl <- c(cover_some_best_n_bbssl, mean(cover_some_best_config_bbssl))
  cover_min_n_backward <- c(cover_min_n_backward, mean(cover_min_config_backward))
  cover_min_n_bbssl <- c(cover_min_n_bbssl, mean(cover_min_config_bbssl))
  
  deviations_SP_hybrid_over_n <- c(deviations_SP_hybrid_over_n, mean(deviations_SP_hybrid_over_conf))
  deviations_sslasso_over_n <- c(deviations_sslasso_over_n, mean(deviations_sslasso_over_conf))
  
  cat("Support Accuracy Backward:", mean(support_acc_n_list_backward), "\n")
  cat("Support Accuracy BBSL:", mean(support_acc_n_list_bbssl), "\n")
  
  
} #that's the end of the loop across nobs, so everything here is done for each n



#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################

setwd(path_output_data)

output_support =  data.frame("support_acc_n_backward" = support_acc_n_backward,
                             "support_acc_n_bbssl" = support_acc_n_bbssl,
                             "deviations_SP_hybrid_over_n" = deviations_SP_hybrid_over_n,
                             "deviations_sslasso_over_n" = deviations_sslasso_over_n,
                             "nobs_list" = nobs_list, "ncomb" = ncomb, "nsim" = nsim)

output_inclusion = data.frame("cover_some_best_n_backward" = cover_some_best_n_backward,
                              "cover_some_best_n_bbssl" = cover_some_best_n_bbssl,
                              "cover_min_n_backward" = cover_min_n_backward,
                              "cover_min_n_bbssl" = cover_min_n_bbssl,
                              "nobs_list" = nobs_list, "ncomb" = ncomb, "nsim" = nsim)



write.csv(output_support, "Simulation_Data/5A_SSLASSO_alpha_support_any.csv")
write.csv(output_inclusion, "Simulation_Data/5A_SSLASSO_alpha_inclusion_any.csv")



#################################################################################
#
# III.END
#
#################################################################################

