#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Show that SP estimator is robust to different Regimes
# AUTHOR:		Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	Jan.10.2022
# LAST MODIFIED: 
# STATUS: 	Draft

#######################################################################################

#The goal of this file is to show that the the SP estimator is robust to different Regimes
#We perform for each regime:
#     - Support Accuracy Simulation:    SP vs LASSO on the alpha space
#     - MSE of best policy effect:      SP vs LASSO on the alpha space
#     - Best policy inclusion measures: SP vs LASSO on the alpha space


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

set.seed(25)

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


#Setting up the Regimes

#Different Regimes
regime = "R1"
c = 0.4


# pos_indices_universe <- c()
# 
# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) == 0) & (sum(candidate) != R*M)) { #all treatment arms are on, but don't pick the very highest intensity (intuition: to fail irrepresentability). 1 means off (lowest intensity)
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

pos_indices_universe = seq(1:length(policy_vectors))


#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################

## SETTING UP SIMULATION VARIABLES ##

nobs_list <- round(logspace(3,4,10)) #logarithmically spaced for computational reasons
ncomb = 5 #5 for paper     # nb of support configurations i.e combinations of effective policies, each simulation having a sample size of n
nsim = 20  #20 for paper


#construct the combinations out of pos_indices_universe which is the set of policies where all arms are on but excluding the one with highest intensity
indices_to_sample = switch(regime, "R1" = M, "R2" = M, "R3" = 2*M, "R4" = M, "R5" = M) #For R3 include the medium policies into the support
pol_sets <- list()
for (i in 1:ncomb) {
  pol_sets[[i]] <- sample(pos_indices_universe,indices_to_sample, replace = FALSE)   
}

#----#

#Initialize support accuracy measures over n
support_acc_backward_over_n <- c()
support_acc_Victor_vanilla_over_n <- c()

#Initialize subset / superset measures
lasso_superset_over_n <- c()
lasso_subset_over_n <- c()
lasso_equal_over_n <- c()
backward_superset_over_n <- c()
backward_subset_over_n <- c()
backward_equal_over_n <- c()

#Initialize MSE measures over n
deviations_SP_hybrid_over_n <- c()
deviations_Victor_vanilla_hybrid_over_n <- c()

#Initialize best policy inclusion measures over n
thickness_backward_over_n <- c()
thickness_Victor_vanilla_over_n <- c()
cover_some_best_backward_over_n <- c()
cover_some_best_vanilla_over_n <- c()
cover_min_backward_over_n <- c()
cover_min_vanilla_over_n <- c()



for (nobs in nobs_list) {
  cat("\nN = ", nobs)
  
  #Treatment effects depending on n
  pos_treat_effects_high <- switch(regime,
                                   "R1" = linspace(1,5,M), #few effective policies
                                   "R2" = linspace(1,5,M),
                                   "R3" = c(linspace(5,10,M),linspace(1,2,M)),
                                   "R4" = linspace(1,5,M)/(nobs^(0.5-c)),
                                   "R5" = linspace(1,5,M)/(nobs^(0.5-c)))
  
  pos_treat_effects_low <-  switch(regime,
                                   "R1" = linspace(1,5,length(policy_vectors)-M)/nobs, #imperfect sparsity for other effects
                                   "R2" = linspace(1,5,length(policy_vectors)-M)/sqrt(nobs),
                                   "R3" = linspace(1,5,length(policy_vectors)-2*M)/sqrt(nobs),
                                   "R4" = rep.int(0, length(policy_vectors)-M),
                                   "R5" = linspace(1,5,length(policy_vectors)-M)/nobs)
  
  
  #----#
  
  #Initialize support accuracy measures over configurations
  support_acc_Victor_vanilla_over_conf <- c()
  support_acc_backward_over_conf <- c()
  
  #Initialize subset/superset measures over configurations
  lasso_superset_over_conf <- c()
  lasso_subset_over_conf <- c()
  lasso_equal_over_conf <- c()
  backward_superset_over_conf <- c()
  backward_subset_over_conf <- c()
  backward_equal_over_conf <- c()
  
  #Initialize MSE measures over configurations
  deviations_SP_hybrid_over_conf <- c()
  deviations_Victor_vanilla_hybrid_over_conf <- c()

  #Initialize best policy inclusion measures over configurations
  thickness_backward_over_conf <- c()
  thickness_Victor_vanilla_over_conf <- c()
  cover_some_best_backward_over_conf <- c()
  cover_some_best_vanilla_over_conf <- c()
  cover_min_backward_over_conf <- c()
  cover_min_vanilla_over_conf <- c()
  
  
  for (comb in 1:ncomb) {
    cat("\n\tConfiguration = ", comb)
    
    pos_indices <- pol_sets[[comb]] #just taking a random combination (comb will range from 1 to ncomb. Combinations in pol_sets were randomly computed)
    
    
    sp_treat_effects_pop <- rep.int(0, length(policy_vectors)) #vector of smart pooling effects
    sp_treat_effects_pop[pos_indices] <- pos_treat_effects_high*stdev#zero except at pos_indices (pos_treatment effect is linspace(1,2,M)/sqrt(nobs))
    sp_treat_effects_pop[-pos_indices] <- pos_treat_effects_low*stdev
    
    unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies
    
    support_pop <- policy_vectors[pos_indices] #True support
    support_pop_binned <- retrieve_binned_support(support_pop, R, M)
    collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M)
    collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M) 
    
    #True best pooled policy
    best_pol_ind_pop <- which(collapsed_pop_trueeffects == max(collapsed_pop_trueeffects)) 
    best_pol_pop <- collapsed_pop_info[best_pol_ind_pop] #this may be more than one policy, so keep it in a list
    
    #True best minimum dosage policy
    best_pol_min <- min(do.call(paste0, as.data.frame(best_pol_pop)))
    
    cat("\n\tBest pool composed of: ",length(best_pol_pop[[1]])/M)
    #----#
    
    #Initialize support accuracy measures over simulations
    support_acc_Victor_vanilla_over_sim <- c()
    support_acc_backward_over_sim <- c()
    
    #Initialize subset/superset measures over simulations
    lasso_superset_over_sim <- c()
    lasso_subset_over_sim <- c()
    lasso_equal_over_sim <- c()
    backward_superset_over_sim <- c()
    backward_subset_over_sim <- c()
    backward_equal_over_sim <- c()
    
    #Initialize MSE measures over simulations
    deviations_SP_hybrid_over_sim <- c()
    deviations_Victor_vanilla_hybrid_over_sim <- c()
    
    #Initialize best policy inclusion measures over simulations
    thickness_backward_over_sim <- c()
    thickness_Victor_vanilla_over_sim <- c()
    cover_some_best_backward_over_sim <- c()
    cover_some_best_vanilla_over_sim <- c()
    cover_min_backward_over_sim <- c()
    cover_min_vanilla_over_sim <- c()
    
    for (sim in 1:nsim) {
      cat("\n\t\tSimulation = ", sim)
      
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
        y <- y + sp_treat_effects_pop[i] * sp_matrix[,i] #keep it relative to standard deviation - best policy is 2*stdev
      }
      y  = y + rnorm(nobs,0,epsilon)
      
      
      #--------------------------#
      # I.BACKWARD ELIMINATION
      #--------------------------#
      
      #Get alpha Support
      support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors) 
      
      #Deduce pooled support
      support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)
      collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
      collapsed_data_df_backward   <- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)
      
      #Post LASSO on the pooled support
      model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward)
      
      #Get the best policy
      best_pol_ind_data_backward <- which(model_data_pl_backward$coefficients[-1] == max(model_data_pl_backward$coefficients[-1])) #this will always be a single policy
      best_pol_data_backward <- collapsed_data_info_backward[[best_pol_ind_data_backward]]
     
      #Apply Winner's curse on the best policy
      hybrid_coef_backward <- get_Andrews_estimates_custom(model_data_pl_backward, type = "hybrid", alpha = 0.05 ,beta = 0.005)
      
      
      
      # print("True support")
      # print(support_pop)
      # print("SP support")
      # print(support_data_backward)
      # print("--------------------")
      # print("true pooled policies")
      # print(collapsed_pop_info)
      # print("SP pooled policies")
      # print(collapsed_data_info_backward)
      # print("--------------------")
      # print("True pooled effects")
      # print(collapsed_pop_trueeffects)
      # print("SP regression")
      # print(model_data_pl_backward)
      # print("--------------------")
      
      # -- Store Results -- #
      
      #a. Support accuracy (cardinal of intersection / cardinal of union)
      support_acc_backward <- length(intersect(support_pop, support_data_backward))/length(union(support_pop, support_data_backward)) 
      support_acc_backward_over_sim <- c(support_acc_backward_over_sim, support_acc_backward)
      
      #b. Is Support strict subset, equal or superset ?
       s1 <- do.call(paste0, do.call("rbind", support_pop) %>% as.data.frame())
       b1 <- do.call(paste0, do.call("rbind", support_data_backward) %>% as.data.frame())

      is_backward_subset <- 0
      is_backward_equal <- 0
      is_backward_superset <- 0
      
      if (all(s1 %in% b1) & (length(s1) == length(b1))){
        is_backward_subset <- 0
        is_backward_equal <- 1
        is_backward_superset <- 0
      } else if (all(s1 %in% b1) & (length(s1) < length(b1))){
        is_backward_subset <- 0
        is_backward_equal <- 0
        is_backward_superset <- 1
      } else if (all(b1 %in% s1) & (length(s1) > length(b1))){
        is_backward_subset <- 1
        is_backward_equal <- 0
        is_backward_superset <- 0
      }
      
        
      backward_subset_over_sim <- c(backward_subset_over_sim,is_backward_subset)
      backward_equal_over_sim <- c(backward_equal_over_sim, is_backward_equal)
      backward_superset_over_sim <- c(backward_superset_over_sim,is_backward_superset)
      
      #c. MSE
      deviation_SP_hybrid <- hybrid_coef_backward[1]  - max(collapsed_pop_trueeffects)
      deviations_SP_hybrid_over_sim <- c(deviations_SP_hybrid_over_sim,deviation_SP_hybrid)
      
      #d. Best Policy Inclusion Measures
      
        #Check if best policy is a subset , and thickness captures to what extent is its a subset
        is_subset_backward <- 0
        thickness_backward <- 0
        includes_some_best_backward <- 0 
        includes_min_best_backward  <- 0
        
        for (j in 1:length(best_pol_pop)) {
          if (any(do.call(paste0,as.data.frame(best_pol_pop[[j]])) %in% do.call(paste0, as.data.frame(best_pol_data_backward)))){  #Some best policy was selected
            includes_some_best_backward <- 1
          }
          
          if (best_pol_min %in% do.call(paste0, as.data.frame(best_pol_data_backward))){        #The minimum dosage best policy was selected
            includes_min_best_backward <- 1
          }
          if (all(do.call(paste0, as.data.frame(best_pol_data_backward)) %in% do.call(paste0,as.data.frame(best_pol_pop[[j]])))) { #Best pooled is subset of true best pooled
            is_subset_backward <- 1
            break
          }
        }
        
        if (is_subset_backward == 1) { #In that case compute support inclusion accuracy
          thickness_backward <- thickness_backward + dim(best_pol_data_backward)[1]/dim(best_pol_pop[[j]])[1]
        }
        
        thickness_backward_over_sim  <- c(thickness_backward_over_sim , thickness_backward)
        cover_some_best_backward_over_sim <- c(cover_some_best_backward_over_sim, includes_some_best_backward)
        cover_min_backward_over_sim <-c(cover_min_backward_over_sim, includes_min_best_backward)
      
      
      
      #--------------------------#
      # II.NAIVE LASSO
      #--------------------------#
      
      #Get alpha Support
      support_data_Victor_vanilla <- retrieve_vanilla_support_Victor(y, sp_matrix, policy_vectors)
      
      #Deduce pooled support
      support_data_binned_Victor_vanilla <- retrieve_binned_support(support_data_Victor_vanilla, R, M)
      collapsed_data_info_Victor_vanilla <- create_collapsed_policies(support_data_binned_Victor_vanilla, R, M)
      collapsed_data_df_Victor_vanilla <- get_collapsed_df(collapsed_data_info_Victor_vanilla, y, beta_matrix, R, M)
      
      #Post LASSO on the pooled support
      model_data_Victor_vanilla <- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_Victor_vanilla)
      
      #Get best policy
      best_pol_ind_data_Victor_vanilla <- which(model_data_Victor_vanilla$coefficients[-1] == max(model_data_Victor_vanilla$coefficients[-1])) #this will always be a single policy
      best_pol_data_Victor_vanilla <- collapsed_data_info_Victor_vanilla[[best_pol_ind_data_Victor_vanilla]]
      
      #Apply winner's curse to the best policy
      hybrid_coef_Victor_vanilla<- get_Andrews_estimates_custom(model_data_Victor_vanilla, type = "hybrid", alpha = 0.05 ,beta = 0.005)
      
      # -- Store Results -- #
      
      #a. Support Accuracy 
        #special circumstances: choose nothing in support:
        if (length(support_data_Victor_vanilla) == 0) { #degenerate case
          
          support_acc_Victor_vanilla <- 0
          
        } else { #not degenerate, something collected in support
          
          #Measure of support accuracy for naive OLS
          support_acc_Victor_vanilla <- length(intersect(support_pop, support_data_Victor_vanilla))/length(union(support_pop, support_data_Victor_vanilla))
        }
        
      support_acc_Victor_vanilla_over_sim <- c(support_acc_Victor_vanilla_over_sim, support_acc_Victor_vanilla)
      
      #b. Is Support strict subset, equal or superset ?
        l1 <- do.call(paste0, do.call("rbind", support_data_Victor_vanilla) %>% as.data.frame())
        
        is_lasso_subset <- 0
        is_lasso_equal <- 0
        is_lasso_superset <- 0
        
        if (all(s1 %in% l1) & (length(s1) == length(l1))){
          is_lasso_subset <- 0
          is_lasso_equal <- 1
          is_lasso_superset <- 0
        } else if (all(s1 %in% l1) & (length(s1) < length(l1))){
          is_lasso_subset <- 0
          is_lasso_equal <- 0
          is_lasso_superset <- 1
        } else if (all(l1 %in% s1) & (length(s1) > length(l1))){
          is_lasso_subset <- 1
          is_lasso_equal <- 0
          is_lasso_superset <- 0
        }
        
        lasso_subset_over_sim <- c(lasso_subset_over_sim,is_lasso_subset)
        lasso_equal_over_sim <- c(lasso_equal_over_sim,is_lasso_equal)
        lasso_superset_over_sim = c(lasso_superset_over_sim,is_lasso_superset)
      
      #c. MSE
      deviation_Victor_vanilla_hybrid <- hybrid_coef_Victor_vanilla[1] - max(collapsed_pop_trueeffects)
      deviations_Victor_vanilla_hybrid_over_sim <- c(deviations_Victor_vanilla_hybrid_over_sim, deviation_Victor_vanilla_hybrid)
      
      #d. Best policy inclusion measures
        #Check if best policy is a subset , and thickness captures to what extent is its a subset
        is_subset_Victor_vanilla <- 0
        thickness_Victor_vanilla <- 0
        includes_some_best_vanilla <- 0 
        includes_min_best_vanilla  <- 0
        
        for (j in 1:length(best_pol_pop)) {
          if (any(do.call(paste0,as.data.frame(best_pol_pop[[j]])) %in% do.call(paste0, as.data.frame(best_pol_data_Victor_vanilla)))){  #Some best policy was selected
            includes_some_best_vanilla <- 1
          }
          
          if (best_pol_min %in% do.call(paste0, as.data.frame(best_pol_data_Victor_vanilla))){        #The minimum dosage best policy was selected
            includes_min_best_vanilla <- 1
          }
          if (all(do.call(paste0, as.data.frame(best_pol_data_Victor_vanilla)) %in% do.call(paste0,as.data.frame(best_pol_pop[[j]])))) {
            is_subset_Victor_vanilla <- 1
            break
          }
        }
        
        if (is_subset_Victor_vanilla == 1) {
          thickness_Victor_vanilla <- thickness_Victor_vanilla + dim(best_pol_data_Victor_vanilla)[1]/dim(best_pol_pop[[j]])[1]
        }
        
        #Store Results
        thickness_Victor_vanilla_over_sim <- c(thickness_Victor_vanilla_over_sim, thickness_Victor_vanilla)
        cover_some_best_vanilla_over_sim <- c(cover_some_best_vanilla_over_sim, includes_some_best_vanilla)
        cover_min_vanilla_over_sim <-c(cover_min_vanilla_over_sim, includes_min_best_vanilla)
        
        
        
        # print(support_pop)
        # print("--------")
        # print(support_data_Victor_vanilla)
        # print("--------")
        # print(support_data_backward)
        # print("--------")
        
        # cat("\n\t\t\tlasso subset = ",is_lasso_subset)
        # cat("\n\t\t\tlasso equal = ",is_lasso_equal)
        # cat("\n\t\t\tlasso superset = ",is_lasso_superset)
        # cat("\n\t\t\tbackward subset = ",is_backward_subset)
        # cat("\n\t\t\tbackward equal = ", is_backward_equal)
        # cat("\n\t\t\tbackward superset = ",is_backward_superset)
        # cat("\n\t\t\tlasso overall measure = ", thickness_Victor_vanilla)
        # cat("\n\t\t\tbackward overall measure = ", thickness_backward)
        
        
    }
    
    #Store Support accuracies over configurations
    support_acc_backward_over_conf <- c(support_acc_backward_over_conf, mean(support_acc_backward_over_sim))
    support_acc_Victor_vanilla_over_conf <- c(support_acc_Victor_vanilla_over_conf, mean(support_acc_Victor_vanilla_over_sim))
    
    #Store MSE over configurations
    # deviations_SP_hybrid_over_conf <- c(deviations_SP_hybrid_over_conf, mean(deviations_SP_hybrid_over_sim))
    # deviations_Victor_vanilla_hybrid_over_conf <- c(deviations_Victor_vanilla_hybrid_over_conf, mean(deviations_Victor_vanilla_hybrid_over_sim))
    deviations_SP_hybrid_over_conf <- c(deviations_SP_hybrid_over_conf, mean(deviations_SP_hybrid_over_sim^2))
    deviations_Victor_vanilla_hybrid_over_conf <- c(deviations_Victor_vanilla_hybrid_over_conf, mean(deviations_Victor_vanilla_hybrid_over_sim^2))
    
    #Store best policy inclusion over configurations
    thickness_backward_over_conf <- c(thickness_backward_over_conf,mean(thickness_backward_over_sim))
    thickness_Victor_vanilla_over_conf <- c(thickness_Victor_vanilla_over_conf,mean(thickness_Victor_vanilla_over_sim))
    cover_some_best_backward_over_conf <- c(cover_some_best_backward_over_conf,mean(cover_some_best_backward_over_sim))
    cover_some_best_vanilla_over_conf <- c(cover_some_best_vanilla_over_conf,mean(cover_some_best_vanilla_over_sim))
    cover_min_backward_over_conf <- c(cover_min_backward_over_conf,mean(cover_min_backward_over_sim))
    cover_min_vanilla_over_conf <- c(cover_min_vanilla_over_conf,mean(cover_min_vanilla_over_sim))
    
    
    lasso_superset_over_conf = c(lasso_superset_over_conf,mean(lasso_superset_over_sim))
    lasso_subset_over_conf <- c(lasso_subset_over_conf,mean(lasso_subset_over_sim))
    lasso_equal_over_conf <- c(lasso_equal_over_conf,mean(lasso_equal_over_sim))
    
    backward_superset_over_conf <- c(backward_superset_over_conf,mean(backward_superset_over_sim))
    backward_subset_over_conf <- c(backward_subset_over_conf,mean(backward_subset_over_sim))
    backward_equal_over_conf <- c(backward_equal_over_conf, mean(backward_equal_over_sim))
    
    
  }
  
  #Store Support accuracies over n
  support_acc_backward_over_n <- c(support_acc_backward_over_n, mean(support_acc_backward_over_conf))
  support_acc_Victor_vanilla_over_n <- c(support_acc_Victor_vanilla_over_n, mean(support_acc_Victor_vanilla_over_conf))
  
  #Store MSE over n
  #deviations_SP_hybrid_over_n <- c(deviations_SP_hybrid_over_n, mean(deviations_SP_hybrid_over_conf^2))
  #deviations_Victor_vanilla_hybrid_over_n <- c(deviations_Victor_vanilla_hybrid_over_n, mean(deviations_Victor_vanilla_hybrid_over_conf^2))
  deviations_SP_hybrid_over_n <- c(deviations_SP_hybrid_over_n, mean(deviations_SP_hybrid_over_conf))
  deviations_Victor_vanilla_hybrid_over_n <- c(deviations_Victor_vanilla_hybrid_over_n, mean(deviations_Victor_vanilla_hybrid_over_conf))
  
  #Store best policy inclusion over n
  thickness_backward_over_n <- c(thickness_backward_over_n,mean(thickness_backward_over_conf))
  thickness_Victor_vanilla_over_n <- c(thickness_Victor_vanilla_over_n,mean(thickness_Victor_vanilla_over_conf))
  cover_some_best_backward_over_n <- c(cover_some_best_backward_over_n,mean(cover_some_best_backward_over_conf))
  cover_some_best_vanilla_over_n <- c(cover_some_best_vanilla_over_n,mean(cover_some_best_vanilla_over_conf))
  cover_min_backward_over_n <- c(cover_min_backward_over_n,mean(cover_min_backward_over_conf))
  cover_min_vanilla_over_n <- c(cover_min_vanilla_over_n,mean(cover_min_vanilla_over_conf))
  
  
  lasso_superset_over_n = c(lasso_superset_over_n,mean(lasso_superset_over_conf))
  lasso_subset_over_n <- c(lasso_subset_over_n,mean(lasso_subset_over_conf))
  lasso_equal_over_n <- c(lasso_equal_over_n,mean(lasso_equal_over_conf))
  
  backward_superset_over_n <- c(backward_superset_over_n,mean(backward_superset_over_conf))
  backward_subset_over_n <- c(backward_subset_over_n,mean(backward_subset_over_conf))
  backward_equal_over_n <- c(backward_equal_over_n, mean(backward_equal_over_conf))
  
  
} #that's the end of the loop across nobs, so everything here is done for each n



#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################

setwd(path_output_data)
output <- data.frame("support_acc_Victor_vanilla_over_n" = support_acc_Victor_vanilla_over_n,
                     "support_acc_backward_over_n" = support_acc_backward_over_n,
                     "deviations_Victor_vanilla_hybrid_over_n" = deviations_Victor_vanilla_hybrid_over_n,
                     "deviations_SP_hybrid_over_n" = deviations_SP_hybrid_over_n,
                     "thickness_Victor_vanilla_over_n" = thickness_Victor_vanilla_over_n,
                     "thickness_backward_over_n" = thickness_backward_over_n,
                     "cover_some_best_vanilla_over_n" = cover_some_best_vanilla_over_n,
                     "cover_some_best_backward_over_n" = cover_some_best_backward_over_n,
                     "cover_min_vanilla_over_n"= cover_min_vanilla_over_n,
                     "cover_min_backward_over_n" = cover_min_backward_over_n,
                     "backward_superset_over_n" = backward_superset_over_n,
                     "backward_subset_over_n" = backward_subset_over_n,
                     "backward_equal_over_n" = backward_equal_over_n,
                     "lasso_superset_over_n" = lasso_superset_over_n,
                     "lasso_subset_over_n" = lasso_subset_over_n,
                     "lasso_equal_over_n" = lasso_equal_over_n,
                     "nobs_list" = nobs_list,
                     "nsim" = rep(nsim,length(nobs_list)),
                     "ncomb" = rep(ncomb,length(nobs_list)),
                     "c" = rep(c, length(nobs_list)))




write.csv(output, paste0("Simulation_Data/3_Regimes_alpha_",regime,"_any.csv"))

#################################################################################
#
# END
#
#################################################################################
## MISCALLENIOUS
# d1 <- do.call("rbind",collapsed_pop_info) %>% as.data.frame()
# d2 <- do.call("rbind",policy_vectors) %>% as.data.frame()
# d3 <- anti_join(d2,d1) #this is taking d1 - d2

