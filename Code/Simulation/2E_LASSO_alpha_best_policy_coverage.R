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

#Generate policy vectors
policy_vectors <- create_policy_vectors(R,M)
G <- create_sp_to_unique_transformation(R,M)

#Omit the pure control
policy_vectors <- policy_vectors[-1] #omit_control_policy
G <- G[-1,-1]

pos_treat_effects <- linspace(1,5,M) #policy effects

pos_indices_universe = seq(1:length(policy_vectors))

# pos_indices_universe <- c()
# 
# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }

# for (i in 1:length(policy_vectors)) {
#   candidate <- policy_vectors[[i]]
#   if ((sum(candidate == 1) == 0) & (sum(candidate) != R*M)) { #all treatment arms are on, but don't pick the very highest intensity (intuition: to fail irrepresentability)
#     pos_indices_universe <- c(pos_indices_universe,i)
#   }
# }


#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################
#nobs_list <- round(logspace(3,4.5,12))
nobs_list <- round(logspace(3,4,10))
ncomb = 5
nsim<- 20

pol_sets <- list()
for (i in 1:ncomb) {
  pol_sets[[i]] <- sample(pos_indices_universe,M, replace = FALSE)
}


#Initialize selection accuracies over n. There are four measures computed per n: 
    # cover_best measures the share of times the method selected is an exact subset of the true support
    # thickness measures the policy inclusion performance using the measure defined in the paper
    # cover_some_best measures the share of times the method selected at least one of the true best
    # cover_min measures the share of times the method selected the minimum dosage best policy
    
cover_best_n_backward <- c()
cover_best_n_Victor_vanilla <- c()
thickness_n_backward <- c()
thickness_n_Victor_vanilla <- c()

cover_some_best_n_backward <- c()
cover_some_best_n_vanilla <- c()
cover_min_n_backward <- c()
cover_min_n_vanilla <- c()

for (nobs in nobs_list) {
  cat("\nN = ", nobs)
  
  #Initialize accuracies over configurations
  cover_best_config_backward <- c()
  cover_best_config_Victor_vanilla <-c()
  thickness_config_backward <- c()
  thickness_config_Victor_vanilla <- c()
  
  cover_some_best_config_backward <- c()
  cover_some_best_config_vanilla <- c()
  cover_min_config_backward <- c()
  cover_min_config_vanilla <- c()
  
  for (comb in 1:ncomb) {
    cat("\n\t Configuration = ", comb)
    
    pos_indices <- pol_sets[[comb]] #just taking a random combination
    #pos_indices <- c(pos_indices, Position(function(x) identical(x, c(2,2,2)), policy_vectors))
    
    #derived quantities
    sp_treat_effects_pop <- rep.int(0, length(policy_vectors)) #vector of smart pooling effects
    sp_treat_effects_pop[pos_indices] <- pos_treat_effects * stdev #zero except at pos_indices - effects scaled by stdev
    unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies
    
    support_pop <- policy_vectors[pos_indices] 
    support_pop_binned <- retrieve_binned_support(support_pop, R, M)
    collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M)
    
    collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M) 
    
    #True best pooled policy
    best_pol_ind_pop <- which(collapsed_pop_trueeffects == max(collapsed_pop_trueeffects)) 
    best_pol_pop <- collapsed_pop_info[best_pol_ind_pop] #this may be more than one policy, so keep it in a list
    
    #True best minimum dosage policy
    best_pol_min <- min(do.call(paste0, as.data.frame(best_pol_pop)))
    
    #Initialize accuracies over simulations
    cover_best_sim_backward <- c()
    cover_best_sim_Victor_vanilla <- c()
    thickness_sim_backward <- c()
    thickness_sim_Victor_vanilla <- c()
    
    cover_some_best_sim_backward  <-c()
    cover_min_sim_backward <-c()
    cover_some_best_sim_vanilla <- c()
    cover_min_sim_vanilla <-c()
    
    for (sim in 1:nsim) {
      cat("\n\t\t Simulation = ", sim)

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
        y <- y + sp_treat_effects_pop[i] * sp_matrix[,i]  #generating outcomes using treatment effects
      }
      
      y  = y + rnorm(nobs,0,epsilon) #adding noise
      
      #--------------------------#
      # I.BACKWARD ELIMINATION
      #--------------------------#
      
      support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors) #retrieve support from BE
      support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)  
      
      collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
      collapsed_data_df_backward<- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)
      
      #Post LASSO on the unique pooled policies
      model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward) 
      
      #get best policy
      best_pol_ind_data_backward <- which(model_data_pl_backward$coefficients[-1] == max(model_data_pl_backward$coefficients[-1])) #this will always be a single policy
      best_pol_data_backward <- collapsed_data_info_backward[[best_pol_ind_data_backward]]
      
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
      
      
      #Store Results
      cover_best_sim_backward <- c(cover_best_sim_backward, is_subset_backward)
      thickness_sim_backward <- c(thickness_sim_backward, thickness_backward)
      cover_some_best_sim_backward <- c(cover_some_best_sim_backward, includes_some_best_backward)
      cover_min_sim_backward <-c(cover_min_sim_backward, includes_min_best_backward)
      
      #--------------------------#
      # II.NAIVE LASSO
      #--------------------------#
      
      support_data_Victor_vanilla <- retrieve_vanilla_support_Victor(y, sp_matrix, policy_vectors)
      support_data_binned_Victor_vanilla <- retrieve_binned_support(support_data_Victor_vanilla, R, M)
      collapsed_data_info_Victor_vanilla <- create_collapsed_policies(support_data_binned_Victor_vanilla, R, M)
      collapsed_data_df_Victor_vanilla <- get_collapsed_df(collapsed_data_info_Victor_vanilla, y, beta_matrix, R, M)

      #Post LASSO on the unique pooled policies
      model_data_pl_Victor_vanilla <- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_Victor_vanilla)

      #get best policy
      best_pol_ind_data_Victor_vanilla <- which(model_data_pl_Victor_vanilla$coefficients[-1] == max(model_data_pl_Victor_vanilla$coefficients[-1])) #this will always be a single policy
      best_pol_data_Victor_vanilla <- collapsed_data_info_Victor_vanilla[[best_pol_ind_data_Victor_vanilla]]

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
      cover_best_sim_Victor_vanilla <-c(cover_best_sim_Victor_vanilla, is_subset_Victor_vanilla)
      thickness_sim_Victor_vanilla <- c(thickness_sim_Victor_vanilla, thickness_Victor_vanilla)
      cover_some_best_sim_vanilla <- c(cover_some_best_sim_vanilla, includes_some_best_vanilla)
      cover_min_sim_vanilla <-c(cover_min_sim_vanilla, includes_min_best_vanilla)
      
      # print("------------------------------------")
      # print("best_pol_pop")
      # print(as.data.frame(best_pol_pop))
      # print(paste0("best_pol_min  ",best_pol_min))
      # print("best_pol_data_backward")
      # print(as.data.frame(best_pol_data_backward))
      # print("best_pol_data_Victor_vanilla")
      # print(as.data.frame(best_pol_data_Victor_vanilla))
      # print(paste0("is_subset_backward= ",is_subset_backward))
      # print(paste0("is_subset_Victor_vanilla= ",is_subset_Victor_vanilla))
      # print(paste0("includes_some_best_backward ",includes_some_best_backward))
      # print(paste0("includes_min_best_backward= ",includes_min_best_backward))
      # print(paste0("includes_some_best_vanilla= ",includes_some_best_vanilla))
      # print(paste0("includes_min_best_vanilla= ",includes_min_best_vanilla))
      # print(paste0("overall_measure_backward= ",thickness_backward))
      # print(paste0("oveall_measure_vanilla= ",thickness_Victor_vanilla))
      
      # print("------------------------------------")
      
    }
    
    #Store Results over configurations
    thickness_config_backward <- c(thickness_config_backward, mean(thickness_sim_backward))
    thickness_config_Victor_vanilla <- c(thickness_config_Victor_vanilla, mean(thickness_sim_Victor_vanilla))
    cover_best_config_backward <- c(cover_best_config_backward, mean(cover_best_sim_backward))
    cover_best_config_Victor_vanilla <- c(cover_best_config_Victor_vanilla, mean(cover_best_sim_Victor_vanilla))
    
    cover_some_best_config_backward <- c(cover_some_best_config_backward, mean(cover_some_best_sim_backward))
    cover_some_best_config_vanilla <- c(cover_some_best_config_vanilla, mean(cover_some_best_sim_vanilla))
    cover_min_config_backward <- c(cover_min_config_backward, mean(cover_min_sim_backward ))
    cover_min_config_vanilla <- c(cover_min_config_vanilla, mean(cover_min_sim_vanilla))
    
    
    
    
    #Print results
    # print("Config Thickness Backward")
    # print(mean(thickness_sim_backward))
    # 
    # print("Config Thickness Victor vanilla")
    # print(mean(thickness_sim_Victor_vanilla))
    # 
    # print("Cover Best Backward")
    # print(mean(cover_best_n_list_backward))
    # 
    # print("Cover Best Victor")
    # print(mean(cover_best_n_list_Victor_vanilla))
    
    # print(mean(cover_some_best_sim_backward))
    # print(mean(cover_some_best_sim_vanilla))
    # print(mean(cover_min_sim_backward ))
    # print(mean(cover_min_sim_vanilla))
    
    
  }
  
  #Store Results over n
  thickness_n_backward <- c(thickness_n_backward, mean(thickness_config_backward))
  thickness_n_Victor_vanilla <- c(thickness_n_Victor_vanilla, mean(thickness_config_Victor_vanilla))
  cover_best_n_backward <- c(cover_best_n_backward, mean(cover_best_config_backward))
  cover_best_n_Victor_vanilla <- c(cover_best_n_Victor_vanilla, mean(cover_best_config_Victor_vanilla))
  
  cover_some_best_n_backward <- c(cover_some_best_n_backward, mean(cover_some_best_config_backward))
  cover_some_best_n_vanilla <- c(cover_some_best_n_vanilla, mean(cover_some_best_config_vanilla))
  cover_min_n_backward <- c(cover_min_n_backward, mean(cover_min_config_backward))
  cover_min_n_vanilla <- c(cover_min_n_vanilla, mean(cover_min_config_vanilla))

  # print("N Thickness Backward")
  # print(mean(thickness_config_backward))
  # 
  # print("N Thickness Victor vanilla")
  # print(mean(thickness_config_Victor_vanilla))
  
} #end of the big n for loop



#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################

setwd(path_output_data)

output <- data.frame("thickness_n_Victor_vanilla" = thickness_n_Victor_vanilla,
                     "thickness_n_backward" = thickness_n_backward,
                     "cover_some_best_n_backward" = cover_some_best_n_backward,
                     "cover_some_best_n_vanilla" = cover_some_best_n_vanilla,
                     "cover_min_n_vanilla" = cover_min_n_vanilla,
                     "cover_min_n_backward"= cover_min_n_backward,
                     "nobs_list" = nobs_list, "ncomb" =  ncomb, "nsim" = nsim)



write.csv(output, "Simulation_Data/2E_LASSO_alpha_best_policy_coverage.csv")

#################################################################################
#
# END
#
#################################################################################
