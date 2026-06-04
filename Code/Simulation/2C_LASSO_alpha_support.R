######################################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Show that Puffered LASSO outperforms Naive LASSO on support selection (alpha space)
# AUTHOR:		Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	02/11/2021
# LAST MODIFIED: 
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

pos_treat_effects <- linspace(1,5,M) #policy effects from 1 to 5 (in the end its going to be a multiple of standard deviation)

pos_indices_universe <- seq(1:length(policy_vectors))

# pos_indices_universe <- c()
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


#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################

## SETTING UP SIMULATION VARIABLES ##


nobs_list <- round(logspace(3,4,10))
nobs = nobs_list[3]
ncomb = 5;     # nb of support configurations i.e combinations of effective policies, each simulation having a sample size of n
nsim = 20


#construct the combinations out of pos_indices_universe which is the set of policies where all arms are on but excluding the one with highest intensity
pol_sets <- list()
for (i in 1:ncomb) {
  pol_sets[[i]] <- sample(pos_indices_universe,M, replace = FALSE)
}


#Initialise support accuracy measures over n
support_acc_n_backward <- c()
support_acc_n_Victor_vanilla <- c()
support_acc_subset_n_backward <- c()
support_acc_subset_n_Victor_vanilla <- c()

for (nobs in nobs_list) {
  cat("\nN = ", nobs)

  #Initialise accuracy measures over sims
  support_acc_n_list_Victor_vanilla <- c()
  support_acc_subset_n_list_Victor_vanilla <- c()
  support_acc_n_list_backward <- c()
  support_acc_subset_n_list_backward <- c()

  for (comb in 1:ncomb) {
    cat("\n\n\tConfiguration =  ", comb)

    pos_indices <- pol_sets[[comb]]#just taking a random combination (comb will range from 1 to ncomb. Combinations in pol_sets were randomly computed)
    #pos_indices <- c(pos_indices, Position(function(x) identical(x, c(2,2,2)), policy_vectors))#just taking a randomc combination


    sp_treat_effects_pop <- rep.int(0, length(policy_vectors)) #vector of smart pooling effects
    sp_treat_effects_pop[pos_indices] <- pos_treat_effects * stdev #zero except at pos_indices / effects are scaled by stdev

    support_pop <- policy_vectors[pos_indices]  #True support

    unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies

    support_pop_binned <- retrieve_binned_support(support_pop, R, M)
    collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M)

    collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M)

    support_acc_comb_list_backward <- c()
    support_acc_comb_subset_list_backward <- c()
    support_acc_comb_list_Victor_vanilla <- c()
    support_acc_comb_subset_list_Victor_vanilla <- c()


    for (sim in 1:nsim) {
      cat("\n\t\tSimulation =  ", sim)

      #DATA _ ORIGINAL ASSIGNMENTS
      original_treatments <- create_original_treatment_assignments(R,M,nobs)

      #DATA - SP MATRIX
      sp_matrix <- create_sp_matrix(R,M, original_treatments)

      # OMIT THE pure control
      sp_matrix <- sp_matrix[,-1]

      #outcome generation
      y <- rep.int(0,nobs)
      for (i in 1:length(policy_vectors)) {
        y <- y + sp_treat_effects_pop[i] * sp_matrix[,i] #keep it relative to standard deviation - best policy is 2*stdev
      }
      y  = y + rnorm(nobs,0,epsilon)


      #--------------------------#
      # I.BACKWARD ELIMINATION
      #--------------------------#

      #Estimated Support from polling/pruning
      support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors)

      #Measure of accuracy of the estimated support (cardinal of intersection / cardinal of union)
      support_acc_backward <- length(intersect(support_pop, support_data_backward))/length(union(support_pop, support_data_backward))

      #--------------------------#
      # II.NAIVE LASSO
      #--------------------------#

      #Estimated support from naive OLS
      support_data_Victor_vanilla <- retrieve_vanilla_support_Victor(y, sp_matrix, policy_vectors)

      #special circumstances: choose nothing in support:
      if (length(support_data_Victor_vanilla) == 0) { #degenerate case

        support_acc_Victor_vanilla <- 0

      } else { #not degenerate, something collected in support

        #Measure of support accuracy for naive OLS
        support_acc_Victor_vanilla <- length(intersect(support_pop, support_data_Victor_vanilla))/length(union(support_pop, support_data_Victor_vanilla))
      }

      support_acc_comb_list_Victor_vanilla <- c(support_acc_comb_list_Victor_vanilla, support_acc_Victor_vanilla)
      support_acc_comb_list_backward <- c(support_acc_comb_list_backward, support_acc_backward)

      cat("\n\t\t\tAcc Backward = ",support_acc_backward)
      cat("\n\t\t\tAcc Victor Vanilla = ", support_acc_Victor_vanilla)
    }

    support_acc_n_list_backward <- c(support_acc_n_list_backward, mean(support_acc_comb_list_backward))
    support_acc_n_list_Victor_vanilla <- c(support_acc_n_list_Victor_vanilla, mean(support_acc_comb_list_Victor_vanilla))
  }

  support_acc_n_backward <- c(support_acc_n_backward, mean(support_acc_n_list_backward))
  support_acc_n_Victor_vanilla <- c(support_acc_n_Victor_vanilla, mean(support_acc_n_list_Victor_vanilla))

  print("Support Accuracy Backward")
  print(mean(support_acc_n_list_backward))

  print("Support Accuracy Victor vanilla")
  print(mean(support_acc_n_list_Victor_vanilla))

} #that's the end of the loop across nobs, so everything here is done for each n



#################################################################################
#
# III.EXPORT RESULTS
#
#################################################################################

setwd(path_output_data)

output <- data.frame("support_acc_n_Victor_vanilla"= support_acc_n_Victor_vanilla,
                     "support_acc_n_backward" = support_acc_n_backward,
                     "nobs_list" = nobs_list, "ncomb" = ncomb, "nsim" = nsim)



write.csv(output, "Simulation_Data/2C_LASSO_alpha_support_any.csv")


#################################################################################
#
# III.END
#
#################################################################################

