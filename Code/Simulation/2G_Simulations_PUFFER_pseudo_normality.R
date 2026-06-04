#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Check Pseudo-true asymptotic normality
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

rm(list = ls()) #clean environment

library('tidyverse')
library('pracma')
library('gtools')
library('numbers')
library('glmnet')
library('purrr')
library('hdm')
library('broom')
library('pracma')
library("stringr")
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
R = c(5,5,3)      #intensity levels
M = 3             #number of treatment arms

nobs_list = c(1000,1500,3000,5000,15000) #in the paper
#nobs_list = c(10000)

#given M treatment arms and R intensities create all possible policy combinations and omit pure control
policy_vectors <- create_policy_vectors(R,M)[-1]
G <- create_sp_to_unique_transformation(R,M)[-1,-1]


#set policy effects from 1 to 2 (in the end its going to be a multiple of standard deviation)
pos_treat_effects <- linspace(1,5,M) 
stdev <- 1
epsilon = 2.3*stdev

pos_indices_universe <- seq(1:length(policy_vectors))
pos_indices <- sample(pos_indices_universe,M, replace = FALSE) #randomly choose M such indices hence cardinality of the support is M


#oracle information
sp_treat_effects_pop <- rep.int(0, length(policy_vectors))    #vector of smart pooling effects
sp_treat_effects_pop[pos_indices] <- pos_treat_effects*stdev  #zero except at pos_indices
unique_treat_effects_pop <- G %*% sp_treat_effects_pop        #note this is for the slim policies

support_pop <- policy_vectors[pos_indices]                                #true support: list of effective policies  
support_pop_binned <- retrieve_binned_support(support_pop, R, M)          #this function is organizing the policies into profiles
collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M) #this function is creating a disjoint union of profiles (like the intersections on the Hasse diagram)
#all policies within the region have the same effect

initial_collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M)  #this function gives the effect of every pooled policy

unique_policy_support_pop <- get_unique_policy_support_pop(collapsed_pop_info)



#################################################################################
#
# II. RUN SIMULATIONS
#
#################################################################################

## SETTING UP SIMULATION VARIABLES ##
nsim = 200 #200 in paper

normalized_discrepancy_list <- c()

for (nobs in nobs_list){
  cat("\nN = ", nobs)
  normalized_discrepancy = list() #to store normalized discrepancies for each simulation
  true_support_detected = 0 #count number of times the true support is detected
  
  for (sim in 1:nsim) {
    cat("\n\tSimulation: ",sim)
    
    collapsed_pop_trueeffects = initial_collapsed_pop_trueeffects   #this is to avoid recomputing them every time
    
    #DATA - original assignments
    original_treatments <- create_original_treatment_assignments(R,M,nobs) #construct treatment assignments given R,M and n
    
    #DATA - sp matrix
    sp_matrix <- create_sp_matrix(R,M, original_treatments) #construct the sp_matrix (X) corresponding to these treatment assignments
    
    #Omit the pure control
    sp_matrix <- sp_matrix[,-1]
    beta_matrix <- sp_matrix %*% solve(G) # Y= X alpha + e becomes Y = X G^-1 G alpha + e - THIS IS THE UNIQUE SLIM POLICIES
    
    #Generate Outcome
    y <- rep.int(0,nobs)
    
    for (i in 1:length(policy_vectors)) {
      y <- y + sp_treat_effects_pop[i] * sp_matrix[,i] 
    }
    y  = y + rnorm(nobs,0,epsilon) #add noise to the outcome

    #Backward Elimination Procedure
    support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors) 
    support_acc_backward <- length(intersect(support_pop, support_data_backward))/length(union(support_pop, support_data_backward)) 
    cat("\n\t\tSupp acc ",support_acc_backward)
    
    support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)  #get the pooled policies
    collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
    
    collapsed_data_df_backward<- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)
    
    model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward)
    
    
    #Store Results
    if (support_acc_backward == 1){
      true_support_detected <- true_support_detected + 1
    } else{
      #re-define pseudo true treatment effects -> TE if the pooling was the one selected smart pooling and pruning
      collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_data_info_backward, unique_treat_effects_pop, R, M)
    }
    
    normalized_discrepancy[[sim]] = (model_data_pl_backward$coefficients[-1] - collapsed_pop_trueeffects)/model_data_pl_backward$std.error[-1]
  
  }
  
  normalized_discrepancy_df <- do.call("bind_rows", normalized_discrepancy) %>% as.data.frame()
  
  cat("\n\t\tShare of true supports detected: ",true_support_detected/nsim)
  
  
  ## EXPORT RESULTS ##
  setwd(path_output_data)
  
  names(normalized_discrepancy_df) = names(normalized_discrepancy_df) %>% str_remove("`") %>% str_remove("`")
  output = normalized_discrepancy_df %>% mutate("support_detection_rate" = true_support_detected/nsim)
  write.csv(output, paste0("Simulation_Data/2G_Puffer_normality_n=",nobs,".csv"))
  
}

#################################################################################
#
# END
#
#################################################################################
# 