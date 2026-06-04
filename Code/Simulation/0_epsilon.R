#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE:  Find epsilon such as the simulation has the same R^2 as the paper 
# AUTHOR:		Anirudh Sankar & Louis-Mael Jean & Elsa Trezeguet
# CREATED:	08/12/2021
# LAST MODIFIED: 08/12/2021
# STATUS: 	Draft

#######################################################################################

rm(list = ls())

library('data.table')
library('dplyr')
library('clubSandwich')
library('stargazer')
library('Hmisc')
library('readstata13')
library('car')
library('miceadds')
library('multiwayvcov')
library('estimatr')
library('hdm')
library('lme4')
library('stringr')
library('hash')

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



###### Selection of the Outcome: 
outcome = "shot_Measles1"
#outcome = "shots_per_dollar"

###########################################################
## R^2 in the real data set : 
#############################################################


#--READ VILLAGEXMONTH dataset
villagexmonth_level <- fread(paste0(path_prepared_data,"/Tablet_VillageXMonth_Costs.csv"),header = TRUE, sep = ",", data.table = FALSE)

#--SELECT SEEDS RISK EXPT
villagexmonth_level <- villagexmonth_level %>% filter(seedsrisk == 1)

#--SELECT ONLY FIRST IMPLEMENTATION
villagexmonth_level <- villagexmonth_level %>% filter(first_implementation == 1)



#--CREATE ALL POLICIES -- 
source('create_sp_variables.R')


#District-Time FE
villagexmonth_level$fes <- group_indices(villagexmonth_level, id_district, created_year, created_month)

#create FES dummies

fes_dummies <- data.frame(lme4::dummy(villagexmonth_level$fes))

villagexmonth_level <- cbind(villagexmonth_level, fes_dummies)

#-- LASSO SELECTION

variables <- grep("^SP", names(villagexmonth_level), value = TRUE)
variables <- variables[variables != "SP_noSeedXnoIncentiveXnoReminder"]

variables_expanded <- c(variables, colnames(fes_dummies))


current_variables <- variables_expanded

deselect_list <- c()

deselect_pval <- c()

pval_cutoff <- 5 * 10^(-13)


current_sp_formula <- as.formula(paste0(outcome,"~",paste0(current_variables,collapse = "+")))

current_model_ols <- estimatr::lm_robust(formula = current_sp_formula, data = villagexmonth_level, weights = village_population, se_type = "classical")

current_max_pval <- max(current_model_ols$p.value)

print(current_max_pval)


while (current_max_pval > pval_cutoff) {
  
  deselect_name <- names(current_model_ols$p.value[which.max(current_model_ols$p.value)])
  
  deselect_list <- c(deselect_list, deselect_name)
  
  deselect_pval <- c(deselect_pval, current_model_ols$p.value[which.max(current_model_ols$p.value)])
  
  current_variables <- current_variables[current_variables != deselect_name]
  
  current_sp_formula <- as.formula(paste0(outcome,"~",paste0(current_variables,collapse = "+")))
  print(current_sp_formula)
  
  current_model_ols <- estimatr::lm_robust(formula = current_sp_formula, data = villagexmonth_level, weights = village_population, se_type = "classical")
  
  current_max_pval <- max(current_model_ols$p.value)
  print(current_max_pval)
}

support_SP <- variables_expanded[!(variables_expanded %in% deselect_list)]

toRemove <- grep(pattern = "ntercept",x = support_SP) #don't want intercept to be lasso selected!
if (length(toRemove) > 0) {
  support_SP <- support_SP[-toRemove]
}


fes_chosen <- grep("^X", support_SP, value = TRUE)
support_SP_policies <- grep("^X", support_SP, value = TRUE, invert = TRUE)

#
# 
# # -- AUTOMATED POOLING
# 

source("pooling_functions.R")

treatment_profiles <- list(c("random", "noIncentive", "noReminder"), #seeds only
                           c("trusted", "noIncentive", "noReminder"),
                           c("gossip", "noIncentive", "noReminder"),
                           
                           c("noSeed", "slope", "noReminder"), #incentives only
                           c("noSeed", "flat", "noReminder"),
                           
                           c("noSeed", "noIncentive", "SMS"), #SMS only
                           
                           c("random", "slope", "noReminder"), #seeds and incentives (slopes) only
                           c("trusted", "slope", "noReminder"),
                           c("gossip", "slope", "noReminder"),
                           
                           c("random", "flat", "noReminder"), #seeds and incentives (flats) only
                           c("trusted", "flat", "noReminder"),
                           c("gossip", "flat", "noReminder"),
                           
                           c("random", "noIncentive", "SMS"), #seeds and SMS only
                           c("trusted", "noIncentive", "SMS"),
                           c("gossip", "noIncentive", "SMS"),
                           
                           c("noSeed", "slope", "SMS"), #incentives (slope) and SMS only
                           c("noSeed", "flat", "SMS"), #incentives (flats) and SMS only
                           
                           c("random", "slope", "SMS"), #seeds and incentives (slopes) and SMS
                           c("trusted", "slope", "SMS"),
                           c("gossip", "slope", "SMS"),
                           
                           c("random", "flat", "SMS"), #seeds and incnetives (flats) and SMS
                           c("trusted", "flat", "SMS"),
                           c("gossip", "flat", "SMS"))

final_data <- villagexmonth_level

for (tp_name in treatment_profiles) {
  tp_support <- get_relevant_sp_in_tp(support_SP, tp_name)
  final_data <- add_pooled_policies_tp(final_data, tp_support, tp_name)
}

pooled_policies <- grep("^POOLED_", colnames(final_data), value = TRUE)


source("map_key_policy_names.R")

#
# POSTLASSO 
#

variables_pooled_expanded <- c(pooled_policies,fes_chosen)


formula_pl <- as.formula(paste0(outcome,"~",paste0(variables_pooled_expanded ,collapse = "+")))
model_pl <- estimatr::lm_robust(formula = formula_pl, data = final_data, clusters = id_sc, weights = village_population, se_type = "CR0")

real_R_squared <- model_pl$r.squared


############################################################################################
### Reproduction of R^2 in simulations: 
############################################################################################

#path_LM2 = "~/Dropbox (MIT)/Smart Pooling and Pruning/Code/Simulation"
#setwd(path_LM2)

set.seed(NULL)
set.seed(100)

source('crossprod_functions.R')

source('simulation_functions.R')

R = c(5,5,3)
M = 3

epsilon_list <- linspace(0.5,1.5,11)


policy_vectors <- create_policy_vectors(R,M)
G <- create_sp_to_unique_transformation(R,M, policy_vectors)


#omit the pure control

policy_vectors <- policy_vectors[-1] #omit_control_policy
G <- G[-1,-1]


pos_indices_universe <- c()


for (i in 1:length(policy_vectors)) {
  candidate <- policy_vectors[[i]]
  if ((sum(candidate == 1) != 0)) { #first M-1 not where all treatments are on
    pos_indices_universe <- c(pos_indices_universe,i)
  }
}

#ncomb = 5;

pol_sets <- list()

# for (i in 1:ncomb) {
#   pol_sets[[i]] <- sample(pos_indices_universe,M-1, replace = FALSE)
#   pol_sets[[i]] <- c(pol_sets[[i]], Position(function(x) identical(x, c(2,2,2,2)), policy_vectors))
# }

pos_treat_effects <- linspace(1,2,M) #policy effects from 1 to 2 (in the end its going to be a mulitple of standard deviation)

nsim<- 20
nobs <- 1000
stdev <- 1

pos_indices <- sample(pos_indices_universe,M-1, replace = FALSE)
pos_indices <- c(pos_indices, Position(function(x) identical(x, replicate(M,2)), policy_vectors)) #just taking a random combination



sp_treat_effects_pop <- rep.int(0, length(policy_vectors)) #vector of smart pooling effects
sp_treat_effects_pop[pos_indices] <- pos_treat_effects*stdev#zero except at pos_indices
unique_treat_effects_pop <- G %*% sp_treat_effects_pop #note this is for the slim policies

## The true best policies of the model:
best_policy <- names(unique_treat_effects_pop[,1])[which(unique_treat_effects_pop == max(unique_treat_effects_pop))]
efficient_best_policy <- best_policy[1]
efficient_best_policy_coef <- unique_treat_effects_pop[which(names(unique_treat_effects_pop[,1]) == efficient_best_policy)]


support_pop <- policy_vectors[pos_indices] 
support_pop_binned <- retrieve_binned_support(support_pop, R, M)
collapsed_pop_info <- create_collapsed_policies(support_pop_binned, R, M)

collapsed_pop_trueeffects <- get_pop_treat_effects(collapsed_pop_info, unique_treat_effects_pop, R, M) 


epsilon_R_squared <- list()
mean_epsilon_R_squared <- c()

for(i in 1:length(epsilon_list)){
  epsilon <- epsilon_list[i]
  
  epsilon_R_squared_i <- c()
  mean_epsilon_R_squared[i] <- 0
  
  for(sim in 1:nsim){
    print("sim")
    print(sim)
    
    #DATA _ ORIGINAL ASSIGNMENTS
    original_treatments <- create_original_treatment_assignments(R,M,nobs)
    
    #DATA - SP MATRIX
    sp_matrix <- create_sp_matrix(R,M, original_treatments, policy_vectors)
    
    # OMIT THE pure control
    #sp_matrix <- sp_matrix[,-1]
    
    beta_matrix <- sp_matrix %*% solve(G) # Y= X alpha + e becomes Y = X G^-1 G alpha + e - THIS IS THE UNIQUE SLIM POLICIES
    
    stdev = 1
    epsilon = 2.3*stdev
    
    #outcome generation
    y <- rep.int(0,nobs)
    
    for (j in 1:length(policy_vectors)) {
      y <- y + sp_treat_effects_pop[j] * stdev * sp_matrix[,j] #keep it releative to standard deviation - best policy is 2*stdev
    }

    y  = y + rnorm(nobs,0,epsilon)
    
    #### Puffer estimate: 
    support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors)
    support_acc_backward <- length(intersect(support_pop, support_data_backward))/length(union(support_pop, support_data_backward))
    
    support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)
    collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
    
    collapsed_data_df_backward<- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)
    
    model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward)
  
    R_squared_i <- model_data_pl_backward$r.squared
    
    epsilon_R_squared_i[sim] <- R_squared_i
    mean_epsilon_R_squared[i] <-  mean_epsilon_R_squared[i]+ R_squared_i
  }
  mean_epsilon_R_squared[i] <- mean_epsilon_R_squared[i]/nsim
  epsilon_R_squared[[i]] <- epsilon_R_squared_i
}
