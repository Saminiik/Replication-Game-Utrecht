#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Appendix C: Failure of irrepresentability when TVA is not pre-conditionned with Puffer
# AUTHOR:		Anirudh Sankar
# CREATED:	XX
# LAST MODIFIED: 07 August 2023 (LMJ)

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
library('purrr')
library(xtable)
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

source("crossprod_functions.R")
source("simulation_functions.R")


#################################################################################
#
# I.SIMULATION
#
#################################################################################


##--- Simulation Setup
set.seed(25)
nobs = 10000

R_grid = c(3,4,5)
M_grid <- c(2,3,4)

irrep_table <- data.frame()


##--- Run Simulation
for (R in R_grid) {
  for(M in M_grid) {
    cat("R =",R, ", M =",M,"\n")
    
    policy_vectors <- create_policy_vectors(R,M)
    
    #DATA _ ORIGINAL ASSIGNMENTS
    original_treatments <- create_original_treatment_assignments(R,M,nobs)
    
    #DATA - SP MATRIX
    sp_matrix <- create_sp_matrix(R,M, original_treatments)
  
    G <- create_sp_to_unique_transformation(R,M)
    
    # OMIT THE pure control
    sp_matrix <- sp_matrix[,-1]
    policy_vectors <- policy_vectors[-1] #omit_control_policy
    G <- G[-1,-1]
    beta_matrix <- sp_matrix %*% solve(G) # Y= X alpha + e becomes Y = X G^-1 G alpha + e - THIS IS THE UNIQUE SLIM POLICIES
    
    sp_matrix_centered <- apply(sp_matrix, 2, function(x) (x - mean(x)))
    sp_matrix_norm2 <- apply(sp_matrix_centered, 2, function(x) (x/sqrt(sum(x^2/nobs))) )
    
    C = (1/nobs) * t(sp_matrix) %*% sp_matrix
    C_norm = (1/nobs) * t(sp_matrix_norm2) %*% sp_matrix_norm2
    
    beta_matrix_centered <- apply(beta_matrix, 2, function(x) (x - mean(x)))
    beta_matrix_norm2 <- apply(beta_matrix_centered, 2, function(x) (x/sqrt(sum(x^2)/nobs)) )
    
    C_norm_beta = (1/nobs) * t(beta_matrix_norm2) %*% beta_matrix_norm2
    C_norm_centered = (1/nobs) * t(beta_matrix_centered) %*% beta_matrix_centered
    
    zero <- dim(sp_matrix_norm2)[2] #R-1 x R-1...R-1 highest intensity SP variable
    
    design_1 <- sp_matrix_norm2[, -zero] #X_n(1)
    design_2 <- sp_matrix_norm2[, zero] #X_n(2)
    
    C_norm_21 <- (1/nobs) * t(design_2) %*% design_1
    C_norm_11 <- (1/nobs) * t(design_1) %*% design_1
    C_norm_12 <- (1/nobs) * t(design_1) %*% design_2
    C_norm_22 <- (1/nobs) * t(design_2) %*% design_2
    
    design_1_unscaled <- sp_matrix_centered[,-zero]
    design_2_unscaled <- sp_matrix_centered[,zero]
    
    C_norm_21_unscaled <- (1/nobs) * t(design_2_unscaled) %*% design_1_unscaled
    C_norm_11_unscaled <- (1/nobs) * t(design_1_unscaled) %*% design_1_unscaled
    C_norm_12_unscaled <- (1/nobs) * t(design_1_unscaled) %*% design_2_unscaled
    C_norm_22_unscaled <- (1/nobs) * t(design_2_unscaled) %*% design_2_unscaled
    
    C_condition <-  C_norm_21 %*% solve(C_norm_11) 
    C_condition_filtered <- ifelse(abs(C_condition)<1e-10,0, C_condition)
    
    C_condition_unscaled <-  C_norm_21_unscaled %*% solve(C_norm_11_unscaled) 
    C_condition_filtered_unscaled <- ifelse(abs(C_condition_unscaled)<1e-10,0, C_condition_unscaled)
    
    key_l1 <- sum(abs(C_condition_filtered))
    key_l1_unscaled <- sum(abs(C_condition_filtered_unscaled))
    
    print(key_l1_unscaled)
    
    irrep_table <- rbind(irrep_table,c(R,M,key_l1,key_l1_unscaled))
  }
}


##--- Format and print table
colnames(irrep_table) <- c("R", "M", "L1 norm (standardized vars)", "L1 norm (unstandardized vars)")
print(xtable(irrep_table, type = "latex", digits = c()), file = paste0(path_tables,"/Simulation/Sim-Irrepresentability-Failure.tex"), include.rownames=FALSE)


#################################################################################
#
# END
#
#################################################################################


