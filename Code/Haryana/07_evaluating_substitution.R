#Author: Laure Heidmann
#Name: 4_Evaluating_Substitution
#Description: Evaluating if the treatment effect is a true effect or a substitution effect
#Created: LH 17/07/2018
#Updated: LH 17/07/2018

#Name in old directory: 4_Evaluating_Substitution.R 

rm(list=ls())

library('this.path')
library('clubSandwich')

setwd(this.path::here())
source("ecma_directory.R")

datadir <- paste0(ecmadir,"Data/")
path_input_data <- paste0(datadir,"Input Data/")
path_prepared_data <- paste0(datadir,"Prepared Data/")
path_output_data <- paste0(datadir, "Output Data/")
path_figures <- paste0(ecmadir,"Figures/")
path_tables <- paste0(ecmadir,"Tables/")
path_functions <- paste0(ecmadir,"Code/Helper_functions/")


source(file = paste0(path_functions,"/AppendixI_Regression_Functions.R"))
source(file = paste0(path_functions,"/AppendixI_2_Treatments_Effects_New_Outcomes.R"))

#1. Simple linear regression
#2. PROBIT regression


# using the observations from endline survey
endline_data <- fread(paste0(path_input_data,"/Prepared_Endline.csv"),header = TRUE, sep = ",", data.table = FALSE)

# selecting outcomes
outcomes <- c(paste0("atleast",c(2,3,4,5,6,7)), "atleastMeasles")
outcomes_label <- c(paste0("At Least ",c(2,3,4,5,6,7)), "Measles 1")
regression_title <- "Immunization Outcomes Restricted to Unmatched Children"

# taking into account only children that are old enough to get the number of vaccines in consideration
for (l in c(1:7,"Measles")){
  endline_data[which(endline_data[,paste0("age_vacc_",l)] != 1), paste0("atleast",l)] <- NA
}

#### 1. ALL SUB-TREATMENTS
directory <- paste0(path_tables,"Substitution/")
invisible(running_new_outcomes(endline_data[which(endline_data$matched==0),], outcomes, outcomes_label, regression_title, directory, 1))
