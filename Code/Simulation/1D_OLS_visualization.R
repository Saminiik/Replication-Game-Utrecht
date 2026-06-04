#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Show that even when support accuracy is 100%, OLS visualization is bad
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
library('pracma')
library("stringr")
library("broom")
library('plyr')
library('ggpubr')
library('scales')
library('car')
library('ggnewscale')
library('RColorBrewer')
library('colorspace')
library('reshape2')

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

set.seed(75)

#################################################################################
#
# I.GENERATE ORACLE INFORMATION FOR SIMULATIONS
#
#################################################################################


R = c(5,5,3)
M = 3

stdev <- 1
epsilon = stdev

#Generate policy vectors
policy_vectors <- create_policy_vectors(R,M)
G <- create_sp_to_unique_transformation(R,M)

#Omit the pure control
policy_vectors <- policy_vectors[-1] #omit_control_policy
G <- G[-1,-1]

#pos_treat_effects <- linspace(1,1.2,M) #policy effects
pos_treat_effects <- c(0.3,0.8,1.2,1.5)

pos_indices <- c(4,6,48,65)

#derived quantities
sp_treat_effects_pop <- rep.int(0, length(policy_vectors)) #vector of smart pooling effects
sp_treat_effects_pop[pos_indices] <- pos_treat_effects * stdev #zero except at pos_indices - effects scaled by stdev
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
#nobs_list <- round(logspace(3,4,10))
nobs = 4000
nsim = 300
#nsim = 10

ols_coefficients_list <- list()
puffer_coefficients_list <- list()
support_accuracies <- c()

ols_models <- list()
puffer_models <- list()

k = 1
for (sim in 1:nsim) {
  cat("\n Simulation ", sim)
  
  #DATA _ ORIGINAL ASSIGNMENTS
  original_treatments <- create_original_treatment_assignments(R,M,nobs)
  
  #DATA - SP MATRIX
  sp_matrix <- create_sp_matrix(R,M, original_treatments, policy_vectors) #no need to omit control since we feed policy_vectors
  
  beta_matrix <- sp_matrix %*% solve(G) # Y= X alpha + e becomes Y = X G^-1 G alpha + e - THIS IS THE UNIQUE SLIM POLICIES
  
  #outcome generation
  y <- rep.int(0,nobs)
  
  for (j in 1:length(policy_vectors)) {
    y <- y + sp_treat_effects_pop[j]* sp_matrix[,j]
  }
  
  y  = y + rnorm(nobs,0,epsilon)
  
  support_data_backward <- retrieve_backward_elimination_support(y,sp_matrix, policy_vectors) 
  support_acc_backward <- length(intersect(support_pop, support_data_backward))/length(union(support_pop, support_data_backward)) 
  
  cat("\n \tSupport acc = ",support_acc_backward)
  support_accuracies <- c(support_accuracies,support_acc_backward)
  
 
  if (support_acc_backward == 1){
    
    
    support_data_binned_backward <- retrieve_binned_support(support_data_backward, R, M)
    collapsed_data_info_backward <- create_collapsed_policies(support_data_binned_backward, R, M)
    collapsed_data_df_backward   <- get_collapsed_df(collapsed_data_info_backward, y, beta_matrix, R, M)
    
    model_data_pl_backward<- estimatr::lm_robust(formula = as.formula("y~ ."), data = collapsed_data_df_backward)
    Puffer_coefficients <- model_data_pl_backward$coefficients[-1]
    names(Puffer_coefficients) = c("Pool 1", "Pool 2", "Pool 3", "Pool 4")
    
    puffer_coefficients_list[[k]] <- Puffer_coefficients
    puffer_models[[k]] <- model_data_pl_backward
    
    #------------------------------------------------#
    # OLS ON BETA MATRIX (unique policy regression)
    #------------------------------------------------#
    
    beta_df_pop <- as.data.frame(cbind(y, beta_matrix))

    #OLS model
    model_beta_ols <- estimatr::lm_robust(formula = as.formula("y~ ."), data = beta_df_pop)
    model_terms <- model_beta_ols$term[-1] %>% str_remove("`") %>% str_remove("`")
    
    #Get best coefficients / best policies
    ols_coefficients <- model_beta_ols$coefficients[-1]
    names(ols_coefficients) = model_terms
    
    ols_coefficients_list[[k]] <-  ols_coefficients
    ols_models[[k]] <- model_beta_ols
    
    k <- k + 1
  }
  
}

ols_coefficients_df <- do.call("rbind", ols_coefficients_list) %>% as.data.frame()
ols_coefficients_long <- reshape2::melt(ols_coefficients_df, value.name = "estimate") %>% rename(replace = c("variable" = "Policies"))

puffer_coefficients_df <- do.call("rbind", puffer_coefficients_list) %>% as.data.frame()
puffer_coefficients_long <- reshape2::melt(puffer_coefficients_df, value.name = "estimate") %>% rename(replace = c("variable" = "Policies"))

for (q in 2:dim(collapsed_data_df_backward)[2]){
  assign(paste0("policy",q-1), str_split(names(collapsed_data_df_backward), "X")[[q]])
}

ols_coefficients_long <- ols_coefficients_long %>% mutate("Pooling" = case_when(
  Policies %in% policy1 ~ "Pool1",
  Policies %in% policy2 ~ "Pool2",
  Policies %in% policy3 ~ "Pool3",
  Policies %in% policy4 ~ "Pool4",
  TRUE ~ "Pruned"))


write_csv(ols_coefficients_long,paste0(path_output_data,"Simulation_Data/1D_OLS_density.csv"))
write_csv(puffer_coefficients_long,paste0(path_output_data,"Simulation_Data/1D_Puffer_density.csv"))


#Only one needs to be run for this one
m = 1 #pick a striking example

model_beta_ols = ols_models[[m]]
model_beta_ols$term = model_beta_ols$term %>% str_remove("`") %>% str_remove("`")
results_ols <- model_beta_ols %>% tidy()
results_ols <- results_ols[-1,]



results_ols = results_ols %>% mutate("Pooled policies" = case_when(
  term %in% policy1 ~ "Pool 1",
  term %in% policy2 ~ "Pool 2",
  term %in% policy3 ~ "Pool 3",
  term %in% policy4 ~ "Pool 4",
  TRUE ~ "Pruned"))

model_data_pl_backward = puffer_models[[m]]
model_data_pl_backward$term = model_data_pl_backward$term %>% str_remove("`") %>% str_remove("`") 
results_backward <- model_data_pl_backward %>% tidy() %>% mutate("term" = c("Intercept", "Pool 1", "Pool 2", "Pool 3", "Pool 4"))
results_backward<- results_backward[-1,]

#Pool 3 vs Pool 1 (two lowest)
t1 = linearHypothesis(model_data_pl_backward, "`c(1, 2, 2)Xc(1, 3, 2)Xc(1, 4, 2)Xc(1, 5, 2)Xc(1, 2, 3)Xc(1, 3, 3)Xc(1, 4, 3)Xc(1, 5, 3)` = `c(1, 3, 1)Xc(1, 4, 1)Xc(1, 5, 1)`")
p1 = format(round(t1$`Pr(>Chisq)`[2],5), nsmall = 5)

#Pool 1 vs Pool 2 (middle ones)
t2 = linearHypothesis(model_data_pl_backward, "`c(1, 3, 1)Xc(1, 4, 1)Xc(1, 5, 1)`  = `c(4, 2, 1)Xc(5, 2, 1)Xc(4, 3, 1)Xc(5, 3, 1)Xc(4, 4, 1)Xc(5, 4, 1)Xc(4, 5, 1)Xc(5, 5, 1)`")
p2 = format(round(t2$`Pr(>Chisq)`[2],5), nsmall = 5)

#Pool 2 vs Pool 4 (strong ones)
t3 = linearHypothesis(model_data_pl_backward, "`c(4, 2, 1)Xc(5, 2, 1)Xc(4, 3, 1)Xc(5, 3, 1)Xc(4, 4, 1)Xc(5, 4, 1)Xc(4, 5, 1)Xc(5, 5, 1)` = `c(5, 2, 3)Xc(5, 3, 3)Xc(5, 4, 3)Xc(5, 5, 3)`")
p3 = round(t3$`Pr(>Chisq)`[2],5)

cols = hue_pal()(5)

setwd(path_figures)

p_ols_coefplot <- ggplot(results_ols, aes(x=term, y=estimate, colour = `Pooled policies`)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high), width=.1, lwd = 1) + 
  scale_colour_manual(name = "Legend",values = c("Pruned" = cols[1], "Pool 1" = cols[2],"Pool 2" = cols[3],"Pool 3" = cols[4],"Pool 4" = cols[5])) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        plot.margin = margin(t=5, r = 5, b = 5, l=20)) + 
  xlab("Policy combinations") +
  ylab("Estimated effect") +
  ggtitle("OLS estimates on unique policy regression") +
  theme(plot.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15, face = "bold"))

ggsave("Simulation_Figures/FigureE1/1D_ols_coefplots.pdf", plot=p_ols_coefplot, width = 18, height = 6)


p_backward_coefplot <- ggplot(results_backward, aes(x=fct_reorder(term, c(2,3,1,4)), y=estimate, colour = term)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high), width=.1, lwd = 1) + 
  geom_bracket(data = results_backward, xmin = "Pool 3", xmax = "Pool 1", y.position = 1, label = paste0(" p = ",p1), color = "black", size = 0.11, tip.length = 0.02, label.size = 5) + 
  geom_bracket(data = results_backward, xmin = "Pool 1", xmax = "Pool 2", y.position = 1.4, label = paste0("p = ",p2), color = "black", size = 0.11, tip.length = 0.02,label.size = 5) + 
  geom_bracket(data = results_backward, xmin = "Pool 2", xmax = "Pool 4", y.position = 1.7, label = paste0("p = ",p3), color = "black", size = 0.11, tip.length = 0.02,label.size = 5) + 
  scale_colour_manual(name = "Legend",values = c("Pool 1" = cols[2],"Pool 2" = cols[3],"Pool 3" = cols[4],"Pool 4" = cols[5])) + 
  scale_y_continuous(limits = c(0, 2)) +
  xlab("Policy combinations") +
  ylab("Estimated effect") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.margin = margin(t=5, r = 5, b = 5, l=20),
        plot.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  ggtitle("TVA estimates on unique policy regression")

ggsave("Simulation_Figures/FigureE1/1D_puffer_coefplots.pdf", plot=p_backward_coefplot, width = 18, height = 6)

detach("package:reshape2",unload=TRUE)
#p_coefs <- ggarrange(p_ols_coefplot,p_backward_coefplot,labels = c("A", "B"), nrow = 2, common.legend  = TRUE, legend = "right")

#ggsave("1D_combined_coefplots.pdf", plot=p_coefs, width = 14, height = 8.5)


# setwd(path_data)
# write.csv(results_ols,"Simulation_Data/1D_OLS_coef.csv")


#################################################################################
#
# END
#
#################################################################################

