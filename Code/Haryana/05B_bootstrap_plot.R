#######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE:  Plots for Bootstrapping 
# AUTHOR:		Louis-Mael Jean & Elsa Trezeguet
# CREATED:	08/11/2021
# MODIFIED: 
# STATUS: 	Draft

#Name in old directory: 05B_bootstrap_plot.R 

#######################################################################################



#################################################################################
#
# 0.ENVIRONMENT SET UP
#
#################################################################################

library('dplyr')
library('ggnewscale')
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
setwd(path_output_data)

set.seed(NULL)
set.seed(534)


#################################################################################
#
# I. DRAW PLOT
#
#################################################################################


for (outcome in c("shot_Measles1","shots_per_dollar")){
    
  simulations_data <- read.csv(paste0("Bootstrap_Results/Bootstrap_simulations_data_",outcome,".csv")) %>%
    mutate(pl_names = sub(",", ",\n",pl_names))
  best_policies_data <- read.csv(paste0("Bootstrap_Results/_best_policies_data_",outcome,".csv")) %>%
    mutate(bootstrapped_best_policy_name = sub(",", ",\n",bootstrapped_best_policy_name))
  
  #Get some variables back based on the data
  initial_best_policy = filter(best_policies_data, Is_true == 1)$bootstrapped_best_policy_name
  nsamples = dim(best_policies_data)[1] - 1
  best_policy_accuracy = best_policies_data$best_policy_accuracy %>% unique()
  
  
  #Recomputing distinct supports counting initial as distinct
  distinct_supports = simulations_data$support_category %>% unique() %>% length()
  
  
  #Legend depends on the number of distinct supports
  scales  <- switch(toString(distinct_supports),
                    "2" = data.frame(cols = c("0" = "black","1" = "#0F425CFF"), shapes = c("0" = 22, "1" = 16), sizes = c("0" = 5, "1" = 5)),
                    "3" = data.frame(cols = c("0" = "black","1" = "#0F425CFF","2" = "#CC8214FF"), shapes = c("0" = 22, "1" = 16, "2" = 16), sizes = c("0" = 5, "1" = 5, "2" = 5)),
                    "4" = data.frame(cols = c("0" = "black","1" = "#0F425CFF","2" = "#CC8214FF", "3" = "forestgreen"), shapes = c("0" = 22, "1" = 16, "2" = 16, "3" = 16), sizes = c("0" = 5, "1" = 5, "2" = 5, "3" = 5)))
  
  axis_cols = ifelse(levels(factor(simulations_data$pl_names)) == initial_best_policy,"#800000FF","black")
  
  total_graph <- ggplot() +
    geom_point(data = simulations_data ,aes(x = pl_names, y= pl_effects, color = factor(support_category), shape = factor(support_category), size = factor(support_category)), alpha = 0.8) +
    scale_color_manual(name = "Selected pooled policies", values = scales$cols, labels = c("Data Support", "Bootstrap Support 1", "Bootstrap Support 2", "Bootstrap Support 3"),
                       guide = guide_legend(order = 1)) +
    scale_size_manual(name = "Selected pooled policies", values = scales$sizes,labels = c("Data Support", "Bootstrap Support 1", "Bootstrap Support 2", "Bootstrap Support 3"),
                      guide = guide_legend(order = 1)) + 
    scale_shape_manual(name = "Selected pooled policies", values = scales$shapes,labels = c("Data Support", "Bootstrap Support 1", "Bootstrap Support 2", "Bootstrap Support 3"),
                       guide = guide_legend(order = 1)) + 
    new_scale_color() +
    geom_point(data = best_policies_data, aes(x=bootstrapped_best_policy_name,y=bootstrapped_best_coef, color = factor(Is_true)), size = 1.5) +
    scale_color_manual(name = "Selected best policy",values = c("0" = "firebrick", "1" = "chartreuse3"), 
                       labels = c("Bootstrap winner's curse\nadjusted estimate", "Data winner's curse\nadjusted estimate" ),
                       guide = guide_legend(order = 2)) + 
    geom_point(data = filter(simulations_data, support_category == 0),aes(x = pl_names, y= pl_effects), color = "black", shape = 22, size = 7, fill = "white") +
    geom_point(data = filter(best_policies_data, Is_true == 1), aes(x=bootstrapped_best_policy_name,y=bootstrapped_best_coef), color = "chartreuse3", shape = 16, size = 3) +
    theme_minimal() + 
    theme(axis.text.x = element_text(colour = axis_cols, size = 12, angle = 90),
          axis.title = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
          legend.title = element_text(size=14, face = "bold"),
          legend.text = element_text(size=12, margin = margin(b = 5, t = 5, unit = "pt")),
          plot.title = element_text(size = 16, face = "bold")) + 
    guides(fill = guide_legend(byrow = TRUE)) + 
    ggtitle(paste0("Post-LASSO Estimates for Bootstrapped Samples (", nsamples," simulations).\nBest Policy Selection Accuracy = ",round(best_policy_accuracy,3))) +
    ylab("Treatment Effects") + 
    xlab("Pooled policy name")
  
  ggsave(paste0(path_figures,"Bootstrap/bootstrapping_",outcome,"_n=", nsamples, "_WC.pdf"), plot = total_graph, width = 14, height = 8)
   
}






