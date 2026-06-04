#######################################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE: 	Constructs the tables and figures from all simulations
# AUTHOR:	  Louis-Mael Jean & Elsa Trezeguet
# CREATED:	Feb 14 2022
# LAST MODIFIED: 
# STATUS: 	Draft

#######################################################################################################


#################################################################################
#
# 0.SET UP
#
#################################################################################

rm(list = ls())

library('tidyverse')
library('pracma')
library('gtools')
library('ggsci')
library('numbers')
library('glmnet')
library('purrr')
library('hdm')
library('pracma')
library('ggpubr')
library('scales')
library('ggnewscale')
library('this.path')

rename = dplyr::rename

setwd(this.path::here())
source("ecma_directory.R")

datadir <- paste0(ecmadir, "Data/")
path_input_data <- paste0(datadir,"Input Data/")
path_prepared_data <- paste0(datadir,"Prepared Data/")
path_output_data <- paste0(datadir, "Output Data/")
path_figures <- paste0(ecmadir,"Figures/")
path_tables <- paste0(ecmadir,"Tables/")
path_functions <- paste0(ecmadir,"Code/Helper_functions/")


setwd(path_figures)

#################################################################################
#
# I. OLS PLOTS
#
#################################################################################

output_wide <- read.csv(paste0(path_output_data,"/Simulation_Data/1B_OLS_wide_harmonized_any.csv"))
output_long <- read.csv(paste0(path_output_data,"/Simulation_Data/1B_OLS_long_harmonized_any.csv"))

#Boxplot data
c1 = output_long %>% dplyr::select(c("n","comb","sim","True_best","OLS_best", "OLS_hybrid", "OLS_abs_diff","OLS_prct_diff")) %>% 
  mutate("method" = "Hybrid OLS") %>% 
  rename("estimate" = "OLS_best", "hybrid_estimate" = "OLS_hybrid",
         "absolute_shrunk" = "OLS_abs_diff", "percentage_shrunk" = "OLS_prct_diff")

c2 = output_long %>% dplyr::select(c("n","comb","sim","True_best","SP_best", "SP_hybrid", "SP_abs_diff","SP_prct_diff")) %>% 
  mutate("method" = "Hybrid TVA") %>% 
  rename("estimate" = "SP_best", "hybrid_estimate" = "SP_hybrid",
         "absolute_shrunk" = "SP_abs_diff", "percentage_shrunk" = "SP_prct_diff")


df_boxplot = rbind(c1,c2)

#-------------------------------#
# 1B: OLS MSE
#-------------------------------#

p_ols1 <- ggplot(output_wide) + 
  geom_point(aes(y = mse_ols_over_n , x = nobs, color = "Best policy (OLS)", fill = "Best policy (OLS)", shape = "Best policy (OLS)"), size = 3) +
  geom_point(aes(y = mse_hybrid_ols_over_n, x = nobs, color = "Hybrid best\npolicy (OLS)", fill = "Hybrid best\npolicy (OLS)", shape = "Hybrid best\npolicy (OLS)"), size = 3) +
  geom_point(aes(y = mse_hybrid_puffer_over_n, x = nobs, color = "Hybrid best\npolicy (TVA)", fill = "Hybrid best\npolicy (TVA)", shape = "Hybrid best\npolicy (TVA)"), size = 3) +
  xlab('Sample size n') +
  ylab('MSE') + theme_bw() + 
  ggtitle("MSE") +
  scale_color_manual(name = 'Legend', guide = 'legend', values = c('Best policy (OLS)' = "#800000FF", "Hybrid best\npolicy (OLS)" = '#FFA319FF',  "Hybrid best\npolicy (TVA)" = '#115F83FF')) +
  scale_fill_manual(name = 'Legend', values = c('Best policy (OLS)' = "#800000FF", "Hybrid best\npolicy (OLS)" = '#FFA319FF',  "Hybrid best\npolicy (TVA)" = '#115F83FF'))+
  scale_shape_manual(name = 'Legend', values = c('Best policy (OLS)' = 22, "Hybrid best\npolicy (OLS)" = 24,  "Hybrid best\npolicy (TVA)" = 21)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  scale_y_continuous(breaks = seq(0,6.5,by=0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 14),
        axis.text = element_text(size = 12), axis.title = element_text(size =  14))

ggsave("Simulation_Figures/Figure3/1B_OLSvsPuffer_MSE_any.pdf", plot = p_ols1, width = 8, height = 6)

#----------------------------------#
# 1B: OLS BEST POLICY INCLUSION
#----------------------------------#

p_ols2 <- ggplot(output_wide) + 
  geom_jitter(aes(y = ols_some_best_acc_over_n, x = nobs, color = "Some best (OLS)", fill = "Some best (OLS)", shape = "Some best (OLS)"),size = 3,height=0,width=250) +
  geom_point(aes(y = ols_min_best_acc_over_n, x = nobs, color = "Minimum best (OLS)", fill = "Minimum best (OLS)", shape = "Minimum best (OLS)"),size = 3) +
  geom_jitter(aes(y = puffer_some_best_acc_over_n, x = nobs, color = "Some best (TVA)", fill = "Some best (TVA)", shape = "Some best (TVA)"), size = 3,height=0,width=250) +
  geom_point(aes(y = puffer_min_best_acc_over_n, x = nobs, color = "Minimum best (TVA)", fill = "Minimum best (TVA)", shape = "Minimum best (TVA)"), size = 3) +
  geom_hline(aes(yintercept = theoretical_random_rate[1],linetype = "Theoretical random \n selection rate"), color = "black") +
  xlab('Sample size n') +
  ylab('Accuracy')+ 
  scale_color_manual(name = 'Performance Metric', values = c("Some best (OLS)" = "#800000FF",  "Minimum best (OLS)" = "#FFA319FF",
                                                             "Some best (TVA)" = "#115F83FF",  "Minimum best (TVA)" = "seagreen3")) +
  scale_fill_manual(name = 'Performance Metric', values = c("Some best (OLS)" = "#800000FF",  "Minimum best (OLS)" = "#FFA319FF",
                                                            "Some best (TVA)" = "#115F83FF",  "Minimum best (TVA)" = "seagreen3")) +
  scale_shape_manual(name = 'Performance Metric',values = c("Some best (OLS)" = 21,  "Minimum best (OLS)" = 21,"Some best (TVA)" = 23,  "Minimum best (TVA)" = 23)) + 
  scale_linetype_manual(name = "Legend", values = c("Theoretical random \n selection rate" = 2)) +
  ggtitle("Best policy inclusion") +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) + theme_bw() + 
  scale_y_continuous(limits = c(0, 1), breaks = linspace(0,1,11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 14),
        axis.text = element_text(size = 12), axis.title = element_text(size =  14))

ggsave("Simulation_Figures/Figure3/1B_OLS_inclusion_accuracy_any.pdf", plot = p_ols2, width = 8, height = 6)

#----------------------------------#
# 1B: OLS BEST POLICY DISTRIBUTION
#----------------------------------#

p_ols3 <- ggplot(filter(output_long,comb==1)) +
  geom_boxplot(mapping = aes(y = OLS_best, x = n, group = n, fill = "OLS estimates"), color = "black") +
  geom_point(aes(x = n, y = True_best, color = "True best"),fill = '#115F83FF', shape = 23, size = 3)+
  scale_fill_manual(name = "Boxplot", values = c('OLS estimates' = '#FFA319FF')) + 
  scale_color_manual(name = 'Points', values = c('True best'='#115F83FF')) +
  ggtitle("Best policy estimate (pre WC-adjutment)") +
  labs(colour = "Legend") + 
  xlab('Sample size n') +
  ylab('Estimate') +
  scale_x_continuous(breaks = linspace(0,10000,11)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))



#ggsave("1B_OLS_best_policy_estimates_any.pdf", plot = p_ols3, width = 20, height = 10)

#--------------------------#
# 1B: OLS SHRINKAGE
#--------------------------#

p_ols4 <- ggplot() +
  geom_boxplot(mutate(df_boxplot, n = n-1), mapping = aes(x = factor(n), y = percentage_shrunk*100, fill = method)) +
  scale_fill_manual(name = "Estimation type", values = c('#FFA319FF', '#115F83FF')) + 
  ggtitle("Shrinkage")+
  xlab('Sample size n') +
  ylab("Shrinkage (%)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 14),
        axis.text = element_text(size = 12), axis.title = element_text(size =  14))

ggsave("Simulation_Figures/Figure3/1B_Shrunk_vs_own_beta_prct_BOXPLOT_any.pdf", plot = p_ols4, width = 8, height = 6)

filter(df_boxplot, method == "Hybrid OLS")$percentage_shrunk %>% max()


#--------------------------------------#
# Final OLS combined plot
#--------------------------------------#

p_inclusion <- p_ols2
p_estimates <- p_ols3
p_shrinkage <- p_ols4
p_MSE_ols       <- p_ols1


p_ols_combined <- ggarrange(p_inclusion, p_estimates, p_shrinkage,p_MSE_ols,labels = c("A","B", "C","D"),nrow  = 2, ncol = 2)
#ggsave("1B_OLS_output.pdf", plot = p_ols_combined, width = 15, height = 10)

###Note: Go to "Puffer Asymptotic Normality" section to have the updated combined graph








#--------------------------------------#
# OLS coefplot visualization
#--------------------------------------#

ols_coefficients_long <- read.csv(paste0(path_output_data,"/Simulation_Data/1D_OLS_density.csv"))
puffer_coefficients_long <- read.csv(paste0(path_output_data,"/Simulation_Data/1D_Puffer_density.csv"))

#PLEASE CHECK THE FOLLOWING NUMBERS IN 1D_OLS_visualization
collapsed_pop_trueeffects <- c(0.8,1.2,0.3,1.5)
nsim = 300

cols = hue_pal()(5)
reds = colorRampPalette(c("red","#F8766D"))
yellows = colorRampPalette(c("#D19300","#A3A500"))
greens = colorRampPalette(c("green","#00BF7D"))
blues = colorRampPalette(c("cyan","#00B0F6"))
purples = colorRampPalette(c("lightpink","#E76BF3"))


p_ols_densities <- ggplot() + 
  geom_density(data = filter(ols_coefficients_long,Pooling == "Pruned"), aes(x = estimate, fill = Policies), alpha = 0.7)+
  scale_fill_manual(name = "Pruned", values = reds(51),guide = guide_legend(ncol = 6, order = 5)) + 
  new_scale_fill() +
  geom_density(data = filter(ols_coefficients_long,Pooling == "Pool1"), aes(x = estimate, fill = Policies), alpha = 0.7)+
  scale_fill_manual(name = "Pool 1", values = yellows(8),guide = guide_legend(ncol = 5, order= 1)) + 
  new_scale_fill() +
  geom_density(data = filter(ols_coefficients_long,Pooling == "Pool2"), aes(x = estimate, fill = Policies), alpha = 0.7)+
  scale_fill_manual(name = "Pool 2", values = greens(8),guide = guide_legend(ncol = 5, order= 2)) + 
  new_scale_fill() +
  geom_density(data = filter(ols_coefficients_long,Pooling == "Pool3"), aes(x = estimate, fill = Policies), alpha = 0.7)+
  scale_fill_manual(name = "Pool 3", values = blues(8),guide = guide_legend(ncol = 5, order= 3)) + 
  new_scale_fill() +
  geom_density(data = filter(ols_coefficients_long,Pooling == "Pool4"), aes(x = estimate, fill = Policies), alpha = 0.7)+
  scale_fill_manual(name = "Pool 4", values = purples(5),guide = guide_legend(ncol = 5, order= 4)) + 
  scale_x_continuous(limits = c(-1, 2.5), breaks = c(-1,0,collapsed_pop_trueeffects)) + 
  geom_vline(xintercept = c(collapsed_pop_trueeffects),linetype = "dashed", lwd = 1) +
  xlab("OLS estimates") + 
  ylab("Density") + 
  ggtitle(paste0("Distribution of OLS estimates by policy for ",nsim," simulations")) + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"),
        legend.title=element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none")


p_puffer_densities <- ggplot(puffer_coefficients_long) + 
  geom_density(aes(x = estimate, fill = Policies), alpha = 0.7)+
  geom_vline(xintercept = c(collapsed_pop_trueeffects),linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 2.5), breaks = c(-1,0,collapsed_pop_trueeffects)) + 
  scale_fill_manual(name = "Legend",
                    values = c("Pruned" = "#FB3B36", "Pool 1" = cols[2],"Pool 2" = cols[3],"Pool 3" = cols[4],"Pool 4" = cols[5]), 
                    guide = guide_legend(ncol= 1)) + 
  xlab("TVA estimates") + 
  ylab("Density") + 
  ggtitle(paste0("Distribution of TVA estimates by policy for ",nsim," simulations")) + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"),
      legend.title=element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15, face = "bold"),
      legend.position = "left")

p_densities <- ggarrange(p_ols_densities,p_puffer_densities,labels = c("A", "B"), nrow = 2, align = "v")

ggsave("Simulation_Figures/FigureE1/1D_ols_densities.pdf", plot=p_ols_densities, width = 8, height = 6)
ggsave("Simulation_Figures/FigureE1/1D_puffer_densities.pdf", plot=p_puffer_densities, width = 9, height = 6)
#ggsave("1D_densities_combined.pdf", plot=p_densities, width = 20, height = 15)


#################################################################################
#
# II. LASSO PLOTS
#
#################################################################################


#--------------------------------------#
# 2A: LASSO BETA MSE UNCONDITIONAL
#--------------------------------------#

output2A_uncond <- read.csv(paste0(path_output_data,"/Simulation_Data/2A_LASSO_beta_MSE_unconditional_any.csv"))
#output2A_uncond <- read.csv(paste0(path_data,"/Simulation_Data/Saving_copies_of_big_sims/2A_LASSO_beta_MSE_unconditional copie.csv"))

p_lasso1 <- ggplot(output2A_uncond) + 
  geom_point(aes(y = deviations_beta_hybrid_over_n , x = nobs_list, color = "No pooling, only pruning", fill = "No pooling, only pruning", shape = 'No pooling, only pruning'),size=3) +
  geom_point(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = "TVA (pooling and pruning)", fill = "TVA (pooling and pruning)", shape = 'TVA (pooling and pruning)'),size=3) +
  xlab('Sample size n') +
  ylab('MSE of WC-adjusted estimate') + theme_bw() + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('No pooling, only pruning' = '#FFA319FF', 'TVA (pooling and pruning)'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('No pooling, only pruning' = '#FFA319FF', 'TVA (pooling and pruning)'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('No pooling, only pruning' = 21, 'TVA (pooling and pruning)'=24)) +
  scale_y_continuous(limits = c(0.0, 2.8), breaks = linspace(0,2.5,6)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  ggtitle("MSE (unconditional)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("2A_LASSO_beta_MSE_uncond.pdf", plot=p_lasso1, width = 12, height = 8.5)


#-------------------------------#
# 2A: LASSO BETA MSE CONDITIONAL
#-------------------------------#
output2A_cond <- read.csv(paste0(path_output_data,"/Simulation_Data/2A_LASSO_beta_MSE_conditional_any.csv"))
#output2A_cond <- read.csv(paste0(path_data,"/Simulation_Data/Saving_copies_of_big_sims/2A_LASSO_beta_MSE_conditional copie.csv"))

p_lasso2 <- ggplot(output2A_cond) + 
  geom_point(aes(y = deviations_beta_hybrid_over_n , x = nobs_list, color = "No pooling, only pruning", fill = "No pooling, only pruning", shape = 'No pooling, only pruning'),size=3) +
  geom_point(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = "TVA (pooling and pruning)", fill = "TVA (pooling and pruning)", shape = 'TVA (pooling and pruning)'),size=3) +
  xlab('Sample size n') +
  ylab('MSE of WC-adjusted estimate') + theme_bw() +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('No pooling, only pruning' = '#FFA319FF', 'TVA (pooling and pruning)'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('No pooling, only pruning' = '#FFA319FF', 'TVA (pooling and pruning)'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('No pooling, only pruning' = 21, 'TVA (pooling and pruning)'=24)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  ggtitle("MSE (conditional)") + 
  scale_y_continuous(limits = c(0.0, 2.8), breaks = linspace(0,2.5,6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("2A_LASSO_beta_MSE_cond.pdf", plot=p_lasso2, width = 12, height = 8.5)




#-------------------------------#
# 2B: LASSO BETA SHRINKAGE
#-------------------------------#

#coef_shrinkage_df <- read.csv(paste0(path_data,"/Simulation_Data/2B_LASSO_beta_Shrinkage.csv"))
coef_shrinkage_df <- read.csv(paste0(path_output_data,"/Simulation_Data/2B_LASSO_beta_Shrinkage_any.csv"))

#Relative terms
p_lasso3 <- ggplot(data = coef_shrinkage_df) + 
  geom_point(aes(y = Lasso_prct_diff , x = n, color = "Hybrid LASSO estimates", fill = "Hybrid LASSO estimates", shape = "Hybrid LASSO estimates")) +
  geom_point(aes(y = SP_prct_diff, x = n, color = "Hybrid TVA estimates", fill = "Hybrid TVA estimates", shape = "Hybrid TVA estimates")) +
  xlab('Sample size n') +
  ylab('% Shrinkage due to WC adjustment') + theme_bw() + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Hybrid LASSO estimates' = '#FFA319FF', 'Hybrid TVA estimates'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('Hybrid LASSO estimates' = '#FFA319FF', 'Hybrid TVA estimates'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('Hybrid LASSO estimates' = 21, 'Hybrid TVA estimates'=24)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  #scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#ggsave("2B_LASSO_beta_shrinkage_relative.pdf", plot=p_lasso3, width = 12, height = 8.5)


#Absolute terms
p_lasso3b <- ggplot(data = coef_shrinkage_df) + 
  geom_point(aes(y = Lasso_abs_diff , x = n, color = "Hybrid LASSO estimate", fill = "Hybrid LASSO estimate", shape = "Hybrid LASSO estimate")) +
  geom_point(aes(y = SP_abs_diff, x = n, color = "Hybrid Puffer estimate", fill = "Hybrid Puffer estimate", shape = "Hybrid Puffer estimate")) +
  xlab('Sample size n') +
  ylab('Shrinkage due to WC adjustment') + theme_bw() + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Hybrid LASSO estimate' = '#FFA319FF', 'Hybrid Puffer estimate'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('Hybrid LASSO estimate' = '#FFA319FF', 'Hybrid Puffer estimate'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('Hybrid LASSO estimate' = 21, 'Hybrid Puffer estimate'=24)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#ggsave("2B_LASSO_beta_shrinkage_absolute.pdf", plot=p_lasso3b, width = 12, height = 8.5)


#Boxplot version
c1 <- coef_shrinkage_df[,c("Lasso_prct_diff", "n")] %>% rename(c("prct_diff" = "Lasso_prct_diff")) %>% mutate("method" = "Hybrid LASSO estimates", n = factor(n))
c2 <- coef_shrinkage_df[,c("SP_prct_diff", "n")] %>% rename(c("prct_diff" = "SP_prct_diff")) %>% mutate("method" = "Hybrid TVA estimates", n = factor(n))
coef_shrinkage_bp_df = rbind(c1,c2)

p3d <- ggplot() +
  geom_boxplot(coef_shrinkage_bp_df, mapping = aes(x = n, y = 100*prct_diff, fill = method)) +
  scale_fill_manual(name = "Legend", values = c('#FFA319FF', '#115F83FF')) + 
  ggtitle("Shrinkage imposed by WC adjustment") +
  xlab('Sample size n') +
  ylab('% Shrinkage due to WC adjustment') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "bold"))

#ggsave("2B_LASSO_beta_shrinkage_relative_boxplot.pdf", plot=p3d, width = 12, height = 8.5)



p_lasso3c <- ggplot(data = coef_shrinkage_df) + 
  geom_point(aes(y = Lasso_best , x = n+50, color = "LASSO pre-WC adjustment", fill = "LASSO pre-WC adjustment", shape = "LASSO pre-WC adjustment")) +
  geom_point(aes(y = Lasso_hybrid , x = n+50, color = 'LASSO post-WC adjustment', fill = 'LASSO post-WC adjustment', shape = 'LASSO post-WC adjustment')) +
  geom_point(aes(y = SP_best, x = n-50, color = 'Puffer pre-WC adjustment', fill = 'Puffer pre-WC adjustment', shape = 'Puffer pre-WC adjustment')) +
  geom_point(aes(y = SP_hybrid, x = n-50, color = 'Puffer post-WC adjustment', fill = 'Puffer post-WC adjustment', shape = 'Puffer post-WC adjustment')) +
  geom_hline(yintercept=5, linetype = "dashed") +
  xlab('Sample size n') +
  ylab('Estimated Best Coefficient') + theme_bw() + 
  scale_color_manual(name = 'Legend', guide = 'legend', values =c('LASSO pre-WC adjustment' = 'blue', 'LASSO post-WC adjustment'='cyan3',
                                                                  'Puffer pre-WC adjustment' = 'red','Puffer post-WC adjustment' = '#FFA319FF')) +
  scale_fill_manual(name = 'Legend', values =c('LASSO pre-WC adjustment' = 'blue', 'LASSO post-WC adjustment'='cyan3',
                                               'Puffer pre-WC adjustment' = 'red','Puffer post-WC adjustment' = '#FFA319FF')) +  
  scale_shape_manual(name = 'Legend', values =c('LASSO pre-WC adjustment' = 21, 'LASSO post-WC adjustment' = 21,
                                                'Puffer pre-WC adjustment' = 24, 'Puffer post-WC adjustment' = 24)) +
  scale_x_continuous(breaks = linspace(1000,10000,11)) +
  #scale_y_continuous(limits = c(2,8), breaks = linspace(0,8.0,9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#ggsave("2B_LASSO_beta_estimated_best_coef.pdf", plot=p_lasso3c, width = 12, height = 8.5)


#Boxplot version
r1 <- coef_shrinkage_df[,c("Lasso_best", "n")] %>% rename(c("estimate" = "Lasso_best")) %>% mutate("Legend" = "LASSO pre-WC adjustment", n = factor(n))
r2 <- coef_shrinkage_df[,c("Lasso_hybrid", "n")] %>% rename(c("estimate" = "Lasso_hybrid")) %>% mutate("Legend" = "LASSO post-WC adjustment", n = factor(n))
r3 <- coef_shrinkage_df[,c("SP_best", "n")] %>% rename(c("estimate" = "SP_best")) %>% mutate("Legend" = "Puffer pre-WC adjustment", n = factor(n))
r4 <- coef_shrinkage_df[,c("SP_hybrid", "n")] %>% rename(c("estimate" = "SP_hybrid")) %>% mutate("Legend" = "Puffer post-WC adjustment", n = factor(n))
coef_shrinkage_bp_df2 = rbind(r1,r2,r3,r4)

p_lasso3_box <- ggplot() +
  geom_boxplot(coef_shrinkage_bp_df2, mapping = aes(x = n, y = estimate, fill = Legend)) +
  scale_fill_manual(values = c("LASSO pre-WC adjustment" = "red", "LASSO post-WC adjustment" = '#FFA319FF', 
                               "Puffer pre-WC adjustment" = '#ADB17DFF', "Puffer post-WC adjustment" = '#155F83FF')) + 
  geom_hline(yintercept=5, linetype = "dashed") +
  ggtitle("Shrinkage imposed by WC adjustment") +
  xlab('Sample size n') +
  ylab('Estimated Best Coefficient') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("2B_LASSO_beta_estimated_best_coef.pdf", plot=p_lasso3_box, width = 12, height = 8.5)



#-------------------------------#
# Final LASSO on (3.1) combined 
#-------------------------------#
panelA <- p_lasso1
panelB <- p_lasso2
#panelC <- p_lasso3_box
panelC <- p3d

p_lasso_combined1 <- ggarrange(ggarrange(panelA, panelB, labels = c("A","B"), ncol = 2, common.legend = TRUE, legend = "right"),
                               panelC,labels = c("","C"), ncol = 1)
ggsave("Simulation_Figures/FigureE2/2_LASSO_output.pdf", plot = p_lasso_combined1, width = 12, height = 10)





#-------------------------------#
# 2C: LASSO ALPHA SUPPORT
#-------------------------------#

output2C <- read.csv(paste0(path_output_data,"/Simulation_Data/2C_LASSO_alpha_support_any.csv"))
  
p4 <- ggplot(output2C) + 
geom_point(aes(y = support_acc_n_Victor_vanilla , x = nobs_list, color = "Naive LASSO", fill = "Naive LASSO", shape = 'Naive LASSO'),size=2) +
geom_point(aes(y = support_acc_n_backward, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size=2) +
xlab('Sample size n') +
ylab('Avg Support Accuracy') + theme_bw() + 
ggtitle("Support Accuracy") + 
scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +
scale_fill_manual(name = 'Model Type', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +  
scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
scale_x_continuous(breaks = round(linspace(0,10000,11))) +
scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"))

#ggsave("2C_LASSO_alpha_Support_accuracy_any.pdf", plot=p4, width = 12, height = 10)


#-------------------------------#
# 2D: LASSO ALPHA MSE
#-------------------------------#

output2D <- read.csv(paste0(path_output_data,"/Simulation_Data/2D_LASSO_alpha_MSE_any.csv"))

##Plotting MSE of WC adjusted coefficients
p5 <- ggplot(output2D) + 
  geom_point(aes(y = deviations_Victor_vanilla_hybrid_over_n  , x = nobs_list, color = "Naive LASSO", fill = "Naive LASSO", shape = 'Naive LASSO')) +
  geom_point(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA')) +
  xlab('Sample size n') +
  ylab('MSE of WC-adjusted estimate') + theme_bw() + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  ggtitle("MSE") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "bold"))


#ggsave("2D_LASSO_alpha_MSE_best_policy_any.pdf", plot=p5, width = 12, height = 8.5)


#-----------------------------------------------#
# 2E: LASSO ALPHA BEST POLICY INCLUSION
#-----------------------------------------------#

#output2E <- read.csv(paste0(path_output_data,"/Simulation_Data/2E_LASSO_alpha_best_policy_coverage_any.csv"))

output2E <- read.csv(paste0(path_output_data,"/Simulation_Data/2E_LASSO_alpha_best_policy_coverage.csv"))

p_overall <- ggplot(output2E) +
  geom_point(aes(y = thickness_n_Victor_vanilla, x = nobs_list, color = "Naive LASSO",fill = "Naive LASSO", shape = 'Naive LASSO'), size = 2, alpha = 1) +
  geom_point(aes(y = thickness_n_backward, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2, alpha = 1) +
  geom_line(aes(y = thickness_n_Victor_vanilla, x = nobs_list, color = "Naive LASSO")) +
  geom_line(aes(y = thickness_n_backward, x = nobs_list, color = "TVA")) +
  xlab('Sample size n') +
  ylab('Overall Best Policy Inclusion') +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_fill_manual(name = 'Model Type',values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.5,1,6)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#ggsave("2E_LASSO_alpha_overall_any.pdf", plot=p_overall, width = 12, height = 8.5)

p_some_best <- ggplot(output2E) +
  #geom_jitter(aes(y = cover_some_best_n_backward, x = nobs, color = "Some best (OLS)", fill = "Some best (OLS)", shape = "Some best (OLS)"),size = 3,height=0,width=250) +
  geom_jitter(aes(y = cover_some_best_n_backward, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2,height=0,width=150) +
  geom_point(aes(y = cover_some_best_n_vanilla, x = nobs_list, color = "Naive LASSO",fill = "Naive LASSO", shape = 'Naive LASSO'), size = 2, alpha = 1) +
  #geom_line(aes(y = cover_some_best_n_backward, x = nobs_list, color = "Smart Pooling and Pruning")) +
  #geom_line(aes(y = cover_some_best_n_vanilla, x = nobs_list, color = "Naive LASSO")) +
  xlab('Sample size n') +
  ylab('Some Best Policy Inclusion') +
  ggtitle("Some Best Inclusion") +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_fill_manual(name = 'Model Type',values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  scale_x_continuous(breaks = linspace(0,10000,11)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.5,1,6)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("2E_LASSO_alpha_some_best_any.pdf", plot=p_some_best, width = 12, height = 8.5)

p_min_best <- ggplot(output2E) +
  geom_point(aes(y = cover_min_n_vanilla, x = nobs_list, color = "Naive LASSO",fill = "Naive LASSO", shape = 'Naive LASSO'), size = 2, alpha = 1) +
  geom_point(aes(y = cover_min_n_backward, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2, alpha = 1) +
  #geom_line(aes(y = cover_min_n_vanilla, x = nobs_list, color = "Naive LASSO")) +
  #geom_line(aes(y = cover_min_n_backward, x = nobs_list, color = "Smart Pooling and Pruning")) +
  xlab('Sample size n') +
  ylab('Min Best Policy Inclusion') +
  ggtitle("Minimum Best Inclusion") + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_fill_manual(name = 'Model Type',values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  scale_x_continuous(breaks = linspace(0,10000,11)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.5,1,6)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("2E_LASSO_alpha_minbest_any.pdf", plot=p_min_best, width = 12, height = 8.5)


#-------------------------------#
# Final LASSO on (3.2) combined 
#-------------------------------#

p_support <- p4
p_MSE <- p5

p_lasso_combined2 <- ggarrange(p_support, p_MSE, p_some_best,p_min_best,labels = c("A", "B","C", "D"), nrow = 2, ncol = 2)

ggsave("Simulation_Figures/FigureE3/2_LASSO_output2.pdf", plot=p_lasso_combined2,  width = 12, height = 8.5)



#-------------------------------#
# PUFFER ASYMPTOTIC NORMALITY
#-------------------------------#

## Pseudo-truth normality using the mixture distributions ## 
normalized_discrepancies_plotdf <- data.frame()

#for (n in c("1000", "3000")) 
#for (n in c("1000", "1500", "3000", "5000", "10000")) #for the paper
  
for (n in c("1000", "1500", "3000", "5000", "15000")){
  path <- paste0(paste0(path_output_data,"/Simulation_Data/2G_Puffer_normality_n=",n,".csv"))
  normalized_discrepancy_df <- read.csv(path)
  print(paste0("Number of unique policies throughout sims = ", dim(normalized_discrepancy_df)[2]-2))
  normalized_discrepancy_list <- normalized_discrepancy_df[,-c(1,dim(normalized_discrepancy_df)[2])] %>% unlist() %>% unname()
  normalized_discrepancy_list <- normalized_discrepancy_list[!is.na(normalized_discrepancy_list)]
  
  
  new_df <- data.frame("discrepancy" = normalized_discrepancy_list,
             "Legend" = paste0("n = ",n," ; r = ",normalized_discrepancy_df$support_detection_rate[1]),
             "n" = n)
  
  normalized_discrepancies_plotdf <- rbind(normalized_discrepancies_plotdf,new_df)
  
  assign(paste0("normalized_discrepancy_df_",n),normalized_discrepancy_df)
  assign(paste0("normalized_discrepancy_list_",n),normalized_discrepancy_list)
}

#Asymptotic Normality
normalized_discrepancies_plotdf = filter(normalized_discrepancies_plotdf, n!= 5000)

legend <- unique(normalized_discrepancies_plotdf$Legend)

l1 = toString(legend[1])
l2 = toString(legend[2])
l3 = toString(legend[3])
l4 = toString(legend[4])

values = c(l1 = "darkolivegreen2", l2 = "orange", l3 = "cyan4",  l4 = "navy", "Theoretical CDF" = "firebrick") %>%
  setNames(c(l1, l2, l3, l4, "Theoretical CDF"))

         
p_puffer_pseudotruth <- ggplot() +
  stat_ecdf(data = filter(normalized_discrepancies_plotdf,n != 5000),aes(discrepancy, color = Legend), size = 1.5)+
  geom_point(aes(x = linspace(-5,5,3000), y = pnorm(linspace(-5,5,3000),0,1), colour = "Theoretical CDF"), size = 0.5) +
  geom_hline(yintercept = c(1,0),linetype = "dashed") +
  scale_color_manual(name = "Legend", values = values) + 
  xlab("Normalized Discrepancy") +
  ylab("Empirical CDF")+
  ggtitle('TVA Normality (mixture distribution)') + 
  scale_x_continuous(limits = c(-5, 5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 14),
        axis.text = element_text(size = 12), axis.title = element_text(size =  14))


ggsave("Simulation_Figures/Figure3/2G_Puffer_pseudo_truth.pdf", plot=p_puffer_pseudotruth, width = 8, height = 6)


#################################################################################
#
# IV. DEBIASED LASSO PLOTS
#
#################################################################################

debiased_lasso_df <- read.csv(paste0(path_output_data,"/Simulation_Data/4A_Debiased_LASSO_any.csv"))
debiased_lasso_mse_df <- read.csv(paste0(path_output_data,"/Simulation_Data/4A_Debiased_LASSO_mse_any.csv"))

#Boxplot version
c1 <- debiased_lasso_df[,c("OLS_prct_diff", "n")] %>% rename(c("prct_diff" = "OLS_prct_diff")) %>% mutate("method" = "Direct OLS", n = factor(n))
c2 <- debiased_lasso_df[,c("SP_prct_diff", "n")] %>% rename(c("prct_diff" = "SP_prct_diff")) %>% mutate("method" = "TVA", n = factor(n))
c3 <- debiased_lasso_df[,c("DL_prct_diff", "n")] %>% rename(c("prct_diff" = "DL_prct_diff")) %>% mutate("method" = "Debiased LASSO", n = factor(n))
debiasd_lasso_df_bp = rbind(c1,c2,c3)

p_debiased1 <- ggplot() +
  geom_boxplot(debiasd_lasso_df_bp, mapping = aes(x = n, y = 100*prct_diff, fill = method)) +
  scale_fill_manual(name = "Legend", values = c("Direct OLS" = '#FFA319FF', "Debiased LASSO" = 'red',"TVA" = '#115F83FF')) + 
  ggtitle("Shrinkage imposed by WC adjustment") +
  xlab('Sample size n') +
  ylab('% Shrinkage due to WC adjustment') +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "bold"))

#ggsave("4A_DebiasedLASSO_shrinkage_any.pdf", plot=p_debiased1, width = 12, height = 8.5)

p_debiased2 <- ggplot(debiased_lasso_mse_df) + 
  geom_point(aes(y = deviations_ols_hybrid_over_n  , x = nobs_list, color = "Direct OLS", fill = "Direct OLS", shape = 'Direct OLS'), size = 3) +
  geom_point(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'), size = 3) +
  geom_point(aes(y = deviations_debiased_lasso_over_n, x = nobs_list, color = "Debiased LASSO", fill = "Debiased LASSO", shape = 'Debiased LASSO'), size = 3) +
  xlab('Sample size n') +
  ylab('MSE of WC-adjusted estimate') + theme_bw() + 
  ggtitle("MSE") +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Direct OLS' = '#FFA319FF', 'Debiased LASSO' = 'red', 'TVA'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('Direct OLS' = '#FFA319FF', 'Debiased LASSO' = 'red', 'TVA'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('Direct OLS' = 21, 'Debiased LASSO' = 22,'TVA'=24)) +
  scale_x_continuous(limits = c(0, 10500), breaks = linspace(0,10000,11)) +
  #scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "bold"))

#ggsave("4A_DebiasedLASSO_mse_any.pdf", plot=p_debiased2, width = 12, height = 8.5)

p_debiasedLASSO_combined <-  ggarrange(p_debiased1,p_debiased2, labels = c("A", "B"), ncol = 2)

ggsave("Simulation_Figures/FigureE4/4_debiasedLASSO_output.pdf", plot=p_debiasedLASSO_combined ,  width = 15, height = 6)





#################################################################################
#
# V. BAYESIAN BOOTSTRAPPING SPIKE AND SLAB LASSO
#
#################################################################################

#-----------------------------------------------#
# 5A: BBSSL - SUPPORT ACCURACY
#-----------------------------------------------#

bbssl_support <- read.csv(paste0(path_output_data,"/Simulation_Data/5A_SSLASSO_alpha_support_any.csv"))
bbssl_best_inclusion <- read.csv(paste0(path_output_data,"/Simulation_Data/5A_SSLASSO_alpha_inclusion_any.csv"))

p_bbssl3 <- ggplot(bbssl_support) +
  geom_point(aes(y = support_acc_n_backward, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size=2) +
  geom_point(aes(y = support_acc_n_bbssl, x = nobs_list, color = "BBSSL", fill = "BBSSL", shape = 'BBSSL'),size=2) +
  xlab('Sample size n') +
  ylab('Avg Support Accuracy') + theme_bw() +
  ggtitle("Support Accuracy") +
  scale_color_manual(name = 'Model Type', guide = 'legend',
                     values =c('TVA'='#155F83FF', "BBSSL"='#FFA319FF')) +
  scale_fill_manual(name = 'Model Type',
                    values =c('TVA'='#155F83FF', "BBSSL"='#FFA319FF')) +
  scale_shape_manual(name = 'Model Type',
                     values =c('TVA'=24, "BBSSL"=23)) +
  scale_x_continuous(breaks = round(linspace(0,10000,11))) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("5A_BBSSL_support_accuracy_any.pdf", plot=p_bbssl3, width = 12, height = 8.5)


p_bbssl_some_best <- ggplot(bbssl_best_inclusion) +
  geom_jitter(aes(y = cover_some_best_n_backward, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2,height=0,width=150) +
  geom_point(aes(y = cover_some_best_n_bbssl, x = nobs_list, color = "BBSSL",fill = "BBSSL", shape = 'BBSSL'), size = 2, alpha = 1) +
  xlab('Sample size n') +
  ylab('Some Best Policy Inclusion') +
  ggtitle("Some Best Inclusion") +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('TVA' = '#155F83FF','BBSSL' = '#FFA319FF')) +
  scale_fill_manual(name = 'Model Type',values =c('TVA' = '#155F83FF','BBSSL' = '#FFA319FF')) +
  scale_shape_manual(name = 'Model Type', values =c('TVA'=24,'BBSSL' = 21)) +
  scale_x_continuous(breaks = linspace(0,10000,11)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.5,1,6)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("5A_BBSSL_somebest_any.pdf", plot=p_bbssl_some_best, width = 12, height = 8.5)

p_bbssl_min_best <- ggplot(bbssl_best_inclusion) +
  geom_point(aes(y = cover_min_n_bbssl, x = nobs_list, color = "BBSSL",fill = "BBSSL", shape = 'BBSSL'), size = 2, alpha = 1) +
  geom_point(aes(y = cover_min_n_backward, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2, alpha = 1) +
  xlab('Sample size n') +
  ylab('Min Best Policy Inclusion') +
  ggtitle("Minimum Best Inclusion") + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('TVA' = '#155F83FF','BBSSL' = '#FFA319FF')) +
  scale_fill_manual(name = 'Model Type',values =c('TVA' = '#155F83FF','BBSSL' = '#FFA319FF')) +
  scale_shape_manual(name = 'Model Type', values =c('TVA'=24,'BBSSL' = 21)) +
  scale_x_continuous(breaks = linspace(0,10000,11)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.5,1,6)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#ggsave("5A_BBSSL_minbest_any.pdf", plot=p_bbssl_min_best, width = 12, height = 8.5)

p_bbssl_mse <- ggplot(bbssl_support) + 
  geom_point(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'), size = 3) +
  geom_point(aes(y = deviations_sslasso_over_n, x = nobs_list, color = "BBSSL", fill = "BBSSL", shape = 'BBSSL'), size = 3) +
  xlab('sample size n') +
  ylab('MSE of WC-adjusted estimate') + theme_bw() + 
  ggtitle("MSE") +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('TVA' = '#155F83FF','BBSSL' = '#FFA319FF')) +
  scale_fill_manual(name = 'Model Type',values =c('TVA' = '#155F83FF','BBSSL' = '#FFA319FF')) +
  scale_shape_manual(name = 'Model Type', values =c('TVA'=24, 'BBSSL'=23)) +
  scale_x_continuous(breaks = linspace(0,10000,11)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "bold"))

#ggsave("5B_BBSSL_mse_any.pdf", plot=p_bbssl_mse, width = 12, height = 8.5)


p_bbssl_output <- ggarrange(p_bbssl3, p_bbssl_mse, p_bbssl_some_best,p_bbssl_min_best,labels = c("A", "B","C", "D"), nrow = 2, ncol = 2)

ggsave("Simulation_Figures/FigureE5/5_BBSSL_output.pdf", plot=p_bbssl_output, width = 12, height = 8.5)


#################################################################################
#
# IV. REGIME PLOTS
#
#################################################################################

regime = "1"
output_regimes <- read.csv(paste0(paste0(path_output_data,"/Simulation_Data/3_Regimes_alpha_R",regime,"_any.csv")))
#output_regimes <- read.csv(paste0(paste0(path_data,"/Simulation_Data/Saving_copies_of_big_sims/3_Regimes_alpha_",regime," copie.csv"))) #old version

regime_title = switch(regime,
                      "1" = "Figure 17: Regime 1 - Few effective policies & imperfect sparsity",
                      "2" = "Figure 18: Regime 2 - Few effective policies & imperfect sparsity",
                      "3" = "Figure 19: Regime 3 - Few very effective policies, few medium policies & imperfect sparsity",
                      "4" = "Figure 20: Regime 4 - Few slightly effective policies & full sparsity",
                      "5" = "Figure 21: Regime 5 - Few slightly effective policies & partial sparsity")

p_support <- ggplot(output_regimes) + 
  geom_point(aes(y = support_acc_Victor_vanilla_over_n, x = nobs_list, color = "Naive LASSO", fill = "Naive LASSO", shape = 'Naive LASSO')) +
  geom_point(aes(y = support_acc_backward_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA')) +
  xlab('Sample size n') +
  ylab('Avg Support Accuracy') + theme_bw() + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  #scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = margin(1,0.5,1,0.5,"cm"), axis.text.x = element_text(angle = 20))



p_MSE <- ggplot(output_regimes) + 
  geom_point(aes(y = deviations_Victor_vanilla_hybrid_over_n  , x = nobs_list, color = "Naive LASSO", fill = "Naive LASSO", shape = 'Naive LASSO')) +
  geom_point(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA')) +
  xlab('Sample size n') +
  ylab('MSE of WC-adjusted estimate') + theme_bw() + 
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +
  scale_fill_manual(name = 'Model Type', values =c('Naive LASSO' = '#FFA319FF', 'TVA'='#155F83FF')) +  
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  #scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  #scale_y_continuous(limits = c(0.0, 5), breaks = linspace(0,5,11)) +
  scale_y_continuous(limits = c(0.0, 1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = margin(1,0.5,1,0.5,"cm"), axis.text.x = element_text(angle = 10))



p_overall_best <- ggplot(output_regimes) +
  geom_point(aes(y = thickness_Victor_vanilla_over_n, x = nobs_list, color = "Naive LASSO",fill = "Naive LASSO", shape = 'Naive LASSO'), size = 2, alpha = 1) +
  geom_point(aes(y = thickness_backward_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2, alpha = 1) +
  geom_line(aes(y = thickness_Victor_vanilla_over_n, x = nobs_list, color = "Naive LASSO")) +
  geom_line(aes(y = thickness_backward_over_n, x = nobs_list, color = "TVA")) +
  xlab('Sample size n') +
  ylab('Overall Best Policy Inclusion') +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_fill_manual(name = 'Model Type',values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  #scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  #scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.7,1,3)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = margin(1,0.5,1,0.5,"cm"), axis.text.x = element_text(angle = 20))



p_some_best <- ggplot(output_regimes) +
  geom_point(aes(y = cover_some_best_vanilla_over_n, x = nobs_list, color = "Naive LASSO",fill = "Naive LASSO", shape = 'Naive LASSO'), size = 2, alpha = 1) +
  geom_point(aes(y = cover_some_best_backward_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2, alpha = 1) +
  geom_line(aes(y = cover_some_best_vanilla_over_n, x = nobs_list, color = "Naive LASSO")) +
  geom_line(aes(y = cover_some_best_backward_over_n, x = nobs_list, color = "TVA")) +
  xlab('Sample size n') +
  ylab('Some Best Policy Inclusion') +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_fill_manual(name = 'Model Type',values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  #scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.7,1,3)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = margin(1,0.5,1,0.5,"cm"), axis.text.x = element_text(angle = 20))


p_min_best <- ggplot(output_regimes) +
  geom_point(aes(y = cover_min_vanilla_over_n, x = nobs_list, color = "Naive LASSO",fill = "Naive LASSO", shape = 'Naive LASSO'), size = 2, alpha = 1) +
  geom_point(aes(y = cover_min_backward_over_n, x = nobs_list, color = "TVA", fill = "TVA", shape = 'TVA'),size = 2, alpha = 1) +
  geom_line(aes(y = cover_min_vanilla_over_n, x = nobs_list, color = "Naive LASSO")) +
  geom_line(aes(y = cover_min_backward_over_n, x = nobs_list, color = "TVA")) +
  xlab('Sample size n') +
  ylab('Min Best Policy Inclusion') +
  scale_color_manual(name = 'Model Type', guide = 'legend', values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_fill_manual(name = 'Model Type',values =c('Naive LASSO' = '#FFA319FF', 'TVA' = '#155F83FF')) +
  scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'TVA'=24)) +
  #scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  #scale_y_continuous(limits = c(0, 1), breaks = linspace(0,1,11)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = margin(1,0.5,1,0.5,"cm"), axis.text.x = element_text(angle = 20))


p_support_sets_b <- ggplot(output_regimes) +
  geom_point(aes(y = backward_superset_over_n, x = nobs_list, color = "Superset",fill = "Superset")) +
  geom_point(aes(y = backward_subset_over_n, x = nobs_list, color = "Subset", fill = "Subset")) +
  geom_point(aes(y = backward_equal_over_n, x = nobs_list, color = "Correct Support", fill = "Correct Support")) +
  geom_line(aes(y = backward_superset_over_n, x = nobs_list, color = "Superset")) +
  geom_line(aes(y = backward_subset_over_n, x = nobs_list, color = "Subset")) +
  geom_line(aes(y = backward_equal_over_n, x = nobs_list, color = "Correct Support")) +
  xlab('Sample size n') +
  ylab('Share of simulations') +
  scale_color_manual(name = 'Legend', values =c('Superset' = '#FFA319FF', 'Subset' = '#155F83FF',"Correct Support"='firebrick')) +
  scale_fill_manual(name = 'Legend', values =c('Superset' = '#FFA319FF', 'Subset' = '#155F83FF',"Correct Support"='firebrick')) +
  #scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'Smart Pooling and Pruning'=24)) +
  #scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  #scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.7,1,3)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = margin(1,0.5,1,0.5,"cm"), axis.text.x = element_text(angle = 20))


p_support_sets_l <- ggplot(output_regimes) +
  geom_point(aes(y = lasso_superset_over_n, x = nobs_list, color = "Superset",fill = "Superset")) +
  geom_point(aes(y = lasso_subset_over_n, x = nobs_list, color = "Subset", fill = "Subset")) +
  geom_point(aes(y = lasso_equal_over_n, x = nobs_list, color = "Correct Support", fill = "Correct Support")) +
  geom_line(aes(y = lasso_superset_over_n, x = nobs_list, color = "Superset")) +
  geom_line(aes(y = lasso_subset_over_n, x = nobs_list, color = "Subset")) +
  geom_line(aes(y = lasso_equal_over_n, x = nobs_list, color = "Correct Support")) +
  xlab('Sample size n') +
  ylab('Share of simulations') +
  scale_color_manual(name = 'Legend', values =c('Superset' = '#FFA319FF', 'Subset' = '#155F83FF',"Correct Support"='firebrick')) +
  scale_fill_manual(name = 'Legend', values =c('Superset' = '#FFA319FF', 'Subset' = '#155F83FF',"Correct Support"='firebrick')) +
  #scale_shape_manual(name = 'Model Type', values =c('Naive LASSO' = 21, 'Smart Pooling and Pruning'=24)) +
  #scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  #scale_y_continuous(limits = c(0.5, 1), breaks = linspace(0.7,1,3)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin = margin(1,0.5,1,0.5,"cm"), axis.text.x = element_text(angle = 20))


p_subset_combined <- ggarrange(p_support_sets_b, p_support_sets_l, labels = c("TVA","LASSO"), nrow = 2)

#ggsave(paste0("3_Regime_plot_",regime,"_Subsets.pdf"), plot=p_subset_combined, width = 10, height = 8.5)

p_combined <- ggarrange(p_support, p_MSE, p_some_best, p_min_best,labels = c("Support", "MSE","Some Best", "Min Best"),
                        nrow  = 2, ncol = 2, legend = "right",common.legend = TRUE, font.label = list(size = 12))

# p_combined <- annotate_figure(p_combined, top = text_grob(regime_title,
#                                             face = "bold", size = 14, family = "Times"))

ggsave(paste0("Simulation_Figures/FigureE",as.integer(regime)+6,"/3_Regime_plot_",regime,".pdf"), plot=p_combined, width = 10, height = 8.5)

#ggsave(paste0("3_Regime_plot_",regime,"_overall_best.pdf"), plot=p_overall_best, width = 10, height = 8.5)







#################################################################################
#
# V. REGIME PLOTS COMBINED
#
#################################################################################

output_R1 <- read.csv(paste0(path_output_data,"/Simulation_Data/3_Regimes_alpha_R1_any.csv")) %>% mutate(Regime = "R1")
output_R2 <- read.csv(paste0(path_output_data,"/Simulation_Data/3_Regimes_alpha_R2_any.csv")) %>% mutate(Regime = "R2")
output_R3 <- read.csv(paste0(path_output_data,"/Simulation_Data/3_Regimes_alpha_R3_any.csv")) %>% mutate(Regime = "R3")
output_R4 <- read.csv(paste0(path_output_data,"/Simulation_Data/3_Regimes_alpha_R4_any.csv")) %>% mutate(Regime = "R4")
output_R5 <- read.csv(paste0(path_output_data,"/Simulation_Data/3_Regimes_alpha_R5_any.csv")) %>% mutate(Regime = "R5")


combined_output_regimes <- rbind(output_R1,output_R2,output_R3,output_R4,output_R5)



## -- SUPPORT GRAPH -- ##
p_support_TVA <- ggplot(combined_output_regimes) + 
  geom_point(aes(y = support_acc_backward_over_n, x = nobs_list, color = Regime, fill = Regime), size = 4) +
  geom_line(aes(y = support_acc_backward_over_n, x = nobs_list, color = Regime), lwd = 1.5) + 
  ggtitle("TVA Support Accuracies") + 
  xlab('Sample size n') +
  ylab('Avg Support Accuracy') + theme_bw() + 
  scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  scale_colour_manual(values = pal_uchicago("default")(5)) +
  scale_fill_manual(values = pal_uchicago("default")(5)) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(1,0.5,1,0.5,"cm"), 
        axis.text.x = element_text(angle = 20), plot.title = element_text(face = "bold", size=  18),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15,face = "bold"), legend.text = element_text(size = 15))

ggsave("Simulation_Figures/FigureE6/3_Regime_supportTVA.pdf", plot=p_support_TVA, width = 8, height = 6)

## -- MSE GRAPH --#

p_MSE_TVA <- ggplot(combined_output_regimes) + 
  geom_point(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = Regime, fill = Regime), size = 4) +
  geom_line(aes(y = deviations_SP_hybrid_over_n, x = nobs_list, color = Regime), lwd = 1.5) + 
  ggtitle("TVA MSE") + 
  xlab('Sample size n') +
  ylab('MSE of WC-adjusted estimate') + theme_bw() + 
  scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  scale_colour_manual(values = pal_uchicago("default")(5)) +
  scale_fill_manual(values = pal_uchicago("default")(5)) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(1,0.5,1,0.5,"cm"), 
        axis.text.x = element_text(angle = 20), plot.title = element_text(face = "bold", size=  18),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15,face = "bold"), legend.text = element_text(size = 15))

ggsave("Simulation_Figures/FigureE6/3_Regime_MSETVA.pdf", plot=p_MSE_TVA, width = 8, height = 6)

## -- SOME BEST GRAPH --#

p_somebest_TVA <- ggplot(combined_output_regimes) + 
  geom_point(aes(y = cover_some_best_backward_over_n, x = nobs_list, color = Regime, fill = Regime), size = 4) +
  geom_line(aes(y = cover_some_best_backward_over_n, x = nobs_list, color = Regime), lwd = 1.5) + 
  ggtitle("TVA Some Best") + 
  xlab('Sample size n') +
  ylab('Some Best Policy Inclusion') + theme_bw() + 
  scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  scale_colour_manual(values = pal_uchicago("default")(5)) +
  scale_fill_manual(values = pal_uchicago("default")(5)) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(1,0.5,1,0.5,"cm"), 
        axis.text.x = element_text(angle = 20), plot.title = element_text(face = "bold", size=  18),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15,face = "bold"), legend.text = element_text(size = 15))


ggsave("Simulation_Figures/FigureE6/3_Regime_somebestTVA.pdf", plot=p_somebest_TVA, width = 8, height = 6)


## -- MIN BEST GRAPH --#


p_minbest_TVA <- ggplot(combined_output_regimes) + 
  geom_point(aes(y = cover_min_backward_over_n, x = nobs_list, color = Regime, fill = Regime), size = 4) +
  geom_line(aes(y = cover_min_backward_over_n, x = nobs_list, color = Regime), lwd = 1.5) + 
  ggtitle("TVA Min Best") + 
  xlab('Sample size n') +
  ylab('Min Best Policy Inclusion') + theme_bw() + 
  scale_x_continuous(limits = c(1000, 10500), breaks = linspace(1000,10000,10)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = linspace(0,1.0,11)) +
  scale_colour_manual(values = pal_uchicago("default")(5)) +
  scale_fill_manual(values = pal_uchicago("default")(5)) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(1,0.5,1,0.5,"cm"), 
        axis.text.x = element_text(angle = 20), plot.title = element_text(face = "bold", size=  18),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15,face = "bold"), legend.text = element_text(size = 15))

ggsave("Simulation_Figures/FigureE6/3_Regime_minbestTVA.pdf", plot=p_minbest_TVA, width = 8, height = 6)






#################################################################################
#
# END
#
#################################################################################

