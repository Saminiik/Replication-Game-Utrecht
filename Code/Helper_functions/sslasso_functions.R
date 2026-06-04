######################################################################################

# PROJECT: 	Smart Pooling and Pruning
# PURPOSE:  Functions for simulations (bayesian spike and slab lasso)
# AUTHOR:		Anirudh Sankar

#######################################################################################



library("BBSSL")

#------------------------------------------------------------------------------#
retrieve_SSLASSO_support_Victor <- function(y, sp_matrix, policy_vectors, lambda0, lambda1, stdev, get_ind = FALSE) {
    data_matrix <- cbind(y, sp_matrix)
    
    data_matrix <- as.data.frame(data_matrix)
    
    
    colnames(data_matrix)[1] <- "y"
    
    formule <- as.formula("y~ .")
    model_lasso <- rlasso(formula = formule, data = data_matrix,
                          penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE)) #you might adjust the penalty argument here
    initial.beta <- model_lasso$beta
    
    model_sslasso <- SSLASSO_2(sp_matrix, y, initial.beta, penalty = "adaptive", variance = "fixed", lambda0 = lambda0, lambda1 = lambda1, max.iter =1000, warn = TRUE)
    
    support_ind <- which(model_sslasso$select == 1)
    support <- policy_vectors[support_ind] #this retrieves the actual support in a list of vectors (e.g. 2 2, 3 4 etc)
    
    if (get_ind == TRUE) {
        return(support_ind)
    }
    
    return(support)
}

#------------------------------------------------------------------------------#

get_SSLASSO <- function(y, sp_matrix, policy_vectors,lambda0,lambda1, stdev){
  data_matrix <- cbind(y, sp_matrix)
  
  data_matrix <- as.data.frame(data_matrix)
  
  
  colnames(data_matrix)[1] <- "y"
  
  formule <- as.formula("y~ .")
  model_lasso <- rlasso(formula = formule, data = data_matrix,
                        penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE)) #you might adjust the penalty argument here
  initial.beta <- model_lasso$beta
  
  model_sslasso <- SSLASSO_2(sp_matrix, y, initial.beta, penalty = "adaptive",
                             variance = "fixed", lambda0 = lambda0, lambda1 = lambda1,
                             sigma = stdev, max.iter =1000)
  res <- list()
  # Support:
  support_ind <- which(model_sslasso$select == 1)
  support <- policy_vectors[support_ind] #this retrieves the actual support in a list of vectors (e.g. 2 2, 3 4 etc)
  res$support <- support
  # Coefficients of the support:
  coefficients_sslasso <- model_sslasso$beta[support_ind]
  res$coef <- coefficients_sslasso
  
  return(res)
}

######################################################################################
#
# END
#
#######################################################################################

