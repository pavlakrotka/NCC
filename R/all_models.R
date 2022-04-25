#' Wrapper function for simulations analyzing given data with all models
#'
#' @description Analyzes given data using the fixed effect model, pooled and separate analyses, and timemachine and MAP prior approach. Performs inference for all treatment arms in the trial except for the first one
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` or `datasim_cont()` function
#' @param alpha Type I error. Default=0.025
#'
#' @keywords internal
#'
#' @export
#' 
#' @examples 
#' 
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#' 
#' all_models(data = trial_data)
#'
#'
#' @return List containing an indicator whether the null hypothesis was rejected or not for the investigated treatment for all models
#' @author Pavla Krotka


all_models <- function(data, alpha=0.025, ...){
  
  res <- list()
  
  if(all(data$response %in% c(0,1))){ # binary endpoints
    
    for (i in 2:max(data$treatment)) {
      res_arm <- list(fixmodel_bin(data, arm = i, alpha)$reject_h0,
                      sepmodel_bin(data, arm = i, alpha)$reject_h0,
                      poolmodel_bin(data, arm = i, alpha)$reject_h0,
                      timemachine_bin(data, arm = i, alpha, ...)$reject_h0,
                      MAPprior_bin(data, arm = i, alpha, ...)$reject_h0)
      
      names(res_arm) <- c(paste0("reject_h0_fix_", i),
                          paste0("reject_h0_sep_", i),
                          paste0("reject_h0_pool_", i),
                          paste0("reject_h0_timemachine_", i),
                          paste0("reject_h0_MAPprior_", i))
      
      res <- append(res, res_arm)
    }
  } else { # continuous endpoints
    
    for (i in 2:max(data$treatment)) {
      res_arm <- list(fixmodel_cont(data, arm = i, alpha)$reject_h0,
                      sepmodel_cont(data, arm = i, alpha)$reject_h0,
                      poolmodel_cont(data, arm = i, alpha)$reject_h0,
                      timemachine_cont(data, arm = i, alpha, ...)$reject_h0)
      
      names(res_arm) <- c(paste0("reject_h0_fix_", i),
                          paste0("reject_h0_sep_", i),
                          paste0("reject_h0_pool_", i),
                          paste0("reject_h0_timemachine_", i))
      
      res <- append(res, res_arm)
    }
  }
  return(res)
}


