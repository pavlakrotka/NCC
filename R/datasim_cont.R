#' Data simulation for continuous endpoints for a platform trial with an arbitrary number of treatment arms entering sequentially
#'
#' @description Simulates data from a platform trial with continuous endpoints and an arbitrary number of treatment arms entering sequentially. There is the option to either specify the timing of adding arms and the overall sample size in the trial or specify the sample size per arm (assumed equal) and the allocation ratios per period.
#'
#' @param n_total Overall sample size for the trial
#' @param num_arms Number of treatment arms in the trial
#' @param t_arm Timing of adding new arms in terms of number of patients allocated to the control arm
#' @param n_arm Sample size per arm
#' @param alloc_ratios Matrix with allocation ratios in each period
#' @param period_blocks Number to multiply the number of active arms with in order to get the block size per period
#' @param mu0 Response in the control arm. Default=0
#' @param delta Vector with treatment effects for each arm
#' @param lambda Vector with strength of time trend in each arm
#' @param sigma Residual variance
#' @param trend Indicates the time trend pattern ("linear", "stepwise" or "inv_u")
#' @param N_peak Point at which the inverted-u time trend switches direction in terms of overall sample size
#' 
#'
#' @import rlang
#' @import stats
#' 
#' @export
#' @return Data frame: simulated trial data
#' @author Pavla Krotka, Marta Bofill Roig

datasim_cont <- function(n_total, num_arms, t_arm, n_arm, alloc_ratios, period_blocks=2, mu0=0, delta, lambda, sigma, trend, N_peak){
  
  requireNamespace("rlang")
  requireNamespace("stats")
  
  require("rlang")
  require("stats")
  
  if (is.null(n_total)==F & is.null(num_arms)==F & is.null(t_arm)==F){
    
    SS_matrix <- get_ss_matrix(n_total, num_arms, t_arm)
    alloc_ratios <- ifelse(!is.na(SS_matrix), 1, 0)
    
    num_periods <- ncol(alloc_ratios) # total number of periods
    num_arms <- nrow(alloc_ratios)-1 # total number of treatment arms
    
  } else if (is.null(n_arm)==F & is.null(alloc_ratios)==F){
    
    num_periods <- ncol(alloc_ratios) # total number of periods
    num_arms <- nrow(alloc_ratios)-1 # total number of treatment arms
    
    SS_matrix <- matrix(nrow = num_arms+1, ncol = num_periods)
    
    for (i in 1:num_arms) { # get sample sizes for each arm
      SS_matrix[i+1,] <- alloc_ratios[i+1,]/sum(alloc_ratios[i+1,], na.rm = T)*n_arm
    }
    
    SS_matrix[1,] <- na.omit(apply(SS_matrix, 2, unique)) # get sample sizes for control
    
    alloc_ratios[is.na(alloc_ratios)] <- 0
    
  } else {
    stop("Either n_total, num_arms and t_arm or n_arm and alloc_ratios must be specified!")
  }
  
  SS_matrix <- round(SS_matrix)
  
  N_period <- colSums(SS_matrix, na.rm=T) # sample sizes per period
  N_arm <- rowSums(SS_matrix, na.rm=T) # sample sizes per arm
  n_total <- sum(SS_matrix, na.rm = T) # total sample size
  active_arms <- colSums(apply(SS_matrix, 2, is.na)==0) # active arms per period
  block_sizes <- period_blocks*active_arms # block sizes per period
  
  
  t <- c()
  
  for (i in 1:num_periods){
    m_i <- t(replicate(trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i]),
                       sample(rep(rep(c(0:(num_arms)), alloc_ratios[,i]), block_sizes[i]/length(rep(c(0:(num_arms)), alloc_ratios[,i]))))))
    
    t_i <- c(t(m_i), sample(rep(c(0:(num_arms)), alloc_ratios[,i]),
                            size = sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])))
    
    t <- c(t, t_i)
  }
  
  
  for (i in 0:num_arms) {
    assign(paste0("j", i), which(t==i)) # j0, j1, j2 ... position in time (order) of allocated patients in every arm
  }
  
  cj <- rep(1:num_periods, N_period) # period indicator
  
  
  # Simulation of individual trend
  
  if(trend=="linear"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                  lambda = lambda[i+1],
                                                  sample_size = c(0, n_total)))
    }
  }
  
  if(trend=="linear2"){ # trend starts in the second period and is linear
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                  lambda = lambda[i+1],
                                                  sample_size = c(N_period[1], sum(N_period[-1]))))
    }
  }
  
  if(trend=="stepwise"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), sw_trend(cj=cj[eval(sym(paste0("j", i)))],
                                              lambda = lambda[i+1]))
    }
  }
  
  if(trend=="inv_u"){
    
    for (i in 0:num_arms) {
      assign(paste0("j", i, "_1"), which(eval(sym(paste0("j", i))) <= N_peak))
      assign(paste0("j", i, "_2"), which(eval(sym(paste0("j", i))) > N_peak))
      
      assign(paste0("ind_trend", i, "_1"), linear_trend(j = eval(sym(paste0("j", i)))[eval(sym(paste0("j", i, "_1")))],
                                                        lambda = lambda[i+1],
                                                        sample_size = n_total))
      
      assign(paste0("ind_trend", i, "_2"), linear_trend(j = eval(sym(paste0("j", i)))[eval(sym(paste0("j", i, "_2")))]-2*N_peak+1,
                                                        lambda = -lambda[i+1],
                                                        sample_size = n_total))
      
      assign(paste0("ind_trend", i), c(eval(sym(paste0("ind_trend", i, "_1"))), eval(sym(paste0("ind_trend", i, "_2")))))
    }
  }
  
  # Simulation of continuous endpoint
  
  means <- c()
  means[j0] <- ind_trend0
  
  for (i in 1:num_arms) {
    means[eval(sym(paste0("j", i)))] <- eval(sym(paste0("ind_trend", i))) + delta[i]
  }
  
  X <- rnorm(n=n_total, mean=mu0+means, sd=sigma)
  
  Data <- data.frame(response = X,
                     treatment = t,
                     period = rep(1:num_periods, N_period),
                     j = c(1:n_total),
                     means = mu0+means)
  
  for (i in 0:num_arms) {
    Data[ ,paste0("lambda", i)] <- lambda[i+1]
  }
  
  
  return(Data)
}


# test <- datasim_cont(n_total = 1000, num_arms = 3, t_arm = 120, period_blocks = 2, delta=rep(0.25, 3), lambda=rep(0.15, 4), sigma=1, trend="linear")

# test <- datasim_cont(n_total = NULL, num_arms = NULL, t_arm = NULL, 
#                       n_arm = 200, alloc_ratios = matrix(c(1,1,1,
#                                                            1,1,NA,
#                                                            NA,1,1), ncol = 3, byrow = T), period_blocks = 2, delta=rep(0.25, 3), lambda=rep(0.15, 4), sigma=1, trend="linear")





