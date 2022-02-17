#' Data simulation for continuous endpoints for a platform trial with 3 treatment arms entering sequentially
#'
#' @description Simulates data from a platform trial with 3 treatment arms entering sequentially. There is the option to either specify the timing of adding arms and the overall sample size in the trial or specify the sample size per arm (assumed equal) and the allocation ratios per period.
#'
#' @param T_ Timing of adding new arms in terms of number of patients allocated to the control arm
#' @param N Overall sample size for the trial
#' @param SS_arm Sample size per arm
#' @param alloc_ratios Matrix with allocation ratios in each period
#' @param block_sizes Vector with block size in each period
#' @param mu0 Response in the control arm. Default=0
#' @param delta Vector with treatment effects for each arm
#' @param lambda Vector with strength of time trend in each arm
#' @param sigma Residual variance
#' @param N_peak Point at which the inverted-u time trend switches direction in terms of overall sample size
#' @param trend Indicates the time trend pattern ("linear", "stepwise" or "inv_u")
#'
#' @export
#' @return Data frame: simulated trial data
#' @author Pavla Krotka, Marta Bofill Roig

data_sim_cont <- function(T_, N, SS_arm, alloc_ratios, block_sizes, mu0=0, delta, lambda, sigma, N_peak, trend){
  
  require(rlang)
  
  if (is.null(T_)==F & is.null(N)==F){
    
    SS_matrix <- get_ss_matrix(T_, N)
    alloc_ratios <- ifelse(!is.na(SS_matrix), 1, 0)
    
    n_periods <- ncol(alloc_ratios) # total number of periods
    n_arms <- nrow(alloc_ratios) # total number of arms
    
  } else if (is.null(SS_arm)==F & is.null(alloc_ratios)==F){
    
    n_periods <- ncol(alloc_ratios) # total number of periods
    n_arms <- nrow(alloc_ratios) # total number of arms
    
    SS_matrix <- matrix(nrow = n_arms, ncol = n_periods)
    
    for (i in 2:n_arms) { # get sample sizes for each arm
      SS_matrix[i,] <-  alloc_ratios[i,]/sum(alloc_ratios[i,], na.rm = T)*SS_arm
    }
    
    SS_matrix[1,] <- na.omit(apply(SS_matrix, 2, unique)) # get sample sizes for contro
    
    alloc_ratios[is.na(alloc_ratios)] <- 0
    
  } else {
    stop("Either T_ and N or SS_arm and alloc_ratios must be specified!")
  }
  
  
  SS_matrix <- round(SS_matrix)
  
  SS_period <- colSums(SS_matrix, na.rm=T) # sample sizes per period
  SS_arm <- rowSums(SS_matrix, na.rm=T) # sample sizes per arm
  SS_total <- sum(SS_matrix, na.rm = T) # total sample size
  active_arms <- colSums(apply(SS_matrix, 2, is.na)==0) # active arms per period
  
  
  t <- c()
  
  for (i in 1:n_periods){
    m_i <- t(replicate(trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i]),
                       sample(rep(rep(c(0:(n_arms-1)), alloc_ratios[,i]), block_sizes[i]/length(rep(c(0:(n_arms-1)), alloc_ratios[,i]))))))
    
    t_i <- c(t(m_i), sample(rep(c(0:(n_arms-1)), alloc_ratios[,i]),
                            size = sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])))
    
    t <- c(t, t_i)
  }
  
  
  for (i in 0:(n_arms-1)) {
    assign(paste0("j", i), which(t==i)) # j0, j1, j2 ... position in time (order) of allocated patients in every arm
  }
  
  cj <- rep(1:n_periods, SS_period) # period indicator
  
  
  # Simulation of individual trend
  
  if(trend=="linear"){
    for (i in 0:(n_arms-1)) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                  lambda = lambda[i+1],
                                                  sample_size = c(0, SS_total)))
    }
  }
  
  if(trend=="linear2"){ # trend starts in the second period and is linear
    for (i in 0:(n_arms-1)) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                   lambda = lambda[i+1],
                                                   sample_size = c(SS_period[1], sum(SS_period[-1]))))
    }
  }
  
  if(trend=="stepwise"){
    for (i in 0:(n_arms-1)) {
      assign(paste0("ind_trend", i), sw_trend(cj=cj[eval(sym(paste0("j", i)))],
                                              lambda = lambda[i+1]))
    }
  }
  
  if(trend=="inv_u"){
    
    for (i in 0:(n_arms-1)) {
      assign(paste0("j", i, "_1"), which(eval(sym(paste0("j", i))) <= N_peak))
      assign(paste0("j", i, "_2"), which(eval(sym(paste0("j", i))) > N_peak))
      
      assign(paste0("ind_trend", i, "_1"), linear_trend(j = eval(sym(paste0("j", i)))[eval(sym(paste0("j", i, "_1")))],
                                                        lambda = lambda[i+1],
                                                        sample_size = SS_total))
      
      assign(paste0("ind_trend", i, "_2"), linear_trend(j = eval(sym(paste0("j", i)))[eval(sym(paste0("j", i, "_2")))]-2*N_peak+1,
                                                        lambda = -lambda[i+1],
                                                        sample_size = SS_total))
      
      assign(paste0("ind_trend", i), c(eval(sym(paste0("ind_trend", i, "_1"))), eval(sym(paste0("ind_trend", i, "_2")))))
    }
  }
  
  # Simulation of continuous endpoint
  
  means <- c()
  means[j0] <- ind_trend0
  
  for (i in 1:(n_arms-1)) {
    means[eval(sym(paste0("j", i)))] <- eval(sym(paste0("ind_trend", i))) + delta[i]
  }
  
  X <- rnorm(n=SS_total, mean=mu0+means, sd=sigma)
  
  Data <- data.frame(response = X,
                     treatment = t,
                     period = rep(1:n_periods, SS_period),
                     j = c(1:SS_total),
                     means = mu0+means)
  
  for (i in 0:(n_arms-1)) {
    Data[ ,paste0("lambda", i)] <- lambda[i+1]
  }
  
  
  return(Data)
}