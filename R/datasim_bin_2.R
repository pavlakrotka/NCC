#' Simulate binary data from a platform trial with a shared control arm and a given number of experimental treatment arms entering at given time points using a user-specified sample size matrix
#'
#' @description This function simulates data from a platform trial with a given number of experimental treatment arms entering at given time points and a shared control arm. The primary endpoint is a binary endpoint. The user specifies the timing of adding arms in terms of patients recruited to the trial so far and the sample size per experimental treatment arm.
#'
#' @param num_arms Integer. Number of experimental treatment arms in the trial.
#' @param n_arm Integer. Sample size per experimental treatment arm (assumed equal).
#' @param d Integer vector with timings of adding new arms in terms of number of patients recruited to the trial so far. The first entry must be 0, so that the trial starts with at least one experimental treatment arm, and the entries must be non-decreasing. The vector length equals `num_arms`.
#' @param period_blocks Integer. Number to define the size of the blocks for the block randomization. The block size in each period equals `period_blocks`times the number of active arms in the period (see Details). Default=2.
#' @param p0 Double. Response probability in the control arm.
#' @param OR Double vector with treatment effects in terms of odds ratios for each experimental treatment arm compared to control. The elements of the vector (odds ratios) are ordered by the entry time of the experimental treatment arms (e.g., the first entry in the vector corresponds to the odds ratio of the first experimental treatment arm). The vector length equals `num_arms`.
#' @param SS_matrix Matrix with sample sizes per arm (rows) and period (columns).
#' @param lambda Vector containing numerical entries or the string "random", indicating the strength of the time trend in each arm ordered by the entry time of the arms (e.g., the first entry in the vector corresponds to the time trend in the control arm, second entry to the time trend in the first experimental treatment arm). The vector length equals `num_arms`+1, as time trend in the control is also allowed. In case of random time trend, its strenght is generated from a normal distribution.
#' @param trend String indicating the time trend pattern ("linear", "linear_2, "stepwise", "stepwise_2", "inv_u" or "seasonal"). See Details for more information.
#' @param N_peak Integer. Timepoint at which the inverted-u time trend switches direction in terms of overall sample size (i.e. after how many recruited participants the trend direction switches).
#' @param n_wave Integer. Number of cycles (waves) should the seasonal trend have.
#' @param trend_mean Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`). Default: N(0, 0.5).
#' @param trend_var Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`). Default: N(0, 0.5).
#' @param full Logical. Indicates whether the output should be in form of a data frame with variables needed for the analysis only (FALSE) or in form of a list containing more information (TRUE). Default=FALSE.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#'
#' @importFrom stats na.omit
#' @importFrom stats rbinom
#' @importFrom rlang sym
#'
#' @export
#'
#' @details
#'
#' \strong{Design assumptions:}
#'
#' - The simulated platform trial consists of a given number of experimental treatment arms (specified by the argument `num_arms`) and one control arm that is shared across the whole platform.
#' - Participants are indexed by entry order, assuming that at each time unit exactly one participant is recruited and the time of recruitment and observation of the response are equal.
#' - All participants are assumed to be eligible for all arms in the trial, i.e. the same inclusion and exclusion criteria apply to all experimental and control arms.
#' - Equal sample sizes (given by parameter `n_arm`) in all experimental treatment arms are assumed.
#' - The duration of the trial is divided into so-called periods, defined as time intervals bounded by distinct time points of any treatment arm entering or leaving the platform.
#' Hence, multiple treatment arms entering or leaving at the same time point imply the start of only one additional period.
#' - Allocation ratio of 1:1:...:1 in each period. Furthermore, block randomization is used to assign patients to the active arms. Block size in each period = `period_blocks`* (number of active arms in the period).
#' - If the period sample size is not a multiple of the block size, arms for the remaining participants are chosen by sampling without replacement from a vector containing the indices of active arms replicated times `ceiling(remaining sample size/number of active arms)`.
#'
#' \strong{Data generation:}
#'
#' The binary response \eqn{y_j} for patient \eqn{j} is generated according to:
#'
#' \deqn{g(E(y_j)) = \eta_0 + \sum_{k=1}^K \cdot I(k_j=k) + f(j)}
#'
#' where \eqn{g(\cdot)} is the logit link function, and \eqn{\eta_0} (logit function of parameter `p0`) and \eqn{\theta_k} (log of the parameter `OR`) are the log odds in the control arm and the log odds ratio of treatment \eqn{k}. \eqn{K} is the total number of treatment arms in the trial (parameter `num_arms`) and \eqn{k_j} is an indicator of the treatment arm patient \eqn{j} is allocated to.
#'
#' The function \eqn{f(j)} denotes the time trend, whose strength is indicated by \eqn{\lambda_{k_j}} (parameter `lambda`) and which can have the following patterns (parameter `trend`):
#'
#' - \strong{"linear"} - trend starts at the beginning of the trial and the log odds increases or decreases linearly with a slope of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1}}, where \eqn{N} is the total sample size in the trial
#' - \strong{"linear_2"} - trend starts after the first period (i.e. there is no time trend in the first period) and the log odds increases or decreases linearly with a slope of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1}}, where \eqn{N} is the total sample size in the trial
#' - \strong{"stepwise"} - the log odds is constant in each period and increases or decreases by \eqn{\lambda} each time any treatment arm enters or leaves the trial (i.e. in each period), according to the function \eqn{f(j) = \lambda_{k_j} \cdot (c_j - 1)}, where \eqn{c_j} is an index of the period patient \eqn{j} was enrolled in
#' - \strong{"stepwise_2"} - the log odds is constant in each period and increases or decreases by \eqn{\lambda} each time a new treatment arm is added to the trial, according to the function \eqn{f(j) = \lambda_{k_j} \cdot (w_j - 1)}, where \eqn{w_j} is an indicator of how many treatment arms have already entered the ongoing trial, when patient \eqn{j} was enrolled
#' - \strong{"inv_u"} - the log odds increases up to the point \eqn{N_p} (parameter `N_peak`) and decreases afterwards, linearly with a slope of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1} (I(j \leq N_p) - I(j > N_p))}, where \eqn{N_p} indicates the point at which the trend turns from positive to negative in terms of the sample size (note that for negative \eqn{\lambda}, the log odds ratio decreases first and increases after)
#' - \strong{"seasonal"} - the log odds increases and decreases periodically with a magnitude of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \mathrm{sin} \big( \psi \cdot 2\pi \cdot \frac{j-1}{N-1} \big)}, where \eqn{\psi} indicates how many cycles should the time trend have (parameter `n_wave`)
#'
#' Trials with no time trend can be simulated too, by setting all elements of the vector `lambda` to zero and choosing an arbitrary pattern.
#'
#' @examples
#' ss_matrix <- matrix(c(125, 125, 125, 125, NA, 250), nrow = 3, byrow = TRUE)
#' head(datasim_bin_2(SS_matrix = ss_matrix,
#' p0 = 0.7, OR = rep(1.8, 2), lambda = rep(0.15, 3), trend="stepwise_2"))
#'
#'
#' @return Data frame: simulated trial data (if full=FALSE, i.e. default) with the following columns:
#'
#' - `j` - patient recruitment index
#' - `response` - binary response for patient `j`
#' - `treatment`- index of the treatment patient `j` was allocated to
#' - `period` - index of the period patient `j` was recruited in
#'
#' or List (if full=TRUE) containing the following elements:
#'
#' - `Data` - simulated trial data, including an additional column `p` with the probability used for simulating the response for patient `j`
#' - `n_total` - total sample size in the trial
#' - `num_arms` - number of experimental treatment arms in the trial
#' - `SS_matrix` - matrix with the sample sizes per arm and per period
#' - `period_blocks` - number to multiply the number of active arms with, in order to get the block size per period
#' - `p0` - response probability in the control arm
#' - `OR` - odds ratios for each experimental treatment arm
#' - `lambda` - strength of time trend in each arm
#' - `time_dep_effect` - time dependent treatment effects for each experimental treatment arm (for computing the bias)
#' - `trend` - time trend pattern
#'
#'
#' @author Pavla Krotka, Marta Bofill Roig

datasim_bin_2 <- function(num_arms, n_arm, d, period_blocks=2, p0, OR, SS_matrix=NULL, lambda, trend, N_peak, n_wave, trend_mean=0, trend_var=0.5, full=FALSE, check=TRUE){
  
  if (!is.null(SS_matrix)) {
    num_arms <- nrow(SS_matrix)-1 # total number of experimental arms
    n_arm <- NULL
    d <- NULL
  }
  
  if (check) {
    if(!is.null(num_arms) && (!is.numeric(num_arms) | length(num_arms)!=1)){
      stop("Number of experimental treatment arms (`num_arms`) must be one number!")
    }
    
    if(!is.null(n_arm) &&  (!is.numeric(n_arm) | length(n_arm)!=1)){
      stop("Sample size per experimental treatment arm (`n_arm`) must be one number! Equal sample sizes in all experimental arms are assumed!")
    }
    
    if (!is.null(d) && (!is.numeric(d) | length(d)!=num_arms | d[1]!=0)) {
      stop("Timing of adding new experimental arms (`d`) must be a numeric vector of length `num_arms` and the first entry must be 0!")
    }
    
    if(!is.numeric(period_blocks) | length(period_blocks)!=1){
      stop("`period_blocks` must be one number!")
    }
    
    if(!is.numeric(p0) | length(p0)!=1){
      stop("Response in the control arm (`p0`) must be one number!")
    }
    
    if (!is.numeric(OR) | length(OR)!=num_arms) {
      stop("Vector with odds ratios (`OR`) must be a numeric vector of length `num_arms` and must be ordered by the entry time of the treatment arms
           (i.e., the first entry in the vector corresponds to the odds ratio of the first experimental treatment arm etc.)!")
    }
    
    if ((sum(!is.numeric(lambda) | lambda!="random")!=length(lambda)) | length(lambda)!=(num_arms+1)) {
      stop("The entries in the vector with the strength of the time trend in each arm (`lambda`) must be numeric or the string 'random', the length of the vector must be `num_arms`+1 and it must be ordered by the entry time of the arms
           (i.e., the first entry in the vector corresponds to the time trend in the control arm, second entry to the time trend in the first treatment arm etc.)!")
    }
    
    if((trend %in% c("linear", "linear_2", "stepwise", "stepwise_2", "inv_u", "seasonal")==FALSE) | length(trend)!=1){
      stop("Time trend pattern (`trend`) must be one of the following strings: 'linear', 'linear_2', 'stepwise', 'stepwise_2', 'inv_u', 'seasonal'!")
    }
    
    if(trend=="inv_u"){
      if(!is.numeric(N_peak) | length(N_peak)!=1){
        stop("Point at which the inverted-u time trend switches direction (`N_peak`) must be one number!")
      }
    }
    
    if(trend=="seasonal"){
      if(!is.numeric(n_wave) | length(n_wave)!=1){
        stop("Number of cycles in the seasonal trend (`n_wave`) must be one number!")
      }
    }
    
    if("random" %in% lambda){
      if(!is.numeric(trend_mean) | length(trend_mean)!=1){
        stop("The mean of the random time trend (`trend_mean`) must be one number!")
      }
    }
    
    if("random" %in% lambda){
      if(!is.numeric(trend_var) | length(trend_var)!=1){
        stop("The variance of the random time trend (`trend_var`) must be one number!")
      }
    }
    
  }
  
  if (is.null(SS_matrix)) {
    
    SS_matrix <- get_ss_matrix(num_arms, n_arm, d)
    
    SS_matrix <- round(SS_matrix)
    
    alloc_ratios <- ifelse(!is.na(SS_matrix), 1, 0)
    
  } else {
    
    alloc_ratios <- ifelse(is.na(SS_matrix), 0, SS_matrix/SS_matrix[1,])
    
  }
  
  num_periods <- ncol(SS_matrix) # total number of periods
  
  N_period <- colSums(SS_matrix, na.rm=T) # sample sizes per period
  N_arm <- rowSums(SS_matrix, na.rm=T) # sample sizes per arm
  n_total <- sum(SS_matrix, na.rm = T) # total sample size
  active_arms <- colSums(apply(SS_matrix, 2, is.na)==0) # active arms per period
  block_sizes <- period_blocks*colSums(alloc_ratios) # block sizes per period
  
  
  t <- c()
  
  for (i in 1:num_periods){
    m_i <- t(replicate(trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i]),
                       sample(rep(rep(c(0:(num_arms)), alloc_ratios[,i]), block_sizes[i]/length(rep(c(0:(num_arms)), alloc_ratios[,i]))))))
    
    rest <- sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])
    
    t_i <- c(t(m_i), sample(rep(c(0:(num_arms)), alloc_ratios[,i]*ceiling(rest/active_arms[i])),
                            size = sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])))
    
    t <- c(t, t_i)
  }
  
  
  j <- c(1:n_total)
  
  
  for (i in 0:num_arms) {
    assign(paste0("j", i), which(t==i)) # j0, j1, j2 ... position in time (order) of allocated patients in every arm
  }
  
  cj <- rep(1:num_periods, N_period) # period indicator
  
  if(num_periods>1){
    arm_added <- 1
    
    for (i in 2:(num_periods)) {
      prev_per <- max(which(!is.na(SS_matrix[ ,i-1])))
      if (prev_per<nrow(SS_matrix)) {
        arm_added[i] <- ifelse(!is.na(SS_matrix[prev_per+1, i]), arm_added[i-1]+sum(!is.na(SS_matrix[(prev_per+1):nrow(SS_matrix), i])), arm_added[i-1])
      } else {
        arm_added[i] <- arm_added[i-1]
      }
    }
    
    cj_added <- rep(arm_added, times = N_period) # indicator of jumps only if a new treatment enters (double jumps if 2 treatments enter etc...)
    
  } else {
    
    cj_added <- cj
  }
  
  
  # Simulation of individual trend
  
  if(trend=="linear"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                  lambda = lambda[i+1],
                                                  sample_size = c(0, n_total),
                                                  trend_mean = trend_mean,
                                                  trend_var = trend_var))
      
      assign(paste0("all_trend", i), linear_trend(j=j,
                                                  lambda = lambda[i+1],
                                                  sample_size = c(0, n_total),
                                                  trend_mean = trend_mean,
                                                  trend_var = trend_var))
    }
  }
  
  if(trend=="linear_2"){ # trend starts in the second period and is linear
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                  lambda = lambda[i+1],
                                                  sample_size = c(N_period[1], sum(N_period[-1])),
                                                  trend_mean = trend_mean,
                                                  trend_var = trend_var))
      
      assign(paste0("all_trend", i), linear_trend(j=j,
                                                  lambda = lambda[i+1],
                                                  sample_size = c(N_period[1], sum(N_period[-1])),
                                                  trend_mean = trend_mean,
                                                  trend_var = trend_var))
    }
  }
  
  if(trend=="stepwise"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), sw_trend(cj=cj[eval(sym(paste0("j", i)))],
                                              lambda = lambda[i+1],
                                              trend_mean = trend_mean,
                                              trend_var = trend_var))
      
      assign(paste0("all_trend", i), sw_trend(cj=cj,
                                              lambda = lambda[i+1],
                                              trend_mean = trend_mean,
                                              trend_var = trend_var))
    }
  }
  
  if(trend=="stepwise_2"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), sw_trend(cj=cj_added[eval(sym(paste0("j", i)))],
                                              lambda = lambda[i+1],
                                              trend_mean = trend_mean,
                                              trend_var = trend_var))
      
      assign(paste0("all_trend", i), sw_trend(cj=cj_added,
                                              lambda = lambda[i+1],
                                              trend_mean = trend_mean,
                                              trend_var = trend_var))
    }
  }
  
  
  if(trend=="inv_u"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), inv_u_trend(j=eval(sym(paste0("j", i))),
                                                 N_peak = N_peak,
                                                 lambda = lambda[i+1],
                                                 n_total = n_total,
                                                 trend_mean = trend_mean,
                                                 trend_var = trend_var))
      
      assign(paste0("all_trend", i), inv_u_trend(j=j,
                                                 N_peak = N_peak,
                                                 lambda = lambda[i+1],
                                                 n_total = n_total,
                                                 trend_mean = trend_mean,
                                                 trend_var = trend_var))
    }
  }
  
  
  if(trend=="seasonal"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), seasonal_trend(j=eval(sym(paste0("j", i))),
                                                    lambda = lambda[i+1],
                                                    n_wave = n_wave,
                                                    n_total = n_total,
                                                    trend_mean = trend_mean,
                                                    trend_var = trend_var))
      
      assign(paste0("all_trend", i), seasonal_trend(j=j,
                                                    lambda = lambda[i+1],
                                                    n_wave = n_wave,
                                                    n_total = n_total,
                                                    trend_mean = trend_mean,
                                                    trend_var = trend_var))
    }
  }
  
  # Simulation of binary endpoint
  
  p <- c()
  eta0 <- log(p0/(1-p0)) + ind_trend0
  
  for (i in 1:num_arms) {
    assign(paste0("eta", i), log(p0/(1-p0)) + log(OR[i]) + eval(sym(paste0("ind_trend", i))))
  }
  
  for (i in 0:num_arms) {
    p[eval(sym(paste0("j", i)))] <- 1/(1 + exp(-eval(sym(paste0("eta", i)))))
  }
  
  
  if(sum(p<0 | p>1)>0){ # check if all probabilities are between 0 and 1
    stop("p must be between 0 and 1")
  }
  
  X <- rbinom(n = n_total, size = 1, prob = p)
  
  # Compute time dependent treatment effects for bias computation
  
  eta0_all <- log(p0/(1-p0)) + all_trend0
  
  for (i in 1:num_arms) {
    assign(paste0("eta", i, "_all"), log(p0/(1-p0)) + log(OR[i]) + eval(sym(paste0("all_trend", i))))
  }
  
  
  p_time_dep <- data.frame(all_p0 = 1/(1 + exp(-eta0_all)))
  
  for (i in 1:num_arms) {
    p_time_dep[ ,paste0("all_p", i)] <- 1/(1 + exp(-eval(sym(paste0("eta", i, "_all")))))
  }
  
  p_time_dep$period <- cj
  
  time_dep_effect <- rep(NA, num_arms)
  
  for (i in 1:num_arms) {
    active_periods <- which(!is.na(SS_matrix[i+1,]))
    
    p_aux <- p_time_dep[c(1,i+1, ncol(p_time_dep))]
    p_aux <- p_aux[p_aux$period %in% active_periods,]
    
    time_dep_effect[i] <- (mean(p_aux[,2])/(1-mean(p_aux[,2]))) / (mean(p_aux[,1])/(1-mean(p_aux[,1]))) # ODDS_trt/ODDS_ctrl
  }
  
  # Resulting data frame
  
  if (full) {
    Data <- data.frame(j = c(1:n_total),
                       response = X,
                       treatment = t,
                       period = rep(1:num_periods, N_period),
                       p = p)
    
    for (i in 0:(num_arms)) {
      Data[ ,paste0("lambda", i)] <- ifelse(lambda[i+1]=="random", lambda[i+1], as.numeric(lambda[i+1]))
    }
    
    return(list(Data = Data,
                n_total = n_total,
                #n_arm = N_arm[2], # assumed equal
                num_arms = num_arms,
                #d = d,
                SS_matrix = SS_matrix,
                period_blocks = period_blocks,
                p0=p0,
                OR = OR,
                lambda = lambda,
                time_dep_effect = time_dep_effect,
                trend = trend))
    
  } else {
    Data <- data.frame(j = c(1:n_total),
                       response = X,
                       treatment = t,
                       period = rep(1:num_periods, N_period))
    
    return(Data)
  }
}

