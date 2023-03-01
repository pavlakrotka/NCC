#' Simulate continuous data from a platform trial with a shared control arm and a given number of experimental treatment arms entering at given time points
#'
#' @description This function simulates data from a platform trial with a given number of experimental treatment arms entering at given time points and a shared control arm. The primary endpoint is a continuous endpoint. The user specifies the timing of adding arms in terms of patients recruited to the trial so far and the sample size per arm.
#'
#' @param num_arms Integer. Number of experimental treatment arms in the trial.
#' @param n_arm Integer. Sample size per experimental treatment arm (assumed equal).
#' @param d Integer vector with timings of adding new arms in terms of number of patients recruited to the trial so far. The first entry must be 0, so that the trial starts with at least one experimental treatment arm, and the entries must be non-decreasing. The vector length equals `num_arms`.
#' @param period_blocks Integer. Number to define the size of the blocks for the block randomization. The block size in each period equals `period_blocks`times the number of active arms in the period (see Details). Default=2.
#' @param mu0 Double. Response in the control arm. Default=0.
#' @param theta Double vector with treatment effects in terms of difference of means for each experimental treatment arm compared to control. The elements of the vector (treatment effects) are ordered by the entry time of the experimental treatment arms (e.g., the first entry in the vector corresponds to the treatment effect of the first experimental treatment arm). The vector length equals `num_arms`.
#' @param lambda Double vector with strength of time trend in each arm ordered by the entry time of the arms (e.g., the first entry in the vector corresponds to the time trend in the control arm, second entry to the time trend in the first experimental treatment arm). The vector length equals `num_arms`+1, as time trend in the control is also allowed.
#' @param sigma Double. Standard deviation of the responses.
#' @param trend String indicating the time trend pattern ("linear", "linear_2, "stepwise", "stepwise_2", "inv_u" or "seasonal"). See Details for more information.
#' @param N_peak Integer. Timepoint at which the inverted-u time trend switches direction in terms of overall sample size (i.e. after how many recruited participants the trend direction switches).
#' @param n_wave Integer. Number of cycles (waves) should the seasonal trend have.
#' @param full Logical. Indicates whether the output should be in form of a data frame with variables needed for the analysis only (FALSE) or in form of a list containing more information (TRUE). Default=FALSE.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#'
#' @importFrom stats na.omit
#' @importFrom stats rnorm
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
#' - The duration of a platform trial is divided into so-called periods, defined as time intervals bounded by distinct time points of any treatment arm entering or leaving the platform.
#' Hence, multiple treatment arms entering or leaving at the same time point imply the start of only one additional period.
#' - Allocation ratio of 1:1:...:1 in each period. Furthermore, block randomization is used to assign patients to the active arms. Block size in each period = `period_blocks`* (number of active arms in the period).
#' - If the period sample size is not a multiple of the block size, arms for the remaining participants are chosen by sampling without replacement from a vector containing the indices of active arms replicated times `ceiling(remaining sample size/number of active arms)`.
#'
#' \strong{Data generation:}
#'
#' The continuous response \eqn{y_j} for patient \eqn{j} is generated according to:
#'
#' \deqn{E(y_j) = \eta_0 + \sum_{k=1}^K \cdot I(k_j=k) + f(j)}
#'
#' where \eqn{\eta_0} (parameter `mu0`) and \eqn{\theta_k} (parameter `theta`) are the response in the control arm and the effect of treatment \eqn{k}. \eqn{K} is the total number of treatment arms in the trial (parameter `num_arms`) and \eqn{k_j} is an indicator of the treatment arm patient \eqn{j} is allocated to.
#'
#' The function \eqn{f(j)} denotes the time trend, whose strength is indicated by \eqn{\lambda_{k_j}} (parameter `lambda`) and which can have the following patterns (parameter `trend`):
#'
#' - \strong{"linear"} - trend starts at the beginning of the trial and the mean response increases or decreases linearly with a slope of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1}}, where \eqn{N} is the total sample size in the trial
#' - \strong{"linear_2"} - trend starts after the first period (i.e. there is no time trend in the first period) and the mean response increases or decreases linearly with a slope of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1}}, where \eqn{N} is the total sample size in the trial
#' - \strong{"stepwise"} - the mean response is constant in each period and increases or decreases by \eqn{\lambda} each time any treatment arm enters or leaves the trial (i.e. in each period), according to the function \eqn{f(j) = \lambda_{k_j} \cdot (c_j - 1)}, where \eqn{c_j} is an index of the period patient \eqn{j} was enrolled in
#' - \strong{"stepwise_2"} - the mean response is constant in each period and increases or decreases by \eqn{\lambda} each time a new treatment arm is added to the trial, according to the function \eqn{f(j) = \lambda_{k_j} \cdot (w_j - 1)}, where \eqn{w_j} is an indicator of how many treatment arms have already entered the ongoing trial, when patient \eqn{j} was enrolled
#' - \strong{"inv_u"} - the mean response increases up to the point \eqn{N_p} (parameter `N_peak`) and decreases afterwards, linearly with a slope of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1} (I(j \leq N_p) - I(j > N_p))}, where \eqn{N_p} indicates the point at which the trend turns from positive to negative in terms of the sample size (note that for negative \eqn{\lambda}, the mean response decreases first and increases after)
#' - \strong{"seasonal"} - the mean response increases and decreases periodically with a magnitude of \eqn{\lambda}, according to the function \eqn{f(j) = \lambda \cdot \mathrm{sin} \big( \psi \cdot 2\pi \cdot \frac{j-1}{N-1} \big)}, where \eqn{\psi} indicates how many cycles should the time trend have (parameter `n_wave`)
#'
#' Trials with no time trend can be simulated too, by setting all elements of the vector `lambda` to zero and choosing an arbitrary pattern.
#'
#' @examples
#'
#' head(datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear"))
#'
#'
#' @return Data frame: simulated trial data (if full=FALSE, i.e. default) with the following columns:
#'
#' - `j` - patient recruitment index
#' - `response` - continuous response for patient `j`
#' - `treatment`- index of the treatment patient `j` was allocated to
#' - `period` - index of the period patient `j` was recruited in
#'
#' or List (if full=TRUE) containing the following elements:
#'
#' - `Data` - simulated trial data, including an additional column `means` with the theoretical means used for the simulation of the response for patient `j`
#' - `n_total` - total sample size in the trial
#' - `n_arm` - sample size per arm (assumed equal)
#' - `num_arms` - number of experimental treatment arms in the trial
#' - `d` - timings of adding new arms
#' - `SS_matrix` - matrix with the sample sizes per arm and per period
#' - `period_blocks` - number to multiply the number of active arms with, in order to get the block size per period
#' - `mu0` - response in the control arm
#' - `theta` - treatment effects for each experimental treatment arm
#' - `lambda` - strength of time trend in each arm
#' - `time_dep_effect` - time dependent treatment effects for each experimental treatment arm (for computing the bias)
#' - `sigma` - standard deviation of the responses
#' - `trend` - time trend pattern
#'
#'
#'
#' @author Pavla Krotka, Marta Bofill Roig

datasim_cont <- function(num_arms, n_arm, d, period_blocks=2, mu0=0, theta, lambda, sigma, trend, N_peak, n_wave, full=FALSE, check=TRUE){

  if (check) {
    if(!is.numeric(num_arms) | length(num_arms)!=1){
      stop("Number of experimental treatment arms (`num_arms`) must be one number!")
    }

    if(!is.numeric(n_arm) | length(n_arm)!=1){
      stop("Sample size per experimental treatment arm (`n_arm`) must be one number! Equal sample sizes in all experimental arms are assumed!")
    }

    if (!is.numeric(d) | length(d)!=num_arms | d[1]!=0) {
      stop("Timing of adding new experimental arms (`d`) must be a numeric vector of length `num_arms` and the first entry must be 0!")
    }

    if(!is.numeric(period_blocks) | length(period_blocks)!=1){
      stop("`period_blocks` must be one number!")
    }

    if(!is.numeric(mu0) | length(mu0)!=1){
      stop("Response in the control arm (`mu0`) must be one number!")
    }

    if (!is.numeric(theta) | length(theta)!=num_arms) {
      stop("Vector with treatment effects (`theta`) must be a numeric vector of length `num_arms` and must be ordered by the entry time of the treatment arms
           (i.e., the first entry in the vector corresponds to the effect of the first experimental treatment arm etc.)!")
    }

    if (!is.numeric(lambda) | length(lambda)!=(num_arms+1)) {
      stop("Vector with strength of time trend in each arm (`lambda`) must be a numeric vector of length `num_arms`+1 and must be ordered by the entry time of the arms
           (i.e., the first entry in the vector corresponds to the time trend in the control arm, second entry to the time trend in the first treatment arm etc.)!")
    }

    if(!is.numeric(sigma) | length(sigma)!=1){
      stop("Standard deviation of the responses (`sigma`) must be one number!")
    }

    if((trend %in% c("linear", "linear_2", "stepwise", "stepwise_2", "inv_u", "seasonal")==FALSE) | length(trend)!=1){
      stop("Time trend pattern (`trend`) must be one of the following strings: 'linear', 'linear_2',  'stepwise', 'stepwise_2', 'inv_u', 'seasonal'!")
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

  }

  SS_matrix <- get_ss_matrix(num_arms, n_arm, d)

  SS_matrix <- round(SS_matrix)

  alloc_ratios <- ifelse(!is.na(SS_matrix), 1, 0)
  num_periods <- ncol(alloc_ratios) # total number of periods

  N_period <- colSums(SS_matrix, na.rm=T) # sample sizes per period
  N_arm <- rowSums(SS_matrix, na.rm=T) # sample sizes per arm
  n_total <- sum(SS_matrix, na.rm = T) # total sample size
  active_arms <- colSums(apply(SS_matrix, 2, is.na)==0) # active arms per period
  block_sizes <- period_blocks*active_arms # block sizes per period


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
                                                  sample_size = c(0, n_total)))

      assign(paste0("all_trend", i), linear_trend(j=j,
                                                  lambda = lambda[i+1],
                                                  sample_size = c(0, n_total)))
    }
  }

  if(trend=="linear_2"){ # trend starts in the second period and is linear
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), linear_trend(j=eval(sym(paste0("j", i))),
                                                  lambda = lambda[i+1],
                                                  sample_size = c(N_period[1], sum(N_period[-1]))))

      assign(paste0("all_trend", i), linear_trend(j=j,
                                                  lambda = lambda[i+1],
                                                  sample_size = c(N_period[1], sum(N_period[-1]))))
    }
  }

  if(trend=="stepwise"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), sw_trend(cj=cj[eval(sym(paste0("j", i)))],
                                              lambda = lambda[i+1]))

      assign(paste0("all_trend", i), sw_trend(cj=cj,
                                              lambda = lambda[i+1]))
    }
  }

  if(trend=="stepwise_2"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), sw_trend(cj=cj_added[eval(sym(paste0("j", i)))],
                                              lambda = lambda[i+1]))

      assign(paste0("all_trend", i), sw_trend(cj=cj_added,
                                              lambda = lambda[i+1]))
    }
  }


  if(trend=="inv_u"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), inv_u_trend(j=eval(sym(paste0("j", i))),
                                                 N_peak = N_peak,
                                                 lambda = lambda[i+1],
                                                 n_total = n_total))

      assign(paste0("all_trend", i), inv_u_trend(j=j,
                                                 N_peak = N_peak,
                                                 lambda = lambda[i+1],
                                                 n_total = n_total))
    }
  }


  if(trend=="seasonal"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), seasonal_trend(j=eval(sym(paste0("j", i))),
                                                    lambda = lambda[i+1],
                                                    n_wave = n_wave,
                                                    n_total = n_total))

      assign(paste0("all_trend", i), seasonal_trend(j=j,
                                                    lambda = lambda[i+1],
                                                    n_wave = n_wave,
                                                    n_total = n_total))
    }
  }

  # Simulation of continuous endpoint

  means <- c()
  means[j0] <- ind_trend0

  for (i in 1:num_arms) {
    means[eval(sym(paste0("j", i)))] <- eval(sym(paste0("ind_trend", i))) + theta[i]
  }

  X <- rnorm(n=n_total, mean=mu0+means, sd=sigma)

  # Compute time dependent treatment effects for bias computation

  means_time_dep <- data.frame(all_means0 = all_trend0)

  for (i in 1:num_arms) {
    means_time_dep[ ,paste0("all_means", i)] <- eval(sym(paste0("all_trend", i))) + theta[i]
  }

  means_time_dep <- means_time_dep+mu0

  means_time_dep$period <- cj

  time_dep_effect <- rep(NA, num_arms)

  for (i in 1:num_arms) {
    active_periods <- which(!is.na(SS_matrix[i+1,]))

    means_aux <- means_time_dep[c(1,i+1, ncol(means_time_dep))]
    means_aux <- means_aux[means_aux$period %in% active_periods,]

    time_dep_effect[i] <- mean(means_aux[,2]) - mean(means_aux[,1]) # Mean(trt)-Mean(ctrl)
  }

  # Resulting data frame

  if (full) {
    Data <- data.frame(j = c(1:n_total),
                       response = X,
                       treatment = t,
                       period = rep(1:num_periods, N_period),
                       means = mu0+means)

    for (i in 0:num_arms) {
      Data[ ,paste0("lambda", i)] <- lambda[i+1]
    }

    return(list(Data = Data,
                n_total = n_total,
                n_arm = N_arm[2], # assumed equal
                num_arms = num_arms,
                d = d,
                SS_matrix = SS_matrix,
                period_blocks = period_blocks,
                mu0 = mu0,
                theta = theta,
                lambda = lambda,
                time_dep_effect = time_dep_effect,
                sigma = sigma,
                trend = trend))

  } else {
    Data <- data.frame(j = c(1:n_total),
                       response = X,
                       treatment = t,
                       period = rep(1:num_periods, N_period))
    return(Data)
  }
}

