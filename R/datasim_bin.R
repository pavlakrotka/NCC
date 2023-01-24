#' Simulate binary data from a platform trial with a shared control arm and an arbitrary number of treatment arms entering at different time points
#'
#' @description The `datasim_bin()` function simulates data from a platform trial with binary endpoints, an arbitrary number of treatment arms entering at different time points and a shared control arm. The user specifies the timing of adding arms in terms of patients recruited to the trial so far and the the sample size per arm.
#'
#' @param num_arms Number of treatment arms in the trial.
#' @param n_arm Sample size per arm (assumed equal).
#' @param d Vector with timings of adding new arms in terms of number of patients recruited to the trial so far (of length `num_arms`). The first entry must be 0, so that the trial starts with at least one experimental treatment arm.
#' @param period_blocks Number to multiply the number of active arms with, in order to get the block size per period. I.e., block size in each period = `period_blocks`* (number of active arms in the period). Default=2.
#' @param p0 Response in the control arm.
#' @param OR Vector with odds ratios for each treatment arm ordered by the entry time of the treatment arms, i.e., the first entry in the vector corresponds to the odds ratio of the first experimental treatment arm etc. (of length `num_arms`).
#' @param lambda Vector with strength of time trend in each arm ordered by the entry time of the arms (i.e., the first entry in the vector corresponds to the time trend in the control arm, second entry to the time trend in the first treatment arm etc. (of length `num_arms`+1, as time trend in the control is also allowed).
#' @param trend Indicates the time trend pattern ("linear", "stepwise", "stepwise_2", "inv_u" or "seasonal").
#' @param N_peak Point at which the inverted-u time trend switches direction in terms of overall sample size.
#' @param n_wave How many cycles (waves) should a seasonal trend have.
#' @param full Boolean. Indicates whether only variables needed for the analysis should be in the output (FALSE) or also additional information (lambdas, underlying responses, input parameters) should be included (TRUE). Default=FALSE.
#' @param check Boolean. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#'
#' @importFrom stats na.omit
#' @importFrom stats rbinom
#' @importFrom rlang sym
#'
#' @export
#'
#' @details Equal sample sizes in all experimental treatment arms are assumed. Furthermore, allocation ratio of 1:1:...:1 in each period is assumed and block randomization is used to assign patients to the active arms.
#'
#'
#' @examples
#'
#' head(datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise"))
#'
#'
#' @return Data frame: simulated trial data (if full=FALSE, i.e. default) with the following columns:
#'
#' - `j` - patient recruitment index
#' - `response` - continuous response for patient `j`
#' - `treatment`- indicator of the treatment patient `j` was allocated in
#' - `period` - indicator of the periods patient `j` was recruited in
#'
#' or List (if full=TRUE) with the following elements:
#'
#' - `Data` - simulated trial data, including an additional column `p` with the expected value used for the simulation of the response for patient `j`
#' - `n_total` - total sample size in the trial
#' - `n_arm` - sample size per arm (assumed equal)
#' - `num_arms` - number of experimental treatment arms in the trial
#' - `d` - timings of adding new arms
#' - `SS_matrix` - matrix with the sample sizes per arm and per period
#' - `period_blocks` - number to multiply the number of active arms with, in order to get the block size per period
#' - `p0` - response in the control arm
#' - `OR` - odds ratios for each treatment arm
#' - `lambda` - strength of time trend in each arm
#' - `time_dep_effect` - time dependent treatment effects for each treatment arm (for computing the bias)
#' - `trend` - time trend pattern
#'
#'
#' @author Pavla Krotka, Marta Bofill Roig

datasim_bin <- function(num_arms, n_arm, d, period_blocks=2, p0, OR, lambda, trend, N_peak, n_wave, full=FALSE, check=TRUE){

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

    if(!is.numeric(p0) | length(p0)!=1){
      stop("Response in the control arm (`p0`) must be one number!")
    }

    if (!is.numeric(OR) | length(OR)!=num_arms) {
      stop("Vector with odds ratios (`OR`) must be a numeric vector of length `num_arms` and must be ordered by the entry time of the treatment arms
           (i.e., the first entry in the vector corresponds to the odds ratio of the first experimental treatment arm etc.)!")
    }

    if (!is.numeric(lambda) | length(lambda)!=(num_arms+1)) {
      stop("Vector with strength of time trend in each arm (`lambda`) must be a numeric vector of length `num_arms`+1 and must be ordered by the entry time of the arms
           (i.e., the first entry in the vector corresponds to the time trend in the control arm, second entry to the time trend in the first treatment arm etc.)!")
    }

    if((trend %in% c("linear", "stepwise", "stepwise_2", "inv_u", "seasonal")==FALSE) | length(trend)!=1){
      stop("Time trend pattern (`trend`) must be one of the following strings: 'linear', 'stepwise', 'stepwise_2', 'inv_u', 'seasonal'!")
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

    t_i <- c(t(m_i), sample(rep(c(0:(num_arms)), alloc_ratios[,i]),
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

  if(trend=="linear2"){ # trend starts in the second period and is linear
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

  # if(trend=="inv_u"){
  #
  #   for (i in 0:num_arms) {
  #     assign(paste0("j", i, "_1"), which(eval(sym(paste0("j", i))) <= N_peak))
  #     assign(paste0("j", i, "_2"), which(eval(sym(paste0("j", i))) > N_peak))
  #
  #     assign(paste0("j", "_1"), which(j <= N_peak))
  #     assign(paste0("j", "_2"), which(j > N_peak))
  #
  #     assign(paste0("ind_trend", i, "_1"), linear_trend(j = eval(sym(paste0("j", i)))[eval(sym(paste0("j", i, "_1")))],
  #                                                       lambda = lambda[i+1],
  #                                                       sample_size = c(0, n_total)))
  #
  #     assign(paste0("ind_trend", i, "_2"), linear_trend(j = eval(sym(paste0("j", i)))[eval(sym(paste0("j", i, "_2")))]-N_peak+1,
  #                                                       lambda = -lambda[i+1],
  #                                                       sample_size = c(0, n_total)) + ifelse(length(eval(sym(paste0("ind_trend", i, "_1"))))==0, 0,
  #                                                                                             eval(sym(paste0("ind_trend", i, "_1")))[length(eval(sym(paste0("ind_trend", i, "_1"))))]))
  #
  #     assign(paste0("ind_trend", i), c(eval(sym(paste0("ind_trend", i, "_1"))), eval(sym(paste0("ind_trend", i, "_2")))))
  #
  #
  #
  #     assign(paste0("all_trend", i, "_1"), linear_trend(j = j[eval(sym(paste0("j", "_1")))],
  #                                                       lambda = lambda[i+1],
  #                                                       sample_size = c(0, n_total)))
  #
  #     assign(paste0("all_trend", i, "_2"), linear_trend(j = j[eval(sym(paste0("j", "_2")))]-N_peak+1,
  #                                                       lambda = -lambda[i+1],
  #                                                       sample_size = c(0, n_total)) + ifelse(length(eval(sym(paste0("ind_trend", i, "_1"))))==0, 0,
  #                                                                                             eval(sym(paste0("ind_trend", i, "_1")))[length(eval(sym(paste0("ind_trend", i, "_1"))))]))
  #
  #     assign(paste0("all_trend", i), c(eval(sym(paste0("all_trend", i, "_1"))), eval(sym(paste0("all_trend", i, "_2")))))
  #   }
  # }

  if(trend=="inv_u"){
    for (i in 0:num_arms) {
      assign(paste0("ind_trend", i), inv_u_trend(j=eval(sym(paste0("j", i))),
                                                 N_peak = N_peak,
                                                 lambda = lambda[i+1],
                                                 sample_size = n_total))

      assign(paste0("all_trend", i), inv_u_trend(j=j,
                                                 N_peak = N_peak,
                                                 lambda = lambda[i+1],
                                                 sample_size = n_total))
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
      Data[ ,paste0("lambda", i)] <- lambda[i+1]
    }

    return(list(Data = Data,
                n_total = n_total,
                n_arm = N_arm[2], # assumed equal
                num_arms = num_arms,
                d = d,
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

