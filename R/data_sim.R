#' Data simulation for a platform trial with non-concurrent controls with an arbitrary number of treatments and periods
#'
#' @description
#'
#' @param SS_matrix a
#' @param block_sizes a
#' @param alloc_ratios a
#' @param mu0 a
#' @param delta a
#' @param p0 a
#' @param OR a
#' @param lambda a
#' @param sigma a
#' @param N_peak a
#' @param trend a
#' @param trend_param a
#' @param endpoint a
#'
#' @export
#' @return
#' @author
#' @references

data_sim <- function(SS_matrix, block_sizes, alloc_ratios, mu0=0, delta, p0, OR, lambda, sigma, N_peak, trend, trend_param, endpoint){

  require(rlang)

  n_arms <- nrow(SS_matrix) # total number of arms
  n_periods <- ncol(SS_matrix) # total number of periods
  SS_period <- colSums(SS_matrix, na.rm=T) # sample sizes per period
  SS_arm <- rowSums(SS_matrix, na.rm=T) # sample sizes per arm
  SS_total <- sum(SS_matrix, na.rm = T) # total sample size
  active_arms <- colSums(apply(SS_matrix, 2, is.na)==0) # active arms per period

  get_ar <- function(x){ # get allocation ratio
    x/min(x, na.rm = T)
  }

  ar <- apply(SS_matrix, 2, get_ar)

  t <- c()

  for (i in 1:n_periods){
    if(sum(ar[,i]==alloc_ratios[,i], na.rm = T)==length(na.omit(ar[,i]))){
      m_i <- t(replicate(trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i]),
                         sample(rep(rep(c(0:(n_arms-1)), alloc_ratios[,i]), block_sizes[i]/length(rep(c(0:(n_arms-1)), alloc_ratios[,i]))))))

      t_i <- c(t(m_i), sample(rep(c(0:(n_arms-1)), alloc_ratios[,i]),
                              size = sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])))

    }else{stop("wrong allocation ratios!")}

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
                                                  sample_size = SS_total))
    }
  }

  if(trend=="linear2"){ # trend starts in the second period and is linear
    for (i in 0:(n_arms-1)) {
      assign(paste0("ind_trend", i), linear_trend2(j=eval(sym(paste0("j", i))),
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

  if(endpoint=="continuous"){

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
  }

  return(Data)
}