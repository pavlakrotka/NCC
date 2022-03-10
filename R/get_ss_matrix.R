#' Sample size matrix for a platform trial with 3 arms
#'
#' @description Computes the sample size matrix - sample sizes per arm (rows) and per period (columns)
#'
#' @param n_total Overall sample size for the trial
#' @param num_arms Number of treatment arms in the trial
#' @param d Timing of adding new arms in terms of number of patients allocated to the control arm
#'
#' @keywords internal
#'
#' @export
#' @return Sample size matrix
#' @author Pavla Krotka

get_ss_matrix <- function(n_total, num_arms, d) {
  n_arm <- (n_total-(num_arms-1)*d)/(num_arms+1)

  if(n_arm<d | 2*d < n_arm){
    stop("n_arm>=d & 2*d >= n_arm must hold!")
  }

  num_periods <- (2*num_arms)-1

  control_cumsum <- rep(NA, num_periods)
  control_cumsum[length(control_cumsum)] <- n_arm+((num_arms-1)*d) # last entry
  control_cumsum[-length(control_cumsum)][c(TRUE, FALSE)] <- d*seq(1, (num_periods-1)/2) # odd entries
  control_cumsum[-length(control_cumsum)][c(FALSE, TRUE)] <- n_arm + d*seq(0, ((num_periods-1)/2)-1) # even entries

  control <- rep(NA, num_periods)
  control[1] <- control_cumsum[1]
  control[-1] <- diff(control_cumsum)
  control <- round(control)
  control <- control[control!=0]

  SS_matrix <- matrix(nrow = num_arms+1, ncol = length(control))
  SS_matrix[1,] <- control

  for(i in 1:num_arms){
    arm_start <- which(cumsum(control)>(i-1)*round(d))[1]
    control_new <- control
    control_new[if(arm_start==1) NA else c(1:(arm_start-1))] <- 0
    arm_stop <- which(cumsum(control_new)>=round(n_arm))[1]

    SS_matrix[i+1, arm_start:arm_stop] <- control[arm_start:arm_stop]

  }

  return(SS_matrix)
}





# get_ss_matrix(1000, 2, 150)
# get_ss_matrix(1000, 2, 200)
#
# get_ss_matrix(1000, 3, 100)
# get_ss_matrix(1000, 3, 120)
# get_ss_matrix(1000, 3, 130)
# get_ss_matrix(1000, 3, 140)
# get_ss_matrix(1000, 3, 150)
# get_ss_matrix(1000, 3, 166.66)
#
#
# get_ss_matrix(1000, 4, 90)
# get_ss_matrix(1000, 4, 100)
# get_ss_matrix(1000, 4, 120)
# get_ss_matrix(1000, 4, 125)
#
#
# get_ss_matrix(1000, 5, 100)
# get_ss_matrix(1000, 5, 90)














