#' Sample size matrix for a platform trial with 3 arms
#'
#' @description Computes the sample size matrix - sample sizes per arm (rows) and per period (columns)
#'
#' @param num_arms Number of treatment arms in the trial
#' @param n_arm Sample size per arm
#' @param d Vector with timings of adding new arms in terms of number of patients recruited to the trial so far (of length num_arms). The first entry must be 0, as the trial is supposed to start with at least one treatment
#' @param alloc_ratio Vector of length 2. Allocation ratio control vs treatment (for now only 1:1 implemented!)
#'
#'
#' @export
#'
#' @examples
#'
#' get_ss_matrix(num_arms = 3, n_arm = 100, d = c(0, 100, 250))
#'
#' @return Sample size matrix
#' @author Pavla Krotka
#'

get_ss_matrix <- function(num_arms, n_arm, d, alloc_ratio = c(1,1)){

  if(all(d==0)) { # if all arms enter at the beginning
    SS_matrix <- matrix(rep(n_arm, num_arms+1), nrow = num_arms+1)

  } else { # general case

    SS_matrix <- matrix(c(1,1), nrow = 2) # start with 1 patient in control and 1st treatment

    i <- 1
    j <- 1
    l <- 2
    while (sum(SS_matrix[-1,], na.rm = T)<n_arm*num_arms) { # while not all arms are finished

      SS_matrix[,i] <- SS_matrix[,i]+1

      if(sum(SS_matrix, na.rm = T)>=d[j+1] & j<length(d)){ # check if a new arm should enter

        n_add <- sum(d<=sum(SS_matrix, na.rm = T)) - nrow(SS_matrix)+1

        SS_matrix <- rbind(SS_matrix, matrix(NA, ncol = ncol(SS_matrix), nrow = n_add))

        finished_arms <- if(dim(SS_matrix)[2]==1){which(SS_matrix[-1,]>=n_arm)+1} else {which(rowSums(SS_matrix[-1,], na.rm = T)>=n_arm)+1}

        new_period <- rep(1, nrow(SS_matrix))
        new_period[finished_arms] <- NA
        SS_matrix <- cbind(SS_matrix, new_period)
        colnames(SS_matrix) <- NULL

        j <- j+n_add
        i <- i+1

      } else {

        for (k in l:(nrow(SS_matrix))) { # check for all arms if they are finished
          if(sum(SS_matrix[-1,], na.rm = T)<n_arm*num_arms) {
            if(sum(SS_matrix[k,], na.rm = T)>=n_arm){
              finished_arms <- if(dim(SS_matrix)[2]==1){which(SS_matrix[-1,]>=n_arm)+1} else {which(rowSums(SS_matrix[-1,], na.rm = T)>=n_arm)+1}
              new_period <- rep(1, nrow(SS_matrix))
              new_period[finished_arms] <- NA
              SS_matrix <- cbind(SS_matrix, new_period)
              colnames(SS_matrix) <- NULL

              i <- i+1
              l <- l+1
            }
          }
        }
      }
    }

    # correct periods if needed
    del <- NULL
    for (i in 1:(ncol(SS_matrix)-1)) {
      if(sum(is.na(SS_matrix[,i])==is.na(SS_matrix[,i+1]))==nrow(SS_matrix)){
        SS_matrix[,i] <- SS_matrix[,i]+SS_matrix[,i+1]
        del <- c(del, i+1)
      }
    }

    if(!is.null(del)) { SS_matrix <- SS_matrix[,-del] }

  }
  return(SS_matrix)
}


# get_ss_matrix(n_arm = 100, num_arms = 3, d = c(0, 200, 400))
#
# get_ss_matrix(n_arm = 100, num_arms = 3, d = c(0, 24, 60))
#
# get_ss_matrix(n_arm = 100, num_arms = 3, d = c(0, 0, 0))
#
# get_ss_matrix(n_arm = 100, num_arms = 3, d = c(0, 200, 260))
#
#
#
# get_ss_matrix(n_arm = 100, num_arms = 5, d = c(0, 24, 60, 100, 200))
#
# get_ss_matrix(n_arm = 100, num_arms = 4, d = c(0, 200, 200, 250))
#
#
# get_ss_matrix(n_arm = 100, num_arms = 4, d = c(0, 200, 250, 250))
#
# get_ss_matrix(n_arm = 100, num_arms = 3, d = c(0, 100, 250))


###############################################################################################################################

# Old version (equal d for all arms and in terms of control group, not the whole trial)

# get_ss_matrix <- function(n_total, num_arms, d) {
#   n_arm <- (n_total-(num_arms-1)*d)/(num_arms+1)
#
#   if(n_arm<d | 2*d < n_arm){
#     stop("n_arm>=d & 2*d >= n_arm must hold!")
#   }
#
#   num_periods <- (2*num_arms)-1
#
#   control_cumsum <- rep(NA, num_periods)
#   control_cumsum[length(control_cumsum)] <- n_arm+((num_arms-1)*d) # last entry
#   control_cumsum[-length(control_cumsum)][c(TRUE, FALSE)] <- d*seq(1, (num_periods-1)/2) # odd entries
#   control_cumsum[-length(control_cumsum)][c(FALSE, TRUE)] <- n_arm + d*seq(0, ((num_periods-1)/2)-1) # even entries
#
#   control <- rep(NA, num_periods)
#   control[1] <- control_cumsum[1]
#   control[-1] <- diff(control_cumsum)
#   control <- round(control)
#   control <- control[control!=0]
#
#   SS_matrix <- matrix(nrow = num_arms+1, ncol = length(control))
#   SS_matrix[1,] <- control
#
#   for(i in 1:num_arms){
#     arm_start <- which(cumsum(control)>(i-1)*round(d))[1]
#     control_new <- control
#     control_new[if(arm_start==1) NA else c(1:(arm_start-1))] <- 0
#     arm_stop <- which(cumsum(control_new)>=round(n_arm))[1]
#
#     SS_matrix[i+1, arm_start:arm_stop] <- control[arm_start:arm_stop]
#
#   }
#
#   return(SS_matrix)
# }





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




