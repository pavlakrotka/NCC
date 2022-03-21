#' Sample size matrix for a platform trial with 3 arms
#'
#' @description Computes the sample size matrix - sample sizes per arm (rows) and per period (columns)
#'
#' @param n_total Overall sample size for the trial
#' @param n_arm Sample size per arm
#' @param num_arms Number of treatment arms in the trial
#' @param d Timing of adding new arms in terms of number of patients allocated to the control arm
#' @param alloc_ratio Vector of lenght 2. Allocation ratio control vs treatment
#'
#' @keywords internal
#'
#' @export
#' @return Sample size matrix
#' @author Pavla Krotka
#'
#'


get_alloc_matrix <- function(n_total, n_arm, num_arms, d, alloc_ratio = c(1,1)){

  if(length(d)==1){
    d <- rep(d, num_arms-1)
    d <- cumsum(d)
  }

  if(missing(n_arm)){
    n_arm <- (n_total-sum(d))/(num_arms+1)
  }

  if(missing(n_total)){
    n_total <- sum(d) + ((num_arms+1)*n_arm)
  }

  SS_matrix <- matrix(c(1,1), nrow = 2) # start with 1 patient in control and 1st treatment

  i <- 1
  j <- 1
  while (sum(SS_matrix, na.rm = T)<n_total) {

    SS_matrix[,i] <- SS_matrix[,i]+1

    if(sum(SS_matrix, na.rm = T)>=sum(d[j])){
      SS_matrix <- rbind(SS_matrix, rep(NA, ncol(SS_matrix)))
      SS_matrix <- cbind(SS_matrix, rep(1, nrow(SS_matrix)))
      j <- j+1
      i <- i+1
    }

    for (k in 2:(nrow(SS_matrix))) {
      if(sum(SS_matrix[k,], na.rm = T)>=n_arm){
        new_period <- rep(1, nrow(SS_matrix))
        new_period[k] <- NA
        SS_matrix <- cbind(SS_matrix, new_period)
        colnames(SS_matrix) <- NULL

        i <- i+1
      }
    }
  }
  return(SS_matrix)
}

# get_alloc_matrix(n_total = 100, num_arms = 3, d = 12)
















