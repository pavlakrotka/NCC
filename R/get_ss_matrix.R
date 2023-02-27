#' Sample size matrix for a platform trial with a given number of treatment arms
#'
#' @description This function computes the matrix with sample sizes per arm and period. It is used in the functions `datasim_bin()` and `datasim_cont()`.
#'
#' @param num_arms Integer. Number of experimental treatment arms in the trial.
#' @param n_arm Integer. Sample size per experimental treatment arm.
#' @param d Integer vector with timings of adding new arms in terms of number of patients recruited to the trial so far. The first entry must be 0, so that the trial starts with at least one experimental treatment arm, and the entries must be non-decreasing. The vector length equals `num_arms`.
#'
#' @export
#'
#' @examples
#'
#' get_ss_matrix(num_arms = 3, n_arm = 100, d = c(0, 100, 250))
#'
#' @return Sample size matrix, consisting of the sample size per arm and per period, where the arms are represented in the rows (with the control arm in the first row and the experimental arms coming after ordered by entry time) and the periods are represented in the columns.
#' @author Pavla Krotka
#'

get_ss_matrix <- function(num_arms, n_arm, d){
  
  if(all(d==0)) { # if all arms enter at the beginning
    SS_matrix <- matrix(rep(n_arm, num_arms+1), nrow = num_arms+1)
    
  } else { # general case
    
    SS_matrix <- matrix(c(1), nrow = sum(d==0)+1) # start with 1 patient in control and 1st treatment
    
    i <- 1
    j <- sum(d==0)
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
    ncol_old <- 1 # placeholder
    ncol_new <- 2 # placeholder

    while (ncol_old!=ncol_new) {
      ncol_old <- ncol(SS_matrix)
      ncol_new <- ncol(SS_matrix)
      
      for (i in 1:(ncol(SS_matrix)-1)) {
        if(sum(is.na(SS_matrix[,i])==is.na(SS_matrix[,i+1]))==nrow(SS_matrix)){
          SS_matrix[,i] <- SS_matrix[,i]+SS_matrix[,i+1]
          SS_matrix <- SS_matrix[,-(i+1)]
          ncol_new <- ncol(SS_matrix)
          break
        }
      }
    }
    
  }
  return(SS_matrix)
}
