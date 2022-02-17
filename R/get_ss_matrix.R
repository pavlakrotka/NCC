#' Sample size matrix for a platform trial with 3 arms
#'
#' @description Computes the sample size matrix - sample sizes per arm (rows) and per period (columns)
#'
#' @param T_ Timing of adding new arms in terms of number of patients allocated to the control arm
#' @param N Overall sample size for the trial
#' 
#' @keywords internal
#'
#' @export
#' @return Sample size matrix
#' @author Pavla Krotka

get_ss_matrix <- function(T_, N){ 
  n = (N-2*T_)/4
  
  if(2*T_ >= n){
  
  control = round(c(T_, n-T_, 2*T_-n, n-T_, T_))
  control = control[control!=0]
  
  if(!is.na(which(control < 0)[1])){
    control = control[1:(which(control < 0)[1]-1)]
  }
  
  arm_1 = rep(NA, length(control))
  arm_2 = rep(NA, length(control))
  arm_3 = rep(NA, length(control))
  
  i <- 1
  while (sum(arm_1, na.rm = T)<n) {
    arm_1[i] <- control[i]
    i = i+1
  }
  
  if(length(control)>1){
    i <- 2
    while (sum(arm_2, na.rm = T)<n) {
      arm_2[i] <- control[i]
      i = i+1
    }
    
    i <- length(control)
    while (sum(arm_3, na.rm = T)<n) {
      arm_3[i] <- control[i]
      i = i-1
    }
  } else {
    arm_2[1] <- control[1]
    arm_3[1] <- control[1]
  }
  
  } else if(2*T_<n) {
    
    control = round(c(T_, T_, n-2*T_, T_, T_))
    control = control[control!=0]
    
    if(!is.na(which(control < 0)[1])){
      control = control[1:(which(control < 0)[1]-1)]
    }
    
    arm_1 = rep(NA, length(control))
    arm_2 = rep(NA, length(control))
    arm_3 = rep(NA, length(control))
    
    i <- 1
    while (sum(arm_1, na.rm = T)<n) {
      arm_1[i] <- control[i]
      i = i+1
    }
    
    if(length(control)>1){
      i <- 2
      while (sum(arm_2, na.rm = T)<n) {
        arm_2[i] <- control[i]
        i = i+1
      }
      
      i <- length(control)
      while (sum(arm_3, na.rm = T)<n) {
        arm_3[i] <- control[i]
        i = i-1
      }
    } else {
      arm_2[1] <- control[1]
      arm_3[1] <- control[1]
    }
    
    
  }
  
  SS_matrix = matrix(c(control, arm_1, arm_2, arm_3), byrow = T, nrow = 4)
  
  return(SS_matrix)
}



# get_ss_matrix(0, 1000)
# get_ss_matrix(100, 1000)
# get_ss_matrix(120, 1000)
# get_ss_matrix(130, 1000)
# get_ss_matrix(140, 1000)
# get_ss_matrix(150, 1000)
# get_ss_matrix(166.66, 1000)
# 
# get_ss_matrix(90, 1000)
# get_ss_matrix(80, 1000)
# get_ss_matrix(70, 1000)
# get_ss_matrix(60, 1000)
# get_ss_matrix(50, 1000)
# get_ss_matrix(40, 1000)
# get_ss_matrix(30, 1000)
# get_ss_matrix(20, 1000)
# get_ss_matrix(10, 1000)

