#' Indirect treatment-treatment comparison for continuous data
#'
#' @description Performs indirect comparison of two treatments
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arms Indicator of the treatment arms to be compared (vector of length 2)
#' @param alpha Type I error. Default=0.025
#'
#' @importFrom stats lm
#' @importFrom stats pt
#' @importFrom stats coef
#' @importFrom stats confint
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' TT_indirect_cont(data = trial_data, arms = c(1,3))
#'
#' @return List containing the p-value (one-sided), estimated difference between the treatments and an indicator whether the null hypothesis was rejected or not for the investigated comparison
#' @author Pavla Krotka

TT_indirect_cont <- function(data, arms, alpha=0.025){

  # direct estimate for the first treatment
  periods_1 <- unique(data[data$treatment==arms[1],]$period)
  data_new_1 <- data[data$treatment %in% c(0, arms[1]) & data$period %in% periods_1,]

  # fit linear model
  mod_1 <- lm(response ~ as.factor(treatment), data_new_1)
  res_1 <- summary(mod_1)

  # estimate and standard error
  dir_est_1 <- coef(res_1)[paste0("as.factor(treatment)", arms[1]), "Estimate"]
  se_1 <- coef(res_1)[paste0("as.factor(treatment)", arms[1]), "Std. Error"]


  # direct estimate for the second treatment
  periods_2 <- unique(data[data$treatment==arms[2],]$period)
  data_new_2 <- data[data$treatment %in% c(0, arms[2]) & data$period %in% periods_2,]

  # fit linear model
  mod_2 <- lm(response ~ as.factor(treatment), data_new_2)
  res_2 <- summary(mod_2)

  # estimate and standard error
  dir_est_2 <- coef(res_2)[paste0("as.factor(treatment)", arms[2]), "Estimate"]
  se_2 <- coef(res_2)[paste0("as.factor(treatment)", arms[2]), "Std. Error"]

  # t test
  t_stat <- (dir_est_2-dir_est_1)/sqrt(se_1^2+se_2^2)

  # one-sided p-value
  n <- sum(data$treatment==1) # equal sample sizes for all treatments assumed!
  p_val <- 1-pt(t_stat, df=2*n-2)

  # metrics
  #est_diff <- res$coefficients[paste0("as.factor(treatment)", arms[2]), "Estimate"]
  #lower_ci <- confint(mod)[paste0("as.factor(treatment)", arms[2]), 1]
  #upper_ci <- confint(mod)[paste0("as.factor(treatment)", arms[2]), 2]
  reject_h0 <- (p_val < alpha)

  return(list(p_val = p_val,
              est_diff = dir_est_2-dir_est_1,
              #lower_ci = lower_ci,
              #upper_ci = upper_ci,
              reject_h0 = reject_h0))
}
