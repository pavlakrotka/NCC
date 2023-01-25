#' Model-based analysis for treatment-treatment comparison for continuous data adjusting for periods
#'
#' @description This function performs linear regression to compare two treatments taking into account all trial data until the last arm under study leaves the trial and adjusting for periods as factors.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arms Indicator of the treatment arms to be compared (vector of length 2)
#' @param alpha Type I error. Default=0.025
#' @param ... Further arguments for simulation function
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
#' TT_fixmodel_cont(data = trial_data, arms = c(1,3))
#'
#' @return List containing the p-value (one-sided), estimated difference between the treatments, 95% confidence interval, an indicator whether the null hypothesis was rejected or not for the investigated comparison and the fitted model
#' @author Pavla Krotka

TT_fixmodel_cont <- function(data, arms, alpha=0.025, ...){

  data$treatment <- factor(data$treatment, levels = c(arms[1], setdiff(unique(data$treatment), arms[1])))

  max_period <- max(data[data$treatment==max(arms),]$period)
  data_new <- data[data$period %in% c(1:max_period),]

  # fit linear model
  if(max_period==1){ # if only one period in the data, don't use period as covariate
    mod <- lm(response ~ as.factor(treatment), data_new)
  } else {

    mod <- lm(response ~ as.factor(treatment) + as.factor(period), data_new)
  }
  res <- summary(mod)

  # one-sided p-value
  p_val <- pt(coef(res)[paste0("as.factor(treatment)", arms[2]), "t value"], mod$df, lower.tail = FALSE)

  # metrics
  est_diff <- res$coefficients[paste0("as.factor(treatment)", arms[2]), "Estimate"]
  lower_ci <- confint(mod)[paste0("as.factor(treatment)", arms[2]), 1]
  upper_ci <- confint(mod)[paste0("as.factor(treatment)", arms[2]), 2]
  reject_h0 <- (p_val < alpha)

  return(list(p_val = p_val,
              est_diff = est_diff,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0,
              model = mod))
}
