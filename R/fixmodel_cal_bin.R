#' Model-based analysis for binary data adjusting for calendar time units
#'
#' @description Performs logistic regression taking into account all trial data until the arm under study leaves the trial and adjusting for calendar time units as factors
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param unit_size Number of patients per calendar time unit Default=25
#'
#' @importFrom stats glm
#' @importFrom stats pnorm
#' @importFrom stats coef
#' @importFrom stats confint
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' fixmodel_cal_bin(data = trial_data, arm = 3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not for the investigated treatment and the fitted model
#' @author Pavla Krotka

fixmodel_cal_bin <- function(data, arm, alpha=0.025, unit_size=25, ...){

  data$cal_time <- rep(c(1:ceiling((nrow(data)/unit_size))), each=unit_size)[1:nrow(data)]

  max_unit <- max(data[data$treatment==arm,]$cal_time)
  data_new <- data[data$cal_time %in% c(1:max_unit),]

  # fit logistic model
  if(max_unit==1){ # if only one calendar time unit in the data, don't use calendar time unit as covariate
    mod <- glm(response ~ as.factor(treatment), data_new, family = "binomial")
  } else {
    mod <- glm(response ~ as.factor(treatment) + as.factor(cal_time), data_new, family = "binomial")
  }
  res <- summary(mod)

  # one-sided p-value
  p_val <- pnorm(coef(res)[paste0("as.factor(treatment)", arm), "z value"], lower.tail = FALSE)

  # metrics
  treat_effect <- res$coefficients[paste0("as.factor(treatment)", arm), "Estimate"]
  lower_ci <- suppressMessages(confint(mod)[paste0("as.factor(treatment)", arm), 1])
  upper_ci <- suppressMessages(confint(mod)[paste0("as.factor(treatment)", arm), 2])
  reject_h0 <- (p_val < alpha)

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0,
              model = mod))
}
