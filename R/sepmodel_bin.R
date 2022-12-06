#' Separate analysis for binary data
#'
#' @description Performs separate analysis (only taking into account concurrent controls) using a logistic model
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param ... Further arguments for simulation function
#'
#' @importFrom stats glm
#' @importFrom stats pnorm
#' @importFrom stats coef
#' @importFrom stats confint
#' @importFrom stats binomial
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' sepmodel_bin(data = trial_data, arm = 3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not for the investigated treatment and the fitted model
#' @author Pavla Krotka

sepmodel_bin <- function(data, arm, alpha=0.025, ...){

  periods <- unique(data[data$treatment==arm,]$period)
  data_new <- data[data$treatment %in% c(0, arm) & data$period %in% periods,]

  # fit logistic model
  mod <- glm(response ~ as.factor(treatment), data_new, family = binomial)
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







