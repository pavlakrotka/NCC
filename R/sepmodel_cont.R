#' Separate analysis for continuous data
#'
#' @description Performs separate analysis (only taking into account concurrent controls) using a linear model
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
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
#' trial_data <- datasim_cont(n_total = 1000, num_arms = 3, d = 120,
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' sepmodel_cont(data = trial_data, arm = 3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Pavla Krotka

sepmodel_cont <- function(data, arm, alpha=0.025){

  periods <- unique(data[data$treatment==arm,]$period)
  data_new <- data[data$treatment %in% c(0, arm) & data$period %in% periods,]

  # fit linear model
  mod <- lm(response ~ as.factor(treatment), data_new)
  res <- summary(mod)

  # one-sided p-value
  p_val <- pt(coef(res)[paste0("as.factor(treatment)", arm), "t value"], mod$df, lower.tail = FALSE)

  # metrics
  treat_effect <- res$coefficients[paste0("as.factor(treatment)", arm), "Estimate"]
  lower_ci <- confint(mod)[paste0("as.factor(treatment)", arm), 1]
  upper_ci <- confint(mod)[paste0("as.factor(treatment)", arm), 2]
  reject_h0 <- (p_val < alpha)

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0))
}


# mod_test <- lm(response ~ as.factor(treatment), test %>% filter(treatment %in% c(0,2) & period %in% c(1,2,3,4)))
# res_test <- summary(mod_test)
#
# sepmodel_cont(test, arm=2)






