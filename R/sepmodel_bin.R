#' Separate analysis for binary data
#'
#' @description Performs separate analysis (only taking into account concurrent controls) using a logistic model
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function
#' @param arm Indicator of the treatment arm under study to perform inference on
#' @param alpha Type I error. Default=0.025
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
#' trial_data <- datasim_bin(n_total = 1000, num_arms = 3, d = 120,
#' period_blocks = 2, p0 = 0.7, OR = rep(1.4, 3), lambda = rep(0.15, 4), trend = "linear")
#'
#' sepmodel_bin(data=trial_data, arm=3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Pavla Krotka

sepmodel_bin <- function(data, arm, alpha=0.025){

  periods <- unique(data[data$treatment==arm,]$period)
  data_new <- data[data$treatment %in% c(0, arm) & data$period %in% periods,]

  # fit linear model
  mod <- glm(response ~ as.factor(treatment), data_new, family = "binomial")
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
              reject_h0 = reject_h0))
}


# mod_test <- glm(response ~ as.factor(treatment), test %>% filter(treatment %in% c(0,2) & period %in% c(1,2,3,4)), family = "binomial")
# res_test <- summary(mod_test)
#
# sepmodel_bin(test, arm = 2)






