#' Model-based analysis for continuous data adjusting for periods
#'
#' @description Performs linear regression taking into account all trial data until the arm under study leaves the trial and adjusting for periods as factors
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
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' fixmodel_cont(data = trial_data, arm = 3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Pavla Krotka

fixmodel_cont <- function(data, arm, alpha=0.025){

  max_period <- max(data[data$treatment==arm,]$period)
  data_new <- data[data$period %in% c(1:max_period),]

  # fit linear model
  if(max_period==1){ # if only one period in the data, don't use period as covariate
    mod <- lm(response ~ as.factor(treatment), data_new)
  } else {

    mod <- lm(response ~ as.factor(treatment) + as.factor(period), data_new)
  }
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


# mod_test <- lm(response ~ as.factor(treatment), test %>% filter(period %in% c(1,2,3,4)))
# res_test <- summary(mod_test)
#
# fixmodel_cont(test, arm=2)






