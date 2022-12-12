#' Model-based analysis for continuous data using regression splines with knots placed according to calendar times
#'
#' @description Performs linear regression taking into account all trial data until the arm under study leaves the trial and adjusting for time using regression splines with knots placed according to calendar times
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param unit_size Number of patients per calendar time unit Default=25
#' @param ncc Boolean. Whether to include NCC data into the analysis. Default=TRUE
#' @param bs_degree Degree of the piecewise polynomial. Default=3
#' @param ... Further arguments for simulation function
#'
#' @importFrom stats lm
#' @importFrom splines bs
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
#' splines_cal_cont(data = trial_data, arm = 3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not for the investigated treatment and the fitted model
#' @author Pavla Krotka

splines_cal_cont <- function(data, arm, alpha=0.025, unit_size=25, ncc=TRUE, bs_degree=3, ...){

  data$cal_time <- rep(c(1:ceiling((nrow(data)/unit_size))), each=unit_size)[1:nrow(data)]

  min_unit <- min(data[data$treatment==arm,]$cal_time)
  max_unit <- max(data[data$treatment==arm,]$cal_time)

  if (ncc) {
    data_new <- data[data$cal_time %in% c(1:max_unit),]
  } else {
    data_new <- data[data$cal_time %in% c(min_unit:max_unit),]
  }

  bs_knots <- seq(unit_size, (nrow(data_new)-unit_size), by = unit_size) # get knots based on length of cal. time unit

  # fit linear model
  if(length(unique(data_new$cal_time))==1){ # if only one calendar time unit in the data, don't use any knots, just ordinary polynomial regression
    mod <- lm(response ~ as.factor(treatment) + bs(j, degree = bs_degree), data_new)
  } else {

    mod <- lm(response ~ as.factor(treatment) + bs(j, knots = bs_knots, degree = bs_degree), data_new)
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
              reject_h0 = reject_h0,
              knots = bs_knots,
              model = mod))
}
