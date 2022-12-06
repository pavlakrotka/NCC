#' Model-based analysis for continuous data using discontinuous piecewise polynomials per period
#'
#' @description Performs linear regression taking into account all trial data until the arm under study leaves the trial and adjusting for time using discontinuous piecewise polynomials in each period
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param ncc Boolean. Whether to include NCC data into the analysis. Default=TRUE
#' @param poly_degree Degree of the piecewise polynomial. Default=3
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
#' piecewise_cont(data = trial_data, arm = 3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not for the investigated treatment and the fitted model
#' @author Pavla Krotka

piecewise_cont <- function(data, arm, alpha=0.025, ncc=TRUE, poly_degree=3, ...){

  data$j <- 1:nrow(data)

  min_period <- min(data[data$treatment==arm,]$period)
  max_period <- max(data[data$treatment==arm,]$period)

  if (ncc) {
    data_new <- data[data$period %in% c(1:max_period),]
  } else {
    data_new <- data[data$period %in% c(min_period:max_period),]
  }

  # fit linear model
  if(length(unique(data_new$period))==1){ # if only one period in the data, don't use any knots, just ordinary polynomial regression
    mod <- lm(response ~ as.factor(treatment) + poly(j, degree = poly_degree, raw = TRUE), data_new)
  } else {

    mod <- lm(response ~ as.factor(treatment) + poly(j, degree = poly_degree, raw = TRUE)*as.factor(period), data_new)
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
              model = mod))
}
