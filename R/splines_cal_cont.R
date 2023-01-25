#' Spline regression analysis for continuous data with knots placed according to calendar time units
#'
#' @description This function performs linear regression taking into account all trial data until the arm under study leaves the trial and adjusting for time using regression splines with knots placed according to calendar time units.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment', 'response' and 'j'.
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Type I error rate. Default=0.025.
#' @param unit_size Number of patients per calendar time unit. Default=25.
#' @param ncc Boolean. Whether to include NCC data into the analysis. Default=TRUE.
#' @param bs_degree Degree of the polynomial spline. Default=3 for cubic spline.
#' @param check Boolean. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments for simulation function.
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
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not (for the investigated treatment specified in the input), the positions of the knots in terms of patient index, and the fitted model.
#' @author Pavla Krotka

splines_cal_cont <- function(data, arm, alpha=0.025, unit_size=25, ncc=TRUE, bs_degree=3, check=TRUE, ...){

  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response", "j") %in% colnames(data))!=3) {
      stop("The data frame with trial data must contain the columns 'treatment', 'response' and 'j'!")
    }

    if(!is.numeric(arm) | length(arm)!=1){
      stop("The evaluated treatment arm (`arm`) must be one number!")
    }

    if(!is.numeric(alpha) | length(alpha)!=1){
      stop("The significance level (`alpha`) must be one number!")
    }

    if(!is.numeric(unit_size) | length(unit_size)!=1){
      stop("The length of calendar time unit (`unit_size`) must be one number!")
    }

    if(!is.logical(ncc) | length(ncc)!=1){
      stop("The indicator of including NCC data to the analysis (`ncc`) must be TRUE or FALSE!")
    }

    if(!is.numeric(bs_degree) | length(bs_degree)!=1){
      stop("Degree of the piecewise polynomial (`bs_degree`) must be one number!")
    }
  }

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
