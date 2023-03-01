#' Spline regression analysis for continuous data with knots placed according to periods
#'
#' @description This function performs linear regression taking into account all trial data until the arm under study leaves the trial and adjusting for time using regression splines with knots placed according to periods.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment', 'response', 'period' and 'j'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Significance level (one-sided). Default=0.025.
#' @param ncc Logical. Indicates whether to include non-concurrent data into the analysis. Default=TRUE.
#' @param bs_degree Integer. Degree of the polynomial spline. Default=3 for cubic spline.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments passed by wrapper functions when running simulations.
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
#' splines_cont(data = trial_data, arm = 3)
#'
#' @return List containing the following elements regarding the results of comparing `arm` to control:
#'
#' - `p-val` - p-value (one-sided)
#' - `treat_effect` - estimated treatment effect in terms of the difference in means
#' - `lower_ci` - lower limit of the (1-2*`alpha`)*100% confidence interval
#' - `upper_ci` - upper limit of the (1-2*`alpha`)*100% confidence interval
#' - `reject_h0` - indicator of whether the null hypothesis was rejected or not (`p_val` < `alpha`)
#' - `knots` - positions of the knots in terms of patient index
#' - `model` - fitted model
#'
#' @author Pavla Krotka

splines_cont <- function(data, arm, alpha=0.025, ncc=TRUE, bs_degree=3, check=TRUE, ...){

  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response", "period", "j") %in% colnames(data))!=4) {
      stop("The data frame with trial data must contain the columns 'treatment', 'response', 'period' and 'j'!")
    }

    if(!is.numeric(arm) | length(arm)!=1){
      stop("The evaluated treatment arm (`arm`) must be one number!")
    }

    if(!is.numeric(alpha) | length(alpha)!=1 | alpha>=1 | alpha<=0){
      stop("The significance level (`alpha`) must be one number between 0 and 1!")
    }

    if(!is.logical(ncc) | length(ncc)!=1){
      stop("The indicator of including NCC data to the analysis (`ncc`) must be TRUE or FALSE!")
    }

    if(!is.numeric(bs_degree) | length(bs_degree)!=1){
      stop("Degree of the piecewise polynomial (`bs_degree`) must be one number!")
    }
  }

  min_period <- min(data[data$treatment==arm,]$period)
  max_period <- max(data[data$treatment==arm,]$period)

  if (ncc) {
    data_new <- data[data$period %in% c(1:max_period),]
  } else {
    data_new <- data[data$period %in% c(min_period:max_period),]
  }

  bs_knots <- c() # get start of individual periods for knots
  for (i in unique(data_new$period)) {
    bs_knots <- c(bs_knots, max(data_new[data_new$period==i,"j"]))
  }
  bs_knots <- bs_knots[-length(bs_knots)]

  # fit linear model
  if(length(unique(data_new$period))==1){ # if only one period in the data, don't use any knots, just ordinary polynomial regression
    mod <- lm(response ~ as.factor(treatment) + bs(j, degree = bs_degree), data_new)
  } else {

    mod <- lm(response ~ as.factor(treatment) + bs(j, knots = bs_knots, degree = bs_degree), data_new)
  }
  res <- summary(mod)

  # one-sided p-value
  p_val <- pt(coef(res)[paste0("as.factor(treatment)", arm), "t value"], mod$df, lower.tail = FALSE)

  # metrics
  treat_effect <- res$coefficients[paste0("as.factor(treatment)", arm), "Estimate"]
  lower_ci <- confint(mod, level = 1-(2*alpha))[paste0("as.factor(treatment)", arm), 1]
  upper_ci <- confint(mod, level = 1-(2*alpha))[paste0("as.factor(treatment)", arm), 2]
  reject_h0 <- (p_val < alpha)

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0,
              knots = bs_knots,
              model = mod))
}
