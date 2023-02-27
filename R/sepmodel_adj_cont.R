#' Separate analysis for continuous data adjusted for periods
#'
#' @description This function performs separate analysis (only taking into account concurrent controls) using a linear model and adjusting for periods, if the treatment arm stays in the trial for more than one period.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Significance level (one-sided). Default=0.025.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments passed by wrapper functions when running simulations.
#'
#' @importFrom stats lm
#' @importFrom stats pt
#' @importFrom stats coef
#' @importFrom stats confint
#'
#' @export
#' 
#' @details 
#' 
#' The adjusted separate analysis takes into account only the data from the evaluated experimental treatment arm and its concurrent controls and adjusts for the time effect by including the factor period (defined as a time interval bounded by any treatment arm entering or leaving the platform). The time is then modelled as a step-function with jumps at the beginning of each period. 
#' Denoting by \eqn{y_j} the response probability for patient \eqn{j}, by \eqn{k_j} the arm patient \eqn{j} was allocated to, and by \eqn{M} the treatment arm under evaluation, the regression model is given by:
#'
#' \deqn{E(y_j) = \eta_0  + \theta_M \cdot I(k_j=M) + \sum_{s=S_{M_1}+1}^{S_{M_2}} \tau_s \cdot I(t_j \in T_{S_s})}
#'
#' where \eqn{\eta_0} is the response in the concurrent controls;
#' \eqn{\theta_M} represents the treatment effect of treatment \eqn{M} as compared to control;
#' \eqn{\tau_s} indicates the stepwise period effect between periods \eqn{S_{M_1}} and \eqn{s} (\eqn{s = S_{M_1}+1, \ldots, S_{M_2}}), where \eqn{S_{M_1}} and \eqn{S_{M_2}} denote the periods, in which the investigated treatment arm joined and left the trial, respectively.
#' 
#' If the data consists of only one period, the period in not used as covariate.
#'
#' @examples
#'
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' sepmodel_adj_cont(data = trial_data, arm = 3)
#'
#' @return List containing the following elements regarding the results of comparing `arm` to control:
#'
#' - `p-val` - p-value (one-sided)
#' - `treat_effect` - estimated treatment effect in terms of the difference in means
#' - `lower_ci` - lower limit of the (1-2*`alpha`)*100% confidence interval
#' - `upper_ci` - upper limit of the (1-2*`alpha`)*100% confidence interval
#' - `reject_h0` - indicator of whether the null hypothesis was rejected or not (`p_val` < `alpha`)
#' - `model` - fitted model
#'
#' @author Pavla Krotka

sepmodel_adj_cont <- function(data, arm, alpha=0.025, check=TRUE, ...){

  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response", "period") %in% colnames(data))!=3) {
      stop("The data frame with trial data must contain the columns 'treatment', 'response' and 'period'!")
    }

    if(!is.numeric(arm) | length(arm)!=1){
      stop("The evaluated treatment arm (`arm`) must be one number!")
    }

    if(!is.numeric(alpha) | length(alpha)!=1 | alpha>=1 | alpha<=0){
      stop("The significance level (`alpha`) must be one number between 0 and 1!")
    }
  }

  periods <- unique(data[data$treatment==arm,]$period)
  data_new <- data[data$treatment %in% c(0, arm) & data$period %in% periods,]

  # fit linear model
  if(length(periods)==1){ # if only one period in the data, don't use period as covariate
    mod <- lm(response ~ as.factor(treatment), data_new)
  } else {

    mod <- lm(response ~ as.factor(treatment) + as.factor(period), data_new)

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
              model = mod))
}



