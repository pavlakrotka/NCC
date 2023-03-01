#' Frequentist logistic regression model analysis for binary data adjusting for periods
#'
#' @description This function performs logistic regression taking into account all trial data until the arm under study leaves the trial and adjusting for periods as factors.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_bin()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Significance level (one-sided). Default=0.025.
#' @param ncc Logical. Indicates whether to include non-concurrent data into the analysis. Default=TRUE.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments passed by wrapper functions when running simulations.
#'
#' @importFrom stats glm
#' @importFrom stats pnorm
#' @importFrom stats coef
#' @importFrom stats confint
#' @importFrom stats binomial
#'
#' @export
#'
#' @details
#'
#' The model-based analysis adjusts for the time effect by including the factor period (defined as a time interval bounded by any treatment arm entering or leaving the platform). The time is then modelled as a step-function with jumps at the beginning of each period.
#' Denoting by \eqn{y_j} the response probability for patient \eqn{j}, by \eqn{k_j} the arm patient \eqn{j} was allocated to, and by \eqn{M} the treatment arm under evaluation, the regression model is given by:
#'
#' \deqn{g(E(y_j)) = \eta_0  + \sum_{k \in \mathcal{K}_M} \theta_k \cdot I(k_j=k) + \sum_{s=2}^{S_M} \tau_s \cdot I(t_j \in T_{S_s})}
#'
#' where \eqn{g(\cdot)} denotes the logit link function and \eqn{\eta_0} is the log odds in the control arm in the first period;
#' \eqn{\theta_k} represents the log odds ratio of treatment \eqn{k} and control for \eqn{k\in\mathcal{K}_M}, where \eqn{\mathcal{K}_M} is the set of treatments
#' that were active in the trial during periods prior or up to the time when the investigated treatment arm left the trial;
#' \eqn{\tau_s} indicates the stepwise period effect in terms of the log odds ratio between periods 1 and \eqn{s} (\eqn{s = 2, \ldots, S_M}), where \eqn{S_M} denotes the period, in which the investigated treatment arm left the trial.
#'
#' If the data consists of only one period (e.g. in case of a multi-arm trial), the period in not used as covariate.
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' fixmodel_bin(data = trial_data, arm = 3)
#'
#' @return List containing the following elements regarding the results of comparing `arm` to control:
#'
#' - `p-val` - p-value (one-sided)
#' - `treat_effect` - estimated treatment effect in terms of the log-odds ratio
#' - `lower_ci` - lower limit of the (1-2*`alpha`)*100% confidence interval
#' - `upper_ci` - upper limit of the (1-2*`alpha`)*100% confidence interval
#' - `reject_h0` - indicator of whether the null hypothesis was rejected or not (`p_val` < `alpha`)
#' - `model` - fitted model
#'
#' @author Pavla Krotka
#'
#' @references On model-based time trend adjustments in platform trials with non-concurrent controls. Bofill Roig, M., Krotka, P., et al. BMC Medical Research Methodology 22.1 (2022): 1-16.

fixmodel_bin <- function(data, arm, alpha=0.025, ncc=TRUE, check=TRUE, ...){

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

    if(!is.logical(ncc) | length(ncc)!=1){
      stop("The indicator of including NCC data to the analysis (`ncc`) must be TRUE or FALSE!")
    }
  }

  min_period <- min(data[data$treatment==arm,]$period)
  max_period <- max(data[data$treatment==arm,]$period)

  if (ncc) {
    data_new <- data[data$period %in% c(1:max_period),]
  } else {
    data_new <- data[data$period %in% c(min_period:max_period),]
  }

  # fit logistic model
  if(length(unique(data_new$period))==1){ # if only one period in the data, don't use period as covariate
    mod <- glm(response ~ as.factor(treatment), data_new, family = binomial)
  } else {
    mod <- glm(response ~ as.factor(treatment) + as.factor(period), data_new, family = binomial)
  }
  res <- summary(mod)

  # one-sided p-value
  p_val <- pnorm(coef(res)[paste0("as.factor(treatment)", arm), "z value"], lower.tail = FALSE)

  # metrics
  treat_effect <- res$coefficients[paste0("as.factor(treatment)", arm), "Estimate"]
  lower_ci <- suppressMessages(confint(mod, level = 1-(2*alpha))[paste0("as.factor(treatment)", arm), 1])
  upper_ci <- suppressMessages(confint(mod, level = 1-(2*alpha))[paste0("as.factor(treatment)", arm), 2])
  reject_h0 <- (p_val < alpha)

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0,
              model = mod))
}
