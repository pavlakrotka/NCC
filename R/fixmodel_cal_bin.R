#' Frequentist logistic regression model analysis for binary data adjusting for calendar time units
#'
#' @description This function performs logistic regression taking into account all trial data until the arm under study leaves the trial and adjusting for calendar time units as factors.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_bin()` function. Must contain columns named 'treatment' and 'response'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Significance level (one-sided). Default=0.025.
#' @param unit_size Integer. Number of patients per calendar time unit. Default=25.
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
#' The model-based analysis adjusts for the time effect by including the factor calendar time unit (defined as time units of fixed length, defined by `ùnit_size`). The time is then modelled as a step-function with jumps at the beginning of each calendar time unit.
#' Denoting by \eqn{y_j} the response probability for patient \eqn{j}, by \eqn{k_j} the arm patient \eqn{j} was allocated to, and by \eqn{M} the treatment arm under evaluation, the regression model is given by:
#'
#' \deqn{g(E(y_j)) = \eta_0  + \sum_{k \in \mathcal{K}_M} \theta_k \cdot I(k_j=k) + \sum_{c=2}^{C_M} \tau_c \cdot I(t_j \in T_{C_c})}
#'
#' where \eqn{g(\cdot)} denotes the logit link function and \eqn{\eta_0} is the log odds in the control arm in the first calendar time unit;
#' \eqn{\theta_k} represents the log odds ratio of treatment \eqn{k} and control for \eqn{k\in\mathcal{K}_M}, where \eqn{\mathcal{K}_M} is the set of treatments
#' that were active in the trial during calendar time units prior or up to the time when the investigated treatment arm left the trial;
#' \eqn{\tau_c} indicates the stepwise calendar time effect in terms of the log odds ratio between calendar time units 1 and \eqn{c} (\eqn{c = 2, \ldots, C_M}), where \eqn{C_M} denotes the calendar time unit, in which the investigated treatment arm left the trial.
#'
#' If the data consists of only one calendar time unit, the calendar time unit in not used as covariate.
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' fixmodel_cal_bin(data = trial_data, arm = 3)
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

fixmodel_cal_bin <- function(data, arm, alpha=0.025, unit_size=25, ncc=TRUE, check=TRUE, ...){

  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response") %in% colnames(data))!=2) {
      stop("The data frame with trial data must contain the columns 'treatment' and 'response'!")
    }

    if(!is.numeric(arm) | length(arm)!=1){
      stop("The evaluated treatment arm (`arm`) must be one number!")
    }

    if(!is.numeric(alpha) | length(alpha)!=1 | alpha>=1 | alpha<=0){
      stop("The significance level (`alpha`) must be one number between 0 and 1!")
    }

    if(!is.numeric(unit_size) | length(unit_size)!=1){
      stop("The length of calendar time unit (`unit_size`) must be one number!")
    }

    if(!is.logical(ncc) | length(ncc)!=1){
      stop("The indicator of including NCC data to the analysis (`ncc`) must be TRUE or FALSE!")
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


  # fit logistic model
  if(length(unique(data_new$cal_time))==1){ # if only one calendar time unit in the data, don't use calendar time unit as covariate
    mod <- glm(response ~ as.factor(treatment), data_new, family = binomial)
  } else {
    mod <- glm(response ~ as.factor(treatment) + as.factor(cal_time), data_new, family = binomial)
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
