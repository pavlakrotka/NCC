#' Separate analysis for binary data
#'
#' @description This function performs separate analysis (only taking into account concurrent controls) using a logistic model.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_bin()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Significance level (one-sided). Default=0.025.
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
#' The separate analysis takes into account only the data from the evaluated experimental treatment arm and its concurrent controls and uses a logistic regression model to evaluate the given treatment arm. 
#' Denoting by \eqn{y_j} the response probability for patient \eqn{j}, by \eqn{k_j} the arm patient \eqn{j} was allocated to, and by \eqn{M} the treatment arm under evaluation, the regression model is given by:
#'
#' \deqn{g(E(y_j)) = \eta_0  + \theta_M \cdot I(k_j=M)}
#'
#' where \eqn{g(\cdot)} denotes the logit link function and \eqn{\eta_0} is the log odds in the concurrent controls;
#' \eqn{\theta_M} represents the log odds ratio of treatment \eqn{M} and control.
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' sepmodel_bin(data = trial_data, arm = 3)
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

sepmodel_bin <- function(data, arm, alpha=0.025, check=TRUE, ...){

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

  # fit logistic model
  mod <- glm(response ~ as.factor(treatment), data_new, family = binomial)
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







