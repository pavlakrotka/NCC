#' Separate analysis for binary data adjusted for periods
#'
#' @description This function performs separate analysis (only taking into account concurrent controls) using a logistic model and adjusting for periods, if the treatment arm stays in the trial for more than one period.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Type I error rate. Default=0.025.
#' @param check Boolean. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.

#' @param ... Further arguments for simulation function
#'
#' @importFrom stats glm
#' @importFrom stats pnorm
#' @importFrom stats coef
#' @importFrom stats confint
#' @importFrom stats binomial
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' sepmodel_adj_bin(data = trial_data, arm = 3)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not (for the investigated treatment specified in the input), and the fitted model.
#' @author Pavla Krotka

sepmodel_adj_bin <- function(data, arm, alpha=0.025, check=TRUE, ...){

  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response", "period") %in% colnames(data))!=3) {
      stop("The data frame with trial data must contain the columns 'treatment', 'response' and 'period'!")
    }

    if(!is.numeric(arm) | length(arm)!=1){
      stop("The evaluated treatment arm (`arm`) must be one number!")
    }

    if(!is.numeric(alpha) | length(alpha)!=1){
      stop("The significance level (`alpha`) must be one number!")
    }
  }

  periods <- unique(data[data$treatment==arm,]$period)
  data_new <- data[data$treatment %in% c(0, arm) & data$period %in% periods,]

  # fit logistic model
  if(length(periods)==1){ # if only one period in the data, don't use period as covariate
    mod <- glm(response ~ as.factor(treatment), data_new, family = binomial)
  } else {

    mod <- glm(response ~ as.factor(treatment) + as.factor(period), data_new, family = binomial)

  }
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
              reject_h0 = reject_h0,
              model = mod))
}






