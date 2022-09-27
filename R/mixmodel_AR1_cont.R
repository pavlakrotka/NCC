#' Model-based analysis for continuous data adjusting for periods as a random factor with AR1 correlation structure
#'
#' @description Performs linear mixed model regression taking into account all trial data until the arm under study leaves the trial and adjusting for periods as random factors with AR1 correlation structure
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param ci Boolean. Whether confidence intervals should be computed. Default=FALSE
#'
#' @importFrom spaMM fitme
#' @importFrom spaMM get_any_IC
#' @importFrom spaMM summary.HLfit
#' @importFrom spaMM confint.HLfit
#' @importFrom stats pt
#' @importFrom stats coef
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' mixmodel_AR1_cont(data = trial_data, arm = 3, ci = TRUE)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Pavla Krotka

mixmodel_AR1_cont <- function(data, arm, alpha=0.025, ci=FALSE, ...){

  max_period <- max(data[data$treatment==arm,]$period)
  data_new <- data[data$period %in% c(1:max_period),]

  # fit linear mixed model
  if(max_period==1){ # if only one period in the data, don't use period as covariate
    mod <- lm(response ~ as.factor(treatment), data_new)
    res <- summary.HLfit(mod)

    # one-sided p-value
    p_val <- pt(coef(res)[paste0("as.factor(treatment)", arm), "t value"], mod$df, lower.tail = FALSE)

    # treatment effect
    treat_effect <- res$coefficients[paste0("as.factor(treatment)", arm), "Estimate"]

    # confidence intervals
    if (ci) {
      lower_ci <- confint(mod)[paste0("as.factor(treatment)", arm), 1]
      upper_ci <- confint(mod)[paste0("as.factor(treatment)", arm), 2]
    }


  } else {
    mod <- fitme(response ~ as.factor(treatment) + AR1(1 | period), data_new)
    res <- summary(mod, verbose = FALSE)

    # one-sided p-value
    IC <- get_any_IC(mod, verbose = FALSE)
    eff_df <- IC["       effective df:"]
    p_val <- pt(res$beta_table[paste0("as.factor(treatment)", arm), "t-value"], eff_df, lower.tail = FALSE)

    # treatment effect
    treat_effect <- res$beta_table[paste0("as.factor(treatment)", arm), "Estimate"]

    # confidence intervals
    if (ci) {
      lower_ci <- confint.HLfit(mod, paste0("as.factor(treatment)", arm), verbose = FALSE)$interval[1]
      upper_ci <- confint.HLfit(mod, paste0("as.factor(treatment)", arm), verbose = FALSE)$interval[2]
    }
  }

  reject_h0 <- (p_val < alpha)

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = ifelse(exists("lower_ci"), lower_ci, "not computed"),
              upper_ci = ifelse(exists("upper_ci"), upper_ci, "not computed"),
              reject_h0 = reject_h0))
}

