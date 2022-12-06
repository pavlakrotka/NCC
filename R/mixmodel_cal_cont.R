#' Model-based analysis for continuous data adjusting for calendar time units as a random factor
#'
#' @description Performs linear mixed model regression taking into account all trial data until the arm under study leaves the trial and adjusting for calendar time units as random factors
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param ci Boolean. Whether confidence intervals should be computed. Default=FALSE
#' @param unit_size Number of patients per calendar time unit Default=25
#' @param ncc Boolean. Whether to include NCC data into the analysis. Default=TRUE
#' @param ... Further arguments for simulation function
#'
#' @importFrom lmerTest lmer
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
#' mixmodel_cal_cont(data = trial_data, arm = 3, ci = TRUE)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not for the investigated treatment and the fitted model
#' @author Pavla Krotka

mixmodel_cal_cont <- function(data, arm, alpha=0.025, ci=FALSE, unit_size=25, ncc=TRUE, ...){

  data$cal_time <- rep(c(1:ceiling((nrow(data)/unit_size))), each=unit_size)[1:nrow(data)]

  min_unit <- min(data[data$treatment==arm,]$cal_time)
  max_unit <- max(data[data$treatment==arm,]$cal_time)

  if (ncc) {
    data_new <- data[data$cal_time %in% c(1:max_unit),]
  } else {
    data_new <- data[data$cal_time %in% c(min_unit:max_unit),]
  }

  # fit linear mixed model
  if(length(unique(data_new$cal_time))==1){ # if only one calendar time unit in the data, don't use unit as covariate
    mod <- lm(response ~ as.factor(treatment), data_new)
    res <- summary(mod)

    # one-sided p-value
    p_val <- pt(coef(res)[paste0("as.factor(treatment)", arm), "t value"], mod$df, lower.tail = FALSE)

  } else {
    mod <- lmer(response ~ as.factor(treatment) + (1 | cal_time), data_new) # using lmerTest
    res <- summary(mod)

    # one-sided p-value
    p_val <- pt(coef(res)[paste0("as.factor(treatment)", arm), "t value"], coef(res)[paste0("as.factor(treatment)", arm), "df"], lower.tail = FALSE)
  }


  # treatment effect
  treat_effect <- res$coefficients[paste0("as.factor(treatment)", arm), "Estimate"]

  reject_h0 <- (p_val < alpha)

  # confidence intervals
  if (ci) {
    lower_ci <- confint(mod, parallel="no")[paste0("as.factor(treatment)", arm), 1]
    upper_ci <- confint(mod, parallel="no")[paste0("as.factor(treatment)", arm), 2]
  }

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = ifelse(exists("lower_ci"), lower_ci, "not computed"),
              upper_ci = ifelse(exists("upper_ci"), upper_ci, "not computed"),
              reject_h0 = reject_h0,
              model = mod))
}
