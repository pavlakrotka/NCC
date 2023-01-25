#' Mixed regression model analysis for continuous data adjusting for periods as a random factor with AR1 correlation structure
#'
#' @description This function performs linear mixed model regression taking into account all trial data until the arm under study leaves the trial and adjusting for periods as random factors with AR1 correlation structure.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Type I error rate. Default=0.025.
#' @param ci Boolean. Whether confidence intervals should be computed. Default=FALSE.
#' @param ncc Boolean. Whether to include NCC data into the analysis. Default=TRUE.
#' @param check Boolean. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments for simulation function.
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
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not (for the investigated treatment specified in the input), and the fitted model.
#' @author Pavla Krotka

mixmodel_AR1_cont <- function(data, arm, alpha=0.025, ci=FALSE, ncc=TRUE, check=TRUE, ...){

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

    if(!is.logical(ci) | length(ci)!=1){
      stop("The indicator whether confidence intervals should be computed (`ci`) must be TRUE or FALSE!")
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

  # fit linear mixed model
  if(length(unique(data_new$period))==1){ # if only one period in the data, don't use period as covariate
    mod <- lm(response ~ as.factor(treatment), data_new)
    res <- summary(mod)

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
    res <- summary.HLfit(mod, verbose = FALSE)

    # one-sided p-value
    IC <- get_any_IC(mod, verbose = FALSE)
    eff_df <- IC["       effective df:"] # effective degrees of freedom
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
              reject_h0 = reject_h0,
              model = mod))
}

