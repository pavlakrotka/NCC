#' Mixed regression model analysis for continuous data adjusting for calendar time units as a random factor
#'
#' @description This function performs linear mixed model regression taking into account all trial data until the arm under study leaves the trial and adjusting for calendar time units as random factors.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment' and 'response'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Significance level (one-sided). Default=0.025.
#' @param ci Logical. Indicates whether confidence intervals should be computed. Default=FALSE.
#' @param unit_size Integer. Number of patients per calendar time unit. Default=25.
#' @param ncc Logical. Indicates whether to include non-concurrent data into the analysis. Default=TRUE.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments passed by wrapper functions when running simulations.
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

mixmodel_cal_cont <- function(data, arm, alpha=0.025, ci=FALSE, unit_size=25, ncc=TRUE, check=TRUE, ...){

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

    if(!is.logical(ci) | length(ci)!=1){
      stop("The indicator whether confidence intervals should be computed (`ci`) must be TRUE or FALSE!")
    }

    if(!is.numeric(unit_size) | length(unit_size)!=1){
      stop("The length of calendar time unit (`unit_size`) must be one number!")
    }

    if(!is.logical(ncc) | length(ncc)!=1){
      stop("The indicator of including NCC data to the analysis (`ncc`) must be TRUE or FALSE!")
    }
  }

  data$cal_time <- rep(c(1:ceiling((nrow(data)/unit_size))), each=unit_size)[1:nrow(data)]

  # min_unit <- min(data[data$treatment==arm,]$cal_time)
  # max_unit <- max(data[data$treatment==arm,]$cal_time)

  min_id <- min(which(data$treatment==arm))
  max_id <- max(which(data$treatment==arm))

  # if (ncc) {
  #   data_new <- data[data$cal_time %in% c(1:max_unit),]
  # } else {
  #   data_new <- data[data$cal_time %in% c(min_unit:max_unit),]
  # }

  if (ncc) {
    data_new <- data[1:max_id,]
  } else {
    data_new <- data[min_id:max_id,]
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
    lower_ci <- confint(mod, level = 1-(2*alpha), parallel="no")[paste0("as.factor(treatment)", arm), 1]
    upper_ci <- confint(mod, level = 1-(2*alpha), parallel="no")[paste0("as.factor(treatment)", arm), 2]
  }

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = ifelse(exists("lower_ci"), lower_ci, "not computed"),
              upper_ci = ifelse(exists("upper_ci"), upper_ci, "not computed"),
              reject_h0 = reject_h0,
              model = mod))
}
