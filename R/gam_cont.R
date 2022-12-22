#' Model-based analysis for continuous data using generalized additive model
#'
#' @description Performs analysis using a generalized additive model taking into account all trial data until the arm under study leaves the trial
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment', 'response', 'period' and 'j'
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param ci Boolean. Whether confidence intervals should be computed. Default=FALSE
#' @param smoothing_basis String indicating the (penalized) smoothing basis to use. Default="tp"
#' @param basis_dim The dimension of the basis used to represent the smooth term. The default depends on the number of variables that the smooth is a function of. Default=-1.
#' @param gam_method The smoothing parameter estimation method. Default="GCV.Cp"
#' @param check Boolean. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE
#' @param ... Further arguments for simulation function
#'
#' @importFrom mgcv gam
#' @importFrom mgcv s
#' @importFrom mgcv summary.gam
#' @importFrom stats pt
#' @importFrom stats qt
#' @importFrom stats vcov
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' gam_cont(data = trial_data, arm = 3, ci = TRUE)
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, an indicator whether the null hypothesis was rejected or not for the investigated treatment and the fitted model
#' @author Pavla Krotka

gam_cont <- function(data, arm, alpha=0.025, ci = FALSE, smoothing_basis = "tp", basis_dim = -1, gam_method = "GCV.Cp", check=TRUE, ...){

  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response", "period", "j") %in% colnames(data))!=4) {
      stop("The data frame with trial data must contain the columns 'treatment', 'response', 'period' and 'j'!")
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

    if((smoothing_basis %in% c("tp", "ts", "ds", "cr", "cs", "cc", "sos", "ps", "cp", "re", "mrf", "gp", "so")==FALSE) | length(smoothing_basis)!=1){
      stop("Smoothing basis (`smoothing_basis`) must be one of the following strings: 'tp', 'ts', 'ds', 'cr', 'cs', 'cc', 'sos', 'ps', 'cp', 're', 'mrf', 'gp', 'so'!")
    }

    if(!is.numeric(basis_dim) | length(basis_dim)!=1){
      stop("The dimension of the basis used to represent the smooth term (`basis_dim`) must be one number!")
    }

    if((gam_method %in% c("GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", "P-ML")==FALSE) | length(gam_method)!=1){
      stop("The smoothing parameter estimation method (`gam_method`) must be one of the following strings: 'GCV.Cp', 'GACV.Cp', 'REML', 'P-REML', 'ML', 'P-ML'!")
    }
  }

  max_period <- max(data[data$treatment==arm,]$period)

  data_new <- data[data$period %in% c(1:max_period),]

  # fit gam model
  mod <- gam(response ~ as.factor(treatment) + s(j, bs = smoothing_basis, k = basis_dim), data = data_new, method = gam_method)

  res <- summary.gam(mod)

  # one-sided p-value
  rdf <- res$residual.df

  p_val <- unname(pt(res$p.t[paste0("as.factor(treatment)", arm)], rdf, lower.tail = FALSE))

  # metrics
  treat_effect <- unname(res$p.coeff[paste0("as.factor(treatment)", arm)])

  reject_h0 <- (p_val < alpha)

  if (ci) {
    Vb <- vcov(mod, unconditional = TRUE)
    se <- sqrt(diag(Vb))

    lower_ci <- treat_effect - (qt(0.975, df = rdf) * se[paste0("as.factor(treatment)", arm)])
    upper_ci <- treat_effect + (qt(0.975, df = rdf) * se[paste0("as.factor(treatment)", arm)])
  }

  return(list(p_val = p_val,
              treat_effect = treat_effect,
              lower_ci = ifelse(exists("lower_ci"), lower_ci, "not computed"),
              upper_ci = ifelse(exists("upper_ci"), upper_ci, "not computed"),
              reject_h0 = reject_h0,
              model = mod))
}
