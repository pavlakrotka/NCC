#' Wrapper function for simulations performing inference on given treatment arms using given models
#'
#' @description This function analyzes given data using different models as indicated by the user. It performs inference for indicated experimental treatment arms.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_bin()` or `datasim_cont()` function.
#' @param arms Integer vector with treatment arms to perform inference on. These arms are compared to the control group. Default - all arms except the first one.
#' @param models Character vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel"). Available models for continuous endpoints are: 'fixmodel', 'fixmodel_cal', 'gam', 'MAPprior', 'mixmodel', 'mixmodel_cal', 'mixmodel_AR1', 'mixmodel_AR1_cal', 'piecewise', 'piecewise_cal', 'poolmodel', 'sepmodel', 'sepmodel_adj', 'splines', 'splines_cal', 'timemachine'. Available models for binary endpoints are: 'fixmodel', 'fixmodel_cal', 'MAPprior', 'poolmodel', 'sepmodel', 'sepmodel_adj', 'timemachine'.
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints.
#' @param alpha Double. Significance level (one-sided). Default=0.025.
#' @param unit_size Integer. Number of patients per calendar time unit for frequentist models adjusting for calendar time. Default=25.
#' @param ncc Logical. Whether to include NCC data into the analysis using frequentist models. Default=TRUE.
#' @param opt Integer (1 or 2). In the MAP Prior approach, if opt==1, all former periods are used as one source; if opt==2, periods get separately included into the final analysis. Default=2.
#' @param prior_prec_tau Double. Dispersion parameter of the half normal prior, the prior for the between study heterogeneity in the MAP Prior approach. Default=4.
#' @param prior_prec_eta Double. Dispersion parameter of the normal prior, the prior for the control response (log-odds or mean) in the MAP Prior approach. Default=0.001.
#' @param n_samples Integer. Number of how many random samples will get drawn for the calculation of the posterior mean, the p-value and the CI's in the MAP Prior approach. Default=1000.
#' @param n_chains Integer. Number of parallel chains for the rjags model in the MAP Prior approach. Default=4.
#' @param n_iter Integer. Number of iterations to monitor of the jags.model. Needed for coda.samples in the MAP Prior approach. Default=4000.
#' @param n_adapt Integer. Number of iterations for adaptation, an initial sampling phase during which the samplers adapt their behavior to maximize their efficiency. Needed for jags.model in the MAP Prior approach. Default=1000.
#' @param robustify Logical. Indicates whether a robust prior is to be used. If TRUE, a mixture prior is considered combining a MAP prior and a weakly non-informative component prior. Default=TRUE.
#' @param weight Double. Weight given to the non-informative component (0 < weight < 1) for the robustification of the MAP Prior according to Schmidli (2014). Default=0.1.
#' @param ci Logical. Whether confidence intervals for the mixed models should be computed. Default=FALSE.
#' @param prec_theta Double. Precision (\eqn{1/\sigma^2_{\theta}}) of the prior regarding the treatment effect \eqn{\theta}. I.e. \eqn{\theta \sim N(0, \sigma^2_{\theta})} . Default=0.001.
#' @param prec_eta Double. Precision (\eqn{1/\sigma^2_{\eta_0}}) of the prior regarding the control response \eqn{\eta_0}. I.e. \eqn{\eta_0 \sim N(0, \sigma^2_{\eta_0})}. Default=0.001.
#' @param tau_a Double. Parameter \eqn{a} of the Gamma distribution for the precision parameter \eqn{\tau} in the model for the time trend in the Time Machine approach. I.e., \eqn{\tau \sim Gamma(a,b)}. Default=0.1.
#' @param tau_b Double. Parameter \eqn{b} of the Gamma distribution for the precision parameter \eqn{\tau} in the model for the time trend in the Time Machine approach. I.e., \eqn{\tau \sim Gamma(a,b)}. Default=0.01.
#' @param prec_a Double. Parameter \eqn{a} of the Gamma distribution regarding the precision of the responses for continuous endpoints in the Time Machine approach. I.e., \eqn{\sigma \sim Gamma(a,b)}. Default=0.001.
#' @param prec_b Double. Parameter \eqn{b} of the Gamma distribution regarding the precision of the responses for continuous endpoints in the Time Machine approach. I.e., \eqn{\sigma \sim Gamma(a,b)}. Default=0.001.
#' @param bucket_size Integer. Number of patients per time bucket in the Time Machine approach. Default=25.
#' @param smoothing_basis String indicating the (penalized) smoothing basis to use in the GAM models. Default="tp".
#' @param basis_dim Integer. The dimension of the basis used to represent the smooth term in the GAM models. The default depends on the number of variables that the smooth is a function of. Default=-1.
#' @param gam_method String indicating the smoothing parameter estimation method for the GAM models. Default="GCV.Cp".
#' @param bs_degree Integer. Degree of the polynomial splines. Default=3.
#' @param poly_degree Integer. Degree of the discontinuous piecewise polynomials. Default=3.
#'
#' @importFrom rlang try_fetch
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' all_models(data = trial_data, arms = c(2,3), endpoint = "bin")
#'
#'
#' @return List containing an indicator whether the null hypothesis was rejected or not, and the estimated treatment effect for all investigated treatment arms and all models.
#' @author Pavla Krotka


all_models <- function(data, arms, models = c("fixmodel", "sepmodel", "poolmodel"), endpoint, alpha = 0.025,
                       unit_size = 250,
                       ncc = TRUE,
                       opt = 2, prior_prec_tau = 4, prior_prec_eta = 0.001, n_samples = 1000, n_chains = 4, n_iter = 4000, n_adapt = 1000, robustify = TRUE, weight = 0.1, a_0 = 0.9,
                       ci = FALSE,
                       prec_theta = 0.001, prec_eta = 0.001, tau_a = 0.1, tau_b = 0.01, prec_a = 0.001, prec_b = 0.001, bucket_size = 25,
                       smoothing_basis = "tp", basis_dim = -1, gam_method = "GCV.Cp",
                       bs_degree = 3, poly_degree = 3){

  arms <- sort(arms)
  models <- sort(models)

  res <- list()

  for (i in arms) {
    for (j in models) {

      res_i_j <- try_fetch(do.call(paste0(j, "_", endpoint), list(data = data,
                                                                  arm = i,
                                                                  alpha = alpha,
                                                                  unit_size = unit_size,
                                                                  ncc = ncc,
                                                                  opt = opt,
                                                                  prior_prec_tau = prior_prec_tau,
                                                                  prior_prec_eta = prior_prec_eta,
                                                                  n_samples = n_samples,
                                                                  n_chains = n_chains,
                                                                  n_iter = n_iter,
                                                                  n_adapt = n_adapt,
                                                                  robustify = robustify,
                                                                  weight = weight,
                                                                  ci = ci,
                                                                  prec_theta = prec_theta,
                                                                  prec_eta = prec_eta,
                                                                  tau_a = tau_a,
                                                                  tau_b = tau_b,
                                                                  prec_a = prec_a,
                                                                  prec_b = prec_b,
                                                                  bucket_size = bucket_size,
                                                                  smoothing_basis = smoothing_basis,
                                                                  basis_dim = basis_dim,
                                                                  gam_method = gam_method,
                                                                  bs_degree = bs_degree,
                                                                  poly_degree = poly_degree,
                                                                  a_0 = a_0,
                                                                  check = FALSE))[c("reject_h0", "treat_effect")], error = function(cnd) list(reject_h0 = NA, treat_effect = NA))

      names(res_i_j) <- c(paste0("reject_h0_", j, "_", i), paste0("treat_effect_", j, "_", i))

      res <- append(res, res_i_j)
    }
  }
  return(res)
}




