#' Wrapper function for simulations performing inference on given treatment arms using given models
#'
#' @description Analyzes given data using different models as indicated by the user. Performs inference for indicated experimental treatment arms.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` or `datasim_cont()` function.
#' @param arms Vector with treatment arms to perform inference on. Default - all arms except the first one.
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel", "timemachine").
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints.
#' @param alpha Type I error. Default=0.025.
#' @param unit_size Number of patients per calendar time unit for frequentist models adjusting for calendar time. Default=25.
#' @param ncc Boolean. Whether to include NCC data into the analysis using frequentist models. Default=TRUE.
#' @param opt Boolean. In the MAP Prior approach, if opt=1, all former periods are used as one source; if opt=2, periods form different sources get separately included into the final analysis. Default=2.
#' @param prior_prec_tau Dispersion parameter of the half normal prior, the prior for the between study heterogeneity in the MAP Prior approach. Default=4.
#' @param n.samples How many random samples will get drawn for the calculation of the posterior mean and the CIs in the MAP Prior approach. Default=1000.
#' @param n.chains Number of parallel chains for the rjags model in the MAP Prior approach. Default=4.
#' @param n.iter Number of iterations to monitor of the jags.model. Needed for coda.samples in the MAP Prior approach. Default=4000.
#' @param n.adapt Number of iterations for adaptation, an initial sampling phase during which the samplers adapt their behavior to maximize their efficiency. Needed for jags.model in the MAP Prior approach. Default=1000.
#' @param robustify Boolean.
#' @param weight Weight given to the non-informative component (0 < weight < 1) for the robustification of the MAP prior according to Schmidli (2014). Default=0.1.
#' @param ci Boolean. Whether confidence intervals for the mixed models should be computed. Default=FALSE.
#' @param prec_delta Precision of the prior regarding the treatment effect in the Time Machine approach. Default=0.001.
#' @param prec_gamma Precision of the prior regarding the control response in the Time Machine approach. Default=0.001.
#' @param tau_a Parameter \eqn{a} of the Gamma distribution regarding the precision of the drift parameter \eqn{\tau} in the Time Machine approach. I.e., \eqn{\tau \sim Gamma(a,b)}. Default=0.1.
#' @param tau_b Parameter \eqn{b} of the Gamma distribution regarding the precision of the drift parameter \eqn{\tau} in the Time Machine approach. I.e., \eqn{\tau \sim Gamma(a,b)}. Default=0.01.
#' @param prec_a Parameter \eqn{a} of the Gamma distribution regarding the precision of the responses for continuous endpoints in the Time Machine approach. I.e., \eqn{\sigma \sim Gamma(a,b)}. Default=0.001.
#' @param prec_b Parameter \eqn{b} of the Gamma distribution regarding the precision of the responses for continuous endpoints in the Time Machine approach. I.e., \eqn{\sigma \sim Gamma(a,b)}. Default=0.001.
#' @param bucket_size Number of patients per time bucket in the Time Machine approach. Default=25.
#' @param smoothing_basis String indicating the (penalized) smoothing basis to use in the GAM models. Default="tp".
#' @param basis_dim The dimension of the basis used to represent the smooth term in the GAM models. The default depends on the number of variables that the smooth is a function of. Default=-1.
#' @param gam_method The smoothing parameter estimation method for the GAM models. Default="GCV.Cp".
#' @param bs_degree Degree of the polynomial splines. Default=3.
#' @param poly_degree Degree of the discontinuous piecewise polynomials. Default=3.
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


all_models <- function(data, arms, models = c("fixmodel", "sepmodel", "poolmodel", "timemachine"), endpoint, alpha = 0.025,
                       unit_size = 250,
                       ncc = TRUE,
                       opt = 2, prior_prec_tau = 4, n.samples = 1000, n.chains = 4, n.iter = 4000, n.adapt = 1000, robustify = TRUE, weight = 0.1,
                       ci = FALSE,
                       prec_delta = 0.001, prec_gamma = 0.001, tau_a = 0.1, tau_b = 0.01, prec_a = 0.001, prec_b = 0.001, bucket_size = 25,
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
                                                                  n.samples = n.samples,
                                                                  n.chains = n.chains,
                                                                  n.iter = n.iter,
                                                                  n.adapt = n.adapt,
                                                                  robustify = robustify,
                                                                  weight = weight,
                                                                  ci = ci,
                                                                  prec_delta = prec_delta,
                                                                  prec_gamma = prec_gamma,
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
                                                                  check = FALSE))[c("reject_h0", "treat_effect")], error = function(cnd) list(reject_h0 = NA, treat_effect = NA))

      names(res_i_j) <- c(paste0("reject_h0_", j, "_", i), paste0("treat_effect_", j, "_", i))

      res <- append(res, res_i_j)
    }
  }
  return(res)
}




