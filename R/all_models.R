#' Wrapper function for simulations analyzing given data with all models
#'
#' @description Analyzes given data using the fixed effect model, pooled and separate analyses, and timemachine and MAP prior approach. Performs inference for all treatment arms in the trial except for the first one
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` or `datasim_cont()` function
#' @param arms Vector with treatment arms to perform inference on. Default - all arms except the first one
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAP_rjags")
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints
#' @param alpha Type I error. Default=0.025
#' @param ... Further arguments for simulation function
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
#' @return List containing an indicator whether the null hypothesis was rejected or not for the investigated treatment for all models
#' @author Pavla Krotka


all_models <- function(data, arms, models = c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAP_rjags", "mixmodel"), endpoint, alpha = 0.025,
                       unit_size = 250,
                       ncc = TRUE,
                       opt = 1, prior_prec_tau = 1, n.samples = 1000, n.chains = 4, n.iter = 4000, n.adapt = 1000, robustify = TRUE, weight = 0.1,
                       ci = FALSE,
                       prec_delta = 0.001, prec_gamma = 0.001, tau_a = 0.1, tau_b = 0.01, prec_a = 0.001, prec_b = 0.001, bucket_size = 25,
                       smoothing_basis = "tp", basis_dim = -1, gam_method = "GCV.Cp",
                       bs_degree = 3, poly_degree = 3, ...){

  if (endpoint=="bin") {
    models <- models[models!="mixmodel"]
  }

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




