#' Wrapper function for simulations analyzing given data with all models
#'
#' @description Analyzes given data using the fixed effect model, pooled and separate analyses, and timemachine and MAP prior approach. Performs inference for all treatment arms in the trial except for the first one
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` or `datasim_cont()` function
#' @param arms Vector with treatment arms to perform inference on. Default - all arms except the first one
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAPprior")
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints
#' @param alpha Type I error. Default=0.025
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


all_models <- function(data, arms, models = c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAPprior", "mixmodel"), endpoint, alpha=0.025, ...){

  if (endpoint=="cont") {
    models <- models[models!="MAPprior"]
  }

  if (endpoint=="bin") {
    models <- models[models!="mixmodel"]
  }

  arms <- sort(arms)
  models <- sort(models)

  res <- list()

  for (i in arms) {
    for (j in models) {
      res_i_j <- list(do.call(paste0(j, "_", endpoint), list(data, arm = i, alpha))$reject_h0)

      names(res_i_j) <- paste0("reject_h0_", j, "_", i)

      res <- append(res, res_i_j)
    }
  }
  return(res)
}


