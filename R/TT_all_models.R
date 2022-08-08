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
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' TT_all_models(data = trial_data, arms = c(1,3), endpoint = "cont")
#'
#'
#' @return List containing an indicator whether the null hypothesis was rejected or not for the investigated treatment for all models
#' @author Pavla Krotka


TT_all_models <- function(data, arms, models = c("fixmodel", "indirect"), endpoint, alpha=0.025, ...){

  models <- sort(models)

  res <- list()

  for (i in models){
    res_i <- list(try(do.call(paste0("TT_", i, "_", endpoint), list(data, arms, alpha))$reject_h0, silent = TRUE))

    names(res_i) <- paste0("reject_h0_", i)

    res <- append(res, res_i)
  }
  return(res)
}
