#' Wrapper function for simulations performing inference on given treatment arms using given models (treatment-treatment comparisons)
#'
#' @description Analyzes given data using different models as indicated by the user. Compares two indicated treatment arms.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` or `datasim_cont()` function.
#' @param arms Indicator of the treatment arms to be compared (vector of length 2).
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "indirect").
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints.
#' @param alpha Type I error rate. Default=0.025.
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
#' @return List containing an indicator whether the null hypothesis was rejected or not for the investigated treatment for all models.
#' @author Pavla Krotka


TT_all_models <- function(data, arms, models = c("fixmodel", "indirect"), endpoint, alpha=0.025){

  models <- sort(models)

  res <- list()

  for (i in models){
    res_i <- list(try(do.call(paste0("TT_", i, "_", endpoint), list(data = data,
                                                                    arms = arms,
                                                                    alpha = alpha))$reject_h0, silent = TRUE))

    names(res_i) <- paste0("reject_h0_", i)

    res <- append(res, res_i)
  }
  return(res)
}
