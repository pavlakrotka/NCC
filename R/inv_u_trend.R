#' Generation of an inverted-u trend
#'
#' @description Generates inverted-u trend for patient j according to \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1} (I(j \leq N_p) - I(j > N_p))}.
#'
#' @param j Patient recruitment index.
#' @param lambda Strength of time trend.
#' @param N_peak Point at which the time trend switches direction in terms of overall sample size.
#' @param sample_size Total sample size in the trial.
#'
#' @export
#' @return Time trend for patient j.
#' @author Marta Bofill Roig, Pavla Krotka

inv_u_trend <- function(j, lambda, N_peak, sample_size){ifelse(j<=N_peak, lambda*(j-1)/(sample_size-1), (-lambda*(j-N_peak)/(sample_size-1))+(lambda*(N_peak-1)/(sample_size-1)))}
