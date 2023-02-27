#' Generation of an inverted-u trend
#'
#' @description This function generates a time trend for given time points in the trial according to an inverted-u function.
#'
#' @param j Time points for which the trend should be generated.
#' @param lambda Strength of time trend.
#' @param N_peak Point at which the time trend switches direction in terms of overall sample size.
#' @param n_total Total sample size in the trial.
#'
#' @export
#' 
#' @details The time trend is generated according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1} (I(j \leq N_p) - I(j > N_p))}, 
#' where \eqn{N} is the total sample size (parameter `n_total`) and \eqn{N_p} (parameter `N_peak`) indicates the point at which the trend switches direction.
#' 
#' @return Time trend for time points j.
#' @author Marta Bofill Roig, Pavla Krotka

inv_u_trend <- function(j, lambda, N_peak, n_total){ifelse(j<=N_peak, lambda*(j-1)/(n_total-1), (-lambda*(j-N_peak)/(n_total-1))+(lambda*(N_peak-1)/(n_total-1)))}
