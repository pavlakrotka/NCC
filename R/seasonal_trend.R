#' Generation of a seasonal trend
#'
#' @description This function generates a time trend for given time points in the trial according to a periodic function.
#'
#' @param j Time points for which the trend should be generated.
#' @param lambda Strength of time trend.
#' @param n_wave How many cycles (waves) should the time trend have (\eqn{\psi}).
#' @param n_total Total sample size in the trial.
#'
#' @export
#'
#' @details The time trend is generated according to the function \eqn{f(j) = \lambda \cdot \mathrm{sin} \big( \psi \cdot 2\pi \cdot \frac{j-1}{N-1} \big)}, where \eqn{N} is the total sample size  (parameter `n_total`) and the parameter \eqn{\psi} corresponds to the input parameter `n_wave`.
#'
#' @return Time trend for time points j.
#' @author Marta Bofill Roig, Pavla Krotka

seasonal_trend <- function(j, lambda, n_wave, n_total){lambda*sin(n_wave*(2*pi)*(j-1)/(n_total-1))}
