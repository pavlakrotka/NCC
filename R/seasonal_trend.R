#' Generation of a seasonal trend
#'
#' @description Generates seasonal trend for patient j according to \eqn{f(j) = \lambda \cdot \mathrm{sin} \big( \psi \cdot 2\pi \cdot \frac{j-1}{N-1} \big)}.
#'
#' @param j Patient recruitment index.
#' @param lambda Strength of time trend.
#' @param n_wave How many cycles (waves) should the time trend have (\eqn{\psi}).
#' @param n_total Total sample size in the trial.
#'
#' @export
#' @return Time trend for patient j.
#' @author Marta Bofill Roig, Pavla Krotka

seasonal_trend <- function(j, lambda, n_wave, n_total){lambda*sin(n_wave*(2*pi)*(j-1)/(n_total-1))}
