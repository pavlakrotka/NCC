#' Generation of a seasonal trend
#'
#' @description Generates seasonal trend for patient j
#'
#' @param j Patient recruitment index
#' @param lambda Strength of time trend
#' @param n_wave How many cycles (waves) should the time trend have
#' @param n_total Total sample size in the trial
#'
#' @export
#' @return Time trend for patient j
#' @author Marta Bofill Roig, Pavla Krotka

seasonal_trend <- function(j, lambda, n_wave, n_total){lambda*sin(n_wave*(2*pi)*(j-1)/(n_total-1))}
