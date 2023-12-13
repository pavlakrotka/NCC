#' Generation of stepwise trend with equal jumps between periods
#'
#' @description This function generates a stepwise trend for a given period. No time trend is assumed in the first period.
#'
#' @param cj Period indicator.
#' @param lambda Strength of time trend.
#' @param trend_mean Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`).
#' @param trend_var Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`).
#'
#' @export
#'
#' @details The time trend is generated according to the function \eqn{f(j) = \lambda \cdot (c_j - 1)}, where \eqn{c_j} is an index of the period patient \eqn{j} was enrolled in.
#'
#' @return Time trend in period \eqn{c_j}.
#' @author Marta Bofill Roig, Pavla Krotka

sw_trend <- function(cj, lambda, trend_mean, trend_var){
  if (lambda=="random") {
    lambda <- rnorm(1, mean = trend_mean, sd = sqrt(trend_var))
  }
  as.numeric(lambda)*(cj-1)
}
