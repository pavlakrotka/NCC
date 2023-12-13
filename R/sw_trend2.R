#' Generation of stepwise trend with jump sizes adapted to sample size per period
#'
#' @description This function generates a stepwise trend for a given period. Jump sizes are adapted to the period sample size.
#'
#' @param cj Period indicator.
#' @param lambda Strength of time trend.
#' @param ss_period Sample size in a given period.
#' @param ss_total Total sample size.
#' @param trend_mean Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`).
#' @param trend_var Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`).
#'
#' @export
#' @return Time trend in period \eqn{c_j}.
#' @author Marta Bofill Roig, Pavla Krotka

sw_trend2 <- function(cj, lambda, ss_period, ss_total, trend_mean, trend_var){
  if (lambda=="random") {
    lambda <- rnorm(1, mean = trend_mean, sd = sqrt(trend_var))
  }
  (ss_period/ss_total)*as.numeric(lambda)*(cj-1)
}
