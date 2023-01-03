#' Generation of stepwise trend with jump sizes adapted to sample size per period
#'
#' @description Generates stepwise trend for a given period. Jump sizes are adapted to the period sample size.
#'
#' @param cj Period indicator.
#' @param lambda Strength of time trend.
#' @param ss_period Sample size in a given period.
#' @param ss_total Total sample size.
#'
#' @export
#' @return Time trend in period \eqn{c_j}.
#' @author Marta Bofill Roig, Pavla Krotka

sw_trend2 <- function(cj,lambda,ss_period,ss_total){(ss_period/ss_total)*lambda*(cj-1)}
