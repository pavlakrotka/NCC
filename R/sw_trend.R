#' Generation of stepwise trend with equal jumps between periods
#'
#' @description Generates stepwise trend for a given period according to \eqn{f(j) = \lambda \cdot (c_j - 1)}. No time trend is assumed in the first period.
#'
#' @param cj Period indicator.
#' @param lambda Strength of time trend.
#'
#' @export
#' @return Time trend in period \eqn{c_j}.
#' @author Marta Bofill Roig, Pavla Krotka

sw_trend <- function(cj,lambda){lambda*(cj-1)}
