#' Generation of stepwise trend with equal jumps between periods
#'
#' @description This function generates a stepwise trend for a given period. No time trend is assumed in the first period.
#'
#' @param cj Period indicator.
#' @param lambda Strength of time trend.
#'
#' @export
#' 
#' @details The time trend is generated according to the function \eqn{f(j) = \lambda \cdot (c_j - 1)}, where \eqn{c_j} is an index of the period patient \eqn{j} was enrolled in.
#' 
#' @return Time trend in period \eqn{c_j}.
#' @author Marta Bofill Roig, Pavla Krotka

sw_trend <- function(cj,lambda){lambda*(cj-1)}
