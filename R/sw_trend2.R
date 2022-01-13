#' Generation of stepwise trend with jump sizes adapted to sample size per period
#'
#' @description
#'
#' @param cj a
#' @param lambda a
#' @param ss_period a
#' @param ss_total a
#'
#' @export
#' @return
#' @author
#' @references

sw_trend2 <- function(cj,lambda,ss_period,ss_total){(ss_period/ss_total)*lambda*(cj-1)}