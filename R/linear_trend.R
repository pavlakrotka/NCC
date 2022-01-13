#' Generation of a linear trend
#'
#' @description
#'
#' @param j a
#' @param lambda a
#' @param sample_size a
#'
#' @export
#' @return
#' @author
#' @references

linear_trend <- function(j, lambda, sample_size){lambda*(j-1)/(sample_size-1)}