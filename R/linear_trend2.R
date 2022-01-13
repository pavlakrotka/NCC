#' Generation of a linear trend that starts in the second period
#'
#' @description
#'
#' @param j a
#' @param lambda a
#' @param sample_size vector of dimension 2, indicating sample size in the first period and the remaining sample size
#'
#' @export
#' @return
#' @author
#' @references

linear_trend2 <- function(j, lambda, sample_size){ifelse(j<=sample_size[1],0,lambda*(j-1)/(sum(sample_size)-1))}