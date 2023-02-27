#' Generation of a linear trend that starts in a given period
#'
#' @description This function generates a time trend for given time points in the trial according to a linear function.
#'
#' @param j Time points for which the trend should be generated.
#' @param lambda Strength of time trend.
#' @param sample_size Vector of dimension 2, indicating sample size in the trial period until the time trend starts and the remaining sample size.
#'
#' @export
#' 
#' @details The time trend is generated according to the function \eqn{f(j) = \lambda \cdot \frac{j-1}{N-1}}, where \eqn{N} is the total sample size.
#' 
#' @return Time trend for time points j.
#' @author Marta Bofill Roig, Pavla Krotka

linear_trend <- function(j, lambda, sample_size){ifelse(j<=sample_size[1], 0, lambda*(j-1)/(sum(sample_size)-1))}
