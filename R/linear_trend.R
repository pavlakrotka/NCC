#' Generation of a linear trend that starts in a given period
#'
#' @description Generates linear trend for patient j
#'
#' @param j Patient recruitment index
#' @param lambda Strength of time trend
#' @param sample_size Vector of dimension 2, indicating sample size in the trial period until the time trend starts and the remaining sample size.
#' 
#' @export
#' @return Time trend for patient j
#' @author Marta Bofill Roig, Pavla Krotka

linear_trend <- function(j, lambda, sample_size){ifelse(j<=sample_size[1], 0, lambda*(j-1)/(sum(sample_size)-1))}