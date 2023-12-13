#' Generation of an inverted-u trend
#'
#' @description This function generates a time trend for given time points in the trial according to an inverted-u function.
#'
#' @param j Time points for which the trend should be generated.
#' @param lambda Strength of time trend.
#' @param N_peak Point at which the time trend switches direction in terms of overall sample size.
#' @param n_total Total sample size in the trial.
#' @param trend_mean Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`).
#' @param trend_var Integer. In case of random time trends, the strength of the time trend will be generated from N(`trend_mean`, `trend_var`).
#'
#' @export
#'
#' @details The time trend is generated according to the function
#'
#' \deqn{f(j) = \lambda \cdot \frac{j-1}{N-1} \hspace{0.2cm} \mathrm{for} \hspace{0.2cm} j \leq N_p}
#' \deqn{f(j) = -\lambda \cdot \frac{j-N_p}{N-1} + \lambda \cdot \frac{N_p-1}{N-1} \hspace{0.2cm} \mathrm{for} \hspace{0.2cm} j > N_p}
#'
#' where \eqn{N} is the total sample size (parameter `n_total`) and \eqn{N_p} (parameter `N_peak`) indicates the point at which the trend switches direction.
#'
#' @return Time trend for time points j.
#' @author Marta Bofill Roig, Pavla Krotka

inv_u_trend <- function(j, lambda, N_peak, n_total, trend_mean, trend_var){
  if (lambda=="random") {
    lambda <- rnorm(1, mean = trend_mean, sd = sqrt(trend_var))
  }
  ifelse(j<=N_peak, as.numeric(lambda)*(j-1)/(n_total-1), (-as.numeric(lambda)*(j-N_peak)/(n_total-1))+(as.numeric(lambda)*(N_peak-1)/(n_total-1)))
}
