% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inv_u_trend.R
\name{inv_u_trend}
\alias{inv_u_trend}
\title{Generation of an inverted-u trend}
\usage{
inv_u_trend(j, lambda, N_peak, n_total, trend_mean, trend_var)
}
\arguments{
\item{j}{Time points for which the trend should be generated.}

\item{lambda}{Strength of time trend.}

\item{N_peak}{Point at which the time trend switches direction in terms of overall sample size.}

\item{n_total}{Total sample size in the trial.}

\item{trend_mean}{Integer. In case of random time trends, the strength of the time trend will be generated from N(\code{trend_mean}, \code{trend_var}).}

\item{trend_var}{Integer. In case of random time trends, the strength of the time trend will be generated from N(\code{trend_mean}, \code{trend_var}).}
}
\value{
Time trend for time points j.
}
\description{
This function generates a time trend for given time points in the trial according to an inverted-u function.
}
\details{
The time trend is generated according to the function

\deqn{f(j) = \lambda \cdot \frac{j-1}{N-1} \hspace{0.2cm} \mathrm{for} \hspace{0.2cm} j \leq N_p}
\deqn{f(j) = -\lambda \cdot \frac{j-N_p}{N-1} + \lambda \cdot \frac{N_p-1}{N-1} \hspace{0.2cm} \mathrm{for} \hspace{0.2cm} j > N_p}

where \eqn{N} is the total sample size (parameter \code{n_total}) and \eqn{N_p} (parameter \code{N_peak}) indicates the point at which the trend switches direction.
}
\author{
Marta Bofill Roig, Pavla Krotka
}
