#' Compute the autocovariance function of an AR process.
#'
#' This function computes the autocovariance function (ACF) of an autoregressive
#' (AR) process with given AR coefficients \code{phis} and standard deviation
#' \code{sd}. It is an add on to the default R code that only computes the 
#' autocorrelation function.
#'
#' @param n The number of lags for which to compute the ACF.
#' @param phis A numeric vector of AR coefficients.
#' @param sd The standard deviation driving the process.
#' @return A numeric vector representing the ACF of the AR process.
#' @examples
#' ar_acf(5, c(0.5, -0.3), 1)
#'
#' @export
ar_acf <- function(n, phis, sd) {

  variance <- stats::integrate(
    speccy::ar_spectrum, -0.5, 0.5, phis, sd, delta = 1
  )$value

  variance * stats::ARMAacf(ar = phis, lag.max = n - 1)
}