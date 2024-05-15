#' Compute the autoregressive (AR) spectrum.
#'
#' This function computes the autoregressive (AR) spectrum given the frequencies
#' \code{ff}, the AR coefficients \code{phis}, the standard deviation \code{sd},
#' and the sampling interval \code{delta}.
#'
#' @param ff A numeric vector of frequencies.
#' @param phis A numeric vector of AR coefficients.
#' @param sd The standard deviation.
#' @inheritParams get_ff
#' @return A numeric vector representing the AR spectrum.
#' @examples
#' ff <- speccy::get_ff(1000)
#' phis <- c(0.5, -0.3)
#' sd <- 1
#' ar_spectrum(ff, phis, sd)
#'
#' @export
ar_spectrum <- function(ff, phis, sd, delta = 1) {

  d <- 1
  d_conj <- 1
  ff <- ff * delta

  for (j in seq_along(phis)) {
    d <- d - phis[j] * exp(-2i * j * pi * ff)
    d_conj <- d_conj - phis[j] * exp(2i * j * pi * ff)
  }

  return(sd^2 * delta / Re(d * d_conj))
}