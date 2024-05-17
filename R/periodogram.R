#' Compute the periodogram of a time series.
#'
#' This function computes the periodogram of a time series \code{ts}. The
#' periodogram is a plot of the squared magnitude of the discrete Fourier
#' transform (DFT) of the time series against the frequency.
#'
#' @param ts A numeric vector representing the time series.
#' @param h A numeric vector taper weights to apply to the time series (default
#' is NULL, equal weighting is assumed).
#' @param return_ff Logical indicating whether to return the frequency vector
#' (default is TRUE).
#' @inheritParams get_ff
#' @return A list with components \code{ff} (the frequency vector) and \code{I}
#' (the periodogram).
#' @examples
#' ts <- sin(seq(0, 10, by = 0.1))
#' periodogram(ts)
#'
#' @export
periodogram <- function(
  ts,
  h = NULL,
  delta = 1,
  return_ff = TRUE,
  ...
) {

  n <- base::length(ts)

  if (is.null(h)) {
    h <- rep(1, n)
  }

  h <- h / sqrt(sum(h^2))
  ts <- ts * h

  ii <- delta * Mod(stats::fft(ts))^2
  ff <- speccy::get_ff(n, delta, ...)
  ii <- speccy::subset_locations(ii, ...)

  if (return_ff) {
    return(list(ff = ff, estimate = ii))
  } else {
    return(ii)
  }
}