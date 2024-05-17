#' Compute the PSD from an ACF using Bochner's theorem.
#'
#' This function computes the PSD using Bochner's (or Wiener-Khinchin's)
#' theorem. It has an option to bias this calculation to calculate the expected
#' periodogram, as detailed in Sykulski et al. (2019). We name this function
#' Bochner primarily due to the simplicity of the name; although note it may
#' equally be \code{named wiener_khinchin()} or \code{expected_periodogram()}.
#'
#' @param acf A numeric vector representing the autocovariance function.
#' @inheritParams periodogram
#' @inheritParams lag_window
#' @param bias Logical indicating whether to intentionally bias the result to
#' calculate the expected periodogram
#' @return A numeric vector representing the PSD or expected periodogram.
#' @examples
#' acf <- c(1, 0.5, 0.3, 0.1, 0, -0.1, -0.2)
#' bochner(acf)
#'
#' @export
bochner <- function(
  acf,
  h = NULL,
  delta = 1,
  bias = TRUE,
  return_ff = TRUE,
  lag_sequence = NULL,
  ...
) {
  n <- length(acf)
  lag_sequence <- ifelse(is.null(lag_sequence), rep(1, n), lag_sequence)

  if (bias) {
    if (is.null(h)) {
      acf <- (1 - (0:(n - 1)) / n) * acf * lag_sequence
    } else {
      acf <- speccy::convolve_taper(h) * acf * lag_sequence
    }
  }

  acf <- c(acf[1] / 2, acf[2:n])

  psd <- 2 * delta * Re(stats::fft(acf))
  psd <- speccy::subset_locations(psd, ...)

  if (return_ff) {
    ff <- get_ff(n, delta = delta, ...)
    return(list(ff = ff, psd = psd))
  } else {
    return(psd)
  }

}