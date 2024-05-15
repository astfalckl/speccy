#' Compute the multitaper spectral estimate.
#'
#' This function computes the multitaper spectral estimate of a time series
#' \code{ts} using the given set of tapers.
#'
#' @param ts A numeric vector representing the time series.
#' @param tapers A matrix of tapers (each column is a taper).
#' @param acf Logical indicating whether to treat ts as the ACF. Setting to TRUE
#' returns the expected value of the multitaper estimate (default is FALSE).
#' @inheritParams periodogram
#' @return A list with components \code{ff} (the frequency vector) and \code{mt}
#' (the multitaper spectral estimate).
#'
#' @export
multitaper <- function(
  ts,
  tapers,
  delta = 1,
  acf = FALSE,
  return_ff = TRUE,
  incl_boundaries = TRUE,
  one_sided = TRUE,
  positive_freqs = TRUE
) {

  m <- ncol(tapers)
  n <- nrow(tapers)

  ff <- speccy::get_ff(n, delta, incl_boundaries, one_sided, positive_freqs)

  periodograms <- matrix(nrow = length(ff), ncol = m)

  for (i in 1:m) {
    if (acf) {
      periodograms[, i] <- speccy::bochner(
        ts, tapers[, i], delta, TRUE, FALSE,
        incl_boundaries, one_sided, positive_freqs
      )
    } else {
      periodograms[, i] <- speccy::periodogram(
        ts, tapers[, i], delta, FALSE,
        incl_boundaries, one_sided, positive_freqs
      )
    }
  }

  mt <- rowMeans(periodograms)

  if (return_ff) {
    return(list(ff = ff, estimate = mt))
  } else {
    return(mt)
  }

}