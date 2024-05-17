#' Calculate shift index for Welch's estimator
#'
#' This function shifts an index \code{l} by a proportion \code{p}.
#'
#' @param l An integer representing the index to shift.
#' @param p A numeric proportion by which to shift the index.
#' @return The shift index \code{s}.
#' @examples
#' shift_index(10, 0.5)
#'
#' @export
shift_index <- function(l, p) {
    return(round(l * (1 - p), 0))
}

#' Calculate the length of time series from Welch parameters.
#'
#' This function calculates the length of a time series specified by the
#' parameterisation of Welch's method.
#'
#' @param m The number of segments.
#' @param l The length of each segment.
#' @param p The proportion of overlap between segments (default is NULL).
#' @param s The shift value between segments (default is NULL).
#' @return The length of a Welch segment.
#' @examples
#' welch_length(5, 10, p = 0.5)
#'
#' @export
welch_length <- function(m, l, p = NULL, s = NULL) {

  if (!is.null(p)) {
    s <- speccy::shift_index(l, p)
  }

  return((m - 1) * s + l)
}

#' Compute the Welch's spectral estimate.
#'
#' This function computes Welch's spectral estimate of a time series \code{ts}
#' using the specified parameters \code{m}, \code{l}, and \code{s}.
#'
#' @inheritParams periodogram
#' @param m The number of segments to divide the time series into.
#' @param l The length of each segment.
#' @param s The overlap between segments (default is NULL).
#' @param overlap The overlap between segments as a proportion (default is
#' NULL).
#' @return A list with components \code{ff} (the frequency vector) and
#' \code{pwelch} (the Welch's periodogram estimate).
#'
#' @export
welch <- function(
  ts,
  m,
  l,
  s = NULL,
  overlap = NULL,
  delta = 1,
  h = NULL,
  return_ff = TRUE,
  ...
) {

  if (is.null(s) && is.null(overlap)) {
    stop("either s or overlap must be specified")
  }

  if ((!is.null(s)) && (!is.null(overlap))) {
    warning("both s and overlap are specified. Taking s value")
  }

  if (is.null(s) && !is.null(overlap)) {
    if (!0 <= overlap && overlap <= 1) {
      stop("overlap must be between 0 and 1")
    }
    s <- speccy::shift_index(l, overlap)
  }

  n <- speccy::welch_length(m, l, s = s)

  if (length(ts) > n) {
    stop("time-series too long for specified values of m, l, s")
  } else if (length(ts) < n) {
    warning("zero-padding ts")
    ts <- c(ts, rep(0, n - length(ts)))
  }

  if (is.null(h)) {
    h <- rep(1, l)
  }

  ff <- speccy::get_ff(l, delta, ...)
  pxx <- matrix(nrow = length(ff), ncol = m)

  for (i in 1:m){
    start <- (i - 1) * s + 1
    end <- start + l - 1
    pxx[, i] <- speccy::periodogram(
      ts[start:end], h = h, delta = delta, return_ff = FALSE, ...
    )
  }

  pwelch <- rowMeans(pxx)

  if (return_ff) {
    return(list(ff = ff, estimate = pwelch))
  } else {
    return(pwelch)
  }

}