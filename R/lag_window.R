#' Create a design window for lag window estimators.
#'
#' This function creates a design window for lag window estimation given the
#' non-negative values. The first value of \code{x} is the zeroth value, and the
#' remainder is the positive side of the window. The design window is returned
#' for a length \code{n} vector with the zero-frequency listed first.
#'
#' @param x A numeric vector representing the non-negative design kernel.
#' @param n Length of the output design window.
#' @return A normalized design window.
#' @examples
#' make_design_window(rep(1, 32), 1028)
#'
#' @export
make_design_window <- function(x, n) {

  m <- length(x) - 1

  x <- c(x[length(x):2], x)
  x <- x / sum(x)

  c(x[1:(m + 1)], rep(0, n - 2 * m - 1), x[(m + 2):(2 * m + 1)])

}

#' Compute the lag window estimate.
#'
#' This function computes the lag window estimate using either a design kernel
#' or a lag sequence.
#'
#' @inheritParams periodogram
#' @inheritParams multitaper
#' @param design_kernel A numeric vector representing the non-negative side of
#' the design kernel (default is NULL). The first value is taken as the zeroth
#' value, and the remainder is the positive side of the kernel.
#' @param lag_sequence A numeric vector representing the lag sequence (default
#' is NULL).
#' @return A list with components \code{ff} (the frequency vector) and
#' \code{estimate} (the lag window estimate).
#'
#' @export
lag_window <- function(
  ts,
  h = NULL,
  delta = 1,
  design_kernel = NULL,
  lag_sequence = NULL,
  return_ff = TRUE,
  acf = FALSE,
  incl_boundaries = TRUE,
  one_sided = TRUE,
  positive_freqs = TRUE
) {

  if (is.null(design_kernel) && is.null(lag_sequence)) {
    stop("Either design kernel or lag sequence need to be defined")
  }

  if (!is.null(design_kernel) && !is.null(lag_sequence)) {
    stop("Specify only one of design kernel or lag sequence")
  }

  n <- length(ts)

  if (is.null(h)) {
    h <- rep(1, n)
  }

  if (!is.null(design_kernel)) {

    window <- speccy::make_design_window(design_kernel, n)

    if (acf) {
      I <- speccy::bochner(ts, h, delta, TRUE, FALSE, one_sided = FALSE)
    } else {
      I <- speccy::periodogram(ts, h, delta, FALSE, one_sided = FALSE)
    }

    lw <- stats::convolve(I, window, type = "circular")
  }

  if (!is.null(lag_sequence)) {

    if (acf) {
      h_conv <- convolve_taper(h)
      acf <- ts * h_conv
    } else {
      ts <- ts * h / sqrt(sum(h^2))
      acf <- stats::convolve(ts, ts, type = "open")[n:(2 * n - 1)]
    }

    acf <- lag_sequence * acf
    acf[1] <- acf[1] / 2
    lw <- 2 * delta * Re(stats::fft(acf))
  }

  if (return_ff) {
    list(
      ff = speccy::get_ff(
        n, delta, incl_boundaries, one_sided, positive_freqs
      ),
      estimate = speccy::subset_locations(
        lw, incl_boundaries, one_sided, positive_freqs
      )
    )
  } else {
    speccy::subset_locations(
        lw, incl_boundaries, one_sided, positive_freqs
    )
  }

}