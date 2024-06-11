#' Compute the centres of basis bins.
#'
#' This function computes the centres of \code{k} bases, each of width
#' \code{width}.
#'
#' @param k The number of bases.
#' @inheritParams periodogram
#' @return A list containing the centres of the frequency bins and the width of
#' each bin.
#' @examples
#' get_centres(10, 0.1)
#'
#' @export
get_centres <- function(k, delta = 1) {

    half_width <- 0.5 / (2 * k)
    centres <- seq(half_width, (2 * k - 1) * half_width, 2 * half_width) / delta

    return(list(centres = centres, width = 2 * half_width))

}

#' Build biased bases debiased spectral analysis
#'
#' This function biased bases for debiased spectral analysis
#' based on the specified parameters. The bases can be constructed using
#' different methods:
#' - If \code{k} is provided, it generates bases with evenly spaced centres
#'   and widths.
#' - If \code{centres} and \code{widths} are provided, it uses the specified
#'   centres and widths for the bases.
#' - If \code{lowers} and \code{uppers} are provided, it constructs bases with
#'   bounded centres and widths.
#'
#' @inheritParams get_ff
#' @inheritParams periodogram
#' @inheritParams multitaper
#' @inheritParams lag_window
#' @param k The number of bases to build with evenly spaced centres and widths
#'   (default is NULL).
#' @param centres A numeric vector representing the centres of the bases
#'   (default is NULL).
#' @param widths A numeric vector representing the widths of the bases
#'   (default is NULL).
#' @param lowers A numeric vector representing the lower bounds of the bases
#'   (default is NULL).
#' @param uppers A numeric vector representing the upper bounds of the bases
#'   (default is NULL).
#' @return A matrix of bases with dimensions \code{n} by \code{k}.
#' @examples
#' build_bases(100, k = 5)
#'
#' @export
build_bases <- function(
  n,
  delta = 1,
  h = NULL,
  k = NULL,
  centres = NULL,
  widths = NULL,
  lowers = NULL,
  uppers = NULL,
  lag_sequence = NULL,
  weights = NULL,
  ...
) {

  if (is.null(h)) {
    h <- rep(1, n)
  }
  h <- as.matrix(h)
  m <- ncol(h)

  if (is.null(weights)) {
    weights <- rep(1, m)
  }

  weights <- weights / sum(weights)
  weights <- as.matrix(weights, ncol = 1)

  if (!is.null(k)) {

    basis_type <- "even"
    centres <- speccy::get_centres(k, delta)$centres
    widths <- rep(speccy::get_centres(k, delta)$width, k)

    # ARCHIVE: used to call dwelch, think it's fixed, keeping for now.
    # centres <- dwelch::get_centres(n, k, delta)$centres
    # widths <- rep(dwelch::get_centres(n, k, delta)$width, k)

  } else if ((!is.null(centres) && !is.null(widths))) {

    basis_type <- "centres"

  } else if ((!is.null(lowers) && !is.null(uppers))) {

    basis_type <- "bounded"
    widths <- uppers - lowers
    centres <- lowers + widths / 2

  }

  cat(paste("Assuming basis type:"), basis_type)

  k <- length(centres)
  ff <- get_ff(n, delta = delta, ...)
  bases <- matrix(nrow = length(ff), ncol = k)
  multi_bases <- matrix(nrow = length(ff), ncol = m)

  tt <- 0:(n - 1)

  for (i in 1:k) {
    acf_tmp <- 2 * widths[i] * speccy::sinc(pi * tt * widths[i]) *
      cos(2 * pi * centres[i] * tt)

    for (j in 1:m) {
      multi_bases[, j] <- speccy::bochner(
        acf_tmp, h = h[, j], delta, TRUE, FALSE, lag_sequence, ...
      )
    }

    bases[, i] <- multi_bases %*% weights

  }

  colnames(bases) <- seq(1, k, 1)

  bases

}