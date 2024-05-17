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

    is_even <- k %% 2 == 0
    width <- 1 / k

    if (is_even) {
        centres <- seq(-0.5, 0.5 - width, length = k) / delta
    } else {
        centres <- seq(-0.5 + width / 2, 0.5 - width / 2, length = k) / delta
    }

    list(centres = centres, width = width / delta)

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
  ...
) {

  if (is.null(h)) {
    h <- rep(1, n)
  }

  lag_sequence <- ifelse(is.null(lag_sequence), 1, lag_sequence)

  if (!is.null(k)) {

    basis_type <- "even"
    centres <- speccy::get_centres(k, delta)$centres
    widths <- rep(speccy::get_centres(k, delta)$width, k)

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

  tt <- 0:(n - 1)

  for (i in 1:k) {
    acf_tmp <- 2 * widths[i] * speccy::sinc(pi * tt * widths[i]) *
      cos(2 * pi * centres[i] * tt)

    bases[, i] <- speccy::bochner(
      acf_tmp, h = h, delta = delta, return_ff = FALSE, ...
    )
  }

  colnames(bases) <- seq(1, k, 1)

  bases

}