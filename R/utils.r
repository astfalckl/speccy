#' Shifts the zero-frequency component to the center of the spectrum.
#'
#' This function takes a vector \code{x} as input and performs a circular
#' shift on its elements to center the zero-frequency component. If the
#' vector is already zero-frequency centred it will shift back to the
#' zero-frequency component listed first, as is default in \code{R}.
#'
#' @param x A numeric vector.
#' @return A numeric vector with the zero-frequency component centered.
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' fftshift(x)
#'
#' @export
fftshift <- function(x) {
    len <- length(x)
    half <- floor(len / 2)
    c(x[(half + 1):len], x[1:half + (len %% 2 == 0)])
}

#' Computes an open convolution of a taper with itself.
#'
#' This function computes open convolution of a taper \code{h} with itself
#' and normalises.
#'
#' @param h A numeric vector.
#' @return A numeric vector representing the tapered convolution result.
#' @examples
#' h <- c(1, 2, 3, 2, 1)
#' convolve_taper(h)
#'
#' @importFrom stats convolve
#' @export
convolve_taper <- function(h) {
    n <- length(h)
    norm <- stats::convolve(h, h, type = "filter")
    h_conv <- stats::convolve(h, h, type = "open")[n:(2 * n - 1)]
    return(h_conv / norm)
}

#' Generate a frequency vector for FFT.
#'
#' This function generates a frequency vector for FFT with given parameters.
#'
#' @param n The length of the input signal.
#' @param delta The sampling interval (default is 1).
#' @param incl_boundaries Logical indicating whether to include the boundary
#' frequencies at zero and Nyquist (default is TRUE).
#' @param one_sided Logical indicating whether to generate a one-sided frequency
#' vector (default is TRUE).
#' @param positive_freqs Logical indicating whether to include only positive
#' frequencies (default is TRUE).
#' @return A numeric vector representing the frequency vector for FFT.
#' @examples
#' get_ff(10)
#'
#' @export
get_ff <- function(
    n, delta = 1,
    incl_boundaries = TRUE,
    one_sided = TRUE,
    positive_freqs = TRUE
) {

    ff <- seq(0, n - 1) / n

    if (!one_sided && !positive_freqs) {
        ff[(ceiling(n / 2) + 1):n] <- ff[(ceiling(n / 2) + 1):n] - 1
        ff <- speccy::fftshift(ff)
    } else if (one_sided) {
        ff <- ff[1:ceiling(n / 2)]
    }

    # NOTE: could index this, but it's fast as is
    if (!incl_boundaries) {
        ff <- ff[!(ff %in% c(0, 0.5, -0.5, 1))]
    }

    return(ff / delta)
}