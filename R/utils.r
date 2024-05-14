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
    is_even <- len %% 2 == 0
    half <- floor(len / 2)

    if (is_even) {
        x <- c(x[(half + 1):len], x[1:half])
    } else {
        x <- c(x[(half + 2):len], x[1:(half + 1)])
    }

    return(x)
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
#' @param positive_freqs Logical indicating whether to include results as
#' positive frequencies or zero-centred (default is TRUE).
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
    }

    ff <- speccy::subset_locations(
        ff, incl_boundaries, one_sided, positive_freqs
    )

    return(ff / delta)
}

#' Subset frequency locations based on specified criteria.
#'
#' This function subsets frequency locations based on the specified criteria,
#' such as whether to include boundaries, whether to consider only positive
#' frequencies, and whether to return a one-sided frequency vector.
#'
#' @param x A numeric vector with default R frequency orderings.
#' @inheritParams get_ff
#' @return A numeric vector of subsetted frequency orderings.
#' @examples
#' x <- seq(-0.5, 0.5, by = 0.1)
#' subset_locations(x)
#'
#' @export
subset_locations <- function(
    x,
    incl_boundaries = TRUE,
    one_sided = TRUE,
    positive_freqs = TRUE
) {

    n <- length(x)
    is_even <- n %% 2 == 0

    if (!one_sided && !positive_freqs) {
        x <- speccy::fftshift(x)
    } else if (one_sided) {
        x <- x[1:ceiling(n / 2)]
    }

    if (one_sided && !incl_boundaries) {
        x <- x[-1]
    } else if (!incl_boundaries) {
        if (is_even) {
            x <- x[-c(1, n / 2 + 1)]
        } else if (!is_even && positive_freqs) {
            x <- x[-1]
        } else if (!is_even && !positive_freqs) {
            x <- x[-ceiling(n / 2)]
        }
    }

    return(x)
}