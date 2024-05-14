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
    incl_boundaries = TRUE,
    one_sided = TRUE,
    positive_freqs = TRUE
) {

    n <- base::length(ts)

    if (is.null(h)) {
        h <- rep(1, n)
    }

    h <- h / sqrt(sum(h^2))
    ts <- ts * h

    I <- Mod(stats::fft(ts))^2
    ff <- speccy::get_ff(n, delta, one_sided = FALSE)

    if (!one_sided && !positive_freqs) {
        ff[(ceiling(n / 2) + 1):n] <- ff[(ceiling(n / 2) + 1):n] - 1
        ff <- speccy::fftshift(ff)
        I <- speccy::fftshift(I)
    } else if (one_sided) {
        ff <- ff[1:ceiling(n / 2)]
        I <- I[1:ceiling(n / 2)]
    }

    # NOTE: could index this, but it's fast as is
    if (!incl_boundaries) {
        I <- I[!(ff %in% c(0, 0.5, -0.5, 1))]
        ff <- ff[!(ff %in% c(0, 0.5, -0.5, 1))]
    }

    if (return_ff) {
        return(list(ff = ff, I = I))
    } else {
        return(I)
    }
}