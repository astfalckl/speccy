#' Compute the PSD from an ACF using Bochner's theorem.
#'
#' This function computes the PSD using Bochner's (or Wiener-Khinchin's)
#' theorem. 
#' , which involves modifying the autocovariance function (ACF) and then taking the Fourier transform of the modified ACF.
#'
#' @param acf A numeric vector representing the autocovariance function.
#' @param bias 
#' @param full Logical indicating whether to return the full periodogram or only the positive frequencies (default is FALSE).
#' @param zero Logical indicating whether to return the zero-frequency component (default is FALSE).
#' @return A numeric vector representing the periodogram.
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
    incl_boundaries = TRUE,
    one_sided = TRUE,
    positive_freqs = TRUE
) {
    n <- length(acf)
    ff <- get_ff(n, delta, incl_boundaries, one_sided, positive_freqs)
    # nfreq <- dwelch::get_nfreq(n)

    if (bias) {
        if (is.null(h)) {
            acf <- (1 - (0:(n - 1)) / n) * acf
        } else {
            acf <- speccy::convolve_taper(h) * acf
        }
    }

    acf <- c(acf[1] / 2, acf[2:n])

    if (full) {
        dwelch::fftshift(2 * delta * Re(stats::fft(acf)))
    } else if (zero) {
        2 * delta * Re(stats::fft(acf))[1:(nfreq + 1)]
    }else {
        2 * delta * Re(stats::fft(acf))[2:(nfreq + 1)]
    }

}



# bochner <- function(acf, delta = 1, h = NULL, full = FALSE, zero = FALSE) {
#   n <- length(acf)
#   nfreq <- dwelch::get_nfreq(n)

#   if (is.null(h)) {
#     acf <- (1 - (0:(n - 1)) / n) * acf
#   } else {
#     h <- h / sqrt(sum(h^2)) # normalise h
#     h_conv <- stats::convolve(h, h, type = "open")[n:(2 * n - 1)]
#     acf <- h_conv * acf
#   }

#   is_even <- n %% 2 == 0

#   if (is_even) {
#     acf <- c(acf[1] / 2, acf[2:(n - 1)], acf[n] / 2)
#   } else {
#     acf <- c(acf[1] / 2, acf[2:n])
#   }

#   if (full) {
#     dwelch::fftshift(2 * delta * Re(stats::fft(acf)))
#   } else if (zero) {
#     2 * delta * Re(stats::fft(acf))[1:(nfreq + 1)]
#   } else {
#     2 * delta * Re(stats::fft(acf))[2:(nfreq + 1)]
#   }
# }
