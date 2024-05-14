#' Compute the multitaper spectral estimate.
#'
#' This function computes the multitaper spectral estimate of a time series \code{ts}
#' using the given set of tapers. The multitaper method uses multiple tapers to
#' reduce variance in the spectral estimate.
#'
#' @param ts A numeric vector representing the time series.
#' @param tapers A matrix of tapers (each column is a taper).
#' @param delta The spacing between frequency bins (default is 1).
#' @param return_ff Logical indicating whether to return the frequency vector (default is TRUE).
#' @param acf Logical indicating whether to compute the autocovariance function (default is FALSE).
#' @return A list with components \code{ff} (the frequency vector) and \code{mt} (the multitaper spectral estimate).
#' @examples
#' ts <- sin(seq(0, 10, by = 0.1))
#' tapers <- matrix(data = rnorm(100), nrow = 100, ncol = 5)
#' multitaper(ts, tapers)
#'
#' @importFrom dwelch bochner
#' @export
multitaper <- function(
    ts, tapers, delta = 1, return_ff = TRUE, acf = FALSE

) {

    m <- ncol(tapers)
    n <- nrow(tapers)
    periodograms <- matrix(nrow = base::ceiling(n / 2), ncol = m)

    for (i in 1:m) {
        if (acf) {
            periodograms[, i] <- dwelch::bochner(
                ts, h = tapers[, i], zero = TRUE
            )
        } else {
            periodograms[, i] <- speccy::periodogram(
                ts, h = tapers[, i], return_ff = FALSE
            )
        }
    }

    mt <- rowMeans(periodograms)

    if (return_ff) {
        return(list(ff = get_ff(n, delta), mt = mt))
    } else {
        return(mt)
    }

}