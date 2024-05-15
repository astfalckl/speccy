#' Generate a Bartlett sequence.
#'
#' This function generates Bartlett's sequence of length \code{n} with parameter
#' \code{l}.
#'
#' @param n The length of the Bartlett window.
#' @param l The parameter for the Bartlett window.
#' @return A Bartlett window of length \code{n}.
#' @examples
#' w_bartlett(10, 5)
#'
#' @export
w_bartlett <- function(n, l) {
  tau <- 0:(n - 1)
  pmax(1 - tau / (2 * l), rep(0, n))
}

#' Generate a Daniell sequence.
#'
#' This function generates Daniell's sequence of length \code{n} with parameter
#' \code{l}.
#'
#' @param n The length of the Daniell sequence.
#' @param l The parameter for the Daniell sequence.
#' @return A Daniell window of length \code{n}.
#' @examples
#' w_daniell(10, 5)
#'
#' @export
w_daniell <- function(n, l) {
  tau <- 0:(n - 1)
  w <- sin(2 * pi * tau / l) * l / (2 * pi * tau)
  w[1] <- 1
  return(w)
}