#' Sinc function.
#'
#' This function computes the sinc function for a given input \code{x}. The sinc
#' function is defined as sin(x) / x, with special handling for x = 0 where the
#' value is defined as 1.
#'
#' @param x A numeric vector.
#' @param normalise Logical indicating whether to normalize \code{x} by
#' multiplying by pi (default is FALSE).
#' @return A numeric vector representing the sinc function values.
#' @examples
#' x <- seq(-10, 10, length.out = 100)
#' sinc(x)
#'
#' @export
sinc <- function(x, normalise = FALSE) {

  if (normalise) {
    x <- pi * x
  }

  value <- sin(x) / x
  value[1] <- 1
  value
}