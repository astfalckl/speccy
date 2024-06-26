% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fftshift}
\alias{fftshift}
\title{Shifts the zero-frequency component to the center of the spectrum.}
\usage{
fftshift(x)
}
\arguments{
\item{x}{A numeric vector.}
}
\value{
A numeric vector with the zero-frequency component centered.
}
\description{
This function takes a vector \code{x} as input and performs a circular
shift on its elements to center the zero-frequency component. If the
vector is already zero-frequency centred it will shift back to the
zero-frequency component listed first, as is default in \code{R}.
}
\examples{
x <- c(1, 2, 3, 4, 5)
fftshift(x)

}
