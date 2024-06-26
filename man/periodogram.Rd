% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/periodogram.R
\name{periodogram}
\alias{periodogram}
\title{Compute the periodogram of a time series.}
\usage{
periodogram(ts, h = NULL, delta = 1, return_ff = TRUE, ...)
}
\arguments{
\item{ts}{A numeric vector representing the time series.}

\item{h}{A numeric vector taper weights to apply to the time series (default
is NULL, equal weighting is assumed).}

\item{delta}{The sampling interval (default is 1).}

\item{return_ff}{Logical indicating whether to return the frequency vector
(default is TRUE).}

\item{...}{Frequency subsetting options (all default to TRUE).
\code{incl_boundaries} Logical indicating whether to include the boundary
frequencies at zero and Nyquist. \code{one_sided} Logical indicating whether
to generate a one-sided frequency vector. \code{positive_freqs} Logical
indicating whether to include results as positive frequencies or
zero-centred.}
}
\value{
A list with components \code{ff} (the frequency vector) and \code{I}
(the periodogram).
}
\description{
This function computes the periodogram of a time series \code{ts}. The
periodogram is a plot of the squared magnitude of the discrete Fourier
transform (DFT) of the time series against the frequency.
}
\examples{
ts <- sin(seq(0, 10, by = 0.1))
periodogram(ts)

}
