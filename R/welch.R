pwelch <- function(
  ts,
  m,
  l,
  s = NULL,
  overlap = NULL,
  delta = 1,
  h = NULL,
  return_ff = TRUE,
  incl_boundaries = TRUE,
  one_sided = TRUE,
  positive_freqs = TRUE
) {

  if (is.null(s) && is.null(overlap)) {
    stop("either s or overlap must be specified")
  }

  if ((!is.null(s)) && (!is.null(overlap))) {
    warning("both s and overlap are specified. Taking s value")
  }

  if (is.null(s) && !is.null(overlap)) {
    if (!0 <= overlap && overlap <= 1) {
      stop("overlap must be between 0 and 1")
    }
    s <- round(overlap * l, 0)
  }

  n <- (m - 1) * s + l

  if (length(ts) > n) {
    stop("time-series too long for specified values of m, l, s")
  } else if (length(ts) < n) {
    warning("zero-padding ts")
    ts <- c(ts, rep(0, n - length(ts)))
  }

  if (is.null(h)) {
    h <- rep(1, l)
  }

  ff <- speccy::get_ff(l, delta, incl_boundaries, one_sided, positive_freqs)
  pxx <- matrix(nrow = length(ff), ncol = m)

  for (i in 1:m){
    start <- (i - 1) * s + 1
    end <- start + l - 1
    pxx[, i] <- speccy::periodogram(
      ts[start:end], h = h, delta = delta, return_ff = FALSE,
      incl_boundaries = incl_boundaries, one_sided = one_sided,
      positive_freqs = positive_freqs
    )
  }

  pwelch <- rowMeans(pxx)

  if (return_ff) {
    return(list(ff = ff, pwelch = pwelch))
  } else {
    return(pwelch)
  }

}