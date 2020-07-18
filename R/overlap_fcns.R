
#' Rescale a PDF over an interval
#'
#' @param pdf A univariate probability density distribution.
#' @param a,b The new lower and upper limits of the distribution.
#' The density will be rescaled such that the integral over the
#' interval `[a,b]` sums to unity.
#' @return An S3 density object.

dscaler <- function(pdf, a, b){
  x <- pdf$x
  if (a < min(x) | b > max(x)){
    stop('limits outside the range of available data')
  }
  fhat <- stats::approxfun(pdf$x, pdf$y)

  # truncate to the new interval
  inRange <- x >= a & x <= b
  pdf$x <- pdf$x[inRange]
  pdf$y <- pdf$y[inRange]

  # Explicitly add end points at exact a and b values
  # Otherwise PDF will be truncated just shy of the given limits
  # Skip this if x values at a and b exist already
  newX <- pdf$x
  if (min(newX) != a){
    pdf$x <- c(a, pdf$x)
    pdf$y <- c(fhat(a), pdf$y)
  }
  if (max(newX) != b){
    pdf$x <- c(pdf$x, b)
    pdf$y <- c(pdf$y, fhat(b))
  }

  # rescale the density over the new interval
  mu <- 1 / pracma::integral(fhat, a, b)
  pdf$y <- pdf$y * mu
  return(pdf)
}


#' Calculate Hellinger distance
#'
#' The Hellinger distance (H) between two probability measures ranges from 0
#' (identical distributions) to 1 (non-overlapping distributions).
#'
#' Let `f(x)` and `g(x)` be the probability density functions for comparison.
#' The squared Hellinger distance is
#' \eqn{0.5 \int(\sqrt f(x) - \sqrt g(x) ^2) dx}. Conveniently for numeric
#' integration, the expression can be written as the integral of a product
#' instead of a difference: \eqn{1 - \int \sqrt{f(x)g(x)} dx}. H is related
#' to the Bhattacharyya coefficient BC: \eqn{H = \sqrt{1 - BC}}. Elsewhere,
#' Hellinger distances are sometimes reported as the square (\eqn{H^2})
#' [@Warren08; DiCola17], or are not rescaled and so range
#' from 0 to \sqrt(2) [@Nikulin01].
#'
#' @param d1,d2 A density distribution.
#' @param a,b The lower and upper limits of comparison.
#' Each density will be rescaled such that the integral over the
#' interval `[a,b]` sums to unity. If empty, the minimum and maximum
#' of the intersection of `d1$x` and `d2$x` will be used.
#' @return A numeric value between 0 and 1 inclusive.
#' @export

hell <- function (d1, d2, a = NULL, b = NULL) {
  if (min(d1$x) > max(d2$x) | max(d1$x) < min(d2$x)){
    stop('no overlap in ranges over which densities are defined')
  }

  # standardise PDFs over the same x-interval
  if (is.null(a)){
    a <- max(min(d1$x), min(d2$x))
  }
  if (is.null(b)){
    b <- min(max(d1$x), max(d2$x))
  }
  d1scl <- dscaler(d1, a, b)
  d2scl <- dscaler(d2, a, b)

  fun1 <- stats::approxfun(d1scl$x, d1scl$y)
  fun2 <- stats::approxfun(d2scl$x, d2scl$y)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))
  int <- pracma::integral(intgrnd, a, b)
  if (int > 1){
    int <- 1
    # warning('estimated Bhattacharyya distance greater than 1')
  }
  h <- sqrt(1 - int)
  h
}


#' Calculate Schoener's D metric of niche overlap
#'
#' Schoener's D metric quantifies niche overlap between two discretised
#' probability functions. The \code{schoenr} function modifies the metric for
#' continuous probability distributions.
#'
#' D was originally defned as the sum of the absolute difference in resource use
#' frequencies for every category i. The sum is rescaled by 0.5 to give the
#' per-niche value, and then subtracted from 1 so as to be a similarity metric
#' (i.e. D = 1 indicates identical niches). \code{schoenr} replaces the sum
#' with an integral on the absolute difference between two estimated
#' kernel density functions.
#'
#' @inheritParams hell
#' @return A numeric value between 0 and 1 inclusive.
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Schoener68}{kerneval}

schoenr <- function (d1, d2, a = NULL, b = NULL) {
  if (min(d1$x) > max(d2$x) | max(d1$x) < min(d2$x)){
    stop('no overlap in ranges over which densities are defined')
  }

  # standardise PDFs over the same x-interval
  if (is.null(a)){
    a <- max(min(d1$x), min(d2$x))
  }
  if (is.null(b)){
    b <- min(max(d1$x), max(d2$x))
  }
  d1scl <- dscaler(d1, a, b)
  d2scl <- dscaler(d2, a, b)

  fun1 <- stats::approxfun(d1scl$x, d1scl$y)
  fun2 <- stats::approxfun(d2scl$x, d2scl$y)
  difFun <- function(x){
    abs(fun1(x) - fun2(x))
  }
  int <- pracma::integral(difFun, a, b)
  if (int > 2){
    int <- 2
#    warning('estimated integrated overlap greater than 2')
  }
  1 - (0.5 * int)
}
