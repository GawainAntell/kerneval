#' Transform a biased sample to estimate probability density
#'
#' Calculate a kernel density estimate while correcting for selection bias
#' by transforming the data.
#'
#' \code{transdens} implements the strategy of Barmi and Simonoff (2000)
#' to correct for selection bias in kernel density estimation. The method
#' (1) transforms the empirical data based on the cumulative distribution function
#' of the bias function \code{w}, (2) scales the density so it integrates to
#' unity, and then (3) back-transforms the density to the original scale.
#'
#' Depending on the shape of the true probability distribution function and
#' the bias function, analysts would be wise to inspect kernel density plots on
#' the transformed scale, just as one might plot estimates (on the original
#' scale) when selecting a bandwidth. In particular, one should consider whether
#' the transformed distribution has a long tail or otherwise is difficult
#' to estimate. If the density estimation problem seems more straightforward
#' on the original scale, one could weight the kernel density estimate
#' with \code{wdens()} instead of transforming the data.
#'
#' @seealso \code{\link{wdens}}
#'
#' @param x A numeric vector from which the estimate is to be computed.
#' @param w A function that gives the probability of observation
#' at any single value in the range of \code{x}.
#' @param reflect Logical: should boundary reflection be applied?
#' @param a The lower limit for density estimation,
#' on the original, untransformed scale. Default is \code{min(x)}.
#' @param b The upper limit for density estimation,
#' on the original, untransformed scale. Default is \code{max(x)}.
#' @param ... Further arguments passed on to \code{\link[stats]{density}}.
#' @return An S3 density, with the additional attribute \code{$xTrans},
#' which is the location of \code{$x} values along the transformed axis.
#' @export
#' @references
#' \insertRef{Barmi00}{kerneval}

# Barmi & Simonoff 2000 method
# One should NOT provide the arguments to/from or lower/upper.
# Add tweak to uniroot function within inverse, for badly behaved data?
transdens <- function(x, w, reflect = FALSE, a = NULL, b = NULL, ...){
  if (is.null(a)){
    a <- min(x)
  }

  # transform the observations
  cdf <- function(z){
    pracma::integral(w, a, z)
  }
  Y <- sapply(x, cdf)

  # Define upper boundary of estimation on the transformed scale
  # Lower boundary will always be 0 on transformed scale (CDF from a to a)
  if (is.null(b)){
    b <- max(x)
    bYscale <- max(Y)
  } else {
    bYscale <- cdf(b)
  }

  # estimate kernel density, with boundary reflection if specified
  if (reflect){
    kde <- density.reflected(Y, lower = 0, upper = bYscale, ...)
  } else {
    kde <- stats::density(Y, from = 0, to = bYscale, ...)
  }

  # Back-transform from the y=W(x) argument to x.
  inv <- GoFKernel::inverse(cdf, lower = a, upper = b)
  kde$xTrans <- kde$x # save the estimate on the y-scale too
  kde$x <- sapply(kde$x, inv)

  # scale to integrate to one (i.e. multiply by mu-hat)
  kdeFun <- stats::approxfun(kde$x, kde$y)
  muHat <- 1/pracma::integral(kdeFun, min(kde$x), max(kde$x))
  kde$y <- muHat * kde$y

  return(kde)
}
