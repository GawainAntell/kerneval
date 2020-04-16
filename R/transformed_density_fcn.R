#' Transform a Biased Sample to Estimate Density
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
#' @export

# Barmi & Simonoff 2000 method
# One should NOT provide the arguments to/from or lower/upper.
# Add tweak to uniroot function within inverse, for badly behaved data?
transdens <- function(x, w, reflect = FALSE, a=NULL, b=NULL, ...){
  if (is.null(a)){
    test <- c(is.infinite(w(0)), is.na(w(0)))
    if (any(test)){
      a <- min(x)
    } else {
      a <- 0
    }
  }

  # transform the observations
  cdf <- function(z){
    pracma::integral(w, a, z)
  }
  Y <- sapply(x, cdf)

  # define the boundaries of estimation on the transformed scale
  aYscale <- cdf(a)
  if (is.null(b)){
    b <- max(x)
    bYscale <- max(Y)
  } else {
    bYscale <- cdf(b)
  }

  # estimate kernel density, with boundary reflection if specified
  if (reflect){
    kde <- density.reflected(Y, lower = aYscale, upper = bYscale, ...)
  } else {
    kde <- stats::density(Y, from = aYscale, to = bYscale, ...)
  }

  # scale to integrate to one
  n <- length(x)
  muHat <- n * sum( w(x)^-1 )^-1
  kde$y <- muHat * kde$y

  # Back-transform from the y=W(x) argument to x.
  inv <- GoFKernel::inverse(cdf, lower = a, upper = b)

  kde$xTrans <- kde$x
  kde$x <- sapply(kde$x, inv)

  # if w(x) is undefined at 0, then min and max of X can be off
  # (albeit only slightly) from values at which kernel can estimate
#  f <- stats::approxfun(kde$x, kde$y)
#  lwr <- min(kde$x)
#  upr <- max(kde$x)
#  append(list(f=f, lower=lwr, upper=upr), kde)

  return(kde)
}
