
#' Calculate Hellinger's distance
#'
#' @param d1,d2 A density distribution.
#' @param extrap Logical. If axes limits of \code{d1} and \code{d2} differ,
#' should they be extrapolated to match? The density function
#' will estimate 0 probability in the region of extrapolation.
#' @return A numeric value between 0 and 1 inclusive.
#' @export

# Hellinger's H
hell <- function (d1, d2, extrap = TRUE) {
  if (extrap==TRUE){
    a <- min(d1$x, d2$x)
    b <- max(d1$x, d2$x)

    # ensure the density estimate extends along the full x-axis
    if (min(d1$x) > a){
      d1$x <- c(a, d1$x)
      d1$y <- c(0, d1$y)
    }
    if (min(d2$x) > a){
      d2$x <- c(a, d2$x)
      d2$y <- c(0, d2$y)
    }
    if (max(d1$x) < b){
      d1$x <- c(d1$x, b)
      d1$y <- c(d1$y, 0)
    }
    if (max(d2$x) < b){
      d2$x <- c(d2$x, b)
      d2$y <- c(d2$y, 0)
    }
  } else {
    if (min(d1$x) > max(d2$x) | max(d1$x) < min(d2$x)){
      stop('no overlap in ranges over which densities are defined')
    }
  }

  fun1 <- stats::approxfun(d1$x, d1$y)
  fun2 <- stats::approxfun(d2$x, d2$y)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))

  # re-define a and b after extrapolation
  # or define in the first place if no extrapolation
  a <- max(min(d1$x), min(d2$x))
  b <- min(max(d1$x), max(d2$x))
  int <- pracma::integral(intgrnd, a, b)
  if (int > 1){
    int <- 1
#    warning('estimated Bhattacharyya distance greater than 1')
  }
  h <- sqrt(1 - int)
  h
}

#' Calculate Schoener's D metric of niche overlap
#'
#' This is a continuous version of the method introduced by Schoener 1968.
#'
#' @inheritParams hell
#' @return A numeric value between 0 and 1 inclusive.
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Schoener68}{kerneval}

# Schoener's D, continuous version
schoenr <- function (d1, d2, extrap = TRUE) {

  if (extrap==TRUE){
    a <- min(d1$x, d2$x)
    b <- max(d1$x, d2$x)

    # ensure the density estimate extends along the full x-axis
    if (min(d1$x) > a){
      d1$x <- c(a, d1$x)
      d1$y <- c(0, d1$y)
    }
    if (min(d2$x) > a){
      d2$x <- c(a, d2$x)
      d2$y <- c(0, d2$y)
    }
    if (max(d1$x) < b){
      d1$x <- c(d1$x, b)
      d1$y <- c(d1$y, 0)
    }
    if (max(d2$x) < b){
      d2$x <- c(d2$x, b)
      d2$y <- c(d2$y, 0)
    }
  } else {
    if (min(d1$x) > max(d2$x) | max(d1$x) < min(d2$x)){
      stop('no overlap in ranges over which densities are defined')
    }
  }

  fun1 <- stats::approxfun(d1$x, d1$y)
  fun2 <- stats::approxfun(d2$x, d2$y)
  difFun <- function(x){
    abs(fun1(x) - fun2(x))
  }

  # re-define a and b after extrapolation
  # or define in the first place if no extrapolation
  a <- max(min(d1$x), min(d2$x))
  b <- min(max(d1$x), max(d2$x))
  int <- pracma::integral(difFun, a, b)
  if (int > 2){
    int <- 2
#    warning('estimated integrated overlap greater than 2')
  }
  1 - (0.5 * int)
}
