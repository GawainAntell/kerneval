#' Calculate Hellinger's distance
#'
#' @param d1,d2 A density distribution.
#' @param extrap Logical. If axes limits of `d1` and `d2` differ, should
#' they be extrapolated to match? The extended region will have density of 0.
#' @return A numeric value between 0 and 1 inclusive.

# Hellinger's H
hell <- function (d1, d2, extrap=TRUE) {
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
  }

  fun1 <- approxfun(d1$x, d1$y)
  fun2 <- approxfun(d2$x, d2$y)
  intgrnd <- function(x) sqrt(fun1(x) * fun2(x))
  a <- max(min(d1$x), min(d2$x))
  b <- min(max(d1$x), max(d2$x))
  int <- integral(intgrnd, a, b)
  if (int > 1){
    int <- 1
    warning('estimated Bhattacharyya distance greater than 1')
  }
  h <- sqrt(1 - int)
  h
}

#' Calculate Schoener's D metric of niche overlap
#' @inheritParams hell
#' @return A numeric value between 0 and 1 inclusive.

# Schroener's D, continuous version
schronr <- function (d1, d2) {
  if (identical(d1$x, d2$x)){
    diff <- abs(d1$y - d2$y)
    difFun <- approxfun(d1$x, diff)
  } else {
    # TODO: add option here to use only the x-axis of overlap
  #  a <- max(min(d1$x), min(d2$x))
  #  b <- min(max(d1$x), max(d2$x))

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

    fun1 <- approxfun(d1$x, d1$y)
    fun2 <- approxfun(d2$x, d2$y)
    diffFun <- function(x){
      abs(fun1(x) - fun2(x))
    }
  }

  int <- integral(diffFun, a, b)
  if (int > 2){
    int <- 2
    warning('estimated integrated overlap greater than 2')
  }
  1 - (0.5 * int)
}
