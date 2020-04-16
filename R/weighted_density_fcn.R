# TODO edit wdens to avoid message:
# For infinite domains Gauss integration is applied!

#' Borrajo et al. 2017 bandwidth selection
#' @param Y A numeric vector of observations.
#' @param w A function that gives the probability of observation
#' # at any single value in the range of \code{Y}.
#' @param method One of \code{c('rt', 'brt')}.
#' \code{rt} is the rule-of-thumb estimate; \code{brt} is a bootstrap estimate
#' with the rule-of-thumb as the pilot.
#' @keywords internal

strappy <- function(Y, w, method='brt'){
  if (! method %in% c('brt','rt')){
    stop('method must be brt or rt')
  }

  # estimate of mu sub omega
  n <- length(Y)
  mu <- (sum(1 / w(Y)) / n) ^ -1

  # estimate of sigma-sq sub omega
  term1 <- sum(w(Y)) / n
  sgma <- sqrt(mu * term1 - mu^2)

  # estimate of c sub omega
  cHat <- mu * 1 / n * sum(1 / w(Y)^2)

  # estimate rule-of-thumb bandwidth
  K2u <- function(u){
    stats::dnorm(u)^2
  }
  Rk <- pracma::integral(K2u, -Inf, Inf)
  u2k <- function(u){
    u^2 * stats::dnorm(u)
  }
  mu2K <- pracma::integral(u2k, -Inf, Inf)
  # for a Gaussian kernel, mu2k = 1
  rtNum <- Rk * mu * cHat * 8 * sqrt(pi)
  rtDenom <- n * mu2K * 3
  rt <- (rtNum / rtDenom)^(1 / 5) * sgma

  if (method=='rt'){
    return(rt)
  }

  if (method=='brt'){
    # bootstrap estimate of bandwidth,
    # with rule-of-thumb for pilot

    wtsUnscld <- 1 / w(Y)
    wts <- wtsUnscld / sum(wtsUnscld)
    fgDens <- stats::density(Y, bw = rt, kernel = 'gaussian', weights = wts)
    fg <- stats::splinefun(fgDens$x, fgDens$y)

    # take 2nd deriv of smooth f estimate, square it, and integrate
    f2g_y <- fg(fgDens$x, deriv=2)
    f2g <- stats::approxfun(fgDens$x, f2g_y)
    f2g_sq <- function(u){
      f2g(u)^2
    }
    a <- min(fgDens$x)
    b <- max(fgDens$x)
    Rfg <- pracma::integral(f2g_sq, a, b)

    num <- Rk * mu * cHat
    denom <- n * mu2K * Rfg
    brt <- (num / denom)^(1/5)
    return(brt)
  }
}

#' Weight a Biased Sample to Estimate Density
#'
#' @seealso \code{\link{transdens}}
#'
#' @inheritParams transdens
#' @param bw A method to estimate the kernel bandwidth from Borrajo et al. 2017.
#' \code{rt} is the rule-of-thumb estimate; \code{brt} is a bootstrap estimate
#' with the rule-of-thumb as the pilot.
#' @export

# weighted kernel density estimation after Jones 1991
# bw argument can be rt, brt, or a numeric width to use
wdens <- function(x, w, bw='brt', reflect=FALSE, a=NULL, b=NULL, ...){
  wts <- 1/w(x)
  wts <- wts/sum(wts)

  if (is.null(a)){
    a <- min(x)
  }
  if (is.null(b)){
    b <- max(x)
  }

  # select bandwidth
  if (is.numeric(bw)){
    h <- bw
  } else {
    if (! bw %in% c('brt','rt')){
      stop('bw method must be brt or rt')
    }
    h <- strappy(x, w, method = bw)
  }

  if (reflect){
    kde <- density.reflected(x, bw = h, lower = a, upper = b, weights = wts, ...)
  } else {
    kde <- stats::density(x, bw = h, from = a, to = b, weights = wts, ...)
  }

#  f <- stats::approxfun(kde$x, kde$y)
#  lwr <- min(kde$x)
#  upr <- max(kde$x)
#  append(list(f=f, lower=lwr, upper=upr), kde)
  return(kde)
}
