
#' Borrajo et al. 2017 bandwidth selection
#' @param Y A numeric vector of observations.
#' @param w A function that gives the probability of observation
#' # at any single value in the range of \code{Y}.
#' @param method One of \code{c('rt', 'brt')}.
#' \code{rt} is the rule-of-thumb estimate; \code{brt} is a bootstrap estimate
#' with the rule-of-thumb as the pilot.
#' @keywords internal

selectbw <- function(Y, w, method='brt'){
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
  Rk <- pracma::quadinf(K2u, -Inf, Inf)$Q
  u2k <- function(u){
    u^2 * stats::dnorm(u)
  }
  mu2K <- pracma::quadinf(u2k, -Inf, Inf)$Q
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

#' Weight a biased sample to estimate probability density
#'
#' Calculate a kernel density estimate while correcting for selection bias
#' by weighting the kernels.
#'
#' The kernel on each datum is weighted by the inverse of the observation
#' probability at that point, \code{1/w(X)}. Weighted kernel estimation should
#' not be confused with adaptive kenel estimation. Both approaches
#' modify the individual kernels that contribute to an estimate.
#' However, in adaptive KDE the badwidth of each kernel is adjusted,
#' whereas in weighted KDE the bandwidth is constant while the height
#' (total probability density) of each kernel is adjusted.
#'
#' Weighted KDE is the method data that has received the most attention
#' in the statistics literature for selection-biased data, beginning with
#' the foundational paper of Jones (1991). The method is a slight modification
#' of classical KDE and does not add onerous calculation, making it attractive.
#' However, the choice of bandwidth is complicated. \code{wdens} calls the
#' internal \code{selectbw} function to select a bandwidth following the
#' user-specified bootstrap methods of Borrajo and others (2017).
#'
#' @seealso \code{\link{transdens}}
#'
#' @inheritParams transdens
#' @param bw A method to estimate the kernel bandwidth from Borrajo et al. 2017.
#' \code{rt} is the rule-of-thumb estimate; \code{brt} is a bootstrap estimate
#' with the rule-of-thumb as the pilot. Alternatively, a numeric value to pass
#' to \code{\link[stats]{density}}.
#' @return An S3 density.
#' @export
#' @references
#' \insertRef{Jones91}{kerneval}
#'
#' \insertRef{Borrajo17}{kerneval}

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
    h <- selectbw(x, w, method = bw)
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
