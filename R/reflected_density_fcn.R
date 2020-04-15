
# DO NOT EXPORT THIS FUNCTION.
# density.reflected is a function from GoFKernel, and is modified here
# to handle the case where bw arg is given as non-integer, e.g. 'nrd0'.
density.reflected <- function (x, lower = -Inf, upper = Inf,
                               weights = NULL, ...) {
  mantener <- !is.na(x)
  x <- x[mantener]
  if (upper < max(x))
    warning("There are values in the sample higher than the upper limit")
  if (lower > min(x))
    warning("There are values in the sample smaller than the lower limit")
  if (stats::sd(x) == 0) {
    dx <- stats::density(c(x, x[1] + .Machine$double.eps, x[1] -
                      .Machine$double.eps))
  }
  else {
    if (is.null(weights)) {
      pesos <- rep(1/length(x), length(x))
    }
    else {
      pesos <- weights[mantener]
    }
    argumentos <- list(...)

    # GSA revision:
    if ("bw" %in% names(argumentos) & is.numeric(argumentos$bw)) {
      #    if ("bw" %in% names(argumentos)) {
      #      if (is.numeric(argumentos$bw))
      broad <- 4 * argumentos$bw
    }
    else {
      pesos <- pesos/sum(pesos)
      broad <- 4 * stats::density(x, weights = pesos, ...)$bw
    }
    if (is.infinite(lower) & is.infinite(upper)) {
      dx <- stats::density(x, weights = pesos, ...)
    }
    else if (is.infinite(lower) & is.finite(upper)) {
      reflected <- which(x >= (upper - broad))
      x.reflect <- c(x, 2 * upper - x[reflected])
      p.reflect <- c(pesos, pesos[reflected])
      p.reflect <- p.reflect/sum(p.reflect)
      dx <- stats::density(x.reflect, weights = p.reflect, ...)
      dx$y <- (dx$y[dx$x >= lower & dx$x <= upper])
      dx$x <- (dx$x[dx$x >= lower & dx$x <= upper])
      bw <- dx$x[2] - dx$x[1]
      area.under <- sum(dx$y) * bw
      dx$y <- dx$y/area.under
    }
    else if (is.finite(lower) & is.infinite(upper)) {
      reflected <- which(x <= (lower + broad))
      x.reflect <- c(x, -x[reflected] + 2 * lower)
      p.reflect <- c(pesos, pesos[reflected])
      p.reflect <- p.reflect/sum(p.reflect)
      dx <- stats::density(x.reflect, weights = p.reflect, ...)
      dx$y <- dx$y[dx$x >= lower & dx$x <= upper]
      dx$x <- dx$x[dx$x >= lower & dx$x <= upper]
      bw <- dx$x[2] - dx$x[1]
      area.under <- sum(dx$y) * bw
      dx$y <- dx$y/area.under
    }
    else {
      reflected.inf <- which(x <= (lower + broad))
      reflected.sup <- which(x >= (upper - broad))
      x.reflect <- c(x, -x[reflected.inf] + 2 * lower)
      p.reflect <- c(pesos, pesos[reflected.inf])
      x.reflect <- c(x.reflect, 2 * upper - x[reflected.sup])
      p.reflect <- c(p.reflect, pesos[reflected.sup])
      p.reflect <- p.reflect/sum(p.reflect)
      dx <- stats::density(x.reflect, weights = p.reflect, ...)
      dx$y <- dx$y[dx$x >= lower & dx$x <= upper]
      dx$x <- dx$x[dx$x >= lower & dx$x <= upper]
      bw <- dx$x[2] - dx$x[1]
      area.under <- sum(dx$y) * bw
      dx$y <- dx$y/area.under
    }
  }
  return(dx)
}
