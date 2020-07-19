context('Density output')
library(kerneval)

w <- function(x){ x }
x <- rchisq(n = 100, df = 14)
je <- wdens(x, w)
te <- transdens(x, w)

xnorm <- stats::rnorm(100, sd = 10)
kdeNorm <- stats::density(xnorm)
kdeScl <- dscaler(kdeNorm, a = -10, b = 10)

# tests -------------------------------------------------------------------

test_that('output is S3 density',{
  expect_s3_class(je, 'density')
  expect_s3_class(te, 'density')
})

test_that('cropped density endpoints are correct and unique', {
  xTrunc <- kdeScl$x
  lengthAll <- length(xTrunc)
  lengthUniq <- length(unique(xTrunc))
  expect_identical(lengthAll, lengthUniq)
  expect_identical(min(xTrunc), -10)
  expect_identical(max(xTrunc),  10)
})

test_that('density estimate integrates to unity', {
  fWeight <- stats::approxfun(je$x, je$y)
  fTrans  <- stats::approxfun(te$x, te$y)
  fNorm <-   stats::approxfun(kdeScl$x, kdeScl$y)
  intWeight <- stats::integrate(fWeight, min(je$x), max(je$x))
  intTrans  <- stats::integrate(fTrans,  min(te$x), max(te$x))
  intNorm   <- stats::integrate(fNorm, -10, 10)
  expect_equal(intWeight$value, 1, tolerance = 0.1)
  expect_equal(intTrans$value,  1, tolerance = 0.1)
  expect_equal(intNorm$value,   1, tolerance = 0.1)
})

test_that('dscaler throws errors with improper bounds', {
  expect_error(
    dscaler(kdeNorm, a = -100, b = 100),
    'limits outside the range of available data'
  )
  expect_error(
    dscaler(kdeNorm, a = 10, b = -10),
    'upper limit must be greater than lower limit'
  )
})

test_that('wdens additional arguments work', {
  h <- wdens(x, w, give.Rkern=TRUE)
  expect_is(h, 'numeric')

  expect_error(
    wdens(x, w, kernel = 'epanechnikov')
  )
})
