context('Density output')
library(kerneval)

# nVals <- c(50, 200)
# kVals <- c(2, 12)
# kVals2 <- c(3, 16)
n <- 100
k <- 12
k2 <- 16

w <- function(x){ x }
x <- rchisq(n = n, df = (k + 2))
je <- wdens(x, w)
te <- transdens(x, w)
h <- wdens(x, w, give.Rkern=TRUE)

# w2 <- function(x){ 1 / x }
# x2 <- rchisq(n = n, df = (k2 - 2))
# je2 <- wdens(x2, w2)
# te2 <- transdens(x2, w2)

# tests -------------------------------------------------------------------

test_that('output is S3 density',{
  expect_is(je, 'density')
  expect_is(te, 'density')
})
# might need to import the generic for S3 methods and export the S3 constructor

test_that('wdens can return bandwidth value', {
  expect_is(h, 'numeric')
})

test_that('density estimate integrates to unity', {
  fWeight <- stats::approxfun(je$x, je$y)
  fTrans  <- stats::approxfun(te$x, te$y)
  intWeight <- integrate(fWeight, min(je$x), max(je$x))
  intTrans  <- integrate(fTrans, min(te$x), max(te$x))
  expect_equal(intWeight$value, 1, tolerance = 0.1)
  expect_equal(intTrans$value,  1, tolerance = 0.1)
})

# TODO investigate:

# if w(x) is undefined at 0, then min and max of X can be off
# (albeit only slightly) from values at which kernel can estimate

# need to tweak to uniroot function within inverse
# for badly behaved empirical data?
