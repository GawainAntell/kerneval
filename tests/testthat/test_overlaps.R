context('Overlap values')
library(kerneval)

n <- 500
x1 <- stats::rnorm(n, 0, 2)
x2 <- stats::rnorm(n, 2, 3)
d1 <- density(x1)
d2 <- density(x2)

# tests -------------------------------------------------------------------

test_that('identical PDFs give maximum overlap', {
  expect_equivalent(hell(d1, d1), 0)
  expect_equivalent(schoenr(d1, d1), 1)
})

test_that('function is commutative', {
  expect_equivalent(hell(d1, d2),
                    hell(d2, d1))
  expect_equivalent(schoenr(d1, d2),
                    schoenr(d2, d1))
})

test_that('bad limits give error', {
  expect_error(
    hell(d1, d2, a = 5, b = 0),
    'upper limit must be greater than lower limit'
  )
  expect_error(
    schoenr(d1, d2, a = 5, b = 0),
    'upper limit must be greater than lower limit'
  )
  expect_error(
    hell(d1, d2, a = -100, b = 5),
    'limits outside the range of available data'
  )
  expect_error(
    hell(d1, d2, a = -100, b = 5),
    'limits outside the range of available data'
  )
})

test_that('limits arguments are used', {
  h1 <- hell(d1, d2, a = -5, b = 5)
  h2 <- hell(d1, d2, a = 0,  b = 5)
  expect_false(identical(h1, h2))

  D1 <- schoenr(d1, d2, a = -5, b = 5)
  D2 <- schoenr(d1, d2, a = 0,  b = 5)
  expect_false(identical(D1, D2))
})
