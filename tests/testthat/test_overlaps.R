context('Overlap values')
library(kerneval)

normPdf <- function(x, shft){
  stats::dnorm(x, (25 + shft), sd=9.7)
}
normSim <- function(n, shft){
  stats::rnorm(n, (25 + shft), sd=9.7)
}
# With a mean=25, sd=9.7, 99% of normal distribution mass falls from 0-50
# p <- 0.99 + pnorm(0, mean=25, sd=9.7)
# qnorm(p, mean=25, sd=9.7)

xmn <- 0
xmx <- 100
n <- 100
shft <- 25

# equivalency metric (1 minus Hellinger's H)
intgrnd <- function(x) sqrt(normPdf(x, 0) * normPdf(x, shft))
int <- pracma::integral(intgrnd, xmn, xmx)
# for norm vs norm, shft=10, H is .024
# versus 0.35 if integrated over -Inf to Inf
H <- sqrt(1 - int)

# Schoener's D
difFun <- function(x) abs(normPdf(x, 0) - normPdf(x, shft))
int2 <- pracma::integral(difFun, xmn, xmx)
# for norm vs norm, shft=10, I is 0.69
# versus 0.61 if integrated over -Inf to Inf
D <- 1 - int2 * 0.5

x1 <- normSim(n, 0)
x2 <- normSim(n, shft)
x1 <- x1[x1 >= xmn & x1 <= xmx]
x2 <- x2[x2 >= xmn & x2 <= xmx]

# build density on x1 and x2
d1 <- density(x1)
d2 <- density(x2)

# tests -------------------------------------------------------------------

test_that('identical dists return max overlap', {
  expect_equivalent(hell(d1, d1), 0)
  expect_equivalent(schoenr(d1, d1), 1)
})

test_that('estimated overlap is within tolerance of true', {
  expect_equal(hell(d1, d2), H,    tolerance = 0.2)
  expect_equal(schoenr(d1, d2), D, tolerance = 0.2)
})

test_that('function is commutative', {
  expect_equivalent(hell(d1, d2),
                    hell(d2, d1))
  expect_equivalent(schoenr(d1, d2),
                    schoenr(d2, d1))
})

xFar <- x1 + 2 * max(x1)
d3 <- density(xFar)

test_that('non-overlapping dists return min overlap', {
  expect_equal(hell(d1, d3, extrap = TRUE), 1,
               tolerance = 0.1)
  expect_equal(schoenr(d1, d3, extrap = TRUE), 0,
               tolerance = 0.1)
})

test_that('extrapolation works correctly', {
  expect_error(hell(d1, d3, extrap = FALSE),
               'no overlap in ranges over which densities are defined')
})

