context('Overlap values')
library(kerneval)

normPdf <- function(x, shft){
  dnorm(x, (25 + shft), sd=9.7)
}
normSim <- function(n, shft){
  rnorm(n, (25 + shft), sd=9.7)
}
# With a mean=25, sd=9.7, 99% of normal distribution mass falls from 0-50
# p <- 0.99 + pnorm(0, mean=25, sd=9.7)
# qnorm(p, mean=25, sd=9.7)

xmx <- 100
n <- 100
shft <- 25

# equivalency metric (1 minus Hellinger's H)
intgrnd <- function(x) sqrt(normPdf(x, 0) * normPdf(x, shft))
int <- pracma::integral(intgrnd, 0, xmx)
# for norm vs norm, shft=10, H is .024
# versus 0.35 if integrated over -Inf to Inf
H <- sqrt(1 - int)

# Schoener's D
difFun <- function(x) abs(normPdf(x, 0) - normPdf(x, shft))
int2 <- pracma::integral(difFun, 0, xmx)
# for norm vs norm, shft=10, I is 0.69
# versus 0.61 if integrated over -Inf to Inf
D <- 1 - int2 * 0.5

x1 <- normSim(n, 0)
x2 <- normSim(n, shft)
x1 <- x1[x1 >= 0 & x1 <= xmx]
x2 <- x2[x2 >= 0 & x2 <= xmx]

# TODO build density on x1 and x2
# TODO calculate H and D from KDE

# tests -------------------------------------------------------------------

test_that('overlap is between 0 and 1', {

})

test_that('identical dists return max overlap', {

})

test_that('non-overlapping dists return min overlap', {

})

test_that('estimated overlap is within tolerance of true', {

})
