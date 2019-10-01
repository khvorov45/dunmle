# Matrix of pairwise products of x elements
# Arseniy Khvorov
# Created 2019/10/02
# Last edit 2019/10/02

library(sclr)

test_that("x_coeffs are calculated", {
  n <- 10
  x <- matrix(c(rep(1, n), rnorm(n), rnorm(n)), nrow = n, ncol = 3)
  x_coeffs <- get_x_coeffs(x)
  expect_equal(x[, 1] * x[, 1], x_coeffs[, 1])
  expect_equal(x[, 1] * x[, 2], x_coeffs[, 2])
  expect_equal(x[, 1] * x[, 3], x_coeffs[, 3])
  expect_equal(x[, 2] * x[, 2], x_coeffs[, 4])
  expect_equal(x[, 2] * x[, 3], x_coeffs[, 5])
  expect_equal(x[, 3] * x[, 3], x_coeffs[, 6])
})
