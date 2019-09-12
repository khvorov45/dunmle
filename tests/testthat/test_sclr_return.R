# Tests of sclr returning the expected values.
# 
# These 'correct' values are coming from the verion of sclr that was
# verified to be good through simulations.
#
# Arseniy Khvorov
# Created 2019/09/02
# Last edit 2019/09/13

library(sclr)

test_that("Returns the expected parameter names", {
  fit_my_names <- sclr(status ~ logHI, one_titre_data)
  expect_named(fit_my_names$parameters, c("lambda", "beta_0", "beta_logHI"))
  est_conv_names <- sclr(
    status ~ logHI, one_titre_data, conventional_names = TRUE
  )
  expect_named(
    est_conv_names$parameters, 
    c("(Baseline)", "(Intercept)", "logHI")
  )
})

test_that("Return is stable", {
  pars <- do.call(c, lapply(1:10, function(ind) {
    fit <- sclr(status ~ logHI, one_titre_data)
    par <- fit$parameters[["beta_logHI"]]
    return(par)
  }))
  expect_equal(length(unique(round(pars, 4))), 1)
})
