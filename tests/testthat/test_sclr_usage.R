# Tests of sclr returning/throwing errors properly
# Arseniy Khvorov
# Created 2019/09/02
# Last edit 2019/09/13

library(sclr)

test_that(
  "Error with missing parameters", {
    expect_error(sclr(status ~ logHI))
    expect_error(sclr(one_titre_data))
  }
)

test_that(
  "Error with unexpected outcome", {
    # Factor
    dat <- one_titre_data
    dat$status <- as.factor(dat$status)
    expect_error(sclr(status ~ logHI, dat))
    # 1's and 2's
    dat$status <- as.numeric(dat$status)
    expect_error(sclr(status ~ logHI, dat))
  }
)
