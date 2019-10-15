# Simulating sample data
# Arseniy Khvorov
# Created 2019/08/30
# Last edit 2019/09/13

library(dplyr)

sim_dat <- function(nsam, lambda, beta_0, beta_HI, beta_NI, simseed) {
  set.seed(simseed)
  logHI <- rnorm(nsam, 2, 2)
  logNI <- rnorm(nsam, 2, 2)
  status <- rbinom(
    nsam, 
    1, 
    lambda / (1 + exp(beta_0 + beta_HI * logHI + beta_NI * logNI))
  )
  logHIcens <- case_when(
    logHI < log(10) ~ log(5),
    logHI < log(20) ~ log(10),
    logHI < log(40) ~ log(20),
    logHI < log(80) ~ log(40),
    logHI < log(160) ~ log(80),
    logHI < log(320) ~ log(160),
    logHI < log(640) ~ log(320),
    logHI < log(1280) ~ log(640),
    TRUE ~ log(1280)
  )
  logHImid <- case_when(
    near(logHIcens, log(5)) ~ logHIcens,
    near(logHIcens, log(1280)) ~ logHIcens,
    TRUE ~ logHIcens + (log(2) / 2)
  )
  logNIcens <- if_else(logNI < log(5), log(2.5), logNI)
  pop <- tibble(status, logHI, logHIcens, logHImid, logNI, logNIcens)
  if (beta_NI == 0) pop <- select(pop, -logNI, -logNIcens)
  return(pop)
}

nsam <- 5000

one_titre_data <- sim_dat(nsam, 0.5, -5, 2, 0, 20190913)
two_titre_data <- sim_dat(nsam, 0.5, -7.5, 2, 2, 20190914)

usethis::use_data(one_titre_data, two_titre_data, overwrite = TRUE)
