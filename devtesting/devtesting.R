# Package testing script for use during development
# Arseniy Khvorov
# Created 2019/07/29
# Last edit 2019/08/05

options(scipen = 999)

library(dplyr)

# Data
dat <- read.csv("data/ExOneLSInd/data-1.csv")

# Model fit
fit <- sclr(status ~ HI, dat)
summary(fit)

dat$HIcens_mid <- (dat$HI_lb + dat$HI_ub) / 2
dat$HIcens_mid[near(dat$HIcens, log(5))] <- log(5)
dat$HIcens_mid[near(dat$HIcens, log(1280))] <- log(2560)
fit <- sclr(status ~ HIcens_mid, dat)
summary(fit)


# Likelihood
par_vec <- c("lambda" = 0.225, "beta_0" = -10, "beta_HI" = 3)
sclr_log_likelihood(fit, par_vec)

# Experimentation
f_dl <- function(par_vec) {
  par_vec <- as.matrix(par_vec)
  get_dl(y, x, par_vec)
}

f_db0 <- function(par_vec) {
  par_vec <- as.matrix(par_vec)
  get_db0(y, x, par_vec)
}

f_dbt <- function(par_vec) {
  par_vec <- as.matrix(par_vec)
  get_dbt(y, x, par_vec)
}

start_vec <- c("lambda" = 0.5, "beta_0" = 0, "beta_HI" = 0)

f_dl(start_vec)
f_db0(start_vec)
f_dbt(start_vec)

f_score <- function(par_vec) {
  c(
    f_dl(par_vec),
    f_db0(par_vec),
    f_dbt(par_vec)
  )
}

f_score(start_vec)

rootSolve::multiroot(f_score, start = start_vec, jactype = "bandint")

