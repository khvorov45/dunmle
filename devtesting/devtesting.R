# Package testing script for use during development
# Arseniy Khvorov
# Created 2019/07/29
# Last edit 2019/08/13

options(scipen = 999)

library(dplyr)

# Data
dat <- read.csv("data/ExOneLSInd/data-2.csv")

# Model fit
fit <- sclr(status ~ HI + NI, dat)
summary(fit)

# Protection
prot_lvls <- data.frame(NI = log(c(0.1, 10, 40)))

get_protection_level(fit, prot_lvls, "HI")

# Prediction
preddata <- data.frame(HI = seq(0, 8, length.out = 101), NI = 0)
seq(0, 8, length.out = 101)[44]
predict(fit, preddata)[44, ]

xcoef <- get_x_coeffs(fit$x)

lapply(split(xcoef, seq(nrow(xcoef))), function(x) sum(x))

tempmat <- matrix(c(2, 3, 4, 5), nrow = 2)
tempmat2 <- matrix(c(2, 3, 4, 5), nrow = 2)
tempmat * tempmat2


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

