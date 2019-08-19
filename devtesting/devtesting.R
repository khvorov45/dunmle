# Package testing script for use during development
# Arseniy Khvorov
# Created 2019/07/29
# Last edit 2019/08/14

devtools::install()

options(scipen = 999)

fit1 <- sclr(status ~ HI, sclronetitre)
fit2 <- sclr(status ~ HI + NI, sclrtwotitre)

fit1cn <- sclr(status ~ HIcens, sclronetitre)

summary(fit1)
summary(fit2)

broom::tidy(fit1)

baddat <- read.csv("misc/bad_data.csv")
sclr(status ~ HIcens_mid, baddat)

confint(fit1)[, "2.5 %"]


preddata1 <- data.frame(HI = seq(0, 8, length.out = 101))
pred1 <- predict(fit1, preddata1)
head(pred1[, c("HI", "prot_l", "prot_point", "prot_u")])

preddata2 <- data.frame(HI = seq(0, 8, length.out = 101), NI = 1)
predict(fit2, preddata2)

# Model fit
fit <- sclr(status ~ HI + NI, dat)
summary(fit)

# Protection
var_name <- "HI"
prot_lvls1 <- data.frame(1)
class(prot_lvls1)
get_protection_level(fit1, "HI")

prot_lvls2 <- data.frame(NI = log(c(0.1, 10, 40)))
get_protection_level(fit2, "HI", prot_lvls2)


# Prediction
preddata <- data.frame(HI = seq(0, 8, length.out = 101), NI = 0)
seq(0, 8, length.out = 101)[44]
predict(fit2, preddata)[44, ]

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
sclr_log_likelihood(fit1)

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

