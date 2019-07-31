# Package testing script for use during development
# Arseniy Khvorov
# Created 2019/07/29
# Last edit 2019/07/31

# Model fit
dat <- read.csv("data/H3onlyInd-1.csv")

fit <- sclr(status ~ HI, dat, tol = 10^(-7))

summary(fit)
