# Package testing script for use during development
# Arseniy Khvorov
# Created 2019/07/29
# Last edit 2019/08/01

# Model fit
dat <- read.csv("data/H3onlyIndTime-1.csv")

fit <- sclr(status ~ HI, dat)

summary(fit)
