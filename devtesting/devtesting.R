devtools::install_github("khvorov45/dunmle")

devtools::install("G:/Nexus/dunmle")

library(dunmle)

testing_documentation("test")

dat <- read.csv("data/simulated_data_1.csv")

fit <- sclr(status ~ HI, dat, tol = 10^(-7))

mf <- model.frame(status ~ HI, dat)

y <- model.response(mf)
x <- model.matrix(attr(mf, "terms"), mf)
