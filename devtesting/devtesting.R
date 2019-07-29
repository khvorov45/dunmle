devtools::install_github("khvorov45/dunmle")

devtools::install("G:/Nexus/dunmle")

library(dunmle)

testing_documentation("test")

dat <- read.csv("data/simulated_data_1.csv")

sclr(status ~ NIcens, dat)
