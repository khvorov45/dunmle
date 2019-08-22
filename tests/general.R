# General test
# Arseniy Khvorov
# Created 2019/07/29
# Last edit 2019/08/22

fit1 <- sclr(status ~ HI, sclronetitre)
fit2 <- sclr(status ~ HI + NI, sclrtwotitre)

summary(fit1)
summary(fit2)
