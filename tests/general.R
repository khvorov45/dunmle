# General test
# Arseniy Khvorov
# Created 2019/07/29
# Last edit 2019/08/22

fit1 <- sclr(status ~ logHI, sclr_one_titre_data)
fit2 <- sclr(status ~ logHI + logNI, sclr_two_titre_data)

summary(fit1)
summary(fit2)

get_protection_level(fit1, "logHI")
