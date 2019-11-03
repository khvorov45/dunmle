# Simulating sample data
# Arseniy Khvorov
# Created 2019/08/30
# Last edit 2019/10/21

one_titre_data <- sclr_ideal_data(
  n = 5000, theta = 0, beta_0 = -5,
  covariate_list = list(
    logHI = list(gen_fun = function(n) rnorm(n, 2, 2), true_par = 2)
  ),
  seed = 20191021
)
two_titre_data <- sclr_ideal_data(
  n = 5000, theta = 0, beta_0 = -5,
  covariate_list = list(
    logHI = list(gen_fun = function(n) rnorm(n, 2, 2), true_par = 2),
    logNI = list(gen_fun = function(n) rnorm(n, 2, 2), true_par = 1)
  ),
  seed = 20191021
)

usethis::use_data(one_titre_data, two_titre_data, overwrite = TRUE)
