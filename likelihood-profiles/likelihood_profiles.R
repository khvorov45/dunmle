# Likelihood profiles
# Arseniy Khvorov
# Created 2019/10/21
# Last edit 2019/10/21

library(sclr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(viridis)

this_folder <- "likelihood-profiles"

calc_likelihood <- function(lambda, beta_0, fit) {
  sclr_log_likelihood(fit, c("lambda" = lambda, "beta_0" = beta_0))
}

# The value of lambda to combine with beta_0 to give prob
lambda_fun <- function(beta_0, prob) {
  lambda <- prob * (1 + exp(beta_0))
  lambda[lambda < 0 | lambda > 1] <- NA_real_
  lambda
}

add_variable <- function(dat, name, values) {
  og_length <- nrow(dat)
  dat %>%
    slice(rep(1:n(), each = length(values))) %>%
    mutate(!!sym(name) := rep(values, times = og_length))
}

fit <- sclr(status ~ 1, one_titre_data)

# Proportion infected with the sample
prob <- sum(one_titre_data$status == 1) / nrow(one_titre_data)

n_points <- 101

likelihoods <- tibble(lambda = seq(0.001, 1, length.out = n_points + 3)) %>%
  add_variable("beta_0", seq(-5, 5, length.out = n_points)) %>%
  mutate(likelihood = map2_dbl(lambda, beta_0, calc_likelihood, fit))

pl_no_cov <- likelihoods %>%
  ggplot(aes(beta_0, lambda, fill = likelihood)) +
  dark_theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_tile(stat = "identity") +
  scale_fill_viridis() +
  stat_function(fun = lambda_fun, args = list(prob = prob), col = "black") +
  labs(
    title = "Likelihood profile with no covariates",
    caption = paste0(
      "Black line marks lambda and beta_0 values which\n", 
      "combine to give the proportion infected in the sample"
    )
  )
ggsave_dark(
  file.path(this_folder, "no_covariate.pdf"), pl_no_cov, dark = FALSE,
  width = 11, height = 10, units = "cm"
)

lambda_1 <- sclr_ideal_data(lambda = 1, seed = 20191022)

fit_l1 <- sclr(status ~ logHI, lambda_1) # Error - does not converge
tidy(fit_l1)

fit_l1 <- new_sclr

sclr_log_likelihood()

