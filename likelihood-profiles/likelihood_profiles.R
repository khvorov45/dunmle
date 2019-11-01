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
library(broom)

this_folder <- "likelihood-profiles"

calc_likelihood <- function(lambda, beta_0, 
                            beta_logHI = NULL, fit = NULL,
                            x = NULL, y = NULL) {
  pars <- c("lambda" = lambda, "beta_0" = beta_0, "beta_logHI" = beta_logHI)
  if (!is.null(fit)) sclr_log_likelihood(fit, pars = pars)
  else if (is.null(x) || is.null(y)) rlang::abort("specify fit or x and y")
  else sclr_log_likelihood(x = x, y = y, pars = pars)
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

# Fit with no covariates ======================================================

fit <- sclr(status ~ 1, one_titre_data)

# Proportion infected with the sample
prob <- sum(one_titre_data$status == 1) / nrow(one_titre_data)

n_points <- 101

likelihoods <- tibble(lambda = seq(0.001, 1, length.out = n_points + 3)) %>%
  add_variable("beta_0", seq(-5, 5, length.out = n_points)) %>%
  mutate(likelihood = map2_dbl(lambda, beta_0, calc_likelihood, fit = fit))

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
  scale_fill_viridis("Log likelihood") +
  stat_function(fun = lambda_fun, args = list(prob = prob), col = "black")
ggsave_dark(
  file.path(this_folder, "no_covariate.pdf"), pl_no_cov, dark = FALSE,
  width = 11, height = 9, units = "cm"
)

# Fit to data with a true baseline of 1 =======================================

lambda_1 <- sclr_ideal_data(lambda = 1, seed = 20191022)
y <- model.frame(status ~ logHI, data = lambda_1) %>% model.response()
x <- model.matrix(status ~ logHI, data = lambda_1)

fit_l1 <- sclr(status ~ logHI, lambda_1) # Error - does not converge
fit_l1_lr <- glm(status ~ logHI, lambda_1, family = binomial(link = "logit"))
mles <- tidy(fit_l1_lr) %>%
  mutate(
    estimate = -estimate, 
    term = recode(term, "(Intercept)" = "beta_0", "logHI" = "beta_logHI")
  ) %>%
  select(term, estimate, std_error = std.error)

likelihoods_l1 <- tibble(
  lambda = seq(0.001, 1, length.out = 1e3),
  beta_0 = mles$estimate[mles$term == "beta_0"],
  beta_logHI = mles$estimate[mles$term == "beta_logHI"],
  likelihood = pmap_dbl(
    list(lambda, beta_0, beta_logHI), calc_likelihood, x = x, y = y
  )
)

base1 <- likelihoods_l1 %>%
  ggplot(aes(lambda, likelihood)) +
  dark_theme_bw(verbose = FALSE) +
  ylab("Log likelihood") +
  geom_line()
ggsave_dark(
  file.path(this_folder, "baseline_1.pdf"),
  base1, dark = FALSE,
  width = 10, height = 10, units = "cm"
)
