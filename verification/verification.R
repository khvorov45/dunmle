# Script to show that the model does what it's supposed to given ideal data.
# This will simulate many datasets, fit the model to each and summarise the
# results.
# Raw results will be in "results-*.csv"
# Summary will be in "summary-*.csv"
# Parallel computation controlled by "plan" below.
# Number of simualtions controlled by n_simulations below.
#
# Arseniy Khvorov
# Created 2019/10/17
# Last edit 2019/10/17

library(sclr)
library(broom)
library(dplyr)
library(readr)
library(future)
library(furrr)

plan(multiprocess)

this_folder <- "verification"

simulate_ideal_data <- function(n_sample = 1e4, 
                                lambda = 0.5, beta_0 = -5, beta_logtitre = 2,
                                logtitre_mean = 2, logtitre_sd = 2,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ideal_data <- tibble(.rows = n_sample) %>%
    mutate(
      logtitre = rnorm(n(), logtitre_mean, logtitre_sd),
      status = rbinom(
        n(), 1, lambda / (1 + exp(beta_0 + beta_logtitre * logtitre))
      )
    )
  attr(ideal_data, "true_values") <- tibble(
    term = c("lambda", "beta_0", "beta_logtitre"),
    true_value = c(lambda, beta_0, beta_logtitre)
  )
  attr(ideal_data, "seed") <- seed
  ideal_data
}

fit_sclr_model <- function(data, formula = status ~ logtitre) {
  tidy_fit <- sclr(formula, data) %>% tidy()
  if (!is.null(attr(data, "seed"))) tidy_fit$seed <- attr(data, "seed")
  inner_join(tidy_fit, attr(data, "true_values"), by = "term")
}

simulate_one <- function(index = NULL, 
                           init_seed = NULL, 
                           formula = status ~ logtitre, 
                           ...) {
  if (is.null(init_seed)) seed <- NULL
  else if (is.null(index)) seed <- init_seed
  else seed <- init_seed + index
  args <- list(...)
  args$seed <- seed
  result <- do.call(simulate_ideal_data, args) %>% fit_sclr_model(formula)
  if (!is.null(index)) result$index <- index
  result
}

simulate_many <- function(n_simulations = 1e3, 
                            formula = status ~ logtitre,
                            init_seed = NULL,
                            ...) {
  future_map_dfr(1:n_simulations, simulate_one, init_seed, formula, ...)
}

summarise_many <- function(simulation_results) {
  simulation_results %>%
    mutate(
      captures_true = (conf_low < true_value) & (conf_high > true_value)
    ) %>%
    group_by(term, true_value) %>%
    summarise(
      estimate_mean = mean(estimate), 
      estimate_sd = sd(estimate),
      se_mean = mean(std_error),
      coverage = sum(captures_true) / n()
    ) %>%
    ungroup()
}

save_csv <- function(data, folder, name, n_simulations) {
  write_csv(
    data, file.path(folder, paste0(name, "-", n_simulations, "sims.csv"))
  )
}

n_simulations <- 10
simulation_results <- simulate_many(
  n_simulations = n_simulations, init_seed = 20191017
)
save_csv(simulation_results, this_folder, "result", n_simulations)

simulation_summary <- summarise_many(simulation_results)
save_csv(simulation_summary, this_folder, "summary", n_simulations)
