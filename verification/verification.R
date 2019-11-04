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
# Last edit 2019/10/21

library(sclr)
library(broom)
library(dplyr)
library(readr)
library(future)
library(furrr)

plan(multiprocess)

ver_folder <- "verification"

simulate_ideal_data <- function(n, theta, beta_0, covariate_list, seed) {
  sclr_ideal_data(
    n, theta, beta_0, covariate_list, 
    outcome_name = "status",
    seed = seed,
    attach_true_vals = TRUE,
    attach_seed = TRUE
  )
}

fit_sclr_model <- function(data, covariate_list) {
  formula <- paste0("status~", paste(names(covariate_list), collapse = "+"))
  tidy_fit <- sclr(as.formula(formula), data) %>% broom::tidy()
  if (!is.null(attr(data, "seed"))) tidy_fit$seed <- attr(data, "seed")
  left_join(tidy_fit, attr(data, "true_values"), by = "term")
}

simulate_one <- function(index = NULL, init_seed = NULL, ...) {
  if (is.null(init_seed)) {
    seed <- NULL
  } else if (is.null(index)) {
    seed <- init_seed
  } else {
    seed <- init_seed + index
  }
  args <- list(...)
  args$seed <- seed
  result <- do.call(simulate_ideal_data, args) %>%
    fit_sclr_model(args$covariate_list)
  if (!is.null(index)) result$index <- index
  result
}

simulate_many <- function(n_simulations = 1e3, init_seed = NULL, ...) {
  future_map_dfr(1:n_simulations, simulate_one, init_seed, ...)
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

cov_list <- list(
  "logHI" = list(gen_fun = function(n) rnorm(n, 2, 2), true_par = 2),
  "logNI" = list(gen_fun = function(n) rnorm(n, 2, 2), true_par = 1)
)

n_simulations <- 2

simulation_results <- simulate_many(
  n_simulations = n_simulations,
  init_seed = 20191017,
  covariate_list = cov_list,
  n = 1e4, theta = 0, beta_0 = -5
)
save_csv(simulation_results, ver_folder, "result", n_simulations)

simulation_summary <- summarise_many(simulation_results)
save_csv(simulation_summary, ver_folder, "summary", n_simulations)
