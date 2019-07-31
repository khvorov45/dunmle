# Likelihood-related functions

#' Log-likelihood
#' 
#' Computes the log-likelihood of the scaled logit model
#' at a given set of parameter values.
#'
#' @param y A vector of observations.
#' @param x A design matrix.
#' @param pars A named vector of parameter values.
#' 
#' @export
sclr_log_likelihood <- function(y, x, pars) {
  lambda <- pars["lambda"]
  beta_0 <- pars["beta_0"]
  beta_titre <- pars["beta_titre"]
  l <- y * log(lambda) + 
    (1 - y) * log(1 + exp(beta_0 + beta_titre * x[, 2]) - lambda) -
    log(1 + exp(beta_0 + beta_titre * x[, 2]))
  l <- sum(l)
  return(l)
}
