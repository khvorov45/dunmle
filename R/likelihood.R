# Likelihood-related functions

#' Log-likelihood
#' 
#' Computes the log-likelihood of the scaled logit model
#' at a given set of parameter values.
#'
#' @param fit An object returned by \code{\link{sclr}}.
#' @param pars A named vector of parameter values. If \code{NULL}
#' then the MLE estimates from fit will be used
#' 
#' @export
sclr_log_likelihood <- function(fit, pars = NULL) {
  if (is.null(pars)) pars <- fit$parameters
  lambda <- pars["lambda"]
  beta_0 <- pars["beta_0"]
  beta_titre <- pars[3]
  y <- fit$model_response
  x <- fit$model_matrix
  l <- y * log(lambda) + 
    (1 - y) * log(1 + exp(beta_0 + beta_titre * x[, 2]) - lambda) -
    log(1 + exp(beta_0 + beta_titre * x[, 2]))
  l <- sum(l)
  return(l)
}
