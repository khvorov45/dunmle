#' Fits the scaled logit model
#' 
#' Used to fit the scaled logit model from Dunning (2006).
#'
#' @param formula an object of class "formula":
#' a symbolic description of the model to be fitted.
#' @param data a data frame.
#' @param tol Tolerance. Used when \code{n_iter} is \code{NULL}.
#' @param n_iter Number of Newton-Raphson iterations. \code{tol} is ignored
#' when this is not \code{NULL}.
#' 
#' @export
sclr <- function(formula, data, tol = 10^(-7), n_iter = NULL) {
  
  cl <- match.call()
  
  # Manipulations to get to response vector and design matrix
  
  ## The call is passed to model.frame
  mf <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf)
  
  ## Design matirix
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)
  
  ## Response vector
  y <- model.response(mf)
  
  # Actual model fit

  fit <- sclr_fit(y, x, tol, n_iter)

  fit[["call"]] <- cl

  return(fit)
}

#' Fitter function for the scaled logit model
#' 
#' Computing engine behind \code{\link{sclr}}.
#'
#' @param y A vector of observations.
#' @param x A design matrix.
#' @param tol Tolerance. Used when \code{n_iter} is \code{NULL}.
#' @param n_iter Number of Newton-Raphson iterations. \code{tol} is ignored
#' when this is not \code{NULL}.
#' 
#' @export
sclr_fit <- function(y, x, tol = 10^(-7), n_iter = NULL) {
  
  # Parameter vector with initial values
  pars <- c("lambda" = mean(y), "beta_0" = 0, "beta_titre" = 0)
  pars_mat <- matrix(pars)
  rownames(pars_mat) <- names(pars)
  
  # Work out the MLEs
  n_iter_cur <- 1
  while (TRUE) {
    pars_mat_prev <- pars_mat
    jacobian_mat <- get_jacobian(y, x, pars_mat)
    inv_jacobian_mat <- matlib::inv(jacobian_mat)
    scores_mat <- get_scores(y, x, pars_mat)
    pars_mat <- pars_mat_prev - inv_jacobian_mat %*% scores_mat
    n_iter_cur <- n_iter_cur + 1
    if (!is.null(n_iter)) {
      if (n_iter_cur > n_iter) break
    } else {
      deltas <- pars_mat - pars_mat_prev
      if (all(deltas <= tol)) break
    }
  }

  # Build the return list
  parameters <- as.vector(pars_mat)
  names(parameters) <- rownames(pars_mat)
  covariance_mat <- -inv_jacobian_mat
  dimnames(covariance_mat) <- list(names(parameters), names(parameters))
  fit <- list(
    parameters = parameters,
    covariance_mat = covariance_mat,
    log_likelihood = sclr_log_likelihood(y, x, parameters),
    n_converge = n_iter_cur - 1
  )
  return(fit)
}

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
