# Fitting functions

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
#' @param max_tol_it Maximum tolerated iterations. If it fails to coverge 
#' within this number of iterations, will return with an error.
#' 
#' @export
sclr <- function(
  formula, data, tol = 10^(-7), n_iter = NULL, max_tol_it = 10^4
) {
  
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

  fit <- sclr_fit(y, x, tol, n_iter, max_tol_it)
  
  class(fit) <- "sclr"
  
  fit[["call"]] <- cl
  
  fit[["confint"]] <- confint(fit)
  
  fit$model_matrix <- x
  fit$model_response <- y
  fit$log_likelihood <- sclr_log_likelihood(fit)

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
#' @param max_tol_it Maximum tolerated iterations. If it fails to coverge 
#' within this number of iterations, will return with an error.
#' 
#' @export
sclr_fit <- function(y, x, tol = 10^(-7), n_iter = NULL, max_tol_it = 10^4) {
  
  # Parameter vector with initial values
  n_par <- ncol(x) + 1
  pars_mat <- matrix(rep(0, n_par))
  rownames(pars_mat) <- c("lambda", get_par_names(x))
  pars_mat["lambda", ] <- mean(y)
  
  pars_mat_init <- pars_mat
  
  # Work out the MLEs
  n_iter_cur <- 1
  while (TRUE) {
    #cat("top loop\n")
    #print(pars_mat)
    if (is_bad(pars_mat)) pars_mat <- guess_again(pars_mat_init)
    #cat("after guess again\n")
    #print(pars_mat)
    pars_mat_prev <- pars_mat
    jacobian_mat <- get_jacobian(y, x, pars_mat_prev)
    #print(jacobian_mat)
    if (is_bad_jac(jacobian_mat)) {
      pars_mat <- guess_again(pars_mat_init)
      next
    }
    #cat("=======\n")
    if (any(is.na(jacobian_mat))) stop("can't calculate second derivatives")
    inv_jacobian_mat <- base::solve(jacobian_mat)
    scores_mat <- get_scores(y, x, pars_mat_prev)
    pars_mat <- pars_mat_prev - inv_jacobian_mat %*% scores_mat
    #print(pars_mat)
    n_iter_cur <- n_iter_cur + 1
    if (!is.null(n_iter)) {
      if (n_iter_cur > n_iter) break
    } else {
      deltas <- abs(pars_mat - pars_mat_prev)
      if (all(deltas <= tol)) {
        break
      } 
    }
    if (n_iter_cur > max_tol_it) 
      stop("did not converge in ", max_tol_it, " iterations\n")
  }

  # Build the return list
  parameters <- as.vector(pars_mat)
  names(parameters) <- rownames(pars_mat)
  covariance_mat <- -inv_jacobian_mat
  dimnames(covariance_mat) <- list(names(parameters), names(parameters))
  fit <- list(
    parameters = parameters,
    covariance_mat = covariance_mat,
    n_converge = n_iter_cur - 1
  )
  return(fit)
}

# Checks if the current parameter guesses are OK for derivative calculations
is_bad <- function(pars_mat) {
  pars_betas <- get_betas_only(pars_mat)
  if (any(abs(pars_betas) > 100)) return(TRUE)
  if (any(pars_betas[-1] < 0)) return(TRUE)
  lambda <- pars_mat["lambda", ]
  if ((lambda < 0) | (lambda > 1)) return(TRUE)
  return(FALSE)
}

# Comes up with a new initial guess
guess_again <- function(pars_mat) {
  delta <- matrix(c(0, abs(rnorm(nrow(pars_mat) - 1))), ncol = 1) # abs(rnorm)?
  pars_mat <- pars_mat + delta
  return(pars_mat)
}

# Checks if the jacobian is bad (positive diagonal)
is_bad_jac <- function(jac) {
  if (any(diag(jac) > 0)) return(TRUE)
  return(FALSE)
}
