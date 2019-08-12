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
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  ## Design matirix
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)

  ## Response vector
  y <- model.response(mf)
  
  # Actual model fit

  fit <- sclr_fit(y, x, tol, n_iter, max_tol_it)
  
  class(fit) <- "sclr"
  
  fit$confint <- confint(fit)
  fit$log_likelihood <- sclr_log_likelihood(fit)
  
  fit$call <- cl
  fit$model <- mf
  fit$terms <- mt
  fit$x <- x
  fit$y <- y

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
  pars_mat <- matrix(rep(1, n_par))
  rownames(pars_mat) <- c("lambda", get_par_names(x))
  pars_mat["lambda", ] <- mean(y)
  
  pars_mat_init <- pars_mat # Initial guess of values

  # Work out the MLEs
  n_iter_cur <- 1
  while (TRUE) {
    
    # Check that the current set is workable
    if (is_bad(pars_mat)) pars_mat <- guess_again(pars_mat_init)
    pars_mat_prev <- pars_mat # Save for later
    
    # Calculate the commonly occuring expessions
    exp_Xb <- get_exp_Xb(y, x, pars_mat)
    
    # Log-likelihood second derivative (negative of information) matrix
    jacobian_mat <- get_jacobian(y, x, pars_mat, exp_Xb)
    if (is_bad_jac(jacobian_mat)) {
      pars_mat <- guess_again(pars_mat_init)
      next
    }
    
    # Invert to get the negative of the covariance matrix
    inv_jacobian_mat <- try(solve(jacobian_mat), silent = TRUE)
    if (inherits(inv_jacobian_mat, "try-error")) {
      pars_mat <- guess_again(pars_mat_init)
      next
    }
    if (is_bad_jac(inv_jacobian_mat)) {
      pars_mat <- guess_again(pars_mat_init)
      next
    }
    
    # Generalised Newton-Raphson
    scores_mat <- get_scores(y, x, pars_mat, exp_Xb)
    pars_mat <- pars_mat_prev - inv_jacobian_mat %*% scores_mat
    
    # See if this should be the last iteration
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
    
    n_iter_cur <- n_iter_cur + 1
  }

  # Build the return list
  parameters <- as.vector(pars_mat)
  names(parameters) <- rownames(pars_mat)
  covariance_mat <- -inv_jacobian_mat
  dimnames(covariance_mat) <- list(names(parameters), names(parameters))
  fit <- list(
    parameters = parameters,
    covariance_mat = covariance_mat,
    n_converge = n_iter_cur
  )
  return(fit)
}

# Checks if the current parameter guesses are OK for derivative calculations
is_bad <- function(pars_mat) {
  pars_betas <- get_betas_only(pars_mat)
  if (any(abs(pars_betas) > 100)) return(TRUE) # Too big for exp(...)
  lambda <- pars_mat["lambda", ]
  if ((lambda < 0) | (lambda > 1)) return(TRUE) # Outside of the defined region
  return(FALSE)
}

# Comes up with a new initial guess
guess_again <- function(pars_mat) {
  delta <- matrix(
    c(
      runif(1, min = 0, max = 1), 
      rnorm(nrow(pars_mat) - 1, mean = 0, sd = 2)
    ),
    ncol = 1
  )
  rownames(delta) <- rownames(pars_mat)
  return(delta)
}

# Checks if the jacobian is bad (positive diagonal or missing values)
is_bad_jac <- function(jac) {
  if (any(is.na(jac))) return(TRUE)
  if (any(diag(jac) > 0)) return(TRUE)
  return(FALSE)
}
