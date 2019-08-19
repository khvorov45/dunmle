# Fitting functions

#' Fits the scaled logit model
#' 
#' Used to fit the scaled logit model from Dunning (2006).
#' 
#' The model is of the form 
#' \deqn{P(Y = 1) = \lambda(1 - logit^{-1}(\beta_0 + \beta_1X_1 + 
#'   \beta_2X_2 + ... + \beta_kX_k))}
#' Where \eqn{Y} is the binary outcome idicator, 
#' (eg. 1 - infected, 0 - not infected). \eqn{X} - covariate. 
#' \eqn{k} - number of covariates.
#' The likelihood maximisation uses the Newton-Raphson algorithm. 
#' Initial values are always 1 for the covariate coefficients 
#' (and the associated intercept) and the proportion of infected for the 
#' baseline risk. If an algorithm's iteration produces values that cannot 
#' be used, a different set of initial values is chosen randomly and the 
#' algorithm starts over. Unusable parameter values are: baseline risk outside
#' of (0, 1) (likelihood undefined) and covariate coefficients (or the 
#' associated intercept) greater than 100 in magnitude (these terms 
#' appear in the exponent, if they are too big then calculations will fail).
#'
#' @param formula an object of class "formula":
#' a symbolic description of the model to be fitted.
#' @param data a data frame.
#' @param calc_ci Whether to calculate confindence intervals. 
#'   \code{TRUE} by default.
#' @param ci_lvl Confidence interval level for the parameter estimates.
#' @param calc_ll Whether to calculate log likelihood at MLEs. 
#'   \code{TRUE} by default.
#' @param tol Tolerance. Used when \code{n_iter} is \code{NULL}.
#' @param n_iter Number of Newton-Raphson iterations. \code{tol} is ignored
#' when this is not \code{NULL}.
#' @param max_tol_it Maximum tolerated iterations. If it fails to coverge 
#' within this number of iterations, will return with an error.
#' 
#' @return An object of class \code{sclr}. This is a list with the following
#' elements:
#' \describe{
#'   \item{paramters}{Maximum likelihood estimates of the parameter values.}
#'   \item{covariance_mat}{The variance-convariance matrix of the parameter
#'     estimates.}
#'   \item{n_converge}{The number of Newton-Raphson iterations 
#'     (including resets) that were required for convergence.}
#'   \item{confint}{Confidence intervals of the parameter estimates.}
#'   \item{x}{Model matrix derived from \code{formula} and \code{data}.}
#'   \item{y}{Response matrix derived from \code{formula} and \code{data}.}
#'   \item{log_likelihood}{Value of log-likelihood calculated at the ML 
#'     estimates of parameters.}
#'   \item{call}{The original call to \code{sclr}.}
#'   \item{model}{Model frame object derived from 
#'     \code{formula} and \code{data}.}
#'   \item{terms}{Terms object derived from model frame.}
#' }
#' Methods supported: \code{print}, \code{vcov}, \code{coef}, 
#' \code{summary}, \code{predict}, \code{tidy} (\code{broom} package)
#' 
#' @references 
#' Dunning AJ (2006). 
#' "A model for immunological correlates of protection." 
#' Statistics in Medicine, 25(9), 1485-1497. doi: 10.1002/sim.2282.
#' 
#' @export
sclr <- function(
  formula, data, calc_ci = TRUE, ci_lvl = 0.95, calc_ll = TRUE,
  tol = 10^(-7), n_iter = NULL, max_tol_it = 10^4
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
  
  if (calc_ci) fit$confint <- confint(fit, level = ci_lvl)
  
  fit$x <- x
  fit$y <- y
  
  if (calc_ll) fit$log_likelihood <- sclr_log_likelihood(fit)
  
  fit$call <- cl
  fit$model <- mf
  fit$terms <- mt

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
