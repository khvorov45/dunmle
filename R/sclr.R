# Fitting functions
# Arseniy Khvorov
# Created 2019/07/31
# Last edit 2019/10/14

#' Fits the scaled logit model
#'
#' Used to fit the scaled logit model from Dunning (2006).
#'
#' The model is of the form \deqn{P(Y = 1) = \lambda(1 - logit^{-1}(\beta_0 +
#' \beta_1X_1 + \beta_2X_2 + ... + \beta_kX_k))} Where \eqn{Y} is the binary
#' outcome indicator, (eg. 1 - infected, 0 - not infected). \eqn{X} - covariate.
#' \eqn{k} - number of covariates. Computing engine behind the fitting is
#' \code{\link{sclr_fit}}.
#'
#' @param formula an object of class "formula": a symbolic description of the
#'   model to be fitted.
#' @param data a data frame.
#' @param ci_lvl Confidence interval level for the parameter estimates.
#' @inheritParams sclr_fit
#'
#' @return An object of class \code{sclr}. This is a list with the following
#'   elements:
#'
#'   \item{parameters}{Maximum likelihood estimates of the parameter values.}
#'
#'   \item{covariance_mat}{The variance-covariance matrix of the parameter
#'   estimates.}
#'
#'   \item{n_converge}{The number of Newton-Raphson iterations (including
#'   resets) that were required for convergence.}
#'
#'   \item{x}{Model matrix derived from \code{formula} and \code{data}.}
#'
#'   \item{y}{Response matrix derived from \code{formula} and \code{data}.}
#'   
#'   \item{call}{The original call to \code{sclr}.}
#'
#'   \item{model}{Model frame object derived from \code{formula} and 
#'   \code{data}.}
#'
#'   \item{terms}{Terms object derived from model frame.}
#'   
#'   \item{ci}{Confidence intervals of the parameter estimates.}
#'   
#'   \item{log_likelihood}{Value of log-likelihood calculated at the ML
#'   estimates of parameters.} 
#'
#'   Methods supported: \code{\link[=print.sclr]{print}},
#'   \code{\link[=vcov.sclr]{vcov}}, \code{\link[=coef.sclr]{coef}},
#'   \code{\link[=summary.sclr]{summary}}, \code{\link[=predict.sclr]{predict}},
#'   \code{\link[=tidy.sclr]{tidy}} (\code{\link{broom}} package).
#'
#' @references Dunning AJ (2006). "A model for immunological correlates of
#'   protection." Statistics in Medicine, 25(9), 1485-1497.
#'   \url{https://doi.org/10.1002/sim.2282}.
#'
#' @examples
#' library(sclr)
#' fit1 <- sclr(status ~ logHI, one_titre_data)
#' summary(fit1)
#' @importFrom stats confint
#' @importFrom rlang abort
#'
#' @export
sclr <- function(formula, data = NULL, 
                 ci_lvl = 0.95, 
                 tol = 10^(-7), 
                 n_iter = NULL, 
                 max_tol_it = 10^4, 
                 n_conv = 3,
                 conventional_names = FALSE) {
  
  if (!inherits(formula, "formula") || missing(formula)) 
    abort("must supply a formula")

  cl <- match.call()
  mf <- stats::model.frame(formula, data, drop.unused.levels = TRUE)

  # Design matirix
  mt <- attr(mf, "terms")
  x <- stats::model.matrix(mt, mf)

  # Response vector
  y <- stats::model.response(mf)
  if (!all(y %in% c(0, 1))) abort("response should be a vector with 0 and 1")
  if (!is.numeric(y)) abort("response should be numeric")

  # Actual model fit
  fit <- sclr_fit(y, x, tol, n_iter, max_tol_it, n_conv, conventional_names)

  # Build the return list
  fit <- new_sclr(fit, x, y, cl, mf, mt)
  fit$ci <- confint(fit, level = ci_lvl)
  fit$log_likelihood <- sclr_log_likelihood(fit)
  fit
}

#' Create a new \code{sclr} object
#' 
#' \code{new_sclr} creates the object \code{\link{sclr}} returns.
#' \code{is_sclr} checks if the object is of class \code{sclr}.
#'
#' @param fit A list returned by \code{\link{sclr_fit}}.
#' @param x Model matrix.
#' @param y Model response.
#' @param cl Call.
#' @param mf Model frame.
#' @param mt Model terms.
#'
#' @return \code{sclr} object
#' @export
new_sclr <- function(fit, x, y, cl, mf, mt) {
  stopifnot(is.list(fit))
  stopifnot(is.matrix(x))
  stopifnot(is.numeric(y))
  stopifnot(is.call(cl))
  stopifnot(is.data.frame(mf))
  stopifnot(inherits(mt, "terms"))
  fit <- c(fit, list(
    x = x, y = y,
    call = cl, model = mf, terms = mt
  ))
  class(fit) <- "sclr"
  fit
}

#' @rdname new_sclr
#' @export
is_sclr <- function(fit) "sclr" %in% class(sclr)

#' Fitter function for the scaled logit model
#'
#' Computing engine behind \code{\link{sclr}}.
#'
#' The likelihood maximisation uses the Newton-Raphson algorithm. Initial values
#' are always 1 for the covariate coefficients (and the associated intercept)
#' and the proportion of infected for the baseline risk.
#'
#' The algorithm will pick a new guess and restart under a set of conditions.
#'
#' 1) Algorithm's iteration produces estimate guesses that cannot be used -
#' baseline risk outside of (0, 1) since likelihood is undefined.
#'
#' 2) The second derivative matrix produced by the current estimates is "bad" -
#' positive diagonal or missing values due to failing large number calculations.
#'
#' @param y A vector of observations.
#' @param x A design matrix.
#' @param tol Tolerance. Used when \code{n_iter} is \code{NULL}.
#' @param n_iter Number of Newton-Raphson iterations. \code{tol} is ignored when
#'   this is not \code{NULL}.
#' @param max_tol_it Maximum tolerated iterations. If it fails to converge
#'   within this number of iterations, will return with an error.
#' @param n_conv Number of times the algorithm has to converge (to work around
#'   local maxima).
#' @param conventional_names If \code{TRUE}, estimated parameter names will be
#'   (Baseline), (Intercept) and the column names in the model matrix. Otherwise
#'   - lambda, beta_0 and beta_ prefix in front of column names in the model
#'   matrix.
#'
#' @export
sclr_fit <- function(y, x, tol = 10^(-7), n_iter = NULL, max_tol_it = 10^4,
                     n_conv = 3, conventional_names = FALSE) {

  # Parameter vector with initial values
  n_par <- ncol(x) + 1
  pars_mat <- matrix(rep(1, n_par))

  rownames(pars_mat) <- get_par_names(x, conventional_names)
  pars_mat[1, ] <- mean(y)

  pars_mat_init <- pars_mat # Initial guess of values

  x_coeffs <- get_x_coeffs(x) # To avoid recalculations

  # Work out the MLEs
  n_iter_cur <- 1
  conv_count <- 0
  lls <- c()
  rets <- list()
  while (TRUE) {

    # Check that the current set is workable
    if (is_bad(pars_mat)) pars_mat <- guess_again(pars_mat_init)
    pars_mat_prev <- pars_mat # Save for later

    # Calculate the commonly occuring expession
    exp_Xb <- get_exp_Xb(y, x, pars_mat)

    # Log-likelihood second derivative (negative of information) matrix
    jacobian_mat <- get_jacobian(y, x, pars_mat, exp_Xb, x_coeffs)
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
        ll_cur <- sclr_log_likelihood(list(x = x, y = y), pars_mat)
        conv_count <- conv_count + 1
        lls[[conv_count]] <- ll_cur
        rets[[conv_count]] <- list(
          "invjac" = inv_jacobian_mat, "pars" = pars_mat
        )
        if (conv_count == n_conv) break
        pars_mat <- guess_again(pars_mat_init)
        next
      }
    }
    if (n_iter_cur > max_tol_it) {
      stop("did not converge in ", max_tol_it, " iterations\n")
    }
    n_iter_cur <- n_iter_cur + 1
  }

  i_best <- which.max(lls)
  pars_mat <- rets[[i_best]][["pars"]]
  inv_jacobian_mat <- rets[[i_best]][["invjac"]]

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

#' Check the parameter matrix.
#'
#' Checks if the current parameter guesses are OK for derivative calculations.
#' Returns \code{TRUE} if they are and \code{FALSE} otherwise.
#'
#' @param pars_mat The current matrix of parameter guesses.
#'
#' @noRd
is_bad <- function(pars_mat) {
  lambda <- pars_mat[1, ]
  if ((lambda < 0) | (lambda > 1)) {
    return(TRUE)
  } # Outside of the defined region
  return(FALSE)
}

#' Check the second derivative matrix
#'
#' Checks if the current parameter guesses are OK for derivative calculations.
#' Returns \code{TRUE} if they are and \code{FALSE} otherwise.
#'
#' @param jac The current second derivative matrix.
#'
#' @noRd
is_bad_jac <- function(jac) {
  if (any(is.na(jac))) {
    return(TRUE)
  }
  if (any(diag(jac) > 0)) {
    return(TRUE)
  }
  return(FALSE)
}

#' Create a new guess
#'
#' Creates a matrix with new parameter guesses.
#' The matrix has the same structure as the matrix with previous guesses.
#'
#' @param pars_mat The current matrix of parameter guesses.
#'
#' @noRd
guess_again <- function(pars_mat) {
  delta <- matrix(
    c(
      stats::runif(1, min = 0, max = 1),
      stats::rnorm(nrow(pars_mat) - 1, mean = 0, sd = 2)
    ),
    ncol = 1
  )
  rownames(delta) <- rownames(pars_mat)
  return(delta)
}

#' Check for local minimum
#'
#' @param jac Second derivative matrix
#'
#' @noRd
is_minimum <- function(jac) {
  eigenvals <- eigen(jac, only.values = TRUE)$values
  eigenvals[dplyr::near(eigenvals, 0)] <- 0
  if (any(eigenvals > 0)) cat("MINIMUM\n")
}
