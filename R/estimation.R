#' Temporary funciton for testing purposes
#'
#' @param x Anything that can be passed to cat. New bit.
#' @export
#' @examples
#' testing_documentation("test")
testing_documentation <- function(x) {
    cat("You passed", x, "\n")
}

#' Fits the scaled logit model
#'
#' @param formula an object of class "formula":
#' a symbolic description of the model to be fitted.
#' @param data a data frame.
#' @export
sclr <- function(formula, data) {
  
  # Manipulations to get to response vector and design matrix
  
  ## The call is passed to model.frame
  mf <- match.call()
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf)
  
  ## Design matirix
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)
  
  ## Response vector
  y <- model.response(mf)
  
  # Actual model fit
  fit <- sclr_fit(y, x)

  return(fit)
}

#' Fits the scaled logit model
#'
#' @param y A vector of observations.
#' @param x A design matrix.
#' @param n_iter Number of Newton-Raphson iterations
sclr_fit <- function(y, x, n_iter = 100) {
  
  # Parameter vector with initial values
  pars <- c("lambda" = mean(y), "beta_0" = 0, "beta_titre" = 0)
  pars_mat <- matrix(pars)
  rownames(pars_mat) <- names(pars)
  # Work out the MLEs
  for (iter in 1:n_iter) {
    jacobian_mat <- get_jacobian(y, x)
    inv_jacobian_mat <- matlib::inv(jacobian_mat)
    pars_mat <- pars_mat - inv_jacobian_mat %*% scores_mat
  }
}

# Computes the jacobian matrix
get_jacobian <- function(y, x, pars) {
  
  # Get the parameters out
  lambda <- pars["lambda"]
  beta_0 <- pars["beta_0"]
  beta_titre <- pars["beta_titre"]
  
  d2l_dl <- get_d2l_dl(y, x, pars)
}

# Computes the second derivative of log-likelihood with respect to lambda
get_d2l_dl <- function(y, x, pars) {

  # Get the parameters out
  lambda <- pars["lambda"]
  beta_0 <- pars["beta_0"]
  beta_titre <- pars["beta_titre"]

  d2l_dl <- -y * 1/lambda^2 - (1 - y) * (1 / get_denom(y, x, pars))
}

# Computes the commonly-occuring denominator term
get_denom <- function(y, x, pars) {

  # Get the parameters out
  lambda <- pars["lambda"]
  beta_0 <- pars["beta_0"]
  beta_titre <- pars["beta_titre"]
}


#' Computes the log-likelihood
#'
#' @param x Anything that can be passed to cat. New bit.
log_likelihood <- function(x) {
  cat("You passed", x, "\n")
}
