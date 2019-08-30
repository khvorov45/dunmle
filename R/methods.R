# Methods for the sclr class
# Arseniy Khvorov
# Created 2019/07/31
# Last edit 2019/08/30

#' Print a \code{sclr} object.
#' 
#' Summarises a \code{sclr} object for printing. For a dataframe summary, see
#' \code{\link[=tidy.sclr]{tidy}}.
#' 
#' @param fit An object returned by \code{\link{sclr}}.
#' 
#' @export
print.sclr <- function(fit) summary(fit)

#' @rdname print.sclr
#' @export
summary.sclr <- function(fit) {
  cat("Call: ")
  print(fit$call[["formula"]])
  
  cat("\nParameter estimates\n")
  print(fit$parameters)
  
  cat("\n95% confidence intervals\n")
  print(fit$confint)
  
  invisible(NULL)
}

#' Variance-covariance matrix
#' 
#' Returns the estimated variance-covariance matrix.
#' 
#' @inheritParams print.sclr
#'
#' @export
vcov.sclr <- function(fit) {
  return(fit$covariance_mat)
}

#' Coefficients
#' 
#' Returns the estimated model coefficients.
#' 
#' @inheritParams print.sclr
#'
#' @export
coef.sclr <- function(fit) {
  return(fit$parameters)
}

#' Predict method for scaled logit model fit.
#' 
#' Returns only the protection estimates. The only supported interval is
#' a confidence interval (i.e. the interval for the estimated expected value).
#' 
#' The model is \deqn{P(Y = 1) = \lambda(1 - logit^{-1}(\beta_0 +
#' \beta_1X_1 + \beta_2X_2 + ... + \beta_kX_k))} Where \eqn{Y} is the binary
#' outcome idicator, (eg. 1 - infected, 0 - not infected). \eqn{X} - covariate.
#' \eqn{k} - number of covariates.
#' This function calculates \deqn{\beta_0 + \beta_1X_1 + \beta_2X_2 + ..
#' . + \beta_kX_k} transformations at the covariate values found in 
#' \code{newdata} as well as the variance-covariance matrices of those
#' transformations. This is used to calculate the confidence intervals at the
#' given parameter values. The inverse logit transformation is then applied
#' to point estimates and interval bounds.
#' 
#' @param fit Object returned by \code{\link{sclr}}.
#' @param newdata A dataframe with all covariates. Names should be as they
#'   appear in the formula in the call to \code{\link{sclr}}.
#' @param ci_lvl Confidence level for the calculated interval.
#'
#' @return A \code{\link[tibble]{tibble}} obtained by adding the following
#' columns to \code{newdata}:
#' \item{prot_point_lin prot_l_lin prot_u_lin}{Point estimate, low and high 
#' bounds of the linear transformation.}
#' \item{prot_sd_lin}{Estimated standard deviation of the linear 
#' transformation.}
#' \item{prot_point prot_l prot_u}{Inverse logit-transformed 
#' point estimate, low and high bounds of the linear transformation.}
#' 
#' @export
predict.sclr <- function(fit, newdata, ci_lvl = 0.95) {
  
  # Covariates
  model_mat <- model.matrix(delete.response(fit$terms), newdata)
  
  # Estimated parameters
  ests <- coef(fit)
  ests_beta_mat <- matrix(ests[-1], ncol = 1) # Estmated betas
  
  # Point estimates
  prot_point_lin <- apply(model_mat, 1, function(x) x %*% ests_beta_mat)

  # Modified beta covariance matrices
  ests_beta_cov <- vcov(fit)[-1, -1] # Beta covariances
  x_coefs <- get_x_coeffs(model_mat) # x modifiers
  cov_mod_mats <- build_symm_mat(x_coefs) # In a list of matrices
  cov_modified <- lapply(cov_mod_mats, function(x) x * ests_beta_cov)
  
  # Standard deviations associated with each of the point estimates
  sds <- lapply(cov_modified, sum)
  sds <- unlist(sds)
  sds <- sqrt(sds)
  
  # Ranges
  lvl <- qnorm((1 + ci_lvl) / 2)
  prot_l_lin <- prot_point_lin - lvl * sds
  prot_u_lin <- prot_point_lin + lvl * sds
  
  # Collate
  predret <- tibble::tibble(
    prot_point_lin = prot_point_lin,
    prot_sd_lin = sds,
    prot_l_lin = prot_l_lin,
    prot_u_lin = prot_u_lin,
    prot_point = invlogit(prot_point_lin),
    prot_l = invlogit(prot_l_lin),
    prot_u = invlogit(prot_u_lin)
  )
  predret <- dplyr::bind_cols(predret, newdata)
  return(predret)
}

#' Tidy a \code{sclr} object.
#' 
#' Summarises the objects returned by \code{\link{sclr}} 
#' into a \code{\link[tibble]{tibble}}.
#'
#' @param fit An object returned by \code{\link{sclr}}.
#' @param ci_level Confidence level for the intervals.
#'
#' @return A \code{\link[tibble]{tibble}} with one row per model parameter. 
#' Columns:
#' \item{term}{Name of model parameter.}
#' \item{estimate}{Point estimate.}
#' \item{std_error}{Standard error.}
#' \item{conf_low}{Lower bound of the confidence interval.}
#' \item{conf_high}{Upper bound of the confidence interval.}
#' 
#' @importFrom broom tidy
#' 
#' @export
tidy.sclr <- function(fit, ci_level = 0.95) {
  pars <- tibble::tibble(
    term = names(fit$parameters),
    estimate = fit$parameters,
    std_error = sqrt(diag(vcov(fit)))
  )
  cis <- confint(fit, level = ci_level)
  cisdf <- tibble::tibble(
    term = rownames(cis), conf_low = cis[, 1], conf_high = cis[, 2]
  )
  fitsum <- dplyr::inner_join(pars, cisdf, by = "term")
  return(fitsum)
}
