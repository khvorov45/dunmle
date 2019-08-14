# Methods for the sclr class

# Prints as a list
#' @export
print.sclr <- function(fit) {
  class(fit) <- "list"
  print(fit)
  invisible(NULL)
}

# Variance-covariance matrix
#' @export
vcov.sclr <- function(fit) {
  return(fit$covariance_mat)
}

# Returns the coefficients like lm
#' @export
coef.sclr <- function(fit) {
  return(fit$parameters)
}

# Summary
#' @export
summary.sclr <- function(fit) {
  cat("Call: ")
  print(fit$call[["formula"]])
  
  cat("\nParameter estimates\n")
  print(fit$parameters)
  
  cat("\n95% confidence intervals\n")
  print(fit$confint)
}

#' Predict method for scaled logit model fit.
#' 
#' Returns only the protection estimates. The only supported interval is
#' a confidence interval (i.e. the interval for the expected value)
#' 
#' @param fit Object returned by \code{\link{sclr}}.
#' @param newdata A dataframe with all covariates. Names should be as they
#' appear in the formula in the call to \code{\link{sclr}}.
#' @param ci_lvl Confidence level for the calculated interval. 
#' Default is 0.95.
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
  predret <- data.frame(
    prot_point_lin = prot_point_lin,
    prot_sd_lin = sds,
    prot_l_lin = prot_l_lin,
    prot_u_lin = prot_u_lin,
    prot_point = invlogit(prot_point_lin),
    prot_l = invlogit(prot_l_lin),
    prot_u = invlogit(prot_u_lin)
  )
  predret <- cbind(predret, newdata)
  return(predret)
}
