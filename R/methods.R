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
