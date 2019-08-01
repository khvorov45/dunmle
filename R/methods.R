# Methods for the sclr class

# Prints as a list
print.sclr <- function(fit) {
  class(fit) <- "list"
  print(fit)
  invisible(NULL)
}

# Variance-covariance matrix
vcov.sclr <- function(fit) {
  return(fit$covariance_mat)
}

# Returns the coefficients like lm
coef.sclr <- function(fit) {
  return(fit$parameters)
}

# Summary
summary.sclr <- function(fit) {
  cat("Call: ")
  print(fit$call[["formula"]])
  
  cat("\nParameter estimates\n")
  print(fit$parameters)
  
  cat("\n95% confidence intervals\n")
  print(fit$confint)
}
