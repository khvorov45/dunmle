# Functions related to protection estimates

#' Protection level titres.
#' 
#' Calculates titres corresponding to a particular protection level.
#'
#' @param fit Object returned by \code{\link{sclr}}.
#' @param lvl Protection level to find titre values for.
#' @param ci_level Confidence level.
#' @param tol Tolerance.
#'
#' @return Named vector.
#' @export
get_protection_level <- function(
  fit, newdata, var_name, lvl = 0.5, ci_level = 0.95, tol = 10^(-7)
) {
  titre_low <- find_prot_titre_val(
    fit, newdata, var_name, "high_bound", lvl, ci_level
  )
  titre_point <- find_prot_titre_val(
    fit, newdata, var_name "prot", lvl, ci_level
  )
  titre_high <- find_prot_titre_val(
    fit, newdata, var_name, "low_bound", lvl, ci_level
  )
  newdata <- newdata[, names(newdata) != var_name]
  titre <- c(
    "low_log" = titre_low, "point_log" = titre_point, "high_log" = titre_high,
    "low" = exp(titre_low), "point" = exp(titre_point), "high" = exp(titre_high)
  )
  return(titre)
}

# Finds titre corresponding to a particular protection data variable
# (e.g. high bound)
find_prot_titre_val <- function(
  fit, var_name, prot_var_name, prot_var_val, ci_level = 0.95, tol = 10^(-7)
) {
  guess_range <- c(-100, 100)
  while (TRUE) {
    midpoint <- median(guess_range)

    prot_sample <- predict(fit, ci_level, midpoint)
    val <- prot_sample[, prot_var_name]
    
    if (abs(val - prot_var_val) < tol) return(midpoint)
    
    if (val < prot_var_val) guess_range[1] <- midpoint
    else guess_range[2] <- midpoint
  }
}
