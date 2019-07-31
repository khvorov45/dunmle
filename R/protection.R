# Functions related to protection estimates

#' Protection curve data
#' 
#' Calculates the data required for the protection curve.
#'
#' @param fit Object returned by \code{\link{sclr}}.
#' @param ci_level Confidence level for the curve.
#' @param log_titre_range Range of log titres over which to calculate the 
#' protection.
#' @param n_points Amount of points in \code{log_titre_range}.
#'
#' @return A dataframe.
#' @export
get_protection_data <- function(
  fit, ci_level = 0.95, log_titre_range = c(0, 8), n_points = 101
) {
  par_names <- c("beta_0", "beta_titre")
  vcov_mat <- fit$covariance_mat[par_names, par_names]
  pars_mat <- fit$parameters[par_names]
  if (length(log_titre_range) == 2) {
    titre_vals <- seq(
      log_titre_range[1], log_titre_range[2], length.out = n_points
    )
  } else if (length(log_titre_range) == 1) {
    titre_vals = log_titre_range
  } else {
    stop(
      "log_titre_range should be of legnth 1 or 2, not ", 
      length(log_titre_range)
    )
  } 
  
  # Work out protection point estimates
  titre_mats_point <- lapply(titre_vals, function(x) matrix(c(1, x), nrow = 2))
  point_est <- lapply(titre_mats_point, function(x) pars_mat %*% x)
  point_est <- unlist(point_est)
  
  # Work out the variances
  titre_mats <- lapply(titre_vals, function(x) matrix(c(1, 0, 0, x), nrow = 2))
  vcov_mat_mod <- lapply(titre_mats, function(x) x %*% vcov_mat %*% x)
  vars <- lapply(vcov_mat_mod, sum)
  vars <- unlist(vars)
  
  # Work out CIs
  low_coef <- qnorm((1 - ci_level) / 2)
  high_coef <- qnorm((1 + ci_level) / 2)
  sds <- sqrt(vars)
  low_bound <- point_est + low_coef * sds
  high_bound <- point_est + high_coef * sds
  
  # Collate
  prot_data <- data.frame(
    logtitre = titre_vals,
    titre = exp(titre_vals),
    prot_logit = point_est,
    prot = invlogit(point_est),
    vars_logit = vars,
    sds_logit = sds,
    low_bound_logit = low_bound,
    low_bound = invlogit(low_bound),
    high_bound_logit = high_bound,
    high_bound = invlogit(high_bound)
  )
  
  return(prot_data)
}

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
  fit, lvl = 0.5, ci_level = 0.95, tol = 10^(-7)
) {
  titre_low <- find_prot_titre_val(fit, "high_bound", lvl, ci_level)
  titre_point <- find_prot_titre_val(fit, "prot", lvl, ci_level)
  titre_high <- find_prot_titre_val(fit, "low_bound", lvl, ci_level)
  titre <- c(
    "low_log" = titre_low, "point_log" = titre_point, "high_log" = titre_high,
    "low" = exp(titre_low), "point" = exp(titre_point), "high" = exp(titre_high)
  )
  return(titre)
}

# Finds titre corresponding to a particular protection data variable
# (e.g. high bound)
find_prot_titre_val <- function(
  fit, prot_var_name, prot_var_val, ci_level = 0.95, tol = 10^(-7)
) {
  guess_range <- c(-100, 100)
  while (TRUE) {
    midpoint <- median(guess_range)
    prot_sample <- get_protection_data(fit, ci_level, midpoint)
    val <- prot_sample[, prot_var_name]
    
    if (abs(val - prot_var_val) < tol) return(midpoint)
    
    if (val < prot_var_val) guess_range[1] <- midpoint
    else guess_range[2] <- midpoint
  }
}
