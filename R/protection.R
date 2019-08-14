# Functions related to protection estimates

#' Protection level calculations
#' 
#' Calculates covariate values corresponding to a particular protection level.
#' Only accepts one covariate at a time, fixed 
#' values of all the others should be provided.
#'
#' @param fit Object returned by \code{\link{sclr}}.
#' @param newdata A dataframe with all covariates except the one for which
#' protection values should be calculated.

#' @param lvl Protection level to find titre values for. Default is 0.5 (50\%) 
#' @param ci_level Confidence level for the calculated interval. 
#' Default is 0.95.
#' @param tol Tolerance. The values will be found numerically,
#' once the algorithm converges within \code{tol} of \code{lvl} 
#' it stops looking.
#'
#' @return A dataframe. Will have the same variables as \code{newdata} with
#' the addition of the \code{var_name} variable.
#' 
#' @export
get_protection_level <- function(
  fit, var_name, newdata = NULL, 
  lvl = 0.5, ci_level = 0.95, tol = 10^(-7)
) {
  titre_low <- find_prot_titre_val(
    fit, newdata, var_name, "prot_u", lvl, ci_level
  )
  titre_point <- find_prot_titre_val(
    fit, newdata, var_name, "prot_point", lvl, ci_level
  )
  titre_high <- find_prot_titre_val(
    fit, newdata, var_name, "prot_l", lvl, ci_level
  )
  titre <- rbind(titre_low, titre_point, titre_high)
  titre$prot_prob <- lvl
  titre$est <- "point"
  titre$est[titre$protvar == "prot_u"] <- "low bound"
  titre$est[titre$protvar == "prot_l"] <- "upper bound"
  titre$protvar <- NULL
  return(titre)
}

#' Search function for scaled logit protection covariate levels
#' 
#' The search engine behind \code{\link{get_protection_level}}. Should not
#' usually be necessary to call this directly.
#'
#' @param fit Object returned by \code{\link{sclr}}.
#' @param newdata A dataframe with all covariates except the one for which
#' protection values should be calculated.
#' @param var_name Name of the covariate for which the protection values should
#' be calculated. This name should appear in the formula of the call to
#' \code{\link{sclr}} which was used to generate \code{fit}.
#' @param prot_var_name A variable name amoung those returned by
#' \code{\link{predict.sclr}} which needs to equal \code{lvl} at the value of
#' \code{var_name} that is supposed to be found.
#' @param lvl Protection level to find titre values for. Default is 0.5 (50\%).
#' @param ci_level Confidence level for the calculated interval. 
#' Default is 0.95.
#' @param tol Tolerance. The values will be found numerically,
#' once the algorithm converges within \code{tol} of \code{lvl} 
#' it stops looking.
#'
#' @return A dataframe. Will have the same variables as \code{newdata} with
#' the addition of the \code{var_name} variable.
#' 
#' @export
find_prot_titre_val <- function(
  fit, newdata, var_name, prot_var_name, lvl = 0.5, 
  ci_level = 0.95, tol = 10^(-7)
) {
  
  # Initial guess interval
  newdata$guess_low <- -100
  newdata$guess_high <- 100
  
  # Check if the variable is protective
  is_protective <- coef(fit)[grepl(var_name, names(coef(fit)))] > 0
  
  # Binary search
  while (TRUE) {
    
    newdata[, var_name] <- (newdata$guess_low + newdata$guess_high) / 2
    prot_sample <- predict(fit, newdata, ci_level)
    
    curvals <- prot_sample[, prot_var_name]
    notfound <- abs(curvals - lvl) > tol
    
    if (sum(notfound) == 0) {
      newdata[, c("guess_low", "guess_high")] <- NULL
      newdata$protvar <- prot_var_name
      return(newdata)
    }
    
    newdata[notfound, "guess_low"] <- ifelse(
      (curvals[notfound] < lvl) & is_protective,
      newdata[notfound, var_name],
      newdata[notfound, "guess_low"]
    )
    newdata[notfound, "guess_high"] <- ifelse(
      (curvals[notfound] > lvl) & is_protective,
      newdata[notfound, var_name],
      newdata[notfound, "guess_high"]
    )
  }
}
