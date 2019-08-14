# Datasets included with sclr

#' Simulated antibody titre data with one relevant titre.
#'
#' A simulated dataset containing 6000 independent observations on 
#' antibody titres and the corresponding infection status. The data
#' was simulated to resemble real influenza infection and 
#' haemagluttinin titre data.
#' 
#' @format A data frame with 6000 observations and 3 variables:
#' \describe{
#'   \item{HI}{haemogluttinin-inhibiting titre. 
#'   True simulated titre on a log scale.}
#'   \item{HIcens}{haemogluttinin-inhibiting censored (observed) titre. 
#'   The titre value on a log scale that would be observed 
#'   in a real dataset with a 
#'   typical HI assay.}
#'   \item{status}{influenza infection status. 1 - infected. 0 - not infected}
#' }
#' 
#' @section Model:
#' 
#' The model behind the simulation was
#' 
#' \deqn{\lambda * (1 - f(\beta_0 + \beta_1 * HI))}
#' 
#' Where 
#' 
#' \itemize{
#'   \item \eqn{f} - Inverse logit function
#'   \item \eqn{\lambda} = 0.225
#'   \item \eqn{\beta_0} = -10
#'   \item \eqn{\beta_1} = 3
#' }
#' 
#'
"sclronetitre"

#' Simulated antibody titre data with two relevant titres.
#'
#' A simulated dataset containing 6000 independent observations on 
#' antibody titres and the corresponding infection status. The data
#' was simulated to resemble real influenza infection and 
#' haemagluttinin + neuraminidase titre data.
#' 
#' @format A data frame with 6000 observations and 5 variables:
#' \describe{
#'   \item{HI}{haemogluttinin-inhibiting titre. 
#'   True simulated titre on a log scale.}
#'   \item{HIcens}{haemogluttinin-inhibiting censored (observed) titre. 
#'   The titre value that would be observed in a real dataset with a 
#'   typical HI assay.}
#'   \item{NI}{neuraminidase-inhibiting titre. 
#'   True simulated titre on a log scale.}
#'   \item{NIcens}{haemogluttinin-inhibiting censored (observed) titre. 
#'   The titre value that would be observed in a real dataset with a 
#'   typical HI assay.}
#'   \item{status}{influenza infection status. 1 - infected. 0 - not infected}
#' }
#' 
#' @section Model:
#' 
#' The model behind the simulation was
#' 
#' \deqn{\lambda * (1 - f(\beta_0 + \beta_1 * HI + \beta_2 * NI))}
#' 
#' Where 
#' 
#' \itemize{
#'   \item \eqn{f} - Inverse logit function
#'   \item \eqn{\lambda} = 0.225
#'   \item \eqn{\beta_0} = -10
#'   \item \eqn{\beta_1} = 3
#'   \item \eqn{\beta_2} = 2.5
#' }
#'
"sclrtwotitre"
