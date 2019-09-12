# Datasets included with sclr
# Arseniy Khvorov
# Created 2019/08/14
# Last edit 2019/09/13

#' Simulated one-titre antibody data
#'
#' A simulated dataset containing 5000 independent observations on antibody
#' titres and the corresponding infection status. The data was simulated to
#' resemble real influenza infection and haemagglutinin titre data.
#'
#' @format A data frame with 5000 observations and 4 variables: 
#' \describe{
#'   \item{logHI}{haemagglutinin-inhibiting (HI) titre. 
#'   True simulated titre on a log scale.} 
#'   \item{logHIcens}{HI censored (observed) titre.
#'   The titre value on a log scale that would be observed in a real dataset
#'   with a typical HI assay.}
#'   \item{logHImid}{Midpoint of the interval (on a log scale)
#'   to which observed HI values are censored.}
#'   \item{status}{influenza infection status. 1 -
#'   infected. 0 - not infected}
#'  }
#'
#' @section Model:
#'
#'   The model behind the simulation was
#'
#'   \deqn{\lambda * (1 - f(\beta_0 + \beta_1 * HI))}
#'
#'   Where
#'
#'   \itemize{ \item \eqn{f} - Inverse logit function \item \eqn{\lambda} = 0.5
#'   \item \eqn{\beta_0} = -5 \item \eqn{\beta_1} = 2 }
#'   
"one_titre_data"

#' Simulated two-titre antibody data
#'
#' A simulated dataset containing 5000 independent observations on antibody
#' titres and the corresponding infection status. The data was simulated to
#' resemble real influenza infection and haemagglutinin + neuraminidase titre
#' data.
#'
#' @format A data frame with 5000 observations and 6 variables: 
#' \describe{
#'   \item{logHI}{haemagglutinin-inhibiting (HI) titre. 
#'   True simulated titre on a log scale.} 
#'   \item{logHIcens}{HI censored (observed) titre.
#'   The titre value on a log scale that would be observed in a real dataset
#'   with a typical HI assay.}
#'   \item{logHImid}{Midpoint of the interval (on a log scale)
#'   to which observed HI values are censored.}
#'   \item{logNI}{neuraminidase-inhibiting titre. True simulated titre on a
#'   log scale.} 
#'   \item{logNIcens}{haemagglutinin-inhibiting censored (observed)
#'   titre. The titre value that would be observed in a real dataset with a
#'   typical NI assay.} 
#'   \item{status}{influenza infection status. 1 - infected.
#'   0 - not infected} 
#'  }
#'
#' @section Model:
#'
#'   The model behind the simulation was
#'
#'   \deqn{\lambda * (1 - f(\beta_0 + \beta_1 * HI + \beta_2 * NI))}
#'
#'   Where
#'
#'   \itemize{ \item \eqn{f} - Inverse logit function \item \eqn{\lambda} = 0.5
#'   \item \eqn{\beta_0} = -7.5 \item \eqn{\beta_1} = 2 \item \eqn{\beta_2} = 2
#'   }
#'   
"two_titre_data"
