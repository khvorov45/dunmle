# General utility functions

# Returns parameter names based on the model matrix
get_par_names <- function(x) {
  par_names <- c("beta_0", paste0("beta_", colnames(x)[-1]))
  return(par_names)
}

get_betas_only <- function(pars_mat) {
  pars_betas <- pars_mat[grepl("beta", rownames(pars_mat)), ]
  return(pars_betas)
}
