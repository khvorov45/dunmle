# Functions to caclulate the jacobian matrix

# Computes the jacobian matrix
get_jacobian <- function(y, x, pars_mat) {
  
  d2l_dl2 <- get_d2l_dl2(y, x, pars_mat)
  d2l_dldb0 <- get_d2l_dldb0(y, x, pars_mat)
  d2l_dldbt <- get_d2l_dldbt(y, x, pars_mat)
  d2l_db02 <- get_d2l_db02(y, x, pars_mat)
  d2l_db0bt <- get_d2l_db0bt(y, x, pars_mat)
  d2l_dbt2 <- get_d2l_dbt2(y, x, pars_mat)

  ord <- c(
    d2l_dl2, d2l_dldb0, d2l_dldbt, 
    d2l_dldb0, d2l_db02, d2l_db0bt, 
    d2l_dldbt, d2l_db0bt, d2l_dbt2
  )

  jacobian_mat <- matrix(
    ord, nrow = 3,
    dimnames = list(c("dl", "db0", "dbt"), c("dl", "db0", "dbt"))
  )

  return(jacobian_mat)
}

# Computes the second derivative of log-likelihood with respect to beta_titre
get_d2l_dbt2 <- function(y, x, pars_mat) {
  d2l_dbt2 <- (1 - y) * get_d2l_dbt2_f1(y, x, pars_mat) - 
    get_d2l_dbt2_f2(y, x, pars_mat)
  d2l_dbt2 <- sum(d2l_dbt2)
  return(d2l_dbt2)
}

# Computes the first fraction term of d2l_dbt2
get_d2l_dbt2_f1 <- function(y, x, pars_mat) {
  lambda <- pars_mat["lambda", ]
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  f1 <- (x[, 2]^2 * (1 - lambda) * exp(beta_0 + beta_titre * x[, 2])) / 
    get_denom(y, x, pars_mat)^2
  return(f1)
}

# Computes the second fraction term of d2l_dbt2
get_d2l_dbt2_f2 <- function(y, x, pars_mat) {
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  f2 <- (x[, 2]^2 * exp(beta_0 + beta_titre * x[, 2])) / 
    (1 + exp(beta_0 + beta_titre * x[, 2]))^2
  return(f2)
}

# Computes the second derivative of log-likelihood with respect to beta_0 and
# beta_titre
get_d2l_db0bt <- function(y, x, pars_mat) {
  d2l_db0bt <- (1 - y) * get_d2l_db0bt_f1(y, x, pars_mat) - 
    get_d2l_db0bt_f2(y, x, pars_mat)
  d2l_db0bt <- sum(d2l_db0bt)
  return(d2l_db0bt)
}

# Computes the first fraction term of d2l_db0bt
get_d2l_db0bt_f1 <- function(y, x, pars_mat) {
  lambda <- pars_mat["lambda", ]
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  f1 <- (x[, 2] * (1 - lambda) * exp(beta_0 + beta_titre * x[, 2])) / 
    get_denom(y, x, pars_mat)^2
  return(f1)
}

# Computes the second fraction term of d2l_db0bt
get_d2l_db0bt_f2 <- function(y, x, pars_mat) {
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  f2 <- (x[, 2] * exp(beta_0 + beta_titre * x[, 2])) / 
    (1 + exp(beta_0 + beta_titre * x[, 2]))^2
  return(f2)
}

# Computes the second derivative of log-likelihood with respect to beta_0
get_d2l_db02 <- function(y, x, pars_mat) {
  d2l_db02 <- (1 - y) * get_d2l_db02_f1(y, x, pars_mat) - 
    get_d2l_db02_f2(y, x, pars_mat)
  d2l_db02 <- sum(d2l_db02)
  return(d2l_db02)
}

# Computes the first fraction term of d2l_db02
get_d2l_db02_f1 <- function(y, x, pars_mat) {
  lambda <- pars_mat["lambda", ]
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  f1 <- (1 - lambda) * exp(beta_0 + beta_titre * x[, 2]) / 
    get_denom(y, x, pars_mat)^2
  return(f1)
}

# Computes the second fraction term of d2l_db02
get_d2l_db02_f2 <- function(y, x, pars_mat) {
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  f2 <- exp(beta_0 + beta_titre * x[, 2]) / 
    (1 + exp(beta_0 + beta_titre * x[, 2]))^2
  return(f2)
}

# Computes the second derivative of log-likelihood with respect to lambda and
# beta_titre
get_d2l_dldbt <- function(y, x, pars_mat) {
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  d2l_dldbt <- (1 - y) * x[, 2] * exp(beta_0 + beta_titre * x[, 2]) / 
    get_denom(y, x, pars_mat)^2
  d2l_dldbt <- sum(d2l_dldbt)
  return(d2l_dldbt)
}

# Computes the second derivative of log-likelihood with respect to lambda and
# beta_0
get_d2l_dldb0 <- function(y, x, pars_mat) {
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  d2l_dldb0 <- (1 - y) * exp(beta_0 + beta_titre * x[, 2]) / 
    get_denom(y, x, pars_mat)^2
  d2l_dldb0 <- sum(d2l_dldb0)
  return(d2l_dldb0)
}

# Computes the second derivative of log-likelihood with respect to lambda
get_d2l_dl2 <- function(y, x, pars_mat) {
  lambda <- pars_mat["lambda", ]
  d2l_dl2 <- -y / lambda^2 - (1 - y) * (1 / get_denom(y, x, pars_mat)^2)
  d2l_dl2 <- sum(d2l_dl2)
  return(d2l_dl2)
}

# Computes the commonly occuring denominator term
get_denom <- function(y, x, pars_mat) {
  lambda <- pars_mat["lambda", ]
  beta_0 <- pars_mat["beta_0", ]
  beta_titre <- pars_mat["beta_titre", ]
  denom <- 1 + exp(beta_0 + beta_titre * x[, 2]) - lambda
  return(denom)
}