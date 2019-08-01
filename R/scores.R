# Functions to compute the score matrix

# Computes the score matrix
get_scores <- function(y, x, pars) {
  dl <- get_dl(y, x, pars)
  db0 <- get_db0(y, x, pars)
  dbt <- get_dbt(y, x, pars)
  scores <- matrix(c(dl, db0, dbt))
  rownames(scores) <- c("dl", "db0", "dbt")
  return(scores)
}

# Computes the first derivative of log-likelihood with respect to lambda
get_dl <- function(y, x, pars) {
  lambda <- pars["lambda", ]
  beta_0 <- pars["beta_0", ]
  beta_titre <- pars[3, ]
  dl <- y / lambda - (1 - y) / (1 + exp(beta_0 + beta_titre * x[, 2]) - lambda)
  dl <- sum(dl)
  return(dl)
}

# Computes the first derivative of log-likelihood with respect to beta_0
get_db0 <- function(y, x, pars) {
  db0 <- (1 - y) * get_db0_f1(y, x, pars) - get_db0_f2(y, x, pars)
  db0 <- sum(db0)
  return(db0)
}

# Computes the first fraction term of db0
get_db0_f1 <- function(y, x, pars) {
  lambda <- pars["lambda", ]
  beta_0 <- pars["beta_0", ]
  beta_titre <- pars[3, ]
  f1 <- exp(beta_0 + beta_titre * x[, 2]) / 
    (1 + exp(beta_0 + beta_titre * x[, 2]) - lambda)
  return(f1)
}

# Computes the second fraction term of db0
get_db0_f2 <- function(y, x, pars) {
  beta_0 <- pars["beta_0", ]
  beta_titre <- pars[3, ]
  f2 <- exp(beta_0 + beta_titre * x[, 2]) / 
    (1 + exp(beta_0 + beta_titre * x[, 2]))
  return(f2)
}

# Computes the first derivative of log-likelihood with respect to beta_titre
get_dbt <- function(y, x, pars) {
  dbt <- (1 - y) * get_dbt_f1(y, x, pars) - get_dbt_f2(y, x, pars)
  dbt <- sum(dbt)
  return(dbt)
}

# Computes the first fraction term of dbt
get_dbt_f1 <- function(y, x, pars) {
  lambda <- pars["lambda", ]
  beta_0 <- pars["beta_0", ]
  beta_titre <- pars[3, ]
  f1 <- x[, 2] * exp(beta_0 + beta_titre * x[, 2]) / 
    (1 + exp(beta_0 + beta_titre * x[, 2]) - lambda)
  return(f1)
}

# Computes the second fraction term of dbt
get_dbt_f2 <- function(y, x, pars) {
  beta_0 <- pars["beta_0", ]
  beta_titre <- pars[3, ]
  f2 <- x[, 2] * exp(beta_0 + beta_titre * x[, 2]) / 
    (1 + exp(beta_0 + beta_titre * x[, 2]))
  return(f2)
}
