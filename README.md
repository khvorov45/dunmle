
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sclr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/sclr)](https://cran.r-project.org/package=sclr)
[![Travis build
status](https://travis-ci.org/khvorov45/sclr.svg?branch=master)](https://travis-ci.org/khvorov45/sclr)
[![codecov](https://codecov.io/gh/khvorov45/sclr/branch/master/graph/badge.svg)](https://codecov.io/gh/khvorov45/sclr)
<!-- badges: end -->

The goal of sclr is to fit the scaled logit model from Dunning (2006)
using the maximum likelihood method.

## Installation

Install the [CRAN](https://CRAN.R-project.org) version with

``` r
install.packages("sclr")
```

Or the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("khvorov45/sclr")
```

## Model

For log likelihood, scores and second derivatives see `vignette("Math",
"sclr")`. Documentation of the main fitting function `?sclr` has details
on how the model is fit.

## Example

Usage is similar to other model fitting functions like `lm`.

``` r
library(sclr)
fit <- sclr(status ~ logHI, one_titre_data)
summary(fit)
#> Call: status ~ logHI
#> 
#> Parameter estimates
#>     lambda     beta_0 beta_logHI 
#>  0.5036743 -5.0798533  2.0467007 
#> 
#> 95% confidence intervals
#>                 2.5 %     97.5 %
#> lambda      0.4790111  0.5283374
#> beta_0     -5.9537676 -4.2059390
#> beta_logHI  1.7561249  2.3372765
```

For more details see `vignette("Usage", "sclr")`.

## References

Dunning AJ (2006). “A model for immunological correlates of protection.”
Statistics in Medicine, 25(9), 1485-1497. doi: 10.1002/sim.2282.
