# sclr 0.1.0-dev

- Moved the parameter values from rhanamkhv to data-raw folder.

- Added jsonlite to Suggests list.

- Added the ability to return parameter names that are more conventional 
(e.g. "(Intercept)" instead of "beta_0"). See `conventional_names` argument
in `?sclr`.

# sclr 0.1.0

First release.

## Main features

- Fits the scaled logit model using the Newton-Raphson method.

- Supports the predict method for the expected value of the linear beta X part of the model.

- Can look for covariate values corresponding to a particular protection level.
