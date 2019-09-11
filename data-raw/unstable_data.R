# Generation of unstable data
# Arseniy Khvorov
# Created 2019/09/11
# Last edit 2019/09/11

bad_seeds <- c(
  20190986, 20191110, 20191406, 20191452, 20191586, 20191610, 20191693
)
load("data-raw/sim_args.Rdata")
unstable_data <- lapply(bad_seeds, function(one_seed) {
  do.call(rsimkhv::sim_scaled_logit, c(sim_args, simseed = one_seed))
})
usethis::use_data(unstable_data, internal = TRUE, overwrite = TRUE)
