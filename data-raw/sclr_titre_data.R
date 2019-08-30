# Data is simulated using functions from rhanahkhv package. 
#   github.com/khvorov45/rhanamkhv
# 
# Arseniy Khvorov
# Created 2019/08/30
# Last edit 2019/08/30

parvals_one <- jsonlite::fromJSON(system.file(
  "parameter_sets", "one_titre.json", package = "rhanamkhv", mustWork = TRUE
))
parvals_two <- jsonlite::fromJSON(system.file(
  "parameter_sets", "two_titre.json", package = "rhanamkhv", mustWork = TRUE
))

clean_data <- function(dat) {
  dat_clean <- dat %>%
    dplyr::select(logHI, logHIcens, logHImid, logNI, logNIcens, status)
  if (sum(is.na(dat_clean) > 0)) stop("unexpected missing values")
  return(dat_clean)
}

sclr_one_titre_data <- do.call(rhanamkhv::sim_hanam, parvals_one) %>% 
  clean_data()
sclr_two_titre_data <- do.call(rhanamkhv::sim_hanam, parvals_two) %>% 
  clean_data()

usethis::use_data(sclr_one_titre_data, sclr_two_titre_data, overwrite = TRUE)
