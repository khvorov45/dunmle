# Data is simulated using functions from rhanamkhv package. 
#   github.com/khvorov45/rhanamkhv
# 
# Arseniy Khvorov
# Created 2019/08/30
# Last edit 2019/09/04

library(dplyr)

parvals_one <- jsonlite::fromJSON("one_titre.json")
parvals_one$seed <- 20190830
parvals_two <- jsonlite::fromJSON("two_titre.json")
parvals_two$seed <- 20190830

clean_data <- function(dat, n_titre) {
  dat_clean <- dat %>%
    dplyr::select(logHI, logHIcens, logHImid, logNI, logNIcens, status)
  if (n_titre == 1) {
    dat_clean <- dat_clean %>%
      dplyr::select(-logNI, -logNIcens)
  }
  if (sum(is.na(dat_clean) > 0)) stop("unexpected missing values")
  return(dat_clean)
}

sclr_one_titre_data <- do.call(rhanamkhv::sim_hanam, parvals_one) %>% 
  clean_data(1)
sclr_two_titre_data <- do.call(rhanamkhv::sim_hanam, parvals_two) %>% 
  clean_data(2)

usethis::use_data(sclr_one_titre_data, sclr_two_titre_data, overwrite = TRUE)
