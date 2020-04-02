## Code for generating CARDIA dataset
# Consists of multivariate normal draws from CARDIA 1 and 4 FC matrices
library(Matrix)
library(mvtnorm)
load("~/Documents/GitHub/Covariance-Harmonization/BLSA CARDIA Example/blsa_cardia_all.Rdata")
roi_labs <- read.csv("~/Documents/GitHub/Covariance-Harmonization/BLSA CARDIA Example/power264rois.csv")[,2]

## get first 15 subjects from CARDIA 1 and 4
CARDIA_raw <- all_arr[,,c(358:372, 604:618)]

set.seed(88888888)
n <- 100
CARDIA <- array(sapply(1:30, function(x)
  as.numeric(nearPD(cov(rmvnorm(1000, sigma = all_arr[,,x])))$mat)),
  dim = dim(CARDIA_raw),
  dimnames = list(roi_labs, roi_labs, NULL)
  )
usethis::use_data(CARDIA, overwrite = TRUE)

# subset of 10 default mode network ROIs
CARDIA_small <- array(sapply(1:30, function(x)
  as.numeric(nearPD(cov(rmvnorm(1000, sigma = all_arr[74:83,74:83,x])))$mat)),
  dim = c(10, 10, 30))
usethis::use_data(CARDIA_small, overwrite = TRUE)

CARDIA_demo <- all_demo[c(358:372, 604:618), c("site", "age_at_scan", "sex")]
CARDIA_demo$age_at_scan <- round(CARDIA_demo$age_at_scan - 50, 0)
rownames(CARDIA_demo) <- NULL

usethis::use_data(CARDIA_demo, overwrite = TRUE)
