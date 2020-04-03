## Code for generating BLSA/CARDIA dataset
# Consists of multivariate normal draws from all matrices in dataset
library(Matrix)
library(mvtnorm)
load("~/Documents/GitHub/Covariance-Harmonization/BLSA CARDIA Example/blsa_cardia_all.Rdata")
roi_labs <- read.csv("~/Documents/GitHub/Covariance-Harmonization/BLSA CARDIA Example/power264rois.csv")[,2]

# subset of default mode network ROIs, 50 subjects from BLSA, CARDIA 1,3,4
set.seed(88888888)
sub <- c(sample(which(all_demo$site == "BLSA"), 50),
         sample(which(all_demo$site == "CARDIA_1"), 50),
         sample(which(all_demo$site == "CARDIA_3"), 50),
         sample(which(all_demo$site == "CARDIA_4"), 50))
N <- length(sub)
rois <- sort(sample(1:264, 50))
roi_labs <- roi_labs[rois]
dmn <- roi_labs == "Default mode"
sub_arr <- all_arr[rois,rois,sub]
BLSA_CARDIA <- sapply(seq_along(sub_arr[1,1,]), function(x)
  as.matrix(nearPD(sub_arr[,,x])$mat), simplify = "array")
BLSA_CARDIA <- sub_arr
dimnames(BLSA_CARDIA) <- list(roi_labs, roi_labs, NULL)
# BLSA_CARDIA <- array(sapply(1:N, function(x)
#   as.matrix(nearPD(cov(rmvnorm(1000, sigma = sub_arr[dmn,dmn,x])))$mat)),
#   dim = c(sum(dmn), sum(dmn), N),
#   dimnames = list(roi_labs[dmn], roi_labs[dmn], NULL))
usethis::use_data(BLSA_CARDIA, overwrite = TRUE)

BLSA_CARDIA_demo <- all_demo[sub, c("site", "age_at_scan", "sex")]
rownames(BLSA_CARDIA_demo) <- NULL

usethis::use_data(BLSA_CARDIA_demo, overwrite = TRUE)
