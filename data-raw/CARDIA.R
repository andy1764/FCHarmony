## Code for generating CARDIA dataset
# Consists of multivariate normal draws from CARDIA 1 and 4 FC matrices

library(Matrix)
library(mvtnorm)
setwd("~/Documents/GitHub/Covariance-Harmonization/BLSA CARDIA Example")
load("blsa_cardia_all.Rdata")

## get first 15 subjects from CARDIA 1 and 4
CARDIA_raw <- all_arr[,,c(358:372, 604:618)]

set.seed(88888888)
n <- 100
CARDIA <- array(sapply(1:30, function(x)
  as.numeric(nearPD(cov(rmvnorm(1000, sigma = all_arr[,,x])))$mat)),
  dim = dim(CARDIA_raw))

setwd("~/Documents/GitHub/FCharmony")
usethis::use_data(CARDIA)
