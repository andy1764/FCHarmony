## code for generating blsa_cardia_arr and blsa_cardia_demo datasets
load("blsa_cardia.Rdata")
sites <- names(demo)[1:5]
m <- 25
n <- (length(sites)-1)*m + 4
blsa_cardia_arr <- array(as.numeric(c(unlist(arrays[[sites[1]]][,,1:m]),
                                      unlist(arrays[[sites[2]]][,,1:m]),
                                      unlist(arrays[[sites[3]]][,,1:4]),
                                      unlist(arrays[[sites[4]]][,,1:m]),
                                      unlist(arrays[[sites[5]]][,,1:m]))),
                         c(264, 264, n))
usethis::use_data(blsa_cardia_arr, overwrite = TRUE)

blsa_cardia_demo <- NULL
for (site in sites) {
  if (site == "CARDIA_3") {
    blsa_cardia_demo <- rbind(blsa_cardia_demo, demo[[site]][1:4, ])
  } else {
    blsa_cardia_demo <- rbind(blsa_cardia_demo, demo[[site]][1:m, ])
  }
}
usethis::use_data(blsa_cardia_demo, overwrite = TRUE)
