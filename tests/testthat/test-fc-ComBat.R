test_that("FC ComBat can harmonize basic matrices", {
  library(CovBat)
  data(CARDIA)
  bat <- c(rep(0, 15), rep(1, 15))
  out <- cpcharmony(CARDIA, bat, cpc.method = "", err.method = "fc-ComBat")

  # comparison to vectorized
  err_vec <- t(apply(CARDIA, 3, function(x) c(x[lower.tri(x)])))
  err_com <- t(combat_modded(t(err_vec), bat)$dat.combat)

  expect_equal(c(out$dat.out[1, 2, 1], out$dat.out[1, 2, 16]),
               as.numeric(c(err_com[1, 1], err_com[16, 1])))
})
