# Method first implemented by Yu et al. (2018), DOI: 10.1002/hbm.24241
# Apply ComBat to Fisher-transformed lower triangular elements
fcComBat =  function(x, bat, mod = NULL, to.corr = TRUE) {
  n <- dim(x)[3] # store number of obs

  if (to.corr) {x <- array(apply(x, 3, cov2cor), dim(x))}
  if(x[,,1] != cov2cor(x[,,1])) {
    stop("fc-ComBat only takes correlation matrix inputs")
  }

  vec <- atanh(t(apply(x, 3, function(m) c(m[lower.tri(m)]))))
  com_out <- combat_modded(t(vec), bat, mod = mod)
  com_dat <- tanh(t(com_out$dat.combat))
  out <- array(0, dim = dim(x))
  for (i in 1:N) {
    out[,,i][lower.tri(out[,,i])] <- com_dat[i,]
    out[,,i] <- out[,,i] + t(out[,,i])
    diag(out[,,i]) <- diag(x[,,i])
  }

  list(dat.combat=out, combat.out = com_out)
}
