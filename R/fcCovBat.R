# Method first implemented by Yu et al. (2018), DOI: 10.1002/hbm.24241
# Apply ComBat to Fisher-transformed lower triangular elements

#' Functional connectivity CovBat
#'
#' @param x
#' @param bat
#' @param mod
#' @param to.corr
#' @param out.pd
#'
#' @return
#' @export
#'
#' @examples
fcCovBat =  function(x, bat, eb = TRUE, mod = NULL, to.corr = TRUE,
                     out.pd = FALSE, fisher = TRUE) {
  N <- dim(x)[3]
  dnames <- dimnames(x)
  bat <- droplevels(bat)

  if (to.corr) {x <- array(apply(x, 3, cov2cor), dim(x))}

  if (fisher) {
    vec <- atanh(t(apply(x, 3, function(m) c(m[lower.tri(m)]))))
    cov_out <- covbat(t(vec), bat, mod = mod, eb = eb)
    cov_out <- tanh(t(cov_out$dat.covbat))
  } else {
    vec <- t(apply(x, 3, function(m) c(m[lower.tri(m)])))
    cov_out <- covbat(t(vec), bat, mod = mod, eb = eb)
    com_dat <- t(cov_out$dat.covbat)
  }

  out <- array(0, dim = dim(x))
  for (i in 1:N) {
    out[,,i][lower.tri(out[,,i])] <- com_dat[i,]
    out[,,i] <- out[,,i] + t(out[,,i])
    diag(out[,,i]) <- diag(x[,,i])
  }

  if (out.pd) {out <- array(apply(out, 3, function(x)
    as.matrix(nearPD(x, corr = TRUE)$mat)), dim(out))}

  dimnames(out) <- dnames

  list(dat.out = out, covbat.out = cov_out)
}
