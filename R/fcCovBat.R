#' Functional connectivity CovBat
#'
#' Applies CovBat to vectorized covariance or correlation matrices.
#'
#' @param x *p x p x n* covariance or correlation matrix where *p* is the number
#'   of ROIs and *n* is the number of features
#' @param bat Factor (or object coercible by \link[base]{as.factor} to a
#'    factor) of length *n* designating batch IDs.
#' @param mod Optional design matrix of covariates to preserve, usually from
#'   the output of \link[stats]{model.matrix}.
#' @param eb If `TRUE``, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization.
#' @param to.corr If `TRUE`, uses \link[stats]{cov2cor} to convert input
#'   matrices into correlation matrices
#' @param out.pd Whether input should be forced to be positive definite using
#'   \link[Matrix]{nearPD}. This step is unnecessary for many downstream
#'   network analyses so defaults to `FALSE`.
#' @param fisher Whether to Fisher-transform the off-diagonal elements before
#'   applying CovBat, highly recommended that this be set to `TRUE`.
#'
#' @return
#' @import CovBat
#' @export
#'
#' @examples
fcCovBat =  function(x, bat, mod = NULL, eb = TRUE, to.corr = TRUE,
                     out.pd = FALSE, fisher = TRUE) {
  N <- dim(x)[3]
  dnames <- dimnames(x)
  bat <- as.factor(bat)
  bat <- droplevels(bat)

  if (to.corr) {x <- array(apply(x, 3, cov2cor), dim(x))}

  if (fisher) {
    vec <- atanh(t(apply(x, 3, function(m) c(m[lower.tri(m)]))))
    cov_out <- covbat(t(vec), bat, mod = mod, eb = eb)
    cov_dat <- tanh(t(cov_out$dat.covbat))
  } else {
    vec <- t(apply(x, 3, function(m) c(m[lower.tri(m)])))
    cov_out <- covbat(t(vec), bat, mod = mod, eb = eb)
    cov_dat <- t(cov_out$dat.covbat)
  }

  out <- array(0, dim = dim(x))
  for (i in 1:N) {
    out[,,i][lower.tri(out[,,i])] <- cov_dat[i,]
    out[,,i] <- out[,,i] + t(out[,,i])
    diag(out[,,i]) <- diag(x[,,i])
  }

  if (out.pd) {out <- array(apply(out, 3, function(x)
    as.matrix(nearPD(x, corr = TRUE)$mat)), dim(out))}

  dimnames(out) <- dnames

  list(dat.out = out, covbat.out = cov_out)
}
