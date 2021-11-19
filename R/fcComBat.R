# Method first implemented by Yu et al. (2018), DOI: 10.1002/hbm.24241
# Apply ComBat to Fisher-transformed lower triangular elements

#' Functional connectivity ComBat
#'
#' @param x *p x p x n* covariance or correlation matrices where *p* is the number
#'   of ROIs and *n* is the number of subjects.
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
#' @export
#'
#' @examples
fcComBat = function(x, bat, eb = TRUE, mod = NULL, to.corr = TRUE,
                     out.pd = FALSE, fisher = TRUE) {
  N <- dim(x)[3] # store number of obs
  dnames <- dimnames(x)
  bat <- droplevels(bat)

  if (to.corr) {x <- array(apply(x, 3, cov2cor), dim(x))}

  if (fisher) {
    vec <- atanh(t(apply(x, 3, function(m) c(m[lower.tri(m)]))))
    com_out <- CovBat::combat(t(vec), bat, mod = mod, eb = eb)
    com_dat <- tanh(t(com_out$dat.combat))
  } else {
    vec <- t(apply(x, 3, function(m) c(m[lower.tri(m)])))
    com_out <- CovBat::combat(t(vec), bat, mod = mod, eb = eb)
    com_dat <- t(com_out$dat.combat)
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

  list(dat.out = out, combat.out = com_out)
}
