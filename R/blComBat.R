#' Method for harmonization while accounting for subnetworks
#'
#' Global ComBat followed by within-block ComBat
#'
#' @param x *p x p x n* covariance or correlation matrices where *p* is the number
#'   of ROIs and *n* is the number of subjects.
#' @param bat Factor (or object coercible by \link[base]{as.factor} to a
#'    factor) of length *n* designating batch IDs.
#' @param mod Optional design matrix of covariates to preserve, usually from
#'   the output of \link[stats]{model.matrix}.
#' @param eb If `TRUE``, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization.
#' @param blocks Blocks indicating submatrices that are harmonized in the
#'   second-stage ComBat. Defaults to input dimension names.
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
blComBat =  function(x, bat, mod = NULL, blocks = dimnames(x)[[1]], eb = TRUE,
                     to.corr = TRUE,
                     out.pd = FALSE, fisher = TRUE) {
  N <- dim(x)[3] # store number of obs
  dnames <- dimnames(x)
  bat <- droplevels(bat)

  if (to.corr) {x[] <- array(apply(x, 3, cov2cor), dim(x))}

  # perform fcComBat globally
  glcom_out <- fcComBat(x, bat, eb = TRUE, mod = mod, to.corr = to.corr,
                        out.pd = FALSE, fisher = fisher)
  gl_out <- glcom_out$dat.out
  out <- gl_out

  subcom_out <- list()
  for (b in blocks) {
    block <- blocks == b
    # if number of ROIs in subnetwork is greater than zero
    if (sum(block) > 1) {
      subcom_out[[b]] <- fcComBat(out[block, block,],
                                  bat, mod = mod, to.corr = FALSE, fisher = fisher)
      out[block, block,] <- subcom_out[[b]]$dat.out
    }
  }

  if (out.pd) {out <- array(apply(out, 3, function(x)
    as.matrix(nearPD(x, corr = TRUE)$mat)), dim(out))}

  dimnames(out) <- dnames

  list(dat.out = out, global.combat.out = glcom_out,
       sub.combat.out = subcom_out)
}
