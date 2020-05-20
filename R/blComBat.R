#' Experimental method for harmonization while accounting for subnetworks
#'
#' Global ComBat followed by within-block ComBat
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
blComBat =  function(x, bat, eb = TRUE, mod = NULL, to.corr = TRUE,
                     out.pd = FALSE, fisher = TRUE) {
  N <- dim(x)[3] # store number of obs
  dnames <- dimnames(x)[[1]]
  bat <- droplevels(bat)

  if (to.corr) {x[] <- array(apply(x, 3, cov2cor), dim(x))}

  # perform fcComBat globally
  glcom_out <- fcComBat(x, bat, eb = TRUE, mod = mod, to.corr = to.corr,
                        out.pd = FALSE, fisher = fisher)
  gl_out <- glcom_out$dat.out
  out <- gl_out

  subcom_out <- list()
  for (d in unique(dnames)) {
    subcom_out[[d]] <- fcComBat(out[dnames == d, dnames == d,], bat,
                                mod = mod, to.corr = FALSE, fisher = fisher)
    out[dnames == d, dnames == d,] <- subcom_out[[d]]$dat.out
  }

  if (out.pd) {out <- array(apply(out, 3, function(x)
    as.matrix(nearPD(x, corr = TRUE)$mat)), dim(out))}

  dimnames(out) <- dnames

  list(dat.out = out, global.combat.out = glcom_out,
       sub.combat.out = subcom_out)
}
