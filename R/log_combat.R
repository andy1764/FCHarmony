#' Harmonization using ComBat on log-transformed matrices
#'
#' @param x
#' @param bat
#' @param mod
#' @param covbat.var
#' @param mean.only
#' @param score.eb
#' @param pc.sym
#' @param input.pd
#'
#' @return
#'
#' @import CovBat
#' @export
#'
#' @examples
log_combat <- function(x, # array of fc matrices, roi x roi x nsubj
                       bat, # vector of batch numbers
                       mod = NULL, # model matrix
                       eb = TRUE,
                       input.pd = FALSE# force positive semi-definiteness in input array
) {
  N <- dim(x)[3]
  dnames <- dimnames(x)
  bat <- droplevels(bat)

  if (input.pd) {
    x <- array(apply(dat, 3, function(x) {as.numeric(nearPD(x)$mat)}), dim(dat))
  }
  logx <- array(apply(x, 3, logm_eig), dim(x))

  # collapse into vectors
  v_vec <- t(apply(logx, 3, function(x) c(x[lower.tri(x, diag = TRUE)])))

  # ComBat to remove site effect in scores
  scores_com <- combat(t(v_vec), bat, mod = mod, eb = eb)
  est_combat <- t(scores_com$dat.combat)

  out <- array(0, dim = dim(x))
  est_mat <- array(0, dim = dim(x))
  for (i in 1:N) {
    est_mat[,,i][lower.tri(est_mat[,,i], diag = TRUE)] <- est_combat[i,]
    est_mat[,,i] <- est_mat[,,i] + t(est_mat[,,i])
    diag(est_mat[,,i]) <- diag(est_mat[,,i])/2
    out[,,i] <- expm_eig(est_mat[,,i])
  }

  dimnames(out) <- dnames

  list(dat.out=out,
       combat.out = scores_com)
}
