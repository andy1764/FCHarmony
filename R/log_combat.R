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
#' @importFrom expm expm
#' @export
#'
#' @examples
log_combat <- function(x, # array of fc matrices, roi x roi x nsubj
                       bat, # vector of batch numbers
                       mod = NULL, # model matrix
                       input.pd = FALSE # force positive semi-definiteness in input array
) {
  n <- dim(x)[3] # store number of obs
  if (input.pd) {
    x <- array(apply(dat, 3, function(x) {as.numeric(nearPD(x)$mat)}), dim(dat))
  }
  logx <- array(apply(x, 3, logm_eig), dim(x))

  # collapse into vectors
  v_vec <- t(apply(logx, 3, function(x) c(x[lower.tri(x, diag = TRUE)])))

  # ComBat to remove site effect in scores
  scores_com <- combat_modded(t(v_vec), bat, mod = mod)
  est_combat <- t(scores_com$dat.combat)

  out_combat <- array(0, dim = dim(x))
  est_mat <- array(0, dim = dim(x))
  for (i in 1:n) {
    est_mat[,,i][lower.tri(est_mat[,,i], diag = TRUE)] <- est_combat[i,]
    est_mat[,,i] <- est_mat[,,i] + t(est_mat[,,i])
    diag(est_mat[,,i]) <- diag(est_mat[,,i])/2
    out_combat[,,i] <- expm(est_mat[,,i], "R_Eigen")
  }

  list(dat.out=out_combat,
       covbat.out = scores_com)
}
