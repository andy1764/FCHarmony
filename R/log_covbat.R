#' Harmonization using CovBat on log-transformed matrices
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
log_covbat <- function(x, # array of fc matrices, roi x roi x nsubj
                       bat, # vector of batch numbers
                       mod = NULL, # model matrix
                       eb = TRUE,
                       covbat.var = 0.90, # percent of variation explained that determines number of scores to correct
                       mean.only = FALSE, # scale parameter in initial ComBat, works better with scaling
                       score.eb = FALSE, # empirical Bayes for scores
                       pc.sym = TRUE, # perform PCA only on lower triangular
                       input.pd = FALSE # force positive semi-definiteness in input array
) {
  n <- dim(x)[3] # store number of obs
  if (input.pd) {
    x <- array(apply(dat, 3, function(x) {as.numeric(nearPD(x)$mat)}), dim(dat))
  }
  logx <- array(apply(x, 3, logm_eig), dim(x))

  # collapse into vectors
  if (pc.sym) {
    v_vec <- t(apply(logx, 3, function(x) c(x[lower.tri(x, diag = TRUE)])))
  } else {
    v_vec <- t(apply(logx, 3, c))
  }

  # ComBat to remove site effect in scores
  scores_com <- covbat(t(v_vec), bat, eb = eb, score.eb = score.eb, mod = mod,
                       percent.var = covbat.var)
  est_covbat <- t(scores_com$dat.covbat)

  out_covbat <- array(0, dim = dim(x))
  est_mat <- array(0, dim = dim(x))
  for (i in 1:n) {
    if (pc.sym) {
      est_mat[,,i][lower.tri(est_mat[,,i], diag = TRUE)] <- est_covbat[i,]
      est_mat[,,i] <- est_mat[,,i] + t(est_mat[,,i])
      diag(est_mat[,,i]) <- diag(est_mat[,,i])/2
    } else {
      est_mat[,,i] <- matrix(est_covbat[i,], dim(x[,,1]))
    }
    out_covbat[,,i] <- expm(est_mat[,,i], "R_Eigen")
  }

  list(dat.out=out_covbat,
       covbat.out = scores_com)
}
