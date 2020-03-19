#' Harmonization using principal components of log-transformed matrices
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
log_covbat <- function(x, # array of fc matrices, roi x roi x nsubj
                       bat, # vector of batch numbers
                       mod = NULL, # model matrix
                       covbat.var = 0.90, # percent of variation explained that determines number of scores to correct
                       mean.only = FALSE, # scale parameter in initial ComBat, works better with scaling
                       score.eb = FALSE, # empirical Bayes for scores
                       pc.sym = TRUE, # perform PCA only on lower triangular
                       input.pd = TRUE # force positive semi-definiteness in input array
) {
  n <- dim(x)[3] # store number of obs
  logx <- array(dim = dim(x))
  for (i in 1:n) {
    if (input.pd) {x[,,i] <- as.numeric(nearPD(x[,,i])$mat)}
    logx[,,i] <- logm(x[,,i])
  }

  # collapse into vectors
  if (pc.sym) {
    v_vec <- t(apply(logx, 3, function(x) c(x[lower.tri(x, diag = TRUE)])))
  } else {
    v_vec <- t(apply(logx, 3, c))
  }

  v_pc <- prcomp(v_vec) # PC on ComBat-adjusted data

  npc <- which(cumsum(v_pc$sdev^2/sum(v_pc$sdev^2)) > covbat.var)[1]
  # print(npc)
  scores <- v_pc$x[,1:npc]

  # ComBat to remove site effect in scores
  scores_com <- combat_modded(t(scores), bat, eb = score.eb, mod = mod)

  # print(apply(t(scores_com$dat.combat)[bat == 1,], 2, mean))
  # print(apply(t(scores_com$dat.combat)[bat == 1,], 2, var))

  full_scores <- v_pc$x
  full_scores[,1:npc] <- t(scores_com$dat.combat)
  est_covbat <- t(t(full_scores %*% t(v_pc$rotation)) + v_pc$center) # get vectorized CovBat corrected

  out_covbat <- array(0, dim = dim(x))
  est_mat <- array(0, dim = dim(x))
  for (i in 1:n) {
    if (pc.sym) {
      est_mat[,,i][lower.tri(est_mat[,,i], diag = TRUE)] <- est_covbat[i,]
      est_mat[,,i] <- est_mat[,,i] + t(est_mat[,,i]) - diag(diag(est_mat[,,i]))
    } else {
      est_mat[,,i] <- matrix(est_covbat[i,], dim(x[,,1]))
    }
    out_covbat[,,i] <- expm(est_mat[,,i])
  }

  return(list(dat.covbat=out_covbat,
              dat.log.covbat=est_mat,
              dat.input=x,
              dat.log.input=logx,
              scores.gamma=scores_com$gamma.hat, scores.delta=scores_com$delta.hat,
              npc=npc,
              x.pc = v_pc, com.scores = scores_com))
}
