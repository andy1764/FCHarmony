#' Harmonization using Common Principal Components
#'
#' @param dat an array of covariance matrices of dimension
#' \eqn{(p \times p \times N)}
#' @param bat a vector of batch labels as length \eqn{n} vector
#' @param mod a matrix of covariates to preserve during harmonization
#' @param cpc.k numeric, number of common principal components to harmonize
#' @param method character, method for harmonization, either "likelihood",
#' "log-Combat", or "ComBat", currently only ComBat supported. Default is
#' "ComBat"
#' @param cpc.method character, estimation method for common principal
#' components, one of "stepwise", "Flury", or "CAP". Default is "stepwise".
#' @param cpc.cap.oc logical, only used if \code{cpc.method = "CAP"}, if
#' \code{TRUE} imposes orthogonal constraint when identifiying higher-order PCs.
#' Default is \code{FALSE}
#'
#' @return \code{cpcharmony} returns a list containing the following components:
#' \item{dat.out}{harmonized data as an array of input dimension}
#' \item{cpcs}{CPC analysis output list}
#' \item{cpc.harmony}{harmonization method output}
#'
#' @import CovBatHarmonization
#' @importFrom cpca cpc
#' @export
#'
#' @examples

cpcharmony <- function(dat, bat, cpc.k, mod = NULL, method = "ComBat",
                       cpc.method = "stepwise", cpc.cap.oc = FALSE) {
  if (!(method %in% c("likelihood", "log-ComBat", "ComBat"))) {
    stop("Method not available")
  }
  if (!(method %in% c("likelihood", "log-ComBat", "ComBat"))) {
    stop("Common principal components method not available")
  }

  p <- dim(dat)[2]
  N <- dim(dat)[3]
  dat_err <- dat # store error matrices
  dat_out <- array(0, dim = dim(dat)) # initialize output matrices

  if (cpc.method == "stepwise") {
    cpcs <- cpc(dat, k = cpc.k) # get common principal components
    for (i in 1:N) {
      for (j in 1:cpc.k) {
        dat_err[,,i] = dat_err[,,i] - cpcs$D[j,i] * tcrossprod(cpcs$CPC[,j])
      }
    }
  } else if (cpc.method == "Flury") {
    cpcs <- cpc_Flury(dat)
    for (i in 1:N) {
      for (j in 1:cpc.k) {
        dat_err[,,i] = dat_err[,,i] - cpcs$D[j,i] * tcrossprod(cpcs$CPC[,j])
      }
    }
  } else if (cpc.method == "CAP") {

  }

  if (method == "ComBat") {
    cpc_harmony <- combat_modded(cpcs$D, bat, mod = mod)
    for (i in 1:N) {
      dat_out[,,i] <- dat_err[,,i] # add back error
      for (j in 1:cpc.k) {
        dat_out[,,i] = dat_out[,,i] +
          cpc_harmony$dat.combat[j,i] * tcrossprod(cpcs$CPC[,j])
      }
    }
  }

  list(dat.out = dat_out,
       dat.err = dat_err,
       cpcs = cpcs,
       cpc.harmony = cpc_harmony)
}
