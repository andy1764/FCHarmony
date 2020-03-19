#' Harmonization using Common Principal Components
#'
#' @param dat An array of covariance matrices of dimension
#' \eqn{(p \times p \times N)}
#' @param bat A vector of batch labels as length \eqn{n} vector
#' @param mod A matrix of covariates to preserve during harmonization
#' @param cpc.k Numeric, number of common principal components to harmonize
#' @param method Method for harmonization. One of "likelihood", "log-ComBat",
#'   "ComBat". Currently only ComBat supported.
#' @param cpc.method Estimation method for common principal components. One of
#'   "stepwise", "Flury", or "CAP". Otherwise, skips step.
#' @param err.method Harmonization method for error terms. One of "logPC".
#' @param cpc.cap.oc logical, only used if \code{cpc.method = "CAP"}, if
#' \code{TRUE} imposes orthogonal constraint when identifiying higher-order PCs.
#'
#' @return \code{cpcharmony} returns a list containing the following components:
#' \item{dat.out}{harmonized data as an array of input dimension}
#' \item{cpcs}{CPC analysis output list}
#' \item{cpc.harmony}{harmonization method output}
#'
#' @import CovBat
#' @importFrom cpca cpc
#' @export
#'
#' @examples

cpcharmony <- function(dat, bat, cpc.k, mod = NULL, method = "ComBat",
                       cpc.method = "stepwise", err.method = "logPC",
                       cpc.cap.oc = FALSE) {
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

  # Common principal components step
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

  } else {
    cpcs <- NULL
    dat_err = dat
  }

  # Harmonization of error matrices
  if (err.method == "logPC") {
    dat_err = log_covbat(dat_err, bat, mod = mod)$dat.covbat
  }

  if (is.null(cpcs)) {
    cpc_harmony = NULL
    dat_out = dat_err
  } else if (method == "ComBat") {
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
