#' Harmonization using Common Principal Components
#'
#' @param dat An array of covariance matrices of dimension
#' \emph{(p x p x N)}
#' @param bat A vector of batch labels as length \emph{n} vector
#' @param mod A matrix of covariates to preserve during harmonization
#' @param cpc.k Numeric, number of common principal components to harmonize
#' @param cpc.method Estimation method for common principal components. One of
#'   "stepwise", "Flury", or "CAP". Otherwise, skips step.
#' @param err.method Harmonization method for error terms. One of "log-CovBat"
#'   or "fc-Combat"
#' @param method Method for harmonization. One of "likelihood", "log-ComBat",
#'   "ComBat". Currently only ComBat supported.
#' @param force.PD Logical vector. First element indicates whether input should
#'   be forced to be positive definite using \link[Matrix]{nearPD}, second for
#'   output
#' @param to.corr Logical vector. First element indicates whether input should
#'   be forced to be a correlation matrix using  \link[stats]{cov2cor}, second
#'   for output
#' @param cpc.cap.oc Logical, only used if \code{cpc.method = "CAP"}, if
#' \code{TRUE} imposes orthogonal constraint when identifiying higher-order PCs.
#'
#' @return \code{cpcharmony} returns a list containing the following components:
#' \item{dat.out}{harmonized data as an array of input dimension}
#' \item{cpcs}{CPC analysis output list}
#' \item{cpc.harmony}{harmonization method output}
#'
#' @import CovBat
#' @importFrom cpca cpc
#' @importFrom Matrix nearPD
#' @importFrom stats cov2cor
#' @export
#'
#' @examples

cpcharmony <- function(dat, bat, mod = NULL, cpc.method = "stepwise",
                       cpc.k = dim(dat)[1], err.method = "log-CovBat",
                       method = "ComBat",
                       force.PD = c(FALSE, FALSE), to.corr = c(FALSE, FALSE),
                       cpc.cap.oc = FALSE) {
  p <- dim(dat)[2]
  N <- dim(dat)[3]

  if (force.PD[1]) {
    for (i in 1:N) {
      dat[,,i] <- as.numeric(nearPD(dat[,,i])$mat)
    }
  }
  if (to.corr[1]) {
    for (i in 1:N) {
      dat[,,i] <- cov2cor(dat[,,i])
    }
  }

  dat_err <- dat # store error matrices
  dat_out <- array(0, dim = dim(dat)) # initialize output matrices

  #### Common principal components step ####
  switch(
    cpc.method,
    "stepwise" = {
      cpcs <- cpc(dat, k = cpc.k) # get common principal components
      for (i in 1:N) {
        for (j in 1:cpc.k) {
          dat_err[,,i] = dat_err[,,i] - cpcs$D[j,i] * tcrossprod(cpcs$CPC[,j])
        }
      }
    },
    "Flury" = {
      cpcs <- cpc_Flury(dat)
      for (i in 1:N) {
        for (j in 1:cpc.k) {
          dat_err[,,i] = dat_err[,,i] - cpcs$D[j,i] * tcrossprod(cpcs$CPC[,j])
        }
      }
    },
    "default" = {
      # placeholder for now
    },
    {
      cpcs <- NULL
    })

  #### Harmonization of error matrices ####
  switch(
    err.method,
    "log-Covbat" = {
      dat_err = log_covbat(dat_err, bat, mod = mod)$dat.covbat
    },
    "fc-ComBat" =  {
      # Method first implemented by Yu et al. (2018), DOI: 10.1002/hbm.24241
      err_vec <- t(apply(dat_err, 3, function(x) c(x[lower.tri(x)])))
      err_com <- t(combat_modded(t(err_vec), bat, mod = mod)$dat.combat)
      err_out <- array(0, dim = dim(dat_err))
      for (i in 1:N) {
        diag(err_out[,,i]) <- diag(dat_err[,,i])
        err_out[,,i][lower.tri(err_out[,,i])] <- err_com[i,]
        err_out[,,i] <- err_out[,,i] + t(err_out[,,i])
      }
      dat_err <- err_out
    },
    {})

  #### Harmonization of CPC eigenvalues ####
  if (is.null(cpcs)) {
    cpc_harmony = NULL
    dat_out = dat_err
  } else {
    switch(
      method,
      "ComBat" = {
        cpc_harmony <- combat_modded(cpcs$D, bat, mod = mod)
        for (i in 1:N) {
          dat_out[,,i] <- dat_err[,,i] # add back error
          for (j in 1:cpc.k) {
            dat_out[,,i] = dat_out[,,i] +
              cpc_harmony$dat.combat[j,i] * tcrossprod(cpcs$CPC[,j])
          }
        }
      },
      "log-ComBat" = {
        cpc_harmony <- combat_modded(log(cpcs$D), bat, mod = mod)
        for (i in 1:N) {
          dat_out[,,i] <- dat_err[,,i] # add back error
          for (j in 1:cpc.k) {
            dat_out[,,i] = dat_out[,,i] +
              exp(cpc_harmony$dat.combat[j,i]) * tcrossprod(cpcs$CPC[,j])
          }
        }
      },
      {
        cpc_harmony <- NULL
        dat_out = dat_err
      })
  }

  if (force.PD[2]) {
    for (i in 1:N) {
      dat_out[,,i] <- as.numeric(nearPD(dat_out[,,i])$mat)
    }
  }
  if (to.corr[2]) {
    for (i in 1:N) {
      dat_out[,,i] <- cov2cor(dat_out[,,i])
    }
  }

  list(dat.out = dat_out,
       dat.err = dat_err,
       dat.in = dat,
       cpcs = cpcs,
       cpc.harmony = cpc_harmony)
}
