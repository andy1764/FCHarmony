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
#'   output.
#' @param to.corr Logical vector. First element indicates whether input should
#'   be forced to be a correlation matrix using  \link[stats]{cov2cor}, second
#'   for output.
#' @param log.input Logical, whether to use the matrix logarithm of the input.
#' @param err.eb Logical, use empirical Bayes in error harmonization.
#' @param cpc.eb Logical, use empirical Bayes in CPC eigenvalue harmonization.
#' @param debug Logical, include internal objects in output.
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
#' @importFrom expm expm
#' @export
#'
#' @examples

cpcharmony <- function(dat, bat, mod = NULL,
                       cpc.method = c("stepwise", "Flury", "CAP", "none"),
                       cpc.k = dim(dat)[1],
                       err.method = c("ComBat", "CovBat", "PC-ComBat",
                                      "block-ComBat", "remove", "none"),
                       method = c("ComBat", "log-ComBat", "none"),
                       log.input = FALSE,
                       force.PD = c(FALSE, FALSE), to.corr = c(FALSE, FALSE),
                       err.eb = TRUE, cpc.eb = TRUE, debug = FALSE) {
  bat <- droplevels(bat)
  try(cpc.method <- match.arg(cpc.method))
  try(err.method <- match.arg(err.method))
  try(method <- match.arg(method))

  p <- dim(dat)[2]
  N <- dim(dat)[3]
  dnames <- dimnames(dat)

  if (force.PD[1]) {
    dat <- array(apply(dat, 3, function(x) {
      as.numeric(nearPD(x)$mat)
    }), dim(dat))
  }
  if (to.corr[1]) {
    dat <- array(apply(dat, 3, cov2cor), dim(dat))
  }
  if (log.input) {
    dat <- array(apply(dat, 3, logm_eig), dim(dat))
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
      # calculate norm percent as % explained by CPCs, = 1 - norm(err)/norm(in)
      cpcs$norm.percent = mean(1 - apply(dat_err, 3, norm, "f")/
        apply(dat, 3, norm, "f"))
    },
    "Flury" = {
      cpcs <- cpc_Flury(dat)
      for (i in 1:N) {
        for (j in 1:cpc.k) {
          dat_err[,,i] = dat_err[,,i] - cpcs$D[j,i] * tcrossprod(cpcs$CPC[,j])
        }
      }
      cpcs$norm.percent = mean(1 - apply(dat_err, 3, norm, "f")/
                                 apply(dat, 3, norm, "f"))
    },
    {
      cpcs <- NULL
    })
  dat_err_raw <- dat_err

  #### Harmonization of error matrices ####
  switch(
    err.method,
    "ComBat" = {
      # err_vec <- t(apply(dat_err, 3, function(x) c(x[lower.tri(x)])))
      # err_harmony <- combat(t(err_vec), bat, mod = mod, eb = err.eb)
      # err_harm_dat <- t(err_harmony$dat.combat)
      # err_out <- array(0, dim = dim(dat_err))
      # for (i in 1:N) {
      #   err_out[,,i][lower.tri(err_out[,,i])] <- err_harm_dat[i,]
      #   err_out[,,i] <- err_out[,,i] + t(err_out[,,i])
      # }
      err_harmony <- fcComBat(dat_err, bat, mod = mod, to.corr = FALSE,
                              fisher = FALSE)
      dat_err <- err_harmony$dat.out
    },
    "CovBat" = {
      # err_vec <- t(apply(dat_err, 3, function(x) c(x[lower.tri(x)])))
      # err_harmony <- covbat(t(err_vec), bat, mod = mod, eb = err.eb)
      # err_harm_dat <- t(err_harmony$dat.covbat)
      # err_out <- array(0, dim = dim(dat_err))
      # for (i in 1:N) {
      #   err_out[,,i][lower.tri(err_out[,,i])] <- err_harm_dat[i,]
      #   err_out[,,i] <- err_out[,,i] + t(err_out[,,i])
      # }
      err_harmony <- fcCovBat(dat_err, bat, mod = mod, to.corr = FALSE,
                              fisher = FALSE)
      dat_err <- err_harmony$dat.out
    },
    "block-ComBat" = {
      err_harmony <- blComBat(dat_err, bat, mod = mod, to.corr = FALSE,
                              fisher = FALSE)
      dat_err <- err_harmony$dat.out
    },
    "PC-ComBat" = {
      err_harmony = pc_combat(dat_err, bat, mod = mod)
      dat_err = err_harmony$dat.out
    },
    "remove" = {
      dat_err[] <- 0
      err_harmony = NULL
    },
    {
      message("No error harmonization")
      err_harmony = NULL
    })

  #### Harmonization of CPC eigenvalues ####
  if (is.null(cpcs)) {
    message("No CPC harmonization")
    cpc_harmony = NULL
    dat_out = dat_err
  } else {
    switch(
      method,
      "ComBat" = {
        cpc_harmony <- combat(cpcs$D, bat, mod = mod, eb = cpc.eb)
        for (i in 1:N) {
          dat_out[,,i] <- dat_err[,,i] # add back error
          for (j in 1:cpc.k) {
            dat_out[,,i] = dat_out[,,i] +
              cpc_harmony$dat.combat[j,i] * tcrossprod(cpcs$CPC[,j])
          }
        }
      },
      "log-ComBat" = {
        cpc_harmony <- combat(log(cpcs$D), bat, mod = mod, eb = cpc.eb)
        for (i in 1:N) {
          dat_out[,,i] <- dat_err[,,i] # add back error
          for (j in 1:cpc.k) {
            dat_out[,,i] = dat_out[,,i] +
              exp(cpc_harmony$dat.combat[j,i]) * tcrossprod(cpcs$CPC[,j])
          }
        }
      },
      {
        message("No CPC harmonization")
        cpc_harmony <- NULL
        dat_out = dat_err
      })
  }

  dimnames(dat_out) <- dnames

  if (log.input) {
    dat_out <- array(apply(dat_out, 3, expm, "R_Eigen"), dim(dat_out))
  }
  if (force.PD[2]) {
    dat_out <- array(apply(dat_out, 3, function(x) {
      as.numeric(nearPD(x)$mat)
    }), dim(dat_out))
  }
  if (to.corr[2]) {
    dat_out <- array(apply(dat_out, 3, cov2cor), dim(dat_out))
  }

  # remove unnecessary arrays if not debug
  if (debug) {
    list(dat.out = dat_out,
         dat.err.raw = dat_err_raw,
         dat.err = dat_err,
         dat.in = dat,
         cpcs = cpcs,
         cpc.harmony = cpc_harmony,
         err.harmony = err_harmony)
  } else {
    list(dat.out = dat_out,
         cpcs = cpcs)
  }
}
