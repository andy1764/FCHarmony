#' Functional Connectivity Harmonization
#'
#' @param dat An array of covariance matrices of dimension
#' \emph{(p x p x N)}
#' @param bat A vector of batch labels as length \emph{n} vector
#' @param mod A matrix of covariates to preserve during harmonization
#' @param method Harmonization method. One of "log-CovBat"
#'   or "fc-Combat"
#' @param force.PD Logical vector. First element indicates whether input should
#'   be forced to be positive definite using \link[Matrix]{nearPD}, second for
#'   output.
#' @param to.corr Logical vector. First element indicates whether input should
#'   be forced to be a correlation matrix using  \link[stats]{cov2cor}, second
#'   for output.
#' @param eb Logical, use empirical Bayes in harmonization.
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

fcharmony <- function(dat, bat, mod = NULL,
                      method = c("FC-ComBat", "FC-CovBat", "Block-ComBat",
                                 "Log-ComBat", "Log-CovBat", "PC-ComBat"),
                      force.PD = c(FALSE, FALSE), to.corr = c(FALSE, FALSE),
                      eb = TRUE, debug = FALSE) {
  bat <- droplevels(bat)
  method <- match.arg(method)

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

  #### Harmonization step ####
  switch(
    method,
    "FC-ComBat" = {
      harmony <- fcComBat(dat, bat, mod = mod, to.corr = FALSE, eb = eb)
      out <- harmony$dat.out
    },
    "FC-CovBat" = {
      harmony <- fcCovBat(dat, bat, mod = mod, to.corr = FALSE, eb = eb)
      out <- harmony$dat.out
    },
    "Block-ComBat" = {
      harmony <- blComBat(dat, bat, mod = mod, to.corr = FALSE, eb = eb)
      out <- harmony$dat.out
    },
    "PC-ComBat" = {
      harmony = pc_combat(dat, bat, mod = mod)
      out = harmony$dat.out
    },
    "Log-ComBat" = {
      harmony = log_combat(dat, bat, mod = mod, eb = eb)
      out = harmony$dat.out
    },
    "Log-CovBat" = {
      harmony = log_covbat(dat, bat, mod = mod, eb = eb)
      out = harmony$dat.out
    })

  dimnames(out) <- dnames

  if (force.PD[2]) {
    out <- array(apply(out, 3, function(x) {
      as.numeric(nearPD(x)$mat)
    }), dim(out))
  }
  if (to.corr[2]) {
    out <- array(apply(out, 3, cov2cor), dim(out))
  }

  # remove unnecessary arrays if not debug
  if (debug) {
    list(dat.out = out,
         method = method,
         harmony.out = harmony)
  } else {
    list(dat.out = out,
         method = method,
         harmony.out = harmony)
  }
}
