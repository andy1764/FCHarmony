#' Title
#'
#' @param raw
#' @param out
#' @param tests
#' @param metric Metric used to calculate distance. Choose from \code{method} in
#'   \link[CovTools]{CovDist}.
#' @param log.eigen
#' @param bat
#'
#' @return
#' @importFrom expm expm
#' @importFrom CovTools CovDist
#' @importFrom MDMR mdmr
#' @importFrom stats dist
#' @importFrom utils combn
#' @export
#'
#' @examples
test_norms <- function(raw, out, ..., bat = NULL, labs = c("Raw", "Out"),
                       tests = c("Original", "Correlation", "Laplacian",
                                 "MDMR"),
                       metric = "E", lap.thr = 0.25, lap.gam = 0.01) {
  if (is.null(bat)) {stop("Need to specify batch")}
  dat <- list(raw, out, ...)
  N <- length(dat)

  if (length(dat) != length(labs)) {
    message("Number of inputs and labels differs: using default labels")
    labs <- c("Raw", "Out", paste0("Out", 2:(N-1)))
  }
  names(dat) <- labs[1:N]

  # Only calculate matrix logs once
  if (metric == "L") {
    if ("Correlation" %in% tests) {
      dat_c <- lapply(dat, function(x) array(apply(x, 3, cov2cor), dim(x)))
    }
    if ("Laplacian" %in% tests) {
      dat_l <- lapply(dat, function(x) array(apply(x, 3, cov2lap,
                                                   lap.thr, lap.gam), dim(x)))
    }
    dat <- lapply(dat, function(x) array(apply(x, 3, logm_eig), dim(x)))
    metric <- "E"
  } else {
    dat_c <- lapply(dat, function(x) array(apply(x, 3, cov2cor), dim(x)))
    dat_l <- lapply(dat, function(x) array(apply(x, 3, cov2lap,
                                                 lap.thr, lap.gam), dim(x)))
  }

  bat <- droplevels(bat)
  bat_lvl <- levels(bat)
  bat_pairs <- combn(bat_lvl, 2)
  P <- ncol(bat_pairs) # number of pairs
  pair_names <- apply(bat_pairs, 2, paste0, collapse = ",")

  all_out <- list()

  corr_out <- matrix(NA, P, N, dimnames = list(pair_names,
                                              c("Raw.Corr", "Out.Corr")))
  lapl_out <- matrix(NA, P, N, dimnames = list(pair_names,
                                              c("Raw.Lapl", "Out.Lapl")))
  mdmr_out <- matrix(NA, P, N, dimnames = list(pair_names,
                                               c("Raw.MDMR", "Out.MDMR")))



  if ("Original" %in% tests) {
    for (p in 1:P) {
      all_pair <- expand.grid(which(bat == bat_pairs[1, p]),
                              which(bat == bat_pairs[2, p]))
      # Frobenius distance of original matrices
      orig_out <- matrix(NA, P, N, dimnames = list(pair_names,
                                                   paste0(labs, ".Orig")))
      if (metric == "E") {
        raw_orig <- apply(all_pair, 1, function(x)
          norm(raw[,,x[1]] - raw[,,x[2]], "f"))
        out_orig <- apply(all_pair, 1, function(x)
          norm(out[,,x[1]] - out[,,x[2]], "f"))
      } else {
        raw_orig <- apply(all_pair, 1, function(x)
          CovDist(raw[,,x], metric)[1,2])
        out_orig <- apply(all_pair, 1, function(x)
          CovDist(out[,,x], metric)[1,2])
      }
      orig_out[p,] <- c(mean(raw_orig), mean(out_orig))
    }


    }
    if ("Correlation" %in% tests) {
      # Frobenius distance of correlation matrices
      if (metric == "E") {
        raw_corr <- apply(all_pair, 1, function(x)
          norm(raw_c[,,x[1]] - raw_c[,,x[2]], "f"))
        out_corr <- apply(all_pair, 1, function(x)
          norm(out_c[,,x[1]] - out_c[,,x[2]], "f"))
      } else {
        raw_corr <- apply(all_pair, 1, function(x)
          CovDist(raw_c[,,x], metric)[1,2])
        out_corr <- apply(all_pair, 1, function(x)
          CovDist(out_c[,,x], metric)[1,2])
      }

      corr_out[p,] <- c(mean(raw_corr), mean(out_corr))
    }
    if ("Laplacian" %in% tests) {
      # Frobenius distance of Laplacian matrices
      if (metric == "E") {
        raw_lapl <- apply(all_pair, 1, function(x)
          norm(raw_l[,,x[1]] - raw_l[,,x[2]], "f"))
        out_lapl <- apply(all_pair, 1, function(x)
          norm(out_l[,,x[1]] - out_l[,,x[2]], "f"))
      } else {
        raw_lapl <- apply(all_pair, 1, function(x)
          CovDist(raw_l[,,x], metric)[1,2])
        out_lapl <- apply(all_pair, 1, function(x)
          CovDist(out_l[,,x], metric)[1,2])
      }

      lapl_out[p,] <- c(mean(raw_lapl), mean(out_lapl))
    }
  }
  if ("MDMR" %in% tests) {
    if (metric == "E") {
      dist_raw <- dist(t(apply(raw, 3, c)))
      dist_out <- dist(t(apply(out, 3, c)))
    } else {
      dist_raw <- CovDist(raw, metric, as.dist = TRUE)
      dist_out <- CovDist(out, metric, as.dist = TRUE)
    }

    mdmr_res_raw <- mdmr(data.frame(bat), D = dist_raw, seed = 8888)
    mdmr_res_out <- mdmr(data.frame(bat), D = dist_out, seed = 8888)

    mdmr_out[1,] <- c(mdmr_res_raw$pv[1,1], mdmr_res_out$pv[1,1])
  }
  cbind(orig_out, corr_out, lapl_out, mdmr_out)
}
