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
                       metric = "E", lap.thr = 0.25, lap.gam = 0.01,
                       mdmr.perm = NULL) {
  if (is.null(bat)) {stop("Need to specify batch")}
  dat <- list(raw, out, ...)
  L <- length(dat)

  if (length(dat) > length(labs)) {
    message("Not enough labels: using default labels")
    labs <- c("Raw", "Out", paste0("Out", 2:(L-1)))
  }
  labs <- labs[1:L]
  names(dat) <- labs

  # Only calculate matrix logs once
  if (metric == "L") {
    if ("Correlation" %in% tests) {
      dat_c <- lapply(dat, function(x) array(apply(x, 3, cov2cor), dim(x)))
      dat_c <- lapply(dat_c, function(x) array(apply(x, 3, logm_eig), dim(x)))
    }
    if ("Laplacian" %in% tests) {
      dat_l <- lapply(dat, function(x) array(apply(x, 3, cov2lap,
                                                   lap.thr, lap.gam), dim(x)))
      dat_l <- lapply(dat_l, function(x) array(apply(x, 3, logm_eig), dim(x)))
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

  if ("Original" %in% tests) {
    # Frobenius distance of original matrices
    orig_out <- matrix(NA, P, L, dimnames = list(pair_names,
                                                 paste0(labs, ".Orig")))
    for (p in 1:P) {
      all_pair <- expand.grid(which(bat == bat_pairs[1, p]),
                              which(bat == bat_pairs[2, p]))
      if (metric == "E") {
        all_norm <- lapply(dat, function(x)
          apply(all_pair, 1, function(pair)
            norm(x[,,pair[1]] - x[,,pair[2]], "f")))
      } else {
        all_norm <- lapply(dat, function(x)
          apply(all_pair, 1, function(pair)
            CovDist(x[,,pair], metric)[1,2]))
      }
      orig_out[p,] <- sapply(all_norm, mean)
    }
    all_out = c(all_out, list(orig_out))
  }
  if ("Correlation" %in% tests) {
    # Frobenius distance of correlation matrices
    corr_out <- matrix(NA, P, L, dimnames = list(pair_names,
                                                 paste0(labs, ".Corr")))
    for (p in 1:P) {
      if (metric == "E") {
        all_norm <- lapply(dat_c, function(x)
          apply(all_pair, 1, function(pair)
            norm(x[,,pair[1]] - x[,,pair[2]], "f")))
      } else {
        all_norm <- lapply(dat_c, function(x)
          apply(all_pair, 1, function(pair)
            CovDist(x[,,pair], metric)[1,2]))
      }
    corr_out[p,] <- sapply(all_norm, mean)
    }
    all_out = c(all_out, list(corr_out))
  }
  if ("Laplacian" %in% tests) {
    # Frobenius distance of Laplacian matrices
    lapl_out <- matrix(NA, P, L, dimnames = list(pair_names,
                                                 paste0(labs, ".Lapl")))
    for (p in 1:P) {
      if (metric == "E") {
        all_norm <- lapply(dat_l, function(x)
          apply(all_pair, 1, function(pair)
            norm(x[,,pair[1]] - x[,,pair[2]], "f")))
      } else {
        all_norm <- lapply(dat_l, function(x)
          apply(all_pair, 1, function(pair)
            CovDist(x[,,pair], metric)[1,2]))
      }
      lapl_out[p,] <- sapply(all_norm, mean)
    }
    all_out = c(all_out, list(lapl_out))
  }
  if ("MDMR" %in% tests) {
    mdmr_out <- matrix(NA, P, L, dimnames = list(pair_names,
                                                 paste0(labs, ".MDMR.p")))
    if (metric == "E") {
      all_dist <- lapply(dat, function(x) dist(t(apply(x, 3, c))))
    } else {
      all_dist <- lapply(dat, function(x) CovDist(x, metric, as.dist = TRUE))
    }

    if (is.null(mdmr.perm)) {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(bat), D = x, seed = 8888))
    } else {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(bat), D = x, perm.p = mdmr.perm, seed = 8888))
    }

    mdmr_out[1,] <- do.call(cbind, sapply(mdmr_res, getElement, "pv"))[1,]
    all_out = c(all_out, list(mdmr_out))
  }
  do.call(cbind, all_out)
}
