#' Title
#'
#'
#' @param raw
#' @param out
#' @param ...
#' @param bat
#' @param mod
#' @param labs
#' @param tests
#' @param metric
#' @param lap.thr
#' @param lap.gam
#' @param mdmr.perm
#'
#' @return
#' @importFrom CovTools CovDist
#' @importFrom MDMR mdmr
#' @importFrom stats dist
#' @importFrom utils combn
#' @export
#'
#' @examples
test_regress <- function(raw, out, ..., bat = NULL, mod = NULL,
                         labs = c("Raw", "Out"),
                         tests = c("MDMR", "MDMR-C", "MDMR-L"),
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
    if ("MDMR-C" %in% tests) {
      dat_c <- lapply(dat, function(x) array(apply(x, 3, cov2cor), dim(x)))
      dat_c <- lapply(dat_c, function(x) array(apply(x, 3, logm_eig), dim(x)))
    }
    if ("MDMR-L" %in% tests) {
      dat_l <- lapply(dat, function(x) array(apply(x, 3, cov2lap,
                                                   lap.thr, lap.gam), dim(x)))
      dat_l <- lapply(dat_l, function(x) array(apply(x, 3, logm_eig), dim(x)))
    }
    dat <- lapply(dat, function(x) array(apply(x, 3, logm_eig), dim(x)))
    metric <- "E"
  } else {
    if ("MDMR-C" %in% tests) {
      dat_c <- lapply(dat, function(x) array(apply(x, 3, cov2cor), dim(x)))
    }
    if ("MDMR-L" %in% tests) {
      dat_l <- lapply(dat, function(x) array(apply(x, 3, cov2lap,
                                                   lap.thr, lap.gam), dim(x)))
    }
  }

  bat <- droplevels(bat)
  all_out <- list()

  if ("MDMR" %in% tests) {
    if (metric == "E") {
      all_dist <- lapply(dat, function(x) dist(t(apply(x, 3, c))))
    } else {
      all_dist <- lapply(dat, function(x) CovDist(x, metric, as.dist = TRUE))
    }

    if (is.null(mdmr.perm)) {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(covt_mod[,-1], bat), D = x, seed = 8888))
    } else {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(covt_mod[,-1], bat), D = x,
             perm.p = mdmr.perm, seed = 8888))
    }

    mdmr_out <- do.call(rbind, sapply(mdmr_res, getElement, "pv"))
    dimnames(mdmr_out) <- list(labs, rownames(mdmr_res[[1]]$pv))
    all_out = c(all_out, list(mdmr_out))
  }
  if ("MDMR-C" %in% tests) {
    if (metric == "E") {
      all_dist <- lapply(dat_c, function(x) dist(t(apply(x, 3, c))))
    } else {
      all_dist <- lapply(dat_c, function(x) CovDist(x, metric, as.dist = TRUE))
    }

    if (is.null(mdmr.perm)) {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(covt_mod[,-1], bat), D = x, seed = 8888))
    } else {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(covt_mod[,-1], bat), D = x,
             perm.p = mdmr.perm, seed = 8888))
    }

    mdmr_out <- do.call(rbind, sapply(mdmr_res, getElement, "pv"))
    dimnames(mdmr_out) <- list(labs, paste0(rownames(mdmr_res[[1]]$pv), ".Corr"))
    all_out = c(all_out, list(mdmr_out))
  }
  if ("MDMR-L" %in% tests) {
    if (metric == "E") {
      all_dist <- lapply(dat_l, function(x) dist(t(apply(x, 3, c))))
    } else {
      all_dist <- lapply(dat_l, function(x) CovDist(x, metric, as.dist = TRUE))
    }

    if (is.null(mdmr.perm)) {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(covt_mod[,-1], bat), D = x, seed = 8888))
    } else {
      mdmr_res <- lapply(all_dist, function(x)
        mdmr(data.frame(covt_mod[,-1], bat), D = x,
             perm.p = mdmr.perm, seed = 8888))
    }

    mdmr_out <- do.call(rbind, sapply(mdmr_res, getElement, "pv"))
    dimnames(mdmr_out) <- list(labs, paste0(rownames(mdmr_res[[1]]$pv), ".Lapl"))
    all_out = c(all_out, list(mdmr_out))
  }
  do.call(cbind, all_out)
}
