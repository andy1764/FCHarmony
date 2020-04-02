#' Title
#'
#' @param raw
#' @param out
#' @param tests
#' @param metric Metric used to calculate distance. Choose from \code{method} in
#'   \link[shapes]{distcov}.
#' @param log.eigen
#' @param bat
#'
#' @return
#' @importFrom expm expm
#' @importFrom shapes distcov
#' @export
#'
#' @examples
test_norms <- function(raw, out, bat,
                       tests = c("Original", "Correlation", "Laplacian"),
                       metric = "Euclidean", lap.thr = 0.25, lap.gam = 0.01) {
  # Only calculate matrix logs once
  if (metric == "LogEuclidean") {
    if ("Correlation" %in% tests) {
      raw_c <- sapply(seq_along(raw[1,1,]),
                      function(x) logm_eig(cov2cor(raw[,,x])),
                      simplify = "array")
      out_c <- sapply(seq_along(out[1,1,]),
                      function(x) logm_eig(cov2cor(out[,,x])),
                      simplify = "array")
    }

    if ("Laplacian" %in% tests) {
      raw_l <- sapply(seq_along(raw[1,1,]),
                      function(x) logm_eig(cov2lap(raw[,,x], lap.thr, lap.gam)),
                      simplify = "array")
      out_l <- sapply(seq_along(out[1,1,]),
                      function(x) logm_eig(cov2lap(out[,,x], lap.thr, lap.gam)),
                      simplify = "array")
    }

    raw <- sapply(seq_along(raw[1,1,]), function(x) logm_eig(raw[,,x]),
                    simplify = "array")
    out <- sapply(seq_along(out[1,1,]), function(x) logm_eig(out[,,x]),
                  simplify = "array")

    metric <- "Euclidean"
  } else {
    raw_c <- sapply(seq_along(raw[1,1,]), function(x) cov2cor(raw[,,x]),
                    simplify = "array")
    out_c <- sapply(seq_along(out[1,1,]), function(x) cov2cor(out[,,x]),
                    simplify = "array")

    raw_l <- sapply(seq_along(raw[1,1,]),
                    function(x) cov2lap(raw[,,x], lap.thr, lap.gam),
                    simplify = "array")
    out_l <- sapply(seq_along(out[1,1,]),
                    function(x) cov2lap(out[,,x], lap.thr, lap.gam),
                    simplify = "array")
  }

  bat <- droplevels(bat)
  bat_lvl <- levels(bat)
  bat_pairs <- combn(bat_lvl, 2)
  P <- ncol(bat_pairs) # number of pairs
  pair_names <- apply(bat_pairs, 2, paste0, collapse = ",")
  orig_out <- matrix(NA, P, 2, dimnames = list(pair_names,
                                              c("Raw.Orig", "Out.Orig")))
  corr_out <- matrix(NA, P, 2, dimnames = list(pair_names,
                                              c("Raw.Corr", "Out.Corr")))
  lapl_out <- matrix(NA, P, 2, dimnames = list(pair_names,
                                              c("Raw.Lapl", "Out.Lapl")))

  for (p in 1:P) {
    all_pair <- expand.grid(which(bat == bat_pairs[1, p]),
                            which(bat == bat_pairs[2, p]))

    if ("Original" %in% tests) {
      # Frobenius distance of original matrices
      raw_orig <- apply(all_pair, 1, function(x)
        distcov(raw[,,x[1]], raw[,,x[2]], metric))
      out_orig <- apply(all_pair, 1, function(x)
        distcov(out[,,x[1]], out[,,x[2]], metric))

      orig_out[p,] <- c(mean(raw_orig), mean(out_orig))
    }
    if ("Correlation" %in% tests) {
      # Frobenius distance of correlation matrices
      raw_corr <- apply(all_pair, 1, function(x)
        distcov(raw_c[,,x[1]], raw_c[,,x[2]], metric))
      out_corr <- apply(all_pair, 1, function(x)
        distcov(out_c[,,x[1]], out_c[,,x[2]], metric))

      corr_out[p,] <- c(mean(raw_corr), mean(out_corr))
    }
    if ("Laplacian" %in% tests) {
      # Frobenius distance of Laplacian matrices
      raw_lapl <- apply(all_pair, 1, function(x)
        distcov(raw_l[,,x[1]], raw_l[,,x[2]], metric))
      out_lapl <- apply(all_pair, 1, function(x)
        distcov(out_l[,,x[1]], out_l[,,x[2]], metric))

      lapl_out[p,] <- c(mean(raw_lapl), mean(out_lapl))
    }
  }
  cbind(orig_out, corr_out, lapl_out)
}
