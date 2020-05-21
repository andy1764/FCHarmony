#' Convert lower triangular entries to symmetric matrix
#' @param vec Lower triangular elements, typically derived from
#'   \link[base]{lower.tri} with `diag = FALSE`.
#' @noRd
vec2corr <- function(vec, names = NULL, corr = TRUE) {
  d <- (1 + sqrt(8*length(vec)+1))/2 # dim of output
  out <- matrix(0, d, d, dimnames = names)
  out[lower.tri(out, diag=FALSE)] <- vec
  out <- out + t(out)
  if (corr) {diag(out) <- 1}
  out
}

#' Convert lower triangular incl. diagonals to symmetric matrix
#' @param vec Lower triangular elements, typically derived from
#'   \link[base]{lower.tri} with `diag = TRUE`.
#' @noRd
dvec2corr <- function(vec, names = NULL) {
  d <- (sqrt(8*length(vec)+1)-1)/2 # dim of output
  out <- matrix(0, d, d, dimnames = names)
  out[lower.tri(out, diag=TRUE)] <- vec
  out <- out + t(out)
  diag(out) <- diag(out)/2
  out
}

#' Convert covariance to Laplacian with a given threshold
#' @param cov Covariance matrix.
#' @param threshold Correlation threshold for edge
#' @param gamma Added value to ensure positive definite Laplacians
#' @noRd
cov2lap <- function(cov, threshold = 0.5, gamma = 0.01) {
  corr <- cov2cor(as.matrix(nearPD(cov)$mat))
  corr[upper.tri(corr)] <- as.numeric(corr[upper.tri(corr)] >= threshold)
  corr[lower.tri(corr)] <- t(corr)[lower.tri(corr)]
  diag(corr) <- 0
  corr <- -corr
  # gamma addition ensures positive definite
  diag(corr) <- -apply(corr, 1, sum) + gamma
  return(corr)
}

#' Convert correlation matrix to within- and between- connectivities
#' @param x Correlation matrix.
#' @param dims Labels to group over.
#' @param fisher z-transform elements of x prior to grouping
#' @noRd
corr2con <- function(corr, dims = dimnames(corr), fisher = TRUE) {
  if (fisher) {corr[] <- atanh(corr)}
  diag(corr) <- 0 # ignore diagonal to get avg connectivities
  rois <- sort(unique(dims[[1]]))
  p <- length(rois)

  out <- matrix(0, p, p, dimnames = list(rois, rois))
  for (i in 1:p) {
    ind_i <- dims[[1]] == rois[i]
    for (j in 1:i) {
      ind_j <- dims[[2]] == rois[j]
      out[i,j] <- mean(corr[ind_i, ind_j]/2) # divided by two due to symmetry
    }
  }
  out <- out + t(out)
  diag(out) <- diag(out)/2
  out
}

#' Convert correlation matrix to partial correlation matrix
#' @param cor Correlation matrix.
#' @noRd
corr2pcor <- function(cor) {
  inv <- -solve(cor)
  diag(inv) <- -diag(inv)
  cov2cor(inv)
}

#' Get matrix log of SPD matrix
#' @param x Symmetric positive definite matrix.
#' @noRd
logm_eig <- function(x) {
  eig <- eigen(x, symmetric = TRUE)
  if (any(eig$values <= 0)) {stop("Input is not positive definite")}
  eig$vectors %*% diag(log(eig$values)) %*% t(eig$vectors)
}

#' Local efficiency, fixed version of \link[brainGraph]{local_eff}
#'
#' Additional option to specify nodes over which to calculate.
#' @param g \link[igraph]{igraph} object.
#' @param ind Numeric, node indicies to calculate local efficiency for.
#' @noRd
local_eff <- function (g, ind, weights = NULL, use.parallel = TRUE) {
  if (is.null(weights)) {
    A <- as_adj(g, names = FALSE)
    weighted <- NULL
  } else {
    A <- as_adj(g, names = FALSE, attr = "weight")
    weighted <- TRUE
  }
  eff <- rep(0, nrow(A))
  nodes <- which(rowSums((A > 0) + 0) > 1)
  if (!is.null(ind)) {
    if (is.logical(ind)) {ind <- which(ind)}
    nodes <- intersect(nodes, ind)
    }
  X <- apply(A, 1, function(x) which(x > 0))
  # handle case when every subgraph has the same number of nodes
  if (class(X) == "matrix") {X <- split(X, col(X))}
  if (length(nodes) > 0) {
    if (isTRUE(use.parallel)) {
      eff[nodes] <- foreach(i = nodes, .combine = "c") %dopar%
        {
          g.sub <- graph_from_adjacency_matrix(A[X[[i]],
                                                 X[[i]]], mode = "undirected", weighted = weighted)
          efficiency(g.sub, "global")
        }
    }
    else {
      for (i in nodes) {
        g.sub <- graph_from_adjacency_matrix(A[X[[i]],
                                               X[[i]]], mode = "undirected", weighted = weighted)
        eff[i] <- efficiency(g.sub, "global")
      }
    }
  }
  eff
}
