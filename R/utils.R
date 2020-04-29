# simple function to convert lower triangular entries to symmetric matrix
vec2corr <- function(vec, names = NULL, corr = FALSE) {
  d <- (1 + sqrt(8*length(vec)+1))/2 # dim of output
  out <- matrix(0, d, d, dimnames = names)
  out[lower.tri(out, diag=FALSE)] <- vec
  out <- out + t(out)
  if (corr) {diag(out) <- 1}
  out
}

# simple function to convert lower triangular (incl diag) to symmetric matrix
dvec2corr <- function(vec, names = NULL) {
  d <- (sqrt(8*length(vec)+1)-1)/2 # dim of output
  out <- matrix(0, d, d, dimnames = names)
  out[lower.tri(out, diag=TRUE)] <- vec
  out <- out + t(out)
  diag(out) <- diag(out)/2
  out
}

# convert covariance to Laplacian with a given threshold,
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

# convert correlation matrix to connectivity values grouped by ROI
corr2con <- function(corr, dims = NULL) {
  if (is.null(dims)) {dims <- dimnames(corr)}
  rois <- sort(unique(dims[[1]]))
  p <- length(rois)
  corr[lower.tri(corr, diag = TRUE)] <- 0

  out <- matrix(0, p, p, dimnames = list(rois, rois))
  for (i in 1:p) {
    ind_i <- dims[[1]] == rois[i]
    for (j in 1:i) {
      ind_j <- dims[[2]] == rois[j]
      out[i,j] <- mean(corr[ind_i, ind_j])
    }
  }
  out <- out + t(out)
  diag(out) <- diag(out)/2
  out
}

# convert correlation matrix to partial correlation matrix
corr2pcor <- function(cor) {
  inv <- -solve(cor)
  diag(inv) <- -diag(inv)
  cov2cor(inv)
}

# get matrix logarithm of SPD matrix
logm_eig <- function(x) {
  eig <- eigen(x, symmetric = TRUE)
  if (any(eig$values <= 0)) {stop("Input is not positive definite")}
  eig$vectors %*% diag(log(eig$values)) %*% t(eig$vectors)
}

# # get matrix exponential of diagonalizable matrix
# expm_eig <- function(x) {
#   eig <- eigen(x, symmetric = TRUE)
#   eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
# }

# logm_eig <- function(x) {
#   n <- dim(x)[1]
#   eig <- eigen(x, symmetric = TRUE)
#   if (any(eig$values <= 0)) {stop("Input is not positive definite")}
#   leig <- log(eig$values)
#   out <- matrix(0, n, n)
#   for (i in 1:n) {
#     out <- out + leig[i] * tcrossprod(eig$vectors[,i])
#   }
#   out
# }

# experimental: GPU-based matrix log
# logm_eig_gpu <- function(x) {
#   n <- dim(x)[1]
#   y <- vclMatrix(x, type = "double")
#   eig <- eigen(y, symmetric = TRUE)
#   # if (any(eig$values <= 0)) {stop("Input is not positive definite")}
#   leig <- log(eig$values)
#   leig_mat <- identity_matrix(n, "double")
#   diag(leig_mat) <- leig
#   out <- eig$vectors %*% leig_mat %*% t(eig$vectors)
#   as.matrix(out)
# }

# get local efficiency, based on code from brainGraph but fixed
# extra option to specify nodes to calculate for
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
