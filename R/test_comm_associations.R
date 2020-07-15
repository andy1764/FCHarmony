#' Community detection tests
#'
#' Examine if grouping affects community detection.
#'
#' @param x *p x p x n* array of functional connectivity matrices or object of
#'   class 'dist', see \link[stats](dist)
#' @param grp Factor specifying groups over which to aggregate.
#' @param atlas Vector of labels specifying original subnetworks.
#' @param grp.method Method for obtaining group-level functional connectivity.
#' @param comm.method Function from `igraph` that finds communities. See
#'   \link[igraph]{membership} for options. Defaults to
#'   \link[igraph]{cluster_louvain}
#' @param comp.method Metric for comparing community detections.
#' @param fisher Whether to z-transform input FC matrices
#' @param debug
#'
#' @return
#' @import igraph
#' @importFrom vegan adonis2
#' @importFrom mcclust vi.dist
#' @export
#'
#' @examples
test_comm_associations <- function(x, grp,
                                   comm.method = cluster_louvain,
                                   metric = "vi",
                                   threshold = "positive",
                                   adonis.params = NULL,
                                   fisher = TRUE,
                                   debug = FALSE) {
  if (is.null(grp)) {stop("Need to specify group")}
  grp <- droplevels(grp)

  n <- dim(x)[3]

  if (class(x) == "dist") {
    clust_dist <- x
  } else {
    # if threshold specified, transform edges using threshold function
    # also check for predefined thresholds for ease of use
    if (class(threshold) == "character") {
      if (threshold == "positive") {
        threshold <- function(x) {ifelse(x>0, x, 0)}
      }
    }
    if (!is.null(threshold)) {
      x_g <- array(apply(x, 3, function(y) {
        thr <- threshold(y[lower.tri(y)])
        out <- matrix(0, dim(x), dim(x))
        out[lower.tri(out)] <- thr
        out + t(out)
      }), dim(x))
    } else {
      x_g
    }

    if (fisher) {
      x_g <- array(apply(x_g, 3, atanh), dim(x))
    }

    gr <- apply(x_g, 3, graph_from_adjacency_matrix,
                mode = "undirected", weighted = TRUE)

    #### Community detection algorithm ####
    clust <- lapply(gr, comm.method)
    clust_mat <- sapply(clust, function(y) y$membership)

    #### Get distance matrix between communities ####
    clust_dist <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:i) {
        if (i == j) {
          clust_dist[i,j] <- 0
        } else {
          clust_dist[i,j] = igraph::compare(clust[[i]], clust[[j]], metric)
        }
      }
    }
    clust_dist <- as.dist(clust_dist)
  }

  #### Perform PERMANOVA to test for differences between grp ####
  clust_assoc <- do.call(adonis, c(formula = as.formula("clust_dist ~ grp"),
                                      adonis.params))
  out <- clust_assoc$aov.tab

  if (debug) {
    list(out = out,
         clust = clust_mat,
         clust.dist = clust_dist,
         clust.assoc = clust_assoc,
         graphs = gr)
  } else {
    out
  }
}
