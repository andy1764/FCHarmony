#' Community detection tests
#'
#' Examine if grouping affects community detection.
#'
#' @param ...
#' @param labs
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
#' @export
#'
#' @examples
test_communities <- function(..., labs = c("Raw", "Out"), grp, atlas,
                             grp.method = c("average", "multilayer"),
                             comm.method = cluster_louvain,
                             metric = c("adjusted.rand", "nmi"),
                             threshold = "positive",
                             fisher = TRUE,
                             debug = FALSE) {
  if (is.null(grp)) {stop("Need to specify group")}
  dat <- list(...)
  L <- length(dat)

  if (length(dat) > length(labs)) {
    message("Not enough labels: assuming first is raw and using default labels")
    labs <- c("Raw", "Out", paste0("Out", 2:(L-1)))
  }
  labs <- labs[1:L]
  names(dat) <- labs

  grp <- droplevels(grp)
  atlas <- as.numeric(as.factor(atlas)) # convert atlas to numeric

  # if threshold specified, transform edges using threshold function
  # also check for predefined thresholds for ease of use
  if (class(threshold) == "character") {
    if (threshold == "positive") {
      threshold <- function(x) {ifelse(x>0, x, 0)}
    }
  }
  if (!is.null(threshold)) {
    dat_g <- lapply(dat, function(x)
      array(apply(x, 3, function(y) {
        thr <- threshold(y[lower.tri(y)])
        out <- matrix(0, dim(x), dim(x))
        out[lower.tri(out)] <- thr
        out + t(out)
      }), dim(x)))
  } else {
    dat_g <- dat
  }

  if (fisher) {
    dat_g <- lapply(dat_g, function(x)
      array(apply(x, 3, atanh), dim(x)))
  }

  #### Obtain group-level FC matrices by site ####
  switch(
    grp.method[1],
    "average" = {
      avg <- lapply(dat_g, function(x) {
        lapply(setNames(levels(grp), levels(grp)), function(b) {
          out <- apply(x[,,which(grp==b)], c(1,2), mean)
          diag(out) <- 1
          out
        })
      })
      gr <- lapply(avg, lapply, graph_from_adjacency_matrix,
                   mode = "undirected", weighted = TRUE)
    })

  #### Community detection algorithm ####
  clust <- lapply(gr, lapply, comm.method)

  #### Compare communities across site ####
  # compare alignment with atlas
  atl_comp <- lapply(clust, function(x) {
    setNames(
      sapply(1:length(x), function(i) igraph::compare(x[[i]], atlas, metric)),
      levels(grp))
  })

  # compare alignment with eachother
  grp_comp <- lapply(clust, function(x) {
    comp <- matrix(0, length(x), length(x),
                   dimnames = list(levels(grp), levels(grp)))
    for (i in 1:length(x)) {
      for (j in 1:length(x)) {
        if (i == j) {
          comp[i,j] <- 0
        } else {
          comp[i,j] = igraph::compare(x[[i]], x[[j]], metric)
        }
      }
    }
    comp
  })

  if (debug) {
    structure(
      list(communities = clust,
           atlas.comp = atl_comp,
           grp.comp = grp_comp,
           graphs = gr,
           atlas = atlas),
      class = "commTest")
  } else {
    structure(
      list(communities = clust,
           atlas.comp = atl_comp,
           grp.comp = grp_comp),
      class = "commTest")
  }
}
