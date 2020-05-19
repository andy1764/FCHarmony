#' Community detection tests
#' Examine if site effects exist in community detection
#'
#' @param ...
#' @param bat
#' @param labs
#' @param atlas Vector of labels specifying original subnetworks.
#' @param grp.method Method for obtaining group-level functional connectivity.
#' @param comm.method Function from `igraph` that finds communities. See
#' \link[igraph]{membership} for options. Defaults to
#' \link[igraph]{cluster_louvain}
#' @param comp.method Metric for comparing community detections.
#' @param fisher Whether to z-transform input FC matrices
#' @param debug
#'
#' @return
#' @import igraph
#' @export
#'
#' @examples
test_communities <- function(..., labs = c("Raw", "Out"), bat, atlas,
                             grp.method = c("average", "multilayer"),
                             comm.method = cluster_louvain,
                             metric = c("adjusted.rand", "nmi"),
                             fisher = TRUE,
                             debug = FALSE) {
  if (is.null(bat)) {stop("Need to specify batch")}
  dat <- list(...)
  L <- length(dat)

  if (length(dat) > length(labs)) {
    message("Not enough labels: assuming first is raw and using default labels")
    labs <- c("Raw", "Out", paste0("Out", 2:(L-1)))
  }
  labs <- labs[1:L]
  names(dat) <- labs

  bat <- droplevels(bat)
  atlas <- as.numeric(as.factor(atlas)) # convert atlas to numeric

  #### Obtain group-level FC matrices by site ####
  switch(
    grp.method[1],
    "average" = {
      avg <- lapply(dat, function(x) {
        lapply(setNames(levels(bat), levels(bat)), function(b) {
          out <- apply(x[,,which(bat==b)], c(1,2), mean)
          out[out < 0] <- 0
          if (fisher) {out <- atanh(out)}
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
      levels(bat))
  })

  # compare alignment with eachother
  site_comp <- lapply(clust, function(x) {
    comp <- matrix(0, length(x), length(x),
                   dimnames = list(levels(bat), levels(bat)))
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
           bat.comp = site_comp,
           graphs = gr,
           atlas = atlas),
      class = "commTest")
  } else {
    structure(
      list(communities = clust,
           atlas.comp = atl_comp,
           bat.comp = site_comp),
      class = "commTest")
  }
}
