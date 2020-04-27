#' Functional Connectivity Evaluations proposed by Meichen Yu
#' Refer to the 2018 paper for further details:
#' \url{https://doi.org/10.1002/hbm.24241}
#'
#' @param ...
#' @param bat
#'
#' @return
#' @import igraph doParallel glassoFast
#' @importFrom brainGraph efficiency
#' @importFrom stats kruskal.test p.adjust prcomp lm
#' @importFrom Matrix rowSums
#' @export
#'
#' @examples
test_yu <- function(..., bat = NULL, net.rois = NULL,
                    net.cov = NULL, threshold = NULL, labs = c("Raw", "Out"),
                    p.method = "BH", net.method = c("corr", "pcorr", "glasso"),
                    net.weighted = TRUE, net.only = TRUE,
                    glasso.rho = NULL, out.graphs = FALSE) {
  if (is.null(bat)) {stop("Need to specify batch")}
  dat <- list(...)
  L <- length(dat)

  if (length(dat) > length(labs)) {
    message("Not enough labels: assuming first is raw and using default labels")
    labs <- c("Raw", "Out", paste0("Out", 2:(L-1)))
  }
  labs <- labs[1:L]
  names(dat) <- labs

  all_vec <- lapply(dat, function(x)
    t(apply(x, 3, function(y) c(y[lower.tri(y)]))))

  #### Entry-wise Kruskal-Wallis ####
  if (!net.only) {
    kw_p <- lapply(all_vec, function(x) {
      apply(x, 2, function(y) kruskal.test(y, bat)$p.value)
    })
    kw_p_adj <- lapply(kw_p, p.adjust, p.method)
  } else {
    kw_p <- NULL; kw_p_adj <- NULL
  }

  ## Network measure associated with known covariate ####
  switch(
    net.method,
    "corr" = {
      dat <- lapply(dat, function(x) array(apply(x, 3, cov2cor), dim(x)))
    },
    "pcorr" = {
      dat <- lapply(dat, function(x) array(apply(x, 3, corr2pcor), dim(x)))
    },
    "glasso" = {
      # TODO: work on tuning parameter selection
      if (is.null(glasso.rho)) {stop("Need to select tuning parameter")}
      dat <- lapply(dat, function(x)
        array(apply(x, 3, function(y) {
          wi <- -glassoFast(y, glasso.rho)$wi
          diag(wi) <- -diag(wi)
          cov2cor(wi)
        }), dim(x)))
    }
  )

  # if threshold specified, remove edges with corr/partial corr below thres
  if (!is.null(threshold)) {
    dat <- lapply(dat, function(x)
      array(apply(x, 3, function(y) {
        ay <- abs(y)
        y[ay <= threshold] <- 0
      }), dim(x)))
  }

  # if unweighted, convert to unweighted adjacency matrices
  if (is.null(net.weighted)) {
    dat <- lapply(dat, function(x)
      array(apply(x, 3, function(y) {
        y <- abs(y)
        y[y > 0] <- 1
      }), dim(x)))
  }

  # using brainGraph
  gr <- lapply(dat, function(x)
    lapply(1:dim(x)[3], function(i) graph_from_adjacency_matrix(
      abs(x[,,i]), mode = "undirected", diag = FALSE, weighted = TRUE)))
  nodal <- lapply(gr, sapply, efficiency, "nodal")
  global <- lapply(nodal, colMeans)

  ## Calculating ROI group metrics
  local <- lapply(gr, sapply, local_eff, net.rois)
  local_sum <- lapply(local, colSums)
  dat_net <- lapply(dat, function(x) x[net.rois, net.rois,])
  # network connectivity, sum of FC values within network ROIs
  net_con <- lapply(dat_net, function(x) apply(x, 3, sum) - length(net.rois))

  #### Network measures associated with site
  net_res_site <- matrix(0, L, 3, dimnames = list(labs,
                                                  c("Network.connectivity.p",
                                                    "Global.Eff.p",
                                                    "Local.Eff.Sum.p")))
  netc_p <- lapply(net_con, function(x) kruskal.test(x, bat)$p.value)
  glob_p <- lapply(global, function(x) kruskal.test(x, bat)$p.value)
  local_p <- lapply(local_sum, function(x) kruskal.test(x, bat)$p.value)
  net_res_site[,1] <- unlist(netc_p)
  net_res_site[,2] <- unlist(glob_p)
  net_res_site[,3] <- unlist(local_p)

  #### Network measures associated with covariate (OPTIONAL)
  if (is.null(net.cov)) {
    net_res <- NULL
  } else {
    net_res <- matrix(0, L, 6, dimnames = list(labs,
                                               c("Network.connectivity.r2",
                                                 "Network.connectivity.p",
                                                 "Global.Eff.r2",
                                                 "Global.Eff.p",
                                                 "Local.Eff.Sum.r2",
                                                 "Local.Eff.Sum.p")))

    lm_netc <- lapply(net_con, function(x) lm(x ~ net.cov))
    lm_glob <- lapply(global, function(x) lm(x ~ net.cov))
    lm_local <- lapply(local_sum, function(x) lm(x ~ net.cov))

    net_summary <- function(x) c(summary(x)$r.squared, summary(x)$coefficients[2,4])
    net_res[,1:2] <- t(sapply(lm_netc, net_summary))
    net_res[,3:4] <- t(sapply(lm_glob, net_summary))
    net_res[,5:6] <- t(sapply(lm_local, net_summary))
  }

  if (!out.graphs) {gr <- NULL}

  list(
    fc.kw.adj = kw_p_adj,
    fc.kw.p = kw_p,
    net.results = net_res,
    net.results.site = net_res_site,
    eff.res = (list(nodal.eff = nodal,
                    local.eff = local,
                    global.eff = global,
                    local.sum = local_sum)),
    graphs = gr # temporary, for testing
  )
}
