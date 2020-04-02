#' Functional Connectivity Evaluations proposed by Meichen Yu
#' Refer to the 2018 paper for further details:
#' \url{https://doi.org/10.1002/hbm.24241}
#'
#' @param raw
#' @param out
#' @param bat
#'
#' @return
#' @import igraph
#' @importFrom brainGraph efficiency
#' @importFrom stats kruskal.test p.adjust prcomp lm
#' @importFrom Matrix rowSums
#' @export
#'
#' @examples
test_yu <- function(raw, out, bat, net.rois, net.cov, p.method = "BH",
                    to.corr = TRUE, net.only = TRUE) {
  p <- dim(raw)[1]
  N <- dim(raw)[3]
  if (to.corr) {
    raw <- sapply(seq_along(raw[1,1,]), function(x) cov2cor(raw[,,x]),
                  simplify = "array")
    out <- sapply(seq_along(out[1,1,]), function(x) cov2cor(out[,,x]),
                  simplify = "array")
  }
  all_vec <- list(raw = t(apply(raw, 3, function(x) c(x[lower.tri(x)]))),
                  out = t(apply(out, 3, function(x) c(x[lower.tri(x)]))))

  #### Entry-wise Kruskal-Wallis ####
  if (!net.only) {
    kw_p <- lapply(all_vec, function(x) {
      apply(x, 2, function(y) kruskal.test(y, bat)$p.value)
    })
    kw_p_adj <- lapply(kw_p, p.adjust, p.method)
  } else {
    kw_p <- NULL; kw_p_adj <- NULL
  }


  #### PC regression using upper-triangular elements ####

  ## Network measure associated with known covariate ####
  # using brainGraph
  raw_net <- raw[net.rois, net.rois,]
  out_net <- out[net.rois, net.rois,]

  # network connectivity, sum of FC values within network ROIs
  net_con <- list(raw = apply(raw_net, 3, sum) - length(net.rois),
                  out = apply(out_net, 3, sum) - length(net.rois))

  gr <- list(raw = NULL, out = NULL)
  for (i in 1:N) {
    gr$raw[[i]] <- graph_from_adjacency_matrix(
      abs(raw_net[,,i]), mode = "undirected", diag = FALSE, weighted = TRUE)
    gr$out[[i]] <- graph_from_adjacency_matrix(
      abs(out_net[,,i]), mode = "undirected", diag = FALSE, weighted = TRUE)
  }
  nodal <- lapply(gr, sapply, efficiency, "nodal", use.parallel = FALSE)
  global <- lapply(nodal, colMeans)

  local <- lapply(gr, sapply, local_eff, use.parallel = FALSE)
  local_sum <- lapply(local, colSums)

  lm_netc <- lapply(net_con, function(x) lm(x ~ net.cov))
  lm_glob <- lapply(global, function(x) lm(x ~ net.cov))
  lm_local <- lapply(local_sum, function(x) lm(x ~ net.cov))

  net_res <- matrix(0, 2, 6, dimnames = list(c("Raw", "Out"),
                                            c("Network.connectivity.r2",
                                              "Network.connectivity.p",
                                              "Global.Eff.r2",
                                              "Global.Eff.p",
                                              "Local.Eff.Sum.r2",
                                              "Local.Eff.Sum.p")))
  net_summary <- function(x) c(summary(x)$r.squared, summary(x)$coefficients[2,4])
  net_res[,1:2] <- t(sapply(lm_netc, net_summary))
  net_res[,3:4] <- t(sapply(lm_glob, net_summary))
  net_res[,5:6] <- t(sapply(lm_local, net_summary))

  netc_p <- lapply(net_con, function(x) kruskal.test(x, bat)$p.value)
  glob_p <- lapply(global, function(x) kruskal.test(x, bat)$p.value)
  local_p <- lapply(local_sum, function(x) kruskal.test(x, bat)$p.value)

  net_res_site <- matrix(0, 2, 3, dimnames = list(c("Raw", "Out"),
                                                  c("Network.connectivity.p",
                                                    "Global.Eff.p",
                                                    "Local.Eff.Sum.p")))
  net_res_site[,1] <- unlist(netc_p)
  net_res_site[,2] <- unlist(glob_p)
  net_res_site[,3] <- unlist(local_p)

  return(list(
    fc.kw.adj = kw_p_adj,
    fc.kw.p = kw_p,
    net.results = net_res,
    net.results.site = net_res_site,
    nodal.eff = nodal,
    local.eff = local,
    global.eff = global,
    local.sum = local_sum
  ))
}
