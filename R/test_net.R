#' Functional Connectivity Network Metric Evaluations
#'
#' Some proposed by Meichen Yu. Refer to the 2018 paper for further details:
#' \url{https://doi.org/10.1002/hbm.24241}
#'
#' @param ...
#' @param bat
#' @param threshold Function that takes adjacency matrix input and returns an
#' adjacency matrix with thresholded values. Examples include `x > 0.5` for
#' right-tail thresholding and `ifelse(x>0, x, 0)` for keeping positive weights
#' @param mod.clust Function from `igraph` that finds communities. See
#' \link[igraph]{membership} for options. Defaults to
#' \link[igraph]{cluster_louvain}
#'
#' @param net.rois
#' @param net.cov
#' @param labs
#' @param net.method
#' @param net.metrics
#' @param add.metrics Additional metrics as named list indicating metric label.
#'   Each list in `add.metrics` should be a list of length `n` vectors, each
#'   corresponding to an input data frame, with the same names as given in
#'   `labs`
#' @param weighted
#' @param glasso.rho
#' @param out.graphs
#'
#' @return
#' @import igraph doParallel glassoFast
#' @importFrom brainGraph efficiency
#' @importFrom NetworkToolbox clustcoeff gateway participation
#' @importFrom stats kruskal.test p.adjust prcomp lm
#' @importFrom Matrix rowSums
#' @export
#'
#' @examples
test_net <- function(..., bat = NULL, net.rois = NULL,
                     net.cov = NULL, labs = c("Raw", "Out"),
                     net.method = c("corr", "pcorr", "glasso"),
                     net.metrics = c("global", "local", "within", "between",
                                     "clust.modularity", "atlas.modularity",
                                     "gateway.coeff", "part.coeff",
                                     "clust.coeff"),
                     add.metrics = NULL,
                     kw.p = FALSE,
                     threshold = NULL, fisher = TRUE,
                     mod.clust = cluster_louvain,
                     weighted = TRUE, glasso.rho = NULL, debug = FALSE) {
  if (is.null(bat)) {stop("Need to specify batch")}
  dat <- list(...)
  L <- length(dat)

  if (length(dat) > length(labs)) {
    message("Not enough labels: assuming first is raw and using default labels")
    labs <- c("Raw", "Out", paste0("Out", 2:(L-1)))
  }
  labs <- labs[1:L]
  names(dat) <- labs

  roi.names <- dimnames(dat[[1]])[[1]] # assumes same ROI labels throughout
  rois <- as.numeric(as.factor(roi.names))

  all_vec <- lapply(dat, function(x)
    t(apply(x, 3, function(y) c(y[lower.tri(y)]))))

  # #### Entry-wise Kruskal-Wallis ####
  # if (!net.only) {
  #   kw_p <- lapply(all_vec, function(x) {
  #     apply(x, 2, function(y) kruskal.test(y, bat)$p.value)
  #   })
  #   kw_p_adj <- lapply(kw_p, p.adjust, p.method)
  # } else {
  #   kw_p <- NULL; kw_p_adj <- NULL
  # }

  #### Choose correlation structure ####
  switch(
    net.method[1],
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

  # remove diagonal elements
  dat <- lapply(dat, function(x)
    array(apply(x, 3, function(y) {
      diag(y) <- 0
      y
    }), dim(x)))

  # if threshold specified, transform edges using threshold function
  # also check for predefined thresholds for ease of use
  if (class(threshold) == "character") {
    if (threshold == "positive") {
      threshold <- function(x) {ifelse(x<0, x, 0)}
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

  # if unweighted, every non-zero element should be 1 to indicate one edge
  if (is.null(weighted)) {
    dat_g <- lapply(dat_g, function(x)
      array(apply(x, 3, function(y) {
        y[y > 0] <- 1
        y[y < 0] <- 1
        y
      }), dim(x)))
  }

  ## Form graphs
  gr <- lapply(dat_g, function(x)
    lapply(1:dim(x)[3], function(i) graph_from_adjacency_matrix(
      x[,,i], mode = "undirected", diag = FALSE, weighted = weighted)))

  #### Calculate network measures ####
  met_all <- NULL
  lm_all <- NULL
  all_site <- NULL
  all_cov <- NULL

  sapply(
    net.metrics, switch,
    "global" = {
      nodal <- lapply(gr, sapply, efficiency, "nodal")
      global <- lapply(nodal, colMeans)
      met_all$nodal <- nodal
      met_all$global <- global

      met_out <- met_regress(global, bat, net.cov, net_lab = "Global",
                             labs = labs, kw.p = kw.p)
      lm_all$global <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "local" = {
      local <- lapply(gr, sapply, local_eff, net.rois, weighted)
      local_sum <- lapply(local, colSums)
      met_all$local <- local
      met_all$local.sum <- local_sum

      met_out <- met_regress(local_sum, bat, net.cov, net_lab = "Local",
                             labs = labs, kw.p = kw.p)
      lm_all$local <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "within" = {
      # use z-transformed correlation values
      dat_net <- lapply(dat, function(x) atanh(x[net.rois, net.rois,])/2)
      # within-network connectivity, average of FC values within
      within <- lapply(dat_net, function(x) apply(x, 3, mean))
      met_all$within <- within

      met_out <- met_regress(within, bat, net.cov, net_lab = "Within",
                             labs = labs, kw.p = kw.p)
      lm_all$within <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "between" = {
      # use z-transformed correlation values
      dat_net <- lapply(dat, function(x) atanh(x[net.rois, -(net.rois),])/2)
      # between-network connectivity, average of FC values between other ROIs
      between <- lapply(dat_net, function(x) apply(x, 3, mean))
      met_all$between <- between

      met_out <- met_regress(between, bat, net.cov, net_lab = "Between",
                             labs = labs, kw.p = kw.p)
      lm_all$between <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "clust.modularity" = {
      cmodul <- lapply(gr, sapply, function(x) modularity(mod.clust(x)))
      met_all$clust.modularity <- cmodul

      met_out <- met_regress(cmodul, bat, net.cov, net_lab = "Cl.Modul",
                             labs = labs, kw.p = kw.p)
      lm_all$cmodul <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "atlas.modularity" = {
      amodul <- lapply(gr, sapply, function(x)
        modularity(x, as.numeric(as.factor(roi.names))))
      met_all$atlas.modularity <- amodul

      met_out <- met_regress(amodul, bat, net.cov, net_lab = "A.Modul",
                             labs = labs, kw.p = kw.p)
      lm_all$amodul <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "part.coeff" = {
      # # old using brainGraph
      # part <- lapply(gr, sapply, function(x)
      #   part_coeff(x, as.numeric(as.factor(roi.names)))[net.rois[1]])
      part <- lapply(dat_g, apply, 3, function(x)
        mean(participation(x, comm = rois)$positive[net.rois]))
      met_all$part.coeff <- part

      met_out <- met_regress(part, bat, net.cov, net_lab = "Part.Coeff",
                             labs = labs, kw.p = kw.p)
      lm_all$part.coeff <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "gateway.coeff" = {
      # # old using brainGraph
      # gate <- lapply(gr, sapply, function(x)
      #   gateway_coeff(x, as.numeric(as.factor(roi.names)))[net.rois[1]])
      gate <- lapply(dat_g, apply, 3, function(x)
        mean(gateway(x, comm = rois)$positive[net.rois]))
      met_all$gate.coeff <- gate

      met_out <- met_regress(gate, bat, net.cov, net_lab = "Gate.Coeff",
                             labs = labs, kw.p = kw.p)
      lm_all$gate.coeff <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "clust.coeff" = {
      # clustc <- lapply(dat_g, apply, 3, function(x)
      #   clustcoeff(x, weighted)$CC)
      clustc <- lapply(gr, sapply, function(x)
        mean(igraph::transitivity(x, "weighted")))
      met_all$clust.coeff <- clustc

      met_out <- met_regress(clustc, bat, net.cov, net_lab = "Clust.Coeff",
                             labs = labs, kw.p = kw.p)
      lm_all$clust.coeff <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    },
    "weight.clust.coeff" = {
      # using tnet package clustering_w
      # https://toreopsahl.com/2009/04/03/article-clustering-in-weighted-networks/
      tnets <- lapply(gr, lapply, function(x) {
        tnet <- cbind(as_edgelist(x), edge_attr(x)$weight)
        colnames(tnet) <- c("i", "j", "w")
        tnet
      })
    }
    )

  if (!is.null(add.metrics)) {
    for (n in names(add.metrics)) {
      met <- add.metrics[[n]]
      met_all[[n]] <- met

      met_out <- met_regress(met, bat, net.cov, net_lab = n,
                             labs = labs, kw.p = kw.p)
      lm_all[[n]] <- met_out$lm.met
      all_cov <- c(all_cov, list(met_out$out.p))
      all_site <- c(all_site, met_out$bat.p)
    }
  }

  net_res_site <- do.call(cbind, all_site)
  rownames(net_res_site) <- labs

  if (debug) {
    list(net.results = do.call(cbind, all_cov),
         net.results.site = net_res_site,
         metrics = met_all,
         lm.res = lm_all,
         graphs = gr,
         graph.dat = dat_g)
  } else {
    list(net.results = do.call(cbind, all_cov),
         net.results.site = net_res_site)
  }
}

# Helper function to regress network metric on input covariates
met_regress <- function(metric, bat, net.cov = NULL, net_lab = NULL,
                        labs = NULL, kw.p = FALSE) {
  if (!is.null(net.cov)) {
    # Linear model for association with covariates
    lm_met <- lapply(metric, function(x) {
      tryCatch(lm(x ~ net.cov), error = function(cond) {
        message(cond)
        message(paste(": Linear model could not be fit for", net_lab))
        return(NULL)
      })
    })
    if (is.null(lm_met[[1]])) {return(list(out.p = NULL,
                                           bat.p = NULL,
                                           lm.met = NULL))}
    lm_met_bat <- lapply(metric, function(x) lm(x ~ net.cov+bat))
    bat_p <- lapply(labs, function(x)
      anova(lm_met[[x]], lm_met_bat[[x]], test = "F")$`Pr(>F)`[2])

    out <- t(sapply(lm_met_bat, net_summary))
    dimnames(out) <- list(labs, paste(net_lab, c("est", "p"), sep = "."))
  }

  # KW test for association with site
  if (kw.p) {
    bat_p <- lapply(metric, function(x) kruskal.test(x, bat)$p.value)
  }

  out_site <- list(unlist(bat_p))
  names(out_site) <- paste(net_lab, "p", sep=".")

  list(out.p = out,
       bat.p = out_site,
       lm.met = lm_met_bat)
}

# Helper function to summarize network linear models
net_summary <- function(x) c(summary(x)$coefficients[2,1],
                             summary(x)$coefficients[2,4])
