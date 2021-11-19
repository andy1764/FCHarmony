#' Regression tests
#'
#' Examine if harmonization affects batch and covariate associations. Options to
#' regress edge weights and mean edge weights within and between subnetworks.
#'
#' @param ... *p x p x n* covariance or correlation matrices where *p* is the number
#'   of ROIs and *n* is the number of subjects. These are generally
#'   unharmonized or harmonized datasets and should have the same dimensions.
#' @param bat Factor (or object coercible by \link[base]{as.factor} to a
#'    factor) of length *n* designating batch IDs.
#' @param mod Optional design matrix of covariates to regress on, usually from
#'   the output of \link[stats]{model.matrix}.
#' @param roi.names Vector of names for regions of interest. For "Group" tests,
#'   `roi.names` defines subnetworks of interest. If `NULL`, the dimension
#'   names of each dataset are used.
#' @param labs Vector of labels for the harmonization methods, in order of the
#'   inputted datasets.
#' @param tests Vector of tests to apply. "Elem" refers to regression on each
#'   edge weight. "Group" refers to regression on mean edge weights within and
#'   between each subnetwork, as defined by `roi.names`.
#' @param fisher Whether to z-transform input FC matrices
#' @param to.corr Logical, whether input should
#'   be forced to be a correlation matrix using  \link[stats]{cov2cor}
#' @param debug Whether to return intermediate objects for debugging
#'
#' @return
#' @export
#'
#' @examples
test_regress <- function(..., bat = NULL, mod = NULL, roi.names = NULL,
                      labs = c("Raw", "Out"),
                      tests = c("Elem", "Group"),
                      fisher = TRUE,
                      to.corr = FALSE,
                      debug = FALSE) {
  dat <- list(...)
  L <- length(dat)

  if (length(dat) > length(labs)) {
    message("Not enough labels: assuming first is raw and using default labels")
    labs <- c("Raw", "Out", paste0("Out", 2:(L-1)))
  }
  labs <- labs[1:L]
  names(dat) <- labs

  if (to.corr) {
    dat <- lapply(dat, function(x) array(apply(x, 3, cov2cor), dim(x)))
  }

  if (is.null(bat)) {
    design <- mod
    mods <- colnames(design)[-1]
  } else {
    design <- cbind(mod, model.matrix(~bat)[,-1])
    mods <- colnames(design)[-1]
  }

  # if no ROI names specified, take from first dataset
  # otherwise, replace dimnames across all datasets with the ROI names
  if (is.null(roi.names)) {
    roi.names <- dimnames(dat[[1]])[1:2]
  } else {
    dat <- lapply(dat, function(x) {dimnames(x)[1:2] <- roi.names; x})
  }

  cov_p_mat <- list()
  grp_p_mat <- list()
  cpc_p_mat <- list()

  elem_reg <- NULL
  cpc_reg <- NULL
  grp_reg <- NULL
  cpc_bat <- NULL

  out <- list()

  sapply(
    tests, switch,
    "Elem" = {
      vec <- lapply(dat, function(x) t(apply(x, 3, function(x)
        c(x[lower.tri(x, diag = FALSE)]))))
      if (fisher) {vec <- lapply(vec, atanh)}

      # get full linear model
      elem_reg <- lapply(vec, function(x) {
        apply(x, 2, function(y) lm(y ~ mod))
      })

      if (!is.null(bat)) {
        elem_reg_bat <- lapply(vec, function(x) {
          apply(x, 2, function(y) lm(y ~ design))
        })
        bat_p <- lapply(unique(names(dat)), function(x)
          sapply(1:length(elem_reg[[x]]), function(y)
            anova(elem_reg[[x]][[y]],
                  elem_reg_bat[[x]][[y]], test = "F")$`Pr(>F)`[2]))
        names(bat_p) <- names(dat)
        elem_reg <- elem_reg_bat
      }

      elem_p <- lapply(elem_reg, sapply, function(x) anova(x)$'Pr(>F)'[1])

      # arrange as matrices
      cov_p_mat$all <- lapply(elem_p, vec2corr, roi.names)
      cov_p_mat$bat <- lapply(bat_p, vec2corr, roi.names)

      for (cov in mods) {
        i <- which(mods == cov)
        cov_p <- lapply(elem_reg, sapply,
                        function(x) summary(x)$coefficients[i+1, 4])
        cov_p_mat[[cov]] <- lapply(cov_p, vec2corr, roi.names)
      }

      out$elem.p <- cov_p_mat
    },
    "Group" = {
      grp_rois <- sort(unique(dimnames(dat[[1]])[[1]]))
      grp_rois <- list(grp_rois, grp_rois)
      dat_grouped <- lapply(dat, function(x) {
        rois <- sort(unique(dimnames(x)[[1]]))
        array(apply(x, 3, corr2con, fisher = fisher), dim = c(length(rois), length(rois),
                                             dim(x)[3]),
              dimnames = list(rois, rois, NULL))
      })
      vec <- lapply(dat_grouped, function(x)
        t(apply(x, 3, function(x) c(x[lower.tri(x, diag = TRUE)]))))

      # get full linear model
      grp_reg <- lapply(vec, function(x) {
        apply(x, 2, function(y) lm(y ~ mod))
      })

      if (!is.null(bat)) {
        grp_reg_bat <- lapply(vec, function(x) {
          apply(x, 2, function(y) lm(y ~ design))
        })
        bat_p <- lapply(unique(names(dat)), function(x)
          sapply(1:length(grp_reg[[x]]), function(y)
            anova(grp_reg[[x]][[y]],
                  grp_reg_bat[[x]][[y]], test = "F")$`Pr(>F)`[2]))
        names(bat_p) <- names(dat)
        grp_reg <- grp_reg_bat
      }

      grp_p <- lapply(grp_reg, sapply, function(x) anova(x)$'Pr(>F)'[1])
      # arrange as matrices
      grp_p_mat$all <- lapply(grp_p, dvec2corr, grp_rois)
      grp_p_mat$bat <- lapply(bat_p, dvec2corr, grp_rois)

      for (cov in mods) {
        i <- which(mods == cov)
        cov_p <- lapply(grp_reg, sapply,
                        function(x) summary(x)$coefficients[i+1, 4])
        grp_p_mat[[cov]] <- lapply(cov_p, dvec2corr, grp_rois)
      }

      out$group.p <- grp_p_mat
    }
  )

  if (debug) {
    c(out, list(elem.reg = elem_reg,
                group.reg = grp_reg,
                cpc.reg = cpc_reg))
  } else {
    out
  }
}
