#' Title
#'
#' @param ...
#' @param bat
#' @param mod
#' @param labs
#' @param tests
#' @param seed.roi
#' @param to.corr
#'
#' @return
#' @export
#'
#' @examples
test_regress <- function(..., bat = NULL, mod = NULL, roi.names = NULL,
                      labs = c("Raw", "Out"),
                      tests = c("Elem", "Group", "CPC"),
                      incl.bat = TRUE, # include batch in the model
                      cpc.k = 15,
                      fisher = TRUE,
                      seed.roi = NULL, to.corr = FALSE, debug = FALSE) {
  if (is.null(bat)) {stop("Need to specify batch")}
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

  bat <- droplevels(bat)
  if (incl.bat) {
    mod_bat <- cbind(mod, model.matrix(~bat)[,-1])
    mods <- colnames(mod_bat)[-1]
  } else {
    mods <- colnames(mod)[-1]
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

      if (incl.bat) {
        elem_reg_bat <- lapply(vec, function(x) {
          apply(x, 2, function(y) lm(y ~ mod_bat))
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
    # TODO: Group test is still very preliminary, boils down to being a
    # within and between connectivity test
    "Group" = {
      grp_rois <- sort(unique(dimnames(dat[[1]])[[1]]))
      grp_rois <- list(grp_rois, grp_rois)
      dat_grouped <- lapply(dat, function(x) {
        rois <- sort(unique(dimnames(x)[[1]]))
        array(apply(x, 3, corr2con), dim = c(length(rois), length(rois),
                                             dim(x)[3]),
              dimnames = list(rois, rois, NULL))
      })
      vec <- lapply(dat_grouped, function(x)
        t(apply(x, 3, function(x) c(x[lower.tri(x, diag = TRUE)]))))
      if (fisher) {vec <- lapply(vec, atanh)}

      # get full linear model
      grp_reg <- lapply(vec, function(x) {
        apply(x, 2, function(y) lm(y ~ mod))
      })

      if (incl.bat) {
        grp_reg_bat <- lapply(vec, function(x) {
          apply(x, 2, function(y) lm(y ~ mod_bat))
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
    },
    "CPC" = {
      all_cpc <- lapply(dat, cpc, k = cpc.k)
      cpc_reg <- lapply(all_cpc, function(x) {
        lapply(1:dim(x$D)[1], function(y) {
          glm(x$D[y,] ~ mod, family = gaussian(link = "log"))
        })
      })
      if (incl.bat) {
        cpc_reg_bat <- lapply(all_cpc, function(x) {
          lapply(1:dim(x$D)[1], function(y) {
            glm(x$D[y,] ~ mod_bat, family = gaussian(link = "log"))
          })
        })
        bat_p <- lapply(unique(names(dat)), function(x)
          sapply(1:length(cpc_reg[[x]]), function(y)
            anova(cpc_reg[[x]][[y]],
                  cpc_reg_bat[[x]][[y]], test = "LRT")$`Pr(>Chi)`[2]))
        names(bat_p) <- names(dat)
        cpc_reg <- cpc_reg_bat
        cpc_bat <- bat_p
      }

      cpc_p <- lapply(cpc_reg, lapply,
                      function(x) summary(x)$coefficients[-1,4])
      cpc_p_mat <- lapply(cpc_p, function(x)
        matrix(do.call(rbind, x), length(x), length(mods),
               dimnames = list(1:length(x), mods)))

      out$cpcs <- all_cpc
      out$cpc.p <- cpc_p_mat
      out$cpc.bat <- cpc_bat
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
