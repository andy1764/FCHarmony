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
test_regress <- function(..., bat = NULL, mod = NULL,
                      labs = c("Raw", "Out"),
                      tests = c("Elem", "Group", "CPC"),
                      cpc.k = 15,
                      fisher = TRUE,
                      seed.roi = NULL, to.corr = FALSE) {
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
  mods <- colnames(mod)[-1]
  roi_names <- dimnames(dat[[1]])[1:2]

  cov_p_mat <- list()
  grp_p_mat <- list()
  cpc_p_mat <-

  if ("Elem" %in% tests) {
    vec <- lapply(dat, function(x) t(apply(x, 3, function(x) c(x[lower.tri(x, diag = FALSE)]))))
    if (fisher) {vec <- lapply(vec, atanh)}
    # get full linear model
    all_p <- lapply(vec, function(x) {
      apply(x, 2, function(y) anova(lm(y ~ mod))$'Pr(>F)'[1])
    })
    # arrange as matrices
    cov_p_mat$all <- lapply(all_p, vec2corr, roi_names)

    for (cov in mods) {
      cov_p <- lapply(vec, function(x) {
        apply(x, 2, function(y) summary(lm(y ~ mod[,cov]))$coefficients[2,4])
      })
      cov_p_mat[[cov]] <- lapply(cov_p, vec2corr, roi_names)
    }
  }
  if ("Group" %in% tests) {
    roi_names <- sort(unique(dimnames(dat[[1]])[[1]]))
    roi_names <- list(roi_names, roi_names)
    dat_grouped <- lapply(dat, function(x) {
      rois <- sort(unique(dimnames(x)[[1]]))
      array(apply(x, 3, corr2con), dim = c(length(rois), length(rois),
                                           dim(x)[3]),
            dimnames = list(rois, rois, NULL))
    })
    vec <- lapply(dat_grouped, function(x)
      t(apply(x, 3, function(x) c(x[lower.tri(x, diag = TRUE)]))))
    # get full linear model
    all_p <- lapply(vec, function(x) {
      apply(x, 2, function(y) anova(lm(y ~ mod))$'Pr(>F)'[1])
    })
    # arrange as matrices
    grp_p_mat$all <- lapply(all_p, dvec2corr)

    for (cov in mods) {
      cov_p <- lapply(vec, function(x) {
        apply(x, 2, function(y) summary(lm(y ~ mod[,cov]))$coefficients[2,4])
      })
      grp_p_mat[[cov]] <- lapply(cov_p, dvec2corr)
    }
  }
  if ("CPC" %in% tests) {
    all_cpc <- lapply(dat, cpc, k = cpc.k)
    cpc_p <- lapply(all_cpc, function(x) {
      lapply(1:dim(x$D)[1], function(y) {
        summary(glm(x$D[y,] ~ mod[,-1], family = Gamma(link = "log")))
      })
    })
  }

  list(
    elem.p = cov_p_mat,
    group.p = grp_p_mat
  )
}
