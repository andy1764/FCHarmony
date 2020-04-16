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
                      tests = c("Elem", "Group", "Seed"),
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

  cov_p_mat <- list()
  if ("Elem" %in% tests) {
    vec <- lapply(dat, function(x) t(apply(x, 3, function(x) c(x[lower.tri(x, diag = FALSE)]))))
    if (fisher) {vec <- lapply(vec, atanh)}
    # get full linear model
    all_p <- lapply(vec, function(x) {
      matrix(apply(x, 2, function(y) anova(lm(y ~ mod))$'Pr(>F)'[1]), dim(x)[1], dim(x)[2])
    })
    # arrange as matrices
    cov_p_mat$all <- lapply(all_p, function(x) {
      d <- (1 + sqrt(8*ncol(x)+1))/2
      out <- array(c(apply(x, 1, vec2corr)), c(d, d, nrow(x)))})

    for (cov in mods) {
      cov_p <- lapply(vec, function(x) {
        matrix(apply(x, 2, function(y) summary(lm(y ~ mod[,cov]))$coefficients[2,4]), dim(x)[1], dim(x)[2])
      })
      cov_p_mat[[cov]] <- lapply(cov_p, function(x) {
        d <- (1 + sqrt(8*ncol(x)+1))/2
        out <- array(c(apply(x, 1, vec2corr)), c(d, d, nrow(x)))})
    }
  }

  list(
    elem.p <- cov_p_mat
  )
}
