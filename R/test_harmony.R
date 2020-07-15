#' Evaluate harmonization methods
#'
#' Apply downstream analysis to multiple datasets, usually outputs of different harmonization methods.
#'
#' @param dats Named list of functional connectivity arrays, not necessarily of
#'   the same dimension
#' @param test Function to apply to all datasets
#' @param to.corr Logical, whether input should
#'   be forced to be a correlation matrix using  \link[stats]{cov2cor}
#' @param ... Arguments passed to `test`
#'
#' @return
#' @export
#'
#' @examples
test_harmony <- function(dats, test, force.PD = FALSE, to.corr = FALSE, ...) {
  if (force.PD) {
    dats <- lapply(dats, function(x) array(apply(x, 3, function(y) {
      as.numeric(nearPD(y)$mat)
    }), dim(x)))
  }

  if (to.corr) {
    dats <- lapply(dats, function(x) array(apply(x, 3, cov2cor), dim(x)))
  }

  params <- list(...)

  out <- lapply(dats, function(x) do.call(test, c(list(x), params)))
  out
}
