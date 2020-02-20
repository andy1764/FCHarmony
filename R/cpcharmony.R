#' Harmonization using Common Principal Components
#'
#' @param dat array of covariance matrices of dimension \eqn{(p \times p \times
#' N)}
#' @param bat batch labels as length \eqn{n} vector
#' @param cpc.k number of common principal components to harmonize
#'
#' @return
#'
#' @importFrom cpca cpc
#' @export
#'
#' @examples

cpcharmony <- function(dat, bat, cpc.k) {
  cpc(dat, k = cpc.k) # get common principal components
}
