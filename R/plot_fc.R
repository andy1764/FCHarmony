#' Plot functional connectivity matrix
#'
#' Function to create ggplots of either functional connectivity values or
#' p-value matrices from functions such as \link[FCharmony]{test_regress}
#'
#' @param cov Matrix of either functional connectivity values or derived
#'   p-values. Both covariance and correlation matrices are accepted.
#' @param lims Lower and upper bounds to plot. Defaults to being
#' @param subgraphs Vector of subgraph labels. Defaults to the input dimension
#'   names.
#' @param p.val Are the elements p-values?
#' @param p.method Method input to \link[stats]{p.adjust}. One of
#'   \link[stats]{p.adjust.methods}.
#' @param alpha p-value threshold. If `NULL`, cutoffs are not plotted.
#' @param binary Threshold the p-values based on specified threshold.
#' @param diag Whether to plot diagonal elements.
#' @param log.p Apply negative log-transformation to p-values, overridden by
#'   `binary = TRUE`.
#' @param bin.param List of graphical parameters passed to
#'   \link[ggplot2]{geom_text}
#' @param rect.param List of graphical parameters passed to
#'   \link[ggplot2]{geom_rect}
#'
#' @return
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
plot_fc <- function(cov, lims = c(-max(abs(cov)), max(abs(cov))),
                    subgraphs = dimnames(cov)[[1]], p.val = FALSE,
                    p.method = "BH", alpha = 0.05, binary = FALSE,
                    diag = FALSE, log.p = FALSE,
                    bin.param = list(color = "red", size = 2.5),
                    rect.param = list(alpha = 0, size = 1.25)) {
  dat <- cov # avoids problems with modification of cov changing the arguments
  dimnames(dat) <- list(NULL, NULL)

  # p-value transformations/binarization
  if (p.val) {
    colors <- c("white", "blue")

    p_mat <- cov
    if (diag) {
      p <- p.adjust(cov[lower.tri(cov, diag = diag)], method = p.method)
      p_mat[] <- 0
      p_mat[lower.tri(p_mat, diag = diag)] <- p
      p_mat <- p_mat + t(p_mat)
      diag(p_mat) <- diag(p_mat)/2
    } else {
      p <- p.adjust(cov[lower.tri(cov)], method = p.method)
      p_mat[] <- 0
      p_mat[lower.tri(p_mat)] <- p
      p_mat <- p_mat + t(p_mat)
      diag(p_mat) <- 1
    }

    if (log.p) {dat <- -log(cov)}
    if (binary) {dat[] <- as.numeric(p_mat < alpha)} # threshold for sig
  } else {
    colors <- c("red", "white", "blue")
  }

  cov_melt <- melt(dat)

  # Calculate p-value cutoffs
  if (is.null(alpha) | !p.val) {
    cov_melt$stars <- ""
  } else {
    dimnames(p_mat) <- list(NULL, NULL)
    p_melt <- melt(p_mat)
    stars <- cut(p_melt$value, breaks = c(0, alpha, 1), label=c("+",""))
    cov_melt$stars <- stars
  }

  if (is.null(subgraphs)) {
    ggplot(data = cov_melt) +
      geom_tile(aes(x=Var1, y=Var2, fill=value)) +
      do.call(geom_text, c(aes(x=Var1, y=Var2, label=stars), bin.param)) +
      scale_fill_gradientn(colours = colors, limits = lims) +
      labs(fill = "") +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank())
  } else {
    sub <- subgraphs
    dat <- dat[order(sub), order(sub)]
    xmin <- sapply(unique(sub), function(x) min(which(sub == x)))-.5
    xmax <- sapply(unique(sub), function(x) max(which(sub == x)))+.5
    Subgraph <- unique(sub)

    # for plotting rectangles
    rect_df <- data.frame(xmin, xmax, Subgraph)

    ggplot(data = cov_melt) +
      geom_tile(aes(x=Var1, y=Var2, fill=value)) +
      do.call(geom_text, args = c(bin.param,
                                  list(aes(x=Var1, y=Var2, label=stars)))) +
      do.call(geom_rect,
              args = c(rect.param, list(data = rect_df, inherit.aes = FALSE,
                          aes(xmin = xmin, xmax = xmax, ymin = xmin,
                              ymax = xmax, color = Subgraph)))) +
      scale_fill_gradientn(colours = colors, limits = lims) +
      labs(fill = "") +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }
}
