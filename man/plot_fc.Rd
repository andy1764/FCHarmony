% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fc.R
\name{plot_fc}
\alias{plot_fc}
\title{Plot functional connectivity matrix}
\usage{
plot_fc(
  cov,
  lims = c(-max(abs(cov)), max(abs(cov))),
  subgraphs = dimnames(cov)[[1]],
  p.val = FALSE,
  p.method = "BH",
  alpha = 0.05,
  binary = FALSE,
  diag = FALSE,
  log.p = FALSE,
  starchar = "+",
  bin.param = list(color = "red", size = 2.5),
  rect.param = list(alpha = 0, size = 1.25)
)
}
\arguments{
\item{cov}{Matrix of either functional connectivity values or derived
p-values. Both covariance and correlation matrices are accepted.}

\item{lims}{Lower and upper bounds to plot. Defaults to being}

\item{subgraphs}{Vector of subgraph labels. Defaults to the input dimension
names.}

\item{p.val}{Are the elements p-values?}

\item{p.method}{Method input to \link[stats]{p.adjust}. One of
\link[stats]{p.adjust.methods}.}

\item{alpha}{p-value threshold. If \code{NULL}, cutoffs are not plotted.}

\item{binary}{Threshold the p-values based on specified threshold.}

\item{diag}{Whether to plot diagonal elements.}

\item{log.p}{Apply negative log-transformation to p-values, overridden by
\code{binary = TRUE}.}

\item{starchar}{Which character to use for annotating significant elements.}

\item{bin.param}{List of graphical parameters passed to
\link[ggplot2]{geom_text}}

\item{rect.param}{List of graphical parameters passed to
\link[ggplot2]{geom_rect}}
}
\value{

}
\description{
Function to create ggplots of either functional connectivity values or
p-value matrices from functions such as \link[FCharmony]{test_regress}
}
