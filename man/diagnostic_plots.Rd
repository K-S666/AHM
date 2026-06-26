\name{compute_Rhat_curve}
\alias{compute_Rhat_curve}
\alias{plot_Rhat_curve}
\alias{plot_G_graph}
\alias{traceplot_AHMQ}
\title{Diagnostic and hierarchy visualization helpers}
\usage{
compute_Rhat_curve(result, parameters = c("s", "g"), step = 1000L,
  min_draws = step, max_draws = NULL, include_burnin = FALSE)

plot_Rhat_curve(curve, parameters = NULL, threshold = 1.1, cutoff = NULL,
  file = NULL, width = 1200, height = 800, res = 130)

plot_G_graph(G, G_posterior = NULL, type = c("graph", "heatmap", "both"),
  cut_value = 0.2, labels = NULL, file = NULL, width = 1000,
  height = 700, res = 130, main = "Estimated hierarchy G",
  node_col = "white", edge_col = "gray30")

traceplot_AHMQ(result, parameter = c("s", "g"), items = NULL,
  chains = NULL, post_idx = NULL, include_burnin = TRUE,
  discard_burnin = FALSE, panels_per_page = 12L, file = NULL,
  width = 1400, height = 900, res = 130)
}
\arguments{
\item{result}{Output from \code{\link{AHMQ}}.}
\item{parameters}{Chain components for Rhat diagnostics or plotting.}
\item{step}{Prefix spacing, in stored draws.}
\item{min_draws}{First prefix length to evaluate.}
\item{max_draws}{Last prefix length; \code{NULL} uses all stored draws.}
\item{include_burnin}{Whether to include burn-in draws when full chains were
saved with \code{keep_burnin = TRUE}. For Rhat curves the default
\code{FALSE} uses post-burn-in draws, matching the final reported Rhat.}
\item{curve}{Output from \code{compute_Rhat_curve}.}
\item{threshold}{Horizontal reference line for Rhat plots. The default 1.1 is
the usual Rhat convergence threshold.}
\item{cutoff}{Deprecated alias for \code{threshold}.}
\item{file}{Optional PNG file path. If \code{NULL}, the active graphics
device is used.}
\item{width, height, res}{PNG device settings.}
\item{G}{Binary hierarchy adjacency matrix, posterior probability matrix, or
an object returned by \code{\link{Est_fun}} or \code{\link{extract_estimates}}.
The matrix is interpreted as \code{G[prerequisite, target] = 1}.}
\item{G_posterior}{Optional posterior edge-probability matrix. If omitted and
\code{G} is an estimate object, the stored posterior mean of G is used.}
\item{type}{Plot type: \code{"graph"}, \code{"heatmap"}, or \code{"both"}.}
\item{cut_value}{Threshold used to binarize \code{G_posterior} for the graph
when a binary G matrix is not supplied.}
\item{labels}{Optional labels for the attributes.}
\item{main}{Plot title for \code{plot_G_graph}.}
\item{node_col}{Node fill color.}
\item{edge_col}{Arrow color.}
\item{items}{Item indices for trace plots. \code{NULL} plots all items.}
\item{chains}{Chain indices for trace plots. \code{NULL} plots all chains.}
\item{post_idx}{Optional stored draw indices for trace plots. If supplied, it
overrides \code{include_burnin}.}
\item{include_burnin}{Whether to include burn-in draws when full chains were
saved with \code{keep_burnin = TRUE}. If \code{TRUE}, a vertical dashed line
marks the burn-in boundary.}
\item{discard_burnin}{Deprecated alias. Use \code{include_burnin = FALSE}.}
\item{panels_per_page}{Maximum number of item panels per trace-plot page.}
}
\value{
\code{compute_Rhat_curve()} returns a list with data frames \code{max} and
\code{item}. \code{plot_Rhat_curve()} invisibly returns the curve object.
\code{plot_G_graph()} invisibly returns node coordinates.
\code{traceplot_AHMQ()} invisibly returns a data frame listing produced pages.
}
\description{
Utilities for studying convergence and presenting the estimated hierarchy.
The Rhat curve functions support both ordinary \code{AHMQ} results that store
only post-burn-in draws and results created with \code{keep_burnin = TRUE}.
Final reported Rhat values usually use post-burn-in draws. Set
\code{include_burnin = TRUE} only when studying the full warm-up trajectory.
}
\examples{
\dontrun{
fit <- AHMQ(Y, K = 4, keep_burnin = TRUE,
            chain_length = 30000, burn_in = 10000)
curve <- compute_Rhat_curve(fit, parameters = c("s", "g"), step = 1000)
curve_full <- compute_Rhat_curve(fit, parameters = c("s", "g"), step = 1000,
                                 include_burnin = TRUE)
plot_Rhat_curve(curve, file = "rhat_curve.png")

est <- extract_estimates(fit)
plot_G_graph(est$est_G, file = "G_graph.png")
plot_G_graph(est, type = "both", file = "G_both.png")
traceplot_AHMQ(fit, parameter = "s", items = 1:6,
               file = "trace_s", panels_per_page = 6)
traceplot_AHMQ(fit, parameter = "s", items = 2,
               include_burnin = FALSE)
}
}
