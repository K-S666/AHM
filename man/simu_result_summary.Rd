\name{simu_result_summary}
\alias{simu_result_summary}
\alias{print.AHMQ_simulation_summary}
\title{Summarize simulation recovery}
\usage{
simu_result_summary(est, truth = NULL, s0 = NULL, g0 = NULL, Q0 = NULL,
  G0 = NULL, alpha0 = NULL, cut_value = NULL, align_labels = TRUE,
  keep_aligned_samples = FALSE, verbose = FALSE)
}
\arguments{
\item{est}{Output from \code{\link{Est_fun}}. For unknown-Q
\code{\link{AHMQ}} output, if \code{align_labels = TRUE}, it must be created
with \code{return_samples = TRUE}. Fixed-Q \code{\link{AHM}} output is not
truth-label aligned because Q anchors the labels.}
\item{truth}{Optional list returned by \code{\link{simulate_ahmq_data}}.
When supplied, \code{s0}, \code{g0}, \code{Q0}, \code{G0}, and
\code{alpha0} are read from this object unless explicitly supplied.}
\item{s0, g0}{True slip and guess parameter vectors.}
\item{Q0}{True Q-matrix. Required for unknown-Q \code{\link{AHMQ}} recovery
metrics; optional and ignored for fixed-Q \code{\link{AHM}} summaries.}
\item{G0}{True attribute hierarchy adjacency matrix.}
\item{alpha0}{True examinee attribute profile matrix.}
\item{cut_value}{Threshold for binarizing posterior mean \code{G}. Defaults
to \code{est$cut_value}, then \code{0.2}.}
\item{align_labels}{Whether to align estimated labels to the true simulation
labels before computing recovery criteria. Ignored for fixed-Q
\code{\link{AHM}} output.}
\item{keep_aligned_samples}{Whether to keep the truth-aligned posterior
samples in \code{aligned_estimates}.}
\item{verbose}{Whether to print the selected truth-label permutation.}
}
\value{
An object of class \code{"AHMQ_simulation_summary"} containing:
\code{aligned_estimates}, item-parameter bias summaries for \code{s} and
\code{g}, \code{AAR}, \code{PAR}, G recovery metrics (\code{TPR},
\code{TFR}, \code{RTPR}, \code{RTFR}), optional Q recovery metrics
(\code{QRR}, \code{AQRR}) for unknown-Q models, Rhat diagnostics, and
truth-label alignment metadata when applicable.
}
\description{
Computes simulation-only evaluation criteria. For unknown-Q \code{\link{AHMQ}}
output, it can align the estimated attribute order to the known true order.
This prevents equivalent label permutations, such as a true chain
\code{1 -> 2 -> 3} estimated as \code{2 -> 3 -> 1}, from being counted as
recovery errors. For fixed-Q \code{\link{AHM}} output, Q is known and anchors
the labels, so this relabeling and Q recovery are skipped.
}
\details{
\code{summary_AHMQ()} remains the generic posterior summary for both real and
simulated data and does not perform truth-label alignment. Use
\code{simu_result_summary()} when true simulation quantities are available
and recovery rates or parameter biases are needed.

For unknown-Q output, the print method reports the selected truth-label
permutation, then the truth-aligned estimated \code{G} and its recovery
metrics, followed by the truth-aligned structured \code{Q} and Q recovery
metrics. For fixed-Q output, it reports G recovery, item-parameter bias, and
attribute recovery only.
}
