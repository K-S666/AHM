\name{align_samples_to_truth}
\alias{align_samples_to_truth}
\alias{align_estimates_to_truth}
\title{Align simulation estimates to known true attribute labels}
\usage{
align_samples_to_truth(Q_sample, G_sample, alpha_sample, pi_sample, Q0,
  G0 = NULL, cut_value = 0.2, K = ncol(Q0), verbose = FALSE)

align_estimates_to_truth(est, Q0, G0 = NULL, cut_value = NULL,
  verbose = FALSE, keep_samples = TRUE)
}
\arguments{
\item{Q_sample}{Item-by-attribute Q-matrix posterior draws, an array of
dimension \eqn{J \times K \times T}.}
\item{G_sample}{Attribute hierarchy posterior draws, an array of dimension
\eqn{K \times K \times T}.}
\item{alpha_sample}{Attribute profile posterior draws, an array of dimension
\eqn{N \times K \times T}.}
\item{pi_sample}{Latent class probability posterior draws, a matrix of
dimension \eqn{2^K \times T}.}
\item{Q0}{True Q-matrix in the desired simulation attribute order.}
\item{G0}{Optional true attribute hierarchy in the desired simulation
attribute order. If supplied, it is used to break ties after Q matching.}
\item{cut_value}{Threshold for converting posterior mean \code{G} to a binary
hierarchy during tie-breaking and point-estimate recomputation.}
\item{K}{Number of attributes.}
\item{verbose}{Whether to print the selected permutation and discrepancies.}
\item{est}{Output from \code{\link{Est_fun}} created with
\code{return_samples = TRUE}.}
\item{keep_samples}{Whether to keep the truth-aligned posterior samples in
the returned estimate object.}
}
\value{
\code{align_samples_to_truth()} returns posterior samples with components
\code{Q}, \code{G}, \code{alpha}, and \code{pi}, plus the selected
\code{permutation} and Q/G discrepancy scores.

\code{align_estimates_to_truth()} returns an updated \code{Est_fun} result in
which \code{Est_Q}, \code{Est_G}, \code{Est_alpha}, \code{Est_pi}, and
\code{Est_GG} have been recomputed after truth-label alignment.
}
\description{
These functions are intended for simulation studies, where the true Q-matrix
and hierarchy are known. The matching starts from the posterior modal
Q-matrix and chooses the attribute permutation closest to \code{Q0}; if
several permutations tie and \code{G0} is supplied, the hierarchy is used as
a tie-breaker. The selected permutation is always applied simultaneously to
\code{Q}, \code{G}, \code{alpha}, and \code{pi}, because all four quantities
depend on the same latent attribute order.
}
\details{
For real data, use \code{\link{summary_AHMQ}} directly. Real data have no
known \code{Q0} or \code{G0}, so truth-label alignment is not identifiable.
}
