\name{simulate_ahmq_data}
\alias{simulate_ahmq_data}
\alias{simu_alpha}
\alias{simu_Q}
\alias{DINA_data}
\alias{Restricted_Q}
\alias{samevec}
\title{Simulation utilities for AHMQ/DINA models}
\usage{
simulate_ahmq_data(N, J, K, Q = NULL, G = NULL, s = NULL, g = NULL,
  pi = NULL, alpha = NULL, alpha_distribution = c("balanced", "unbalanced"),
  sigma = 0.5, s_range = c(0.05, 0.25), g_range = c(0.05, 0.25),
  seed = NULL)

simu_alpha(N, G, distribution = c("balanced", "unbalanced"),
  sigma = 0.5, seed = NULL)

simu_Q(J, G, seed = NULL, max_tries = 1000L)

DINA_data(N, J, K, s, g, Q, alpha)

Restricted_Q(Q_all, R, J, K)

samevec(Q, q)
}
\arguments{
\item{N}{Number of examinees.}
\item{J}{Number of items.}
\item{K}{Number of attributes.}
\item{Q}{Optional true Q-matrix. If omitted, a random Q-matrix is generated.}
\item{G}{Optional true attribute hierarchy adjacency matrix. If omitted, the
no-prerequisite hierarchy is used.}
\item{s,g}{Slip and guess parameters. Scalars are recycled to item length.}
\item{pi}{Optional latent class probabilities over all \eqn{2^K} profiles.}
\item{alpha}{Optional true attribute profile matrix.}
\item{alpha_distribution}{Attribute profile distribution used when
\code{alpha} is omitted: balanced or unbalanced.}
\item{distribution}{Distribution used by \code{simu_alpha}.}
\item{sigma}{Common off-diagonal covariance/correlation for unbalanced
attribute generation.}
\item{s_range,g_range}{Ranges used to sample slip and guess when omitted.}
\item{seed}{Optional random seed.}
\item{max_tries}{Maximum attempts for identifiable Q-matrix generation.}
\item{Q_all}{Candidate Q-matrix to close under the hierarchy.}
\item{R}{Reachability matrix of the hierarchy.}
\item{q}{One Q-vector used by \code{samevec}.}
}
\value{
\code{simulate_ahmq_data()} returns a list with response data \code{Y}, true
\code{Q}, \code{G}, reachability matrix \code{R}, true \code{alpha}, class
probabilities \code{pi}, item parameters \code{s} and \code{g}, and dimensions.
\code{simu_alpha()} returns hierarchy-consistent attribute profiles.
\code{simu_Q()} returns a structured identifiable Q-matrix under \code{G}.
\code{DINA_data()} returns an \eqn{N \times J} response matrix.
\code{Restricted_Q()} returns a hierarchy-consistent Q-matrix and its decimal
row codes. \code{samevec()} returns a row count.
}
\description{
These helpers generate synthetic DINA data for package examples and simulation
studies. The high-level generator respects an attribute hierarchy by generating
hierarchy-consistent attribute profiles and an identifiable structured Q-matrix.
The Q generator enforces the manuscript conditions: identity-equivalent
\code{Q0}, at least three ones per column in the sparsified Q-matrix, and
\code{K} distinct columns in the densified \code{Q_star}.
}
\examples{
dat <- simulate_ahmq_data(N = 30, J = 8, K = 2, seed = 1)
dim(dat$Y)
dat$Q
}
