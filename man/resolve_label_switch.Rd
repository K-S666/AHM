\name{resolve_label_switch}
\alias{resolve_label_switch}
\title{Resolve attribute label switching in posterior samples}
\usage{
resolve_label_switch(Q_sample, G_sample, alpha_sample, pi_sample, K,
  verbose = TRUE)
}
\arguments{
\item{Q_sample}{Item-by-attribute Q-matrix draws, an array of dimension
\eqn{J \times K \times T}.}
\item{G_sample}{Attribute hierarchy draws, an array of dimension
\eqn{K \times K \times T}.}
\item{alpha_sample}{Attribute profile draws, an array of dimension
\eqn{N \times K \times T}.}
\item{pi_sample}{Latent class probability draws, a matrix of dimension
\eqn{2^K \times T}.}
\item{K}{Number of attributes.}
\item{verbose}{Whether to print relabeling progress.}
}
\value{
A list with relabeled posterior draws: \code{Q}, \code{G}, \code{alpha}, and
\code{pi}.
}
\description{
Uses the posterior mean Q-matrix as a reference and exhaustively searches
attribute permutations to align every posterior draw. The selected permutation
is applied consistently to Q columns, G rows and columns, alpha columns, and
latent class probabilities.
}
