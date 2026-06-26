\name{compute_DIC}
\alias{compute_DIC}
\alias{compute_DIC_fixedQ}
\title{Compute DIC for DINA posterior draws}
\usage{
compute_DIC(y, N, J, K, chain, s, g, alpha, Q)

compute_DIC_fixedQ(y, N, J, K, chain, s, g, alpha, Q)
}
\arguments{
\item{y}{Binary response matrix.}
\item{N}{Number of examinees.}
\item{J}{Number of items.}
\item{K}{Number of attributes.}
\item{chain}{Number of posterior draws.}
\item{s,g}{Matrices of slip and guess draws with dimension
\eqn{J \times chain}.}
\item{alpha}{Array of attribute profile draws with dimension
\eqn{N \times K \times chain}.}
\item{Q}{For \code{compute_DIC()}, an array of Q-matrix draws with dimension
\eqn{J \times K \times chain}. For \code{compute_DIC_fixedQ()}, the known
\eqn{J \times K} Q-matrix.}
}
\value{
Numeric DIC value for the supplied posterior draws.
}
\description{
Computes a DINA deviance information criterion from posterior MCMC draws.
\code{compute_DIC()} is used when Q is sampled by \code{\link{AHMQ}}.
\code{compute_DIC_fixedQ()} is used when Q is known in \code{\link{AHM}}.
Both are used internally by \code{\link{select_chain_by_DIC}} to choose the
chain with the smallest DIC before final summaries.
}
