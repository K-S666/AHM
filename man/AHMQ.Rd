\name{AHMQ}
\alias{AHMQ}
\title{Multi-chain Bayesian DINA estimation with unknown Q and hierarchy}
\usage{
AHMQ(Y, K, N1 = NULL, chain_length = 20000L, burn_in = 10000L,
  keep_burnin = FALSE, chain_num = 4L, a_s0 = 1, a_g0 = 1,
  b_s0 = 1, b_g0 = 1, p_add = 0.5,
  parallel = TRUE, n_cores = min(chain_num, parallel::detectCores(logical = FALSE)),
  progress = TRUE, print_every = 1000L)
}
\arguments{
\item{Y}{An \eqn{N \times J} binary response matrix.}
\item{K}{Number of latent attributes.}
\item{N1}{Mini-batch sample size for updating item parameters, Q, and G.
If \code{NULL}, the default is \code{floor(N / 2)}, where \code{N} is the
number of examinees.}
\item{chain_length}{Total number of MCMC iterations per chain.}
\item{burn_in}{Number of initial iterations treated as burn-in.}
\item{keep_burnin}{If \code{TRUE}, retain burn-in draws in
\code{chain_sample}; if \code{FALSE}, store only post-burn-in draws. Use
\code{TRUE} when full-chain diagnostics over time are needed.}
\item{chain_num}{Number of independent chains.}
\item{a_s0, b_s0, a_g0, b_g0}{Beta prior hyperparameters for slipping and guessing.}
\item{p_add}{Proposal probability for adding an edge in the hierarchy
\code{G} update; removing an edge is proposed with probability
\code{1 - p_add}.}
\item{parallel}{Whether to run chains in parallel.}
\item{n_cores}{Number of workers used when \code{parallel = TRUE}.}
\item{progress}{Whether to print per-chain progress messages.}
\item{print_every}{Print progress every this many MCMC iterations.}
}
\value{
A list of chains. Each element contains the input settings in \code{data} and
posterior draws in \code{chain_sample}: \code{G}, \code{s}, \code{g},
\code{Q}, \code{pi}, and \code{alpha}.
}
\description{
Runs independent MCMC chains for the DINA model when the Q-matrix and the
attribute hierarchy are unknown. The sampler estimates item slipping and
guessing parameters, latent class probabilities, examinee attribute profiles,
the Q-matrix, and the transitive-reduced attribute hierarchy.
}
\examples{
\dontrun{
ecpe <- ecpe_ahmq_data()
fit <- AHMQ(ecpe$Y, K = ecpe$K, N1 = 128,
            chain_length = 2000, burn_in = 1000, chain_num = 4)
est <- Est_fun(fit, cut_value = 0.2)
}
}
