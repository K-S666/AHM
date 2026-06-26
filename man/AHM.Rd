\name{AHM}
\alias{AHM}
\title{Multi-chain Bayesian DINA estimation with a known Q-matrix}
\usage{
AHM(Y, Q, N1 = NULL, chain_length = 20000L, burn_in = 10000L,
  keep_burnin = FALSE, chain_num = 4L, a_s0 = 1, a_g0 = 1,
  b_s0 = 1, b_g0 = 1, p_add = 0.5, parallel = TRUE,
  n_cores = min(chain_num, parallel::detectCores(logical = FALSE)),
  progress = TRUE, print_every = 1000L)
}
\arguments{
\item{Y}{Binary response matrix with dimension \eqn{N \times J}.}
\item{Q}{Known binary Q-matrix with dimension \eqn{J \times K}.}
\item{N1}{Mini-batch sample size for updating item parameters and G. Q is
fixed in this model. If \code{NULL}, the default is \code{floor(N / 2)},
where \code{N} is the number of examinees.}
\item{chain_length}{Total number of MCMC iterations per chain.}
\item{burn_in}{Number of initial iterations treated as burn-in.}
\item{keep_burnin}{If \code{TRUE}, retain burn-in draws in
\code{chain_sample}; if \code{FALSE}, store only post-burn-in draws. This is
mainly useful for full-chain trace plots and Rhat trajectories.}
\item{chain_num}{Number of independent MCMC chains.}
\item{a_s0,b_s0,a_g0,b_g0}{Beta prior hyperparameters for slip and guess.}
\item{p_add}{Proposal probability for adding an edge in the hierarchy
\code{G} update; removing an edge is proposed with probability
\code{1 - p_add}.}
\item{parallel}{Whether to run chains in parallel.}
\item{n_cores}{Number of workers used when \code{parallel = TRUE}.}
\item{progress}{Whether to print per-chain progress messages.}
\item{print_every}{Print progress every this many MCMC iterations.}
}
\value{
A list of chain outputs. Each element contains \code{data}, including the
fixed \code{Q_fixed}, and \code{chain_sample}, including posterior draws for
\code{s}, \code{g}, \code{pi}, \code{alpha}, and \code{G}. Q is fixed and is
not sampled.
}
\description{
\code{AHM()} is the fixed-Q counterpart of \code{\link{AHMQ}}. It estimates
item parameters, individual attributes, latent class probabilities, and the
attribute hierarchy under the DINA model while treating the supplied Q-matrix
as known. Because Q anchors attribute labels, \code{\link{Est_fun}} does not
perform Q-based label switching for \code{AHM()} output and summary functions
do not report Q recovery.
}
\examples{
\dontrun{
dat <- simulate_ahmq_data(N = 100, J = 12, K = 3, seed = 1)
fit <- AHM(dat$Y, dat$Q, chain_length = 2000, burn_in = 1000,
           chain_num = 2)
est <- Est_fun(fit)
summary_AHMQ(est)
}
}
