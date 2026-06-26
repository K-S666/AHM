\name{extract_estimates}
\alias{extract_estimates}
\alias{print.AHMQ_estimates}
\title{Extract AHM/AHMQ point estimates}
\usage{
extract_estimates(object, ...)
}
\arguments{
\item{object}{Output from \code{\link{AHMQ}}, \code{\link{AHM}}, or
\code{\link{Est_fun}}.}
\item{...}{Additional arguments passed to \code{\link{Est_fun}} when
\code{object} is raw MCMC output.}
}
\value{
A list with lower-case estimate fields: \code{est_s}, \code{est_g},
\code{est_alpha}, \code{est_Q}, \code{known_Q}, \code{known_Q_matrix},
\code{est_G}, \code{est_G_posterior}, \code{est_pi}, plus \code{DIC},
\code{DIC_all}, \code{best_chain}, \code{Rhat}, and the original
\code{Est_fun} object. For fixed-Q \code{\link{AHM}} output,
\code{est_Q = NULL} and \code{known_Q_matrix} stores the supplied Q-matrix.
}
\description{
Creates a convenient estimate object for downstream analysis. This is useful
when code should explicitly access estimates such as \code{est$est_s},
\code{est$est_alpha}, \code{est$est_Q}, and \code{est$est_G} rather than
working with the longer internal names returned by \code{Est_fun}.
}
\examples{
dat <- simulate_ahmq_data(N = 30, J = 8, K = 2, seed = 1)
\dontrun{
fit <- AHMQ(dat$Y, K = dat$K, chain_length = 1000, burn_in = 500)
est <- extract_estimates(fit)
est$est_s
est$est_alpha
est$est_Q
est$est_G
}
}
