\name{summary_AHMQ}
\alias{summary_AHMQ}
\alias{print.AHMQ_summary}
\title{Summarize AHM/AHMQ posterior estimates}
\usage{
summary_AHMQ(object, digits = 3, rhat_cutoff = 1.1)
}
\arguments{
\item{object}{Output from \code{\link{Est_fun}}.}
\item{digits}{Number of digits used for printed numeric summaries.}
\item{rhat_cutoff}{Threshold used to flag Rhat convergence.}
}
\value{
An object of class \code{"AHMQ_summary"} with fit information, convergence
diagnostics, parameter estimates, and model dimensions.
}
\description{
Builds a compact report from posterior point estimates, DIC values, selected
chain index, Rhat diagnostics, the posterior mean of \code{G}, the cut-off
value used to finalize \code{G}, and the final estimated \code{G}. For
unknown-Q \code{\link{AHMQ}} output it also prints the final structured
\code{Q}-matrix; for fixed-Q \code{\link{AHM}} output it prints the known
Q-matrix and skips Q-estimation language.
}
\examples{
dat <- simulate_ahmq_data(N = 30, J = 8, K = 2, seed = 1)
\dontrun{
fit <- AHMQ(dat$Y, K = dat$K, chain_length = 1000, burn_in = 500,
            chain_num = 2)
est <- Est_fun(fit)
summary_AHMQ(est)
}
}
