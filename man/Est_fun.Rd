\name{Est_fun}
\alias{Est_fun}
\alias{select_chain_by_DIC}
\alias{compute_Rhat}
\title{Post-process AHM/AHMQ MCMC output}
\usage{
Est_fun(result, burn_in = NULL, chain_length = NULL, cut_value = 0.2,
  compute_rhat = TRUE, rhat_parameters = c("s", "g"),
  return_samples = FALSE, verbose = TRUE)

select_chain_by_DIC(result, post_idx, sample_size)

compute_Rhat(result, post_idx = NULL, parameters = c("s", "g"))
}
\arguments{
\item{result}{Output from \code{\link{AHMQ}} or fixed-Q \code{\link{AHM}}.}
\item{burn_in}{Optional number of stored draws to discard before summaries.
If the \code{AHMQ} result was created with \code{keep_burnin = TRUE},
\code{NULL} automatically uses the burn-in value recorded in the result.}
\item{chain_length}{Optional number of stored draws to use.}
\item{cut_value}{Threshold for binarizing posterior mean hierarchy edges. The
default is 0.2, matching the simulation-study cutoff used for final recovery
evaluation.}
\item{compute_rhat}{Whether to compute Gelman-Rubin diagnostics before chain selection.}
\item{rhat_parameters}{Chain sample components passed to \code{\link{compute_Rhat}}.
Only \code{"s"} and \code{"g"} are used.}
\item{return_samples}{Whether to return posterior samples after any required
label-switch processing. For fixed-Q \code{\link{AHM}} output, samples are
returned without relabeling.}
\item{verbose}{Whether to print label-switching progress for unknown-Q
\code{\link{AHMQ}} output.}
\item{post_idx}{Posterior draw indices used by \code{select_chain_by_DIC}.}
\item{sample_size}{Number of posterior draws used for DIC.}
\item{parameters}{Chain sample components for Rhat diagnostics. Only
\code{"s"} and \code{"g"} are used. If
\code{post_idx = NULL}, \code{compute_Rhat()} discards the recorded burn-in
when full chains were saved with \code{keep_burnin = TRUE}; otherwise it uses
all stored draws.}
}
\value{
\code{Est_fun()} returns point estimates for \code{s}, \code{g}, \code{pi},
\code{G}, and \code{alpha}, plus the DIC value, all chain DIC values, the
selected chain index, and optional Rhat diagnostics. For unknown-Q
\code{\link{AHMQ}} output it also returns estimated \code{Q}; for fixed-Q
\code{\link{AHM}} output it returns the supplied matrix as \code{Q_fixed}.
\code{select_chain_by_DIC()} returns the minimum-DIC chain and raw posterior
draws for downstream relabeling. \code{compute_Rhat()} returns univariate Rhat
vectors, maximum univariate Rhat values, and multivariate PSRF values for
\code{s} and \code{g}.
}
\description{
Selects the best chain by DIC and computes posterior point estimates. For
unknown-Q \code{\link{AHMQ}} output it applies Q-matrix based label-switch
correction and synchronizes the same relabeling to \code{G}, \code{alpha}, and
\code{pi}. For fixed-Q \code{\link{AHM}} output, Q anchors attribute labels, so
label-switch correction and Q estimation are skipped.
}
