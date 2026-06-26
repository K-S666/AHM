\name{AAR}
\alias{AAR}
\alias{PAR}
\alias{G_recovery_metrics}
\alias{Q_recovery_metrics}
\alias{single_run_bias}
\title{Agreement rates for simulation studies}
\usage{
AAR(alpha0, alpha)

PAR(alpha0, alpha)

G_recovery_metrics(G0, G)

Q_recovery_metrics(Q0, Q, G0 = NULL, G = NULL, restrict = TRUE)

single_run_bias(true, estimate)
}
\arguments{
\item{alpha0}{True binary attribute profile matrix.}
\item{alpha}{Estimated binary attribute profile matrix.}
\item{G0}{True binary hierarchy adjacency matrix.}
\item{G}{Estimated binary hierarchy adjacency matrix.}
\item{Q0}{True Q-matrix.}
\item{Q}{Estimated Q-matrix.}
\item{restrict}{Whether to structure Q-matrices by their corresponding
hierarchies before comparison.}
\item{true}{True numeric parameter vector.}
\item{estimate}{Estimated numeric parameter vector.}
}
\value{
\code{AAR()} returns per-attribute agreement rates. \code{PAR()} returns the
proportion of examinees whose full profile is exactly correct.
\code{G_recovery_metrics()} returns TPR, TFR, RTPR, RTFR, and exact recovery
indicators for G and the reachability matrix. \code{Q_recovery_metrics()}
returns QRR by attribute, AQRR, and exact Q recovery. \code{single_run_bias()}
returns bias summaries for one replication.
}
\description{
Computes common simulation-study accuracy metrics for estimated attribute
profiles, Q-matrices, hierarchy matrices, and item parameters.
}
