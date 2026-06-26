\name{ecpe_ahmq_data}
\alias{ecpe_ahmq_data}
\title{Prepare bundled ECPE data for AHMQ}
\usage{
ecpe_ahmq_data()
}
\value{
A list with:
\describe{
\item{\code{Y}}{Binary item response matrix.}
\item{\code{Q}}{Original ECPE Q-matrix.}
\item{\code{K}}{Number of attributes.}
}
}
\description{
Converts the bundled \code{\link{data.ecpe}} object into the matrix inputs
expected by \code{\link{AHMQ}}. The response matrix drops the first column of
\code{data.ecpe$data}, matching the original ECPE analysis script.
}
\examples{
ecpe <- ecpe_ahmq_data()
dim(ecpe$Y)
dim(ecpe$Q)
ecpe$K
}
