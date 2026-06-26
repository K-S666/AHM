\name{data.ecpe}
\alias{data.ecpe}
\docType{data}
\title{ECPE response data from the CDM package}
\format{
A list containing at least:
\describe{
\item{\code{data}}{A data frame whose first column is an examinee/id field and
remaining columns are binary item responses.}
\item{\code{q.matrix}}{The original Q-matrix for the ECPE items.}
}
}
\description{
The ECPE data object copied from the \pkg{CDM} package. It is bundled so the
AHMQ workflow has a real-data example without requiring \pkg{CDM} at runtime.
}
\source{
The \pkg{CDM} package data object \code{data.ecpe}.
}
\keyword{datasets}
