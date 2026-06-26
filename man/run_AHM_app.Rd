\name{run_AHM_app}
\alias{run_AHM_app}
\title{Launch the AHM Shiny application}
\usage{
run_AHM_app(...)
}
\arguments{
\item{...}{Additional arguments passed to \code{shiny::runApp}, such as
\code{port}, \code{host}, or \code{launch.browser}.}
}
\value{
The return value of \code{shiny::runApp}.
}
\description{
Opens the interactive AHM analysis workbench bundled with the package. The app
supports simulated data, the bundled ECPE real-data example, and user-uploaded
response/Q-matrix CSV files.
}
\examples{
\dontrun{
library(AHM)
run_AHM_app()
}
}
