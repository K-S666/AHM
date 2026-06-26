#' Launch the AHM Shiny application
#'
#' Opens the interactive analysis workbench bundled with the package. The app
#' supports simulated data, the bundled ECPE real-data example, and user
#' uploaded response/Q-matrix CSV files. The Shiny dependency is optional for
#' the core package and is checked only when the app is launched.
#'
#' @param ... Additional arguments passed to \code{shiny::runApp}, such as
#'   \code{port}, \code{host}, or \code{launch.browser}.
#' @return The return value of \code{shiny::runApp}.
#' @export
run_AHM_app <- function(...)
{
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The Shiny app requires the 'shiny' package. Install it with install.packages('shiny').",
         call. = FALSE)
  }
  app_dir <- system.file("shiny", package = "AHM")
  if (!nzchar(app_dir)) {
    stop("Cannot find the bundled Shiny app directory.", call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
