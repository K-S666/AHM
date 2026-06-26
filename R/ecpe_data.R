#' ECPE response data from the CDM package
#'
#' The \code{data.ecpe} object is copied from the \pkg{CDM} package and is
#' included here so the AHM/DINA workflow has a real-data example without
#' requiring users to install \pkg{CDM}. The object contains the original ECPE
#' response data and the CDM Q-matrix.
#'
#' @format A list with at least two components:
#' \describe{
#'   \item{\code{data}}{A data frame whose first column is an examinee/id field
#'     and whose remaining columns are binary item responses.}
#'   \item{\code{q.matrix}}{The original CDM Q-matrix for the ECPE items.}
#' }
#' @source The \pkg{CDM} package data object \code{data.ecpe}.
"data.ecpe"

#' Prepare ECPE data for AHMQ examples
#'
#' Converts the bundled \code{\link{data.ecpe}} object into the matrix inputs
#' expected by \code{\link{AHMQ}}. The response matrix drops the first column of
#' \code{data.ecpe$data}, matching the usage in the original analysis script.
#'
#' @return A list with \code{Y}, \code{Q}, and \code{K}.
#' @examples
#' ecpe <- ecpe_ahmq_data()
#' dim(ecpe$Y)
#' ecpe$K
#'
#' \dontrun{
#' fit <- AHMQ(ecpe$Y, K = ecpe$K,
#'             N1 = 128, chain_length = 2000, burn_in = 1000,
#'             chain_num = 4)
#' est <- Est_fun(fit, cut_value = 0.2)
#' est$Est_Q
#' est$Est_G
#' }
#' @export
ecpe_ahmq_data <- function()
{
  if (!exists("data.ecpe", inherits = TRUE)) {
    utils::data("data.ecpe", package = "AHM", envir = environment())
  }
  ecpe <- get("data.ecpe", envir = environment(), inherits = TRUE)

  Y <- as.matrix(ecpe$data[, -1L])
  storage.mode(Y) <- "double"
  Q <- as.matrix(ecpe$q.matrix)
  storage.mode(Q) <- "double"

  list(Y = Y, Q = Q, K = ncol(Q))
}
