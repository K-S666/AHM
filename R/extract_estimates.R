#' Extract AHMQ point estimates
#'
#' Provides a simple, stable interface for downstream analysis. The input can
#' be either raw MCMC output from \code{\link{AHMQ}} or a post-processed object
#' returned by \code{\link{Est_fun}}. Raw MCMC output is first passed through
#' \code{Est_fun}; the returned list then uses lower-case names such as
#' \code{est_s}, \code{est_g}, \code{est_alpha}, \code{est_Q}, and
#' \code{est_G}. For fixed-Q \code{\link{AHM}} output, \code{est_Q} is
#' \code{NULL} and the supplied Q-matrix is returned as \code{known_Q_matrix}.
#'
#' @param object Output from \code{\link{AHMQ}} or \code{\link{Est_fun}}.
#' @param ... Additional arguments passed to \code{\link{Est_fun}} when
#'   \code{object} is raw MCMC output. Useful arguments include
#'   \code{burn_in}, \code{chain_length}, \code{cut_value},
#'   \code{compute_rhat}, and \code{return_samples}.
#' @return A list with item parameter estimates, examinee attribute estimates,
#'   Q and G estimates, posterior mean of G, model fit information, and the
#'   original \code{Est_fun} object under \code{est_object}.
#' @examples
#' dat <- simulate_ahmq_data(N = 30, J = 8, K = 2, seed = 1)
#' \dontrun{
#' fit <- AHMQ(dat$Y, K = dat$K, chain_length = 1000, burn_in = 500)
#' est <- extract_estimates(fit)
#' est$est_s
#' est$est_alpha
#' est$est_Q
#' est$est_G
#' }
#' @export
extract_estimates <- function(object, ...)
{
  is_est_fun <- all(c("Est_s", "Est_g", "Est_alpha", "Est_G") %in%
                      names(object))
  est <- if (is_est_fun) object else Est_fun(object, ...)

  out <- list(
    est_s = est$Est_s,
    est_g = est$Est_g,
    est_alpha = est$Est_alpha,
    est_Q = est$Est_Q,
    known_Q = isTRUE(est$known_Q),
    known_Q_matrix = est$Q_fixed,
    est_G = est$Est_G,
    est_G_posterior = est$Est_GG,
    est_pi = est$Est_pi,
    cut_value = est$cut_value,
    DIC = est$DIC,
    DIC_all = est$DIC_all,
    best_chain = est$best_chain,
    Rhat = est$Rhat,
    runtime = est$runtime,
    est_object = est
  )
  class(out) <- "AHMQ_estimates"
  out
}

#' @export
print.AHMQ_estimates <- function(x, digits = 3, ...)
{
  cat("AHMQ estimates\n")
  if (isTRUE(x$known_Q)) {
    cat("Model: AHM with known Q\n")
  }
  K <- if (!is.null(x$est_Q)) ncol(x$est_Q) else ncol(x$known_Q_matrix)
  cat("Items:", length(x$est_s),
      " Attributes:", K,
      " Examinees:", nrow(x$est_alpha), "\n")
  cat("Best chain:", x$best_chain,
      " DIC:", round(x$DIC, digits), "\n")
  cat("\nFirst item parameter estimates:\n")
  item <- data.frame(
    item = seq_along(x$est_s),
    s = round(x$est_s, digits),
    g = round(x$est_g, digits)
  )
  print(utils::head(item, 10L), row.names = FALSE)
  cat("\nEstimated G:\n")
  print(x$est_G)
  if (isTRUE(x$known_Q)) {
    cat("\nKnown Q-matrix:\n")
    print(x$known_Q_matrix)
  } else {
    cat("\nEstimated structured Q:\n")
    print(x$est_Q)
  }
  invisible(x)
}
