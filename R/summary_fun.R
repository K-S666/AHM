#' Summarize posterior estimates
#'
#' Creates a compact, human-readable summary of posterior point estimates,
#' convergence diagnostics, and model fit measures returned by
#' \code{\link{Est_fun}}. The function is intentionally lightweight: it does
#' not recompute MCMC quantities, but organizes the already computed DIC,
#' best-chain index, Rhat values, and estimated parameters.
#'
#' @param object Output from \code{\link{Est_fun}}.
#' @param digits Number of digits used in printed numeric summaries.
#' @param rhat_cutoff Threshold used to flag convergence in the summary table.
#' @return An object of class \code{"est_summary"} with model-fit,
#'   convergence, and parameter-estimate components. The object is printed
#'   automatically when returned at the console.
#' @examples
#' dat <- simulate_ahmq_data(N = 30, J = 8, K = 2, seed = 1)
#' \dontrun{
#' fit <- AHMQ(dat$Y, K = dat$K, chain_length = 1000, burn_in = 500,
#'             chain_num = 2)
#' est <- Est_fun(fit)
#' summary_est(est)
#' }
#' @export
summary_est <- function(object, digits = 3, rhat_cutoff = 1.1)
{
  required <- c("Est_s", "Est_g", "Est_pi", "Est_G",
                "Est_alpha", "DIC", "DIC_all", "best_chain")
  missing <- setdiff(required, names(object))
  if (length(missing) > 0L) {
    stop("object is missing required Est_fun components: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  rhat <- object$Rhat
  convergence <- NULL
  if (!is.null(rhat)) {
    convergence <- data.frame(
      parameter = names(rhat$max),
      max_Rhat = as.numeric(rhat$max),
      mpsrf = as.numeric(rhat$mpsrf),
      converged = as.numeric(rhat$max) <= rhat_cutoff,
      row.names = NULL
    )
  }

  out <- list(
    fit = list(
      DIC = object$DIC,
      DIC_all = object$DIC_all,
      best_chain = object$best_chain,
      runtime = object$runtime
    ),
    convergence = convergence,
    estimates = list(
      slip = object$Est_s,
      guess = object$Est_g,
      pi = object$Est_pi,
      slip_sd = object$Est_s_sd,
      guess_sd = object$Est_g_sd,
      pi_sd = object$Est_pi_sd,
      G_posterior = object$Est_GG,
      G = object$Est_G,
      Q = if (isTRUE(object$known_Q)) NULL else object$Est_Q,
      alpha = object$Est_alpha
    ),
    dimensions = list(
      J = length(object$Est_s),
      K = ncol(object$Est_alpha),
      N = nrow(object$Est_alpha)
    ),
    known_Q = isTRUE(object$known_Q),
    model = object$model %||% if (isTRUE(object$known_Q)) "AHM" else "AHMQ",
    digits = digits,
    rhat_cutoff = rhat_cutoff,
    cut_value = object$cut_value %||% 0.2
  )
  class(out) <- "est_summary"
  out
}

#' @export
print.est_summary <- function(x, ...)
{
  d <- x$digits
  cat(x$model, " posterior summary\n", sep = "")
  if (isTRUE(x$known_Q)) {
    cat("Q-matrix is fixed and treated as known.\n")
  }
  cat("Items:", x$dimensions$J,
      " Attributes:", x$dimensions$K,
      " Examinees:", x$dimensions$N, "\n")
  cat("Best chain by DIC:", x$fit$best_chain,
      " DIC:", round(x$fit$DIC, d), "\n")
  if (!is.null(x$fit$runtime)) {
    cat("Runtime:", round(as.numeric(x$fit$runtime), d), "seconds\n")
  }

  cat("\nDIC by chain:\n")
  print(round(x$fit$DIC_all, d))

  if (!is.null(x$convergence)) {
    cat("\nRhat diagnostics:\n")
    conv <- x$convergence
    conv$max_Rhat <- round(conv$max_Rhat, d)
    conv$mpsrf <- round(conv$mpsrf, d)
    print(conv, row.names = FALSE)
  }

  cat("\nSlip estimates:\n")
  if (!is.null(x$estimates$slip_sd)) {
    print(round(data.frame(mean = x$estimates$slip,
                           SD = x$estimates$slip_sd), d))
  } else {
    print(round(x$estimates$slip, d))
  }
  cat("\nGuess estimates:\n")
  if (!is.null(x$estimates$guess_sd)) {
    print(round(data.frame(mean = x$estimates$guess,
                           SD = x$estimates$guess_sd), d))
  } else {
    print(round(x$estimates$guess, d))
  }
  if (!is.null(x$estimates$pi_sd)) {
    cat("\nClass probability estimates pi:\n")
    print(round(data.frame(mean = x$estimates$pi,
                           SD = x$estimates$pi_sd), d))
  }

  cat("\nPosterior mean of G:\n")
  print(round(x$estimates$G_posterior, d))
  cat("\nG cut_value:", x$cut_value, "\n")
  cat("\nFinal estimated hierarchy G:\n")
  print(x$estimates$G)
  if (isTRUE(x$known_Q)) {
    cat("\nQ-matrix was supplied as known and is not reported as an estimate.\n")
  } else {
    cat("\nFinal estimated structured Q-matrix:\n")
    print(x$estimates$Q)
  }
  invisible(x)
}

`%||%` <- function(x, y)
{
  if (is.null(x)) y else x
}



