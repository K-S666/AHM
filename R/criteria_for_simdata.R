#' Attribute agreement rate (AAR)
#'
#' For each attribute \eqn{k}, computes the proportion of subjects on whom
#' the estimated attribute profile matches the true profile on that column.
#' Equivalently, \code{1 -} the per-attribute misclassification rate.
#'
#' @param alpha0 True attribute profile matrix (\eqn{N \times K}).
#' @param alpha Estimated attribute profile matrix (\eqn{N \times K}).
#' @return A numeric vector of length \code{K} with the agreement rate for each
#'   attribute.
#' @export
AAR <- function(alpha0, alpha)
{
  N <- nrow(alpha)
  # Per-attribute agreement: 1 minus the mean absolute error rate
  1 - colSums(abs(alpha - alpha0)) / N
}

#' Profile agreement rate (PAR)
#'
#' Proportion of subjects whose full \eqn{K}-dimensional attribute profile
#' matches the true profile exactly (all attributes correct simultaneously).
#'
#' @param alpha0 True attribute profile matrix (\eqn{N \times K}).
#' @param alpha Estimated attribute profile matrix (\eqn{N \times K}).
#' @return Scalar in \eqn{[0, 1]}: fraction of subjects with exact profile match.
#' @export
PAR <- function(alpha0, alpha)
{
  N <- nrow(alpha)
  # Count subjects with row-wise exact match
  num <- sum(apply(alpha == alpha0, 1, all))
  num / N
}

ahmq_offdiag <- function(K)
{
  row(diag(K)) != col(diag(K))
}

ahmq_binary_matrix <- function(x)
{
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[x != 0] <- 1
  x
}

ahmq_rate <- function(num, den)
{
  if (den == 0) {
    return(NA_real_)
  }
  num / den
}

#' Recovery metrics for an estimated attribute hierarchy
#'
#' Computes direct-edge recovery on the adjacency matrix and indirect/direct
#' reachability recovery on the reachability matrix. The diagonal is excluded
#' because self-reachability is fixed and not estimated.
#'
#' @param G0 True binary adjacency matrix.
#' @param G Estimated binary adjacency matrix.
#' @return A named list with \code{TPR}, \code{TFR}, \code{RTPR},
#'   \code{RTFR}, \code{G_exact}, and \code{R_exact}. Here \code{TFR} follows
#'   the notation used in the manuscript tables and is the true negative rate
#'   for non-existing edges.
#' @export
G_recovery_metrics <- function(G0, G)
{
  G0 <- ahmq_binary_matrix(G0)
  G <- ahmq_binary_matrix(G)
  if (!identical(dim(G0), dim(G))) {
    stop("G0 and G must have the same dimensions.", call. = FALSE)
  }

  K <- ncol(G0)
  idx <- ahmq_offdiag(K)
  true_edge <- G0[idx] == 1
  est_edge <- G[idx] == 1

  tp <- sum(est_edge & true_edge)
  fn <- sum(!est_edge & true_edge)
  tn <- sum(!est_edge & !true_edge)
  fp <- sum(est_edge & !true_edge)

  R0 <- Reachability(G0, K)
  R <- Reachability(G, K)
  true_reach <- R0[idx] == 1
  est_reach <- R[idx] == 1
  rtp <- sum(est_reach & true_reach)
  rfn <- sum(!est_reach & true_reach)
  rtn <- sum(!est_reach & !true_reach)
  rfp <- sum(est_reach & !true_reach)

  list(
    TPR = ahmq_rate(tp, tp + fn),
    TFR = ahmq_rate(tn, tn + fp),
    RTPR = ahmq_rate(rtp, rtp + rfn),
    RTFR = ahmq_rate(rtn, rtn + rfp),
    G_exact = isTRUE(all(G == G0)),
    R_exact = isTRUE(all(R == R0))
  )
}

ahmq_structured_Q <- function(Q, G)
{
  Q <- ahmq_binary_matrix(Q)
  G <- ahmq_binary_matrix(G)
  Restricted_Q(Q, Reachability(G, ncol(G)), nrow(Q), ncol(Q))$Q_restrict
}

#' Recovery metrics for an estimated Q-matrix
#'
#' Computes the Q-matrix recovery rate for each attribute and the average
#' Q-matrix recovery rate. By default, both the true and estimated Q-matrices
#' are first converted to hierarchy-structured Q-matrices: the true Q by the
#' true hierarchy \code{G0}, and the estimated Q by the estimated hierarchy
#' \code{G}. This matches the simulation evaluation logic where Q recovery is
#' assessed after deciding the hierarchy.
#'
#' @param Q0 True Q-matrix.
#' @param Q Estimated Q-matrix.
#' @param G0 True hierarchy used to structure \code{Q0}.
#' @param G Estimated hierarchy used to structure \code{Q}.
#' @param restrict If \code{TRUE}, apply hierarchy-based Q restriction before
#'   comparison.
#' @return A named list with per-attribute \code{QRR}, scalar \code{AQRR},
#'   logical \code{Q_exact}, and the compared matrices \code{Q0_compare} and
#'   \code{Q_compare}.
#' @export
Q_recovery_metrics <- function(Q0, Q, G0 = NULL, G = NULL, restrict = TRUE)
{
  Q0 <- ahmq_binary_matrix(Q0)
  Q <- ahmq_binary_matrix(Q)
  if (!identical(dim(Q0), dim(Q))) {
    stop("Q0 and Q must have the same dimensions.", call. = FALSE)
  }

  if (isTRUE(restrict)) {
    if (is.null(G0) || is.null(G)) {
      stop("G0 and G are required when restrict = TRUE.", call. = FALSE)
    }
    Q0 <- ahmq_structured_Q(Q0, G0)
    Q <- ahmq_structured_Q(Q, G)
  }

  J <- nrow(Q0)
  QRR <- 1 - colSums(abs(Q - Q0)) / J
  names(QRR) <- paste0("QRR", seq_along(QRR))
  list(
    QRR = QRR,
    AQRR = mean(QRR),
    Q_exact = isTRUE(all(Q == Q0)),
    Q0_compare = Q0,
    Q_compare = Q
  )
}

#' Bias summaries for item parameters in a single replication
#'
#' @param true True parameter vector.
#' @param estimate Estimated parameter vector.
#' @return A list with element-wise \code{bias}, \code{mean_bias}, and
#'   \code{mean_abs_bias}. For one replication, the absolute bias is the
#'   single-run analogue of RMSE.
#' @export
single_run_bias <- function(true, estimate)
{
  if (length(true) != length(estimate)) {
    stop("true and estimate must have the same length.", call. = FALSE)
  }
  bias <- as.numeric(estimate) - as.numeric(true)
  list(
    bias = bias,
    mean_bias = mean(bias),
    mean_abs_bias = mean(abs(bias))
  )
}

#' Summarize one simulation replication
#'
#' Computes simulation-only evaluation criteria after matching estimated
#' attribute labels to the known true labels. The matching step is performed
#' by \code{\link{align_estimates_to_truth}} and therefore applies the same
#' selected permutation simultaneously to \code{Q}, \code{G}, \code{alpha},
#' and \code{pi}. This avoids treating equivalent label permutations, such as
#' \code{1 -> 2 -> 3} versus \code{2 -> 3 -> 1}, as recovery errors.
#'
#' This function is intended for simulated data. For real data summaries use
#' \code{\link{summary_AHMQ}}, because no true \code{Q}, \code{G}, or
#' \code{alpha} is available for real data label alignment.
#'
#' @param est Output from \code{\link{Est_fun}}. If \code{align_labels = TRUE},
#'   \code{est} must be created with \code{return_samples = TRUE}.
#' @param truth Optional data list returned by \code{\link{simulate_ahmq_data}}.
#'   When supplied, \code{s0}, \code{g0}, \code{Q0}, \code{G0}, and
#'   \code{alpha0} are read from \code{truth} unless explicitly supplied.
#' @param s0,g0 True slip and guess vectors.
#' @param Q0 True Q-matrix.
#' @param G0 True attribute hierarchy.
#' @param alpha0 True examinee attribute profiles.
#' @param cut_value Threshold for binarizing posterior mean \code{G}; defaults
#'   to \code{est$cut_value}, then \code{0.2}.
#' @param align_labels If \code{TRUE}, align estimated labels to \code{Q0}
#'   before computing recovery criteria.
#' @param keep_aligned_samples If \code{TRUE}, keep posterior samples in the
#'   aligned estimate object returned as \code{aligned_estimates}; set
#'   \code{FALSE} to reduce object size.
#' @param verbose If \code{TRUE}, print the true-label matching permutation.
#' @return An object of class \code{"AHMQ_simulation_summary"} containing the
#'   truth-aligned estimates, item-parameter bias summaries, AAR/PAR, G/Q
#'   recovery metrics, Rhat table, and label-alignment metadata.
#' @export
simu_result_summary <- function(est,
                                truth = NULL,
                                s0 = NULL,
                                g0 = NULL,
                                Q0 = NULL,
                                G0 = NULL,
                                alpha0 = NULL,
                                cut_value = NULL,
                                align_labels = TRUE,
                                keep_aligned_samples = FALSE,
                                verbose = FALSE)
{
  known_Q <- isTRUE(est$known_Q)
  if (!is.null(truth)) {
    if (is.null(s0)) s0 <- truth$s
    if (is.null(g0)) g0 <- truth$g
    if (is.null(Q0)) Q0 <- truth$Q
    if (is.null(G0)) G0 <- truth$G
    if (is.null(alpha0)) alpha0 <- truth$alpha
  }
  if (known_Q) {
    if (is.null(s0) || is.null(g0) || is.null(G0) || is.null(alpha0)) {
      stop("s0, g0, G0, and alpha0 are required for fixed-Q simulation summaries.",
           call. = FALSE)
    }
  } else {
    if (is.null(s0) || is.null(g0) || is.null(Q0) ||
        is.null(G0) || is.null(alpha0)) {
      stop("s0, g0, Q0, G0, and alpha0 are required, either directly or via truth.",
           call. = FALSE)
    }
  }
  if (is.null(cut_value)) {
    cut_value <- est$cut_value
  }
  if (is.null(cut_value)) {
    cut_value <- 0.2
  }

  aligned_est <- est
  if (known_Q) {
    align_labels <- FALSE
  }
  if (isTRUE(align_labels)) {
    aligned_est <- align_estimates_to_truth(
      est = est,
      Q0 = Q0,
      G0 = G0,
      cut_value = cut_value,
      verbose = verbose,
      keep_samples = keep_aligned_samples
    )
  }

  s_bias_summary <- single_run_bias(s0, aligned_est$Est_s)
  g_bias_summary <- single_run_bias(g0, aligned_est$Est_g)
  s_g_metrics <- data.frame(
    parameter = c("s", "g"),
    mean_true = c(mean(s0), mean(g0)),
    mean_est = c(mean(aligned_est$Est_s), mean(aligned_est$Est_g)),
    mean_bias = c(s_bias_summary$mean_bias, g_bias_summary$mean_bias),
    mean_abs_bias = c(s_bias_summary$mean_abs_bias,
                      g_bias_summary$mean_abs_bias),
    row.names = NULL
  )

  alpha_aar <- AAR(alpha0, aligned_est$Est_alpha)
  alpha_metrics <- list(
    AAR = alpha_aar,
    mean_AAR = mean(alpha_aar),
    PAR = PAR(alpha0, aligned_est$Est_alpha)
  )

  G_metrics <- G_recovery_metrics(G0, aligned_est$Est_G)
  Q_metrics <- NULL
  if (!known_Q) {
    Q_metrics <- Q_recovery_metrics(Q0, aligned_est$Est_Q,
                                    G0 = G0,
                                    G = aligned_est$Est_G,
                                    restrict = TRUE)
  }

  rhat_table <- NULL
  if (!is.null(aligned_est$Rhat)) {
    rhat_table <- data.frame(
      parameter = names(aligned_est$Rhat$max),
      max_Rhat = as.numeric(aligned_est$Rhat$max),
      mpsrf = as.numeric(aligned_est$Rhat$mpsrf),
      row.names = NULL
    )
  }

  out <- list(
    aligned_estimates = aligned_est,
    s_bias = s_bias_summary,
    g_bias = g_bias_summary,
    s_g_metrics = s_g_metrics,
    alpha_metrics = alpha_metrics,
    G_metrics = G_metrics,
    Q_metrics = Q_metrics,
    rhat = rhat_table,
    cut_value = cut_value,
    known_Q = known_Q,
    align_labels = align_labels,
    truth_label_permutation = aligned_est$truth_label_permutation,
    truth_label_Q_distance = aligned_est$truth_label_Q_distance,
    truth_label_G_distance = aligned_est$truth_label_G_distance
  )
  class(out) <- "AHMQ_simulation_summary"
  out
}

#' @export
print.AHMQ_simulation_summary <- function(x, digits = 3, ...)
{
  cat("AHMQ simulation summary\n")
  if (isTRUE(x$known_Q)) {
    cat("Q-matrix is fixed and treated as known; Q recovery is not evaluated.\n")
  }
  if (isTRUE(x$align_labels) && !is.null(x$truth_label_permutation)) {
    cat("Truth-label permutation:",
        paste(x$truth_label_permutation, collapse = ","), "\n")
  }

  if (!is.null(x$aligned_estimates$Est_G)) {
    cat("\nTruth-aligned estimated hierarchy G:\n")
    print(x$aligned_estimates$Est_G)
  }

  cat("\nG recovery:\n")
  print(data.frame(
    TPR = round(x$G_metrics$TPR, digits),
    TFR = round(x$G_metrics$TFR, digits),
    RTPR = round(x$G_metrics$RTPR, digits),
    RTFR = round(x$G_metrics$RTFR, digits),
    G_exact = x$G_metrics$G_exact,
    R_exact = x$G_metrics$R_exact
  ), row.names = FALSE)

  if (!isTRUE(x$known_Q) && !is.null(x$aligned_estimates$Est_Q)) {
    cat("\nTruth-aligned structured Q-matrix:\n")
    print(x$aligned_estimates$Est_Q)
  }

  if (!isTRUE(x$known_Q)) {
    cat("\nQ recovery after hierarchy restriction:\n")
    print(data.frame(
      t(round(x$Q_metrics$QRR, digits)),
      AQRR = round(x$Q_metrics$AQRR, digits),
      Q_exact = x$Q_metrics$Q_exact,
      check.names = FALSE
    ), row.names = FALSE)
  }

  cat("\nBias for s_j and g_j:\n")
  tmp <- x$s_g_metrics[, c("parameter", "mean_bias", "mean_abs_bias"),
                       drop = FALSE]
  num_cols <- vapply(tmp, is.numeric, logical(1L))
  tmp[num_cols] <- lapply(tmp[num_cols], round, digits = digits)
  print(tmp, row.names = FALSE)
  cat("Note: with one replication, mean_abs_bias is the single-run analogue of RMSE.\n")

  cat("\nAttribute recovery:\n")
  cat("AAR:", paste(round(x$alpha_metrics$AAR, digits), collapse = ", "),
      "\n")
  cat("mean AAR:", round(x$alpha_metrics$mean_AAR, digits), "\n")
  cat("PAR:", round(x$alpha_metrics$PAR, digits), "\n")
  invisible(x)
}
