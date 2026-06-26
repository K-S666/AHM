ahmq_permutations <- function(x)
{
  if (length(x) == 1L) {
    return(list(x))
  }
  out <- vector("list", 0L)
  for (i in seq_along(x)) {
    rest <- x[-i]
    tail_perm <- ahmq_permutations(rest)
    out <- c(out, lapply(tail_perm, function(p) c(x[i], p)))
  }
  out
}

ahmq_permute_draws <- function(Q_sample,
                               G_sample,
                               alpha_sample,
                               pi_sample,
                               perm,
                               K = length(perm))
{
  J <- dim(Q_sample)[1L]
  N <- dim(alpha_sample)[1L]
  L <- 2^K
  sample_size <- dim(Q_sample)[3L]
  alpha_all <- Trans_10to2_mat(K, c(0:(L - 1L)))

  Q_re <- array(0, c(J, K, sample_size))
  G_re <- array(0, c(K, K, sample_size))
  alpha_re <- array(0, c(N, K, sample_size))
  pi_re <- matrix(0, L, sample_size)

  alpha_all_re <- alpha_all
  for (k in seq_len(K)) {
    alpha_all_re[, k] <- alpha_all[, perm[k]]
  }
  new_order <- Trans_2to10_mat(alpha_all_re, K)
  pi_order <- integer(L)
  for (l in seq_len(L)) {
    pi_order[l] <- which(new_order == (l - 1L))[1L]
  }

  for (i in seq_len(sample_size)) {
    Q_re[, , i] <- Q_sample[, perm, i, drop = FALSE][, , 1L]
    alpha_re[, , i] <- alpha_sample[, perm, i, drop = FALSE][, , 1L]

    G_tmp <- G_sample[, perm, i, drop = FALSE][, , 1L]
    G_re[, , i] <- G_tmp[perm, , drop = FALSE]

    pi_re[, i] <- pi_sample[pi_order, i]
  }

  list(Q = Q_re, G = G_re, alpha = alpha_re, pi = pi_re)
}

ahmq_permute_matrix <- function(x, perm)
{
  x <- as.matrix(x)
  x[perm, perm, drop = FALSE]
}

#' Resolve attribute label switching in MCMC samples
#'
#' Iteratively relabels Q-matrix columns (and consistently relabels \code{G},
#' \code{alpha}, and \code{pi}) across posterior draws so that attribute
#' indices are aligned with the posterior mean Q-matrix. The relabeling rule
#' follows the Q-matrix mean reference/Hungarian-style exhaustive assignment
#' used in the simulation scripts.
#'
#' @param Q_sample Item-by-attribute Q draws (\eqn{J \times K \times T}).
#' @param G_sample Attribute hierarchy draws (\eqn{K \times K \times T}).
#' @param alpha_sample Attribute profile draws (\eqn{N \times K \times T}).
#' @param pi_sample Class probability draws (\eqn{L \times T}, \code{L = 2^K}).
#' @param K Number of attributes.
#' @param verbose If \code{TRUE}, print iteration progress (default \code{TRUE}).
#' @return A list with aligned draws: \code{Q}, \code{G}, \code{alpha},
#'   \code{pi}.
resolve_label_switch <- function(Q_sample,
                                 G_sample,
                                 alpha_sample,
                                 pi_sample,
                                 K,
                                 verbose = TRUE)
{
  J <- dim(Q_sample)[1]
  N <- dim(alpha_sample)[1]
  L <- 2^K
  sample_size <- dim(Q_sample)[3]
  permutation <- ahmq_permutations(seq_len(K))
  permutation_length <- length(permutation)
  alpha_all <- Trans_10to2_mat(K, c(0:(L - 1)))

  Q_re_est <- Q_sample
  alpha_re_est <- alpha_sample
  G_re_est <- G_sample
  pi_re_est <- pi_sample
  Q_est <- array(0, c(J, K, permutation_length))
  Q_current <- Q_sample
  Q_pre <- Q_sample
  G_current <- G_sample
  alpha_current <- alpha_sample
  d <- c()
  e <- 1
  iter <- 0

  while (e > 0) {
    iter <- iter + 1
    QQ <- apply(Q_current, c(1, 2), mean)
    Q_pre <- Q_current
    for (i in seq_len(sample_size)) {
      for (j in seq_len(permutation_length)) {
        ind <- permutation[[j]]
        for (k in seq_len(K)) {
          Q_est[, k, j] <- Q_current[, ind[k], i]
        }
        d[j] <- sqrt(sum((Q_est[, , j] - QQ)^2))
      }
      min_ind <- which.min(d)
      Q_re_est[, , i] <- Q_est[, , min_ind]
      Q_current[, , i] <- Q_re_est[, , i]
      ind <- permutation[[min_ind]]
      if (sum(abs(ind - seq_len(K))) > 0) {
        one_draw <- ahmq_permute_draws(
          Q_sample = Q_current[, , i, drop = FALSE],
          G_sample = G_current[, , i, drop = FALSE],
          alpha_sample = alpha_current[, , i, drop = FALSE],
          pi_sample = pi_sample[, i, drop = FALSE],
          perm = ind,
          K = K
        )
        G_re_est[, , i] <- one_draw$G[, , 1L]
        alpha_re_est[, , i] <- one_draw$alpha[, , 1L]
        pi_re_est[, i] <- one_draw$pi[, 1L]
      }
    }
    e <- sum(abs(Q_pre - Q_current))
    if (verbose) {
      cat("\rlabel-switch iter=", iter, "; change=", e)
    }
  }
  if (verbose) {
    cat("\n")
  }

  list(Q = Q_re_est, G = G_re_est, alpha = alpha_re_est, pi = pi_re_est)
}

#' Align posterior samples to known simulation labels
#'
#' In simulation studies the true attribute labels are known, but the MCMC
#' output is still invariant to a common permutation of attribute labels. This
#' function applies one additional global permutation after ordinary
#' label-switch resolution. The permutation is chosen by matching the
#' posterior modal Q-matrix to \code{Q0}; if multiple permutations produce the
#' same Q discrepancy, the posterior modal hierarchy is matched to \code{G0}
#' to break the tie.
#'
#' The selected permutation is applied simultaneously to \code{Q}, \code{G},
#' \code{alpha}, and \code{pi}. This is essential because all four objects are
#' indexed by the same latent attribute order.
#'
#' @param Q_sample Item-by-attribute Q draws (\eqn{J \times K \times T}).
#' @param G_sample Attribute hierarchy draws (\eqn{K \times K \times T}).
#' @param alpha_sample Attribute profile draws (\eqn{N \times K \times T}).
#' @param pi_sample Class probability draws (\eqn{2^K \times T}).
#' @param Q0 True Q-matrix in the desired attribute order.
#' @param G0 Optional true hierarchy in the desired attribute order. When
#'   supplied, it is used to break ties after Q matching.
#' @param cut_value Threshold used to convert posterior mean \code{G} to a
#'   binary hierarchy for tie-breaking.
#' @param K Number of attributes. Defaults to \code{ncol(Q0)}.
#' @param verbose If \code{TRUE}, print the chosen permutation and discrepancies.
#' @return A list with relabeled posterior draws \code{Q}, \code{G},
#'   \code{alpha}, \code{pi}, plus \code{permutation}, \code{Q_distance},
#'   and \code{G_distance}.
#' @export
align_samples_to_truth <- function(Q_sample,
                                   G_sample,
                                   alpha_sample,
                                   pi_sample,
                                   Q0,
                                   G0 = NULL,
                                   cut_value = 0.2,
                                   K = ncol(Q0),
                                   verbose = FALSE)
{
  Q0 <- ahmq_binary_matrix(Q0)
  if (ncol(Q0) != K) {
    stop("Q0 must have K columns.", call. = FALSE)
  }
  if (!is.null(G0)) {
    G0 <- ahmq_binary_matrix(G0)
    if (!identical(dim(G0), c(K, K))) {
      stop("G0 must be a K by K matrix.", call. = FALSE)
    }
  }

  Q_mean <- apply(Q_sample, c(1L, 2L), mean)
  Q_bin <- matrix(0, nrow(Q_mean), ncol(Q_mean))
  Q_bin[Q_mean > 0.5] <- 1

  G_bin <- NULL
  if (!is.null(G0)) {
    G_mean <- apply(G_sample, c(1L, 2L), mean)
    G_bin <- matrix(0, K, K)
    G_bin[G_mean > cut_value] <- 1
    G_bin <- Transitive(G_bin, K)
  }

  permutation <- ahmq_permutations(seq_len(K))
  scores <- data.frame(
    index = seq_along(permutation),
    Q_distance = NA_real_,
    G_distance = NA_real_
  )

  for (i in seq_along(permutation)) {
    perm <- permutation[[i]]
    Q_perm <- Q_bin[, perm, drop = FALSE]
    scores$Q_distance[i] <- sum(abs(Q_perm - Q0))
    if (!is.null(G0)) {
      G_perm <- ahmq_permute_matrix(G_bin, perm)
      G_perm <- Transitive(G_perm, K)
      scores$G_distance[i] <- sum(abs(G_perm - G0))
    }
  }

  min_q <- min(scores$Q_distance)
  candidate <- which(scores$Q_distance == min_q)
  if (!is.null(G0) && length(candidate) > 1L) {
    min_g <- min(scores$G_distance[candidate])
    candidate <- candidate[which(scores$G_distance[candidate] == min_g)]
  }
  selected <- candidate[1L]
  perm <- permutation[[selected]]

  aligned <- ahmq_permute_draws(
    Q_sample = Q_sample,
    G_sample = G_sample,
    alpha_sample = alpha_sample,
    pi_sample = pi_sample,
    perm = perm,
    K = K
  )
  aligned$permutation <- perm
  aligned$Q_distance <- scores$Q_distance[selected]
  aligned$G_distance <- scores$G_distance[selected]
  aligned$alignment_scores <- scores

  if (isTRUE(verbose)) {
    cat("truth-label permutation:",
        paste(perm, collapse = ","), "\n")
    cat("truth-label Q distance:", aligned$Q_distance, "\n")
    if (!is.null(G0)) {
      cat("truth-label G distance:", aligned$G_distance, "\n")
    }
  }

  aligned
}

#' Align an \code{Est_fun} result to known simulation labels
#'
#' Reorders the posterior samples stored in an \code{\link{Est_fun}} result so
#' that the estimated labels match the true simulation labels. After applying
#' the selected permutation to \code{Q}, \code{G}, \code{alpha}, and \code{pi},
#' the function recomputes \code{Est_Q}, \code{Est_G}, \code{Est_alpha},
#' \code{Est_pi}, and \code{Est_GG}. This function is intended for simulation
#' evaluation only; real data do not have a known \code{Q0} or \code{G0}.
#'
#' @param est Output from \code{\link{Est_fun}} with
#'   \code{return_samples = TRUE}.
#' @param Q0 True Q-matrix.
#' @param G0 Optional true hierarchy.
#' @param cut_value Threshold for binarizing posterior mean \code{G}; defaults
#'   to \code{est$cut_value}, then \code{0.2}.
#' @param verbose If \code{TRUE}, print the chosen permutation.
#' @param keep_samples If \code{TRUE}, keep the truth-aligned posterior samples
#'   in \code{est$posterior_samples}. If \code{FALSE}, samples are removed to
#'   reduce object size after recomputing estimates.
#' @return Updated \code{Est_fun} result with truth-aligned estimates and
#'   metadata fields \code{truth_label_permutation},
#'   \code{truth_label_Q_distance}, and \code{truth_label_G_distance}.
#' @export
align_estimates_to_truth <- function(est,
                                     Q0,
                                     G0 = NULL,
                                     cut_value = NULL,
                                     verbose = FALSE,
                                     keep_samples = TRUE)
{
  if (is.null(est$posterior_samples)) {
    stop("est must be created by Est_fun(..., return_samples = TRUE).",
         call. = FALSE)
  }
  if (is.null(cut_value)) {
    cut_value <- est$cut_value
  }
  if (is.null(cut_value)) {
    cut_value <- 0.2
  }

  samples <- est$posterior_samples
  aligned <- align_samples_to_truth(
    Q_sample = samples$Q,
    G_sample = samples$G,
    alpha_sample = samples$alpha,
    pi_sample = samples$pi,
    Q0 = Q0,
    G0 = G0,
    cut_value = cut_value,
    K = ncol(Q0),
    verbose = verbose
  )

  J <- dim(aligned$Q)[1L]
  K <- dim(aligned$Q)[2L]
  N <- dim(aligned$alpha)[1L]
  re_est <- binarize_aligned_samples(
    aligned = aligned,
    s = est$Est_s,
    g = est$Est_g,
    cut_value = cut_value,
    J = J,
    N = N,
    K = K
  )

  est$Est_pi <- re_est$pi
  est$Est_G <- re_est$G
  est$Est_Q <- re_est$Est_Q
  est$Est_alpha <- re_est$alpha
  est$Est_GG <- re_est$GG
  est$cut_value <- cut_value
  est$truth_label_permutation <- aligned$permutation
  est$truth_label_Q_distance <- aligned$Q_distance
  est$truth_label_G_distance <- aligned$G_distance
  est$truth_label_alignment_scores <- aligned$alignment_scores

  if (isTRUE(keep_samples)) {
    est$posterior_samples <- aligned[c("Q", "G", "alpha", "pi")]
  } else {
    est$posterior_samples <- NULL
  }

  est
}
