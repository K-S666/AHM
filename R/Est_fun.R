#' Posterior column indices for \code{AHMQ} chain samples
#'
#' @param result Output of \code{AHMQ}.
#' @param burn_in Burn-in columns to drop; \code{NULL} means none.
#' @param chain_length Total stored iterations; \code{NULL} means infer from data.
#' @return List with \code{post_idx}, \code{sample_size}, \code{burn_in},
#'   \code{chain_length}, and dimensions \code{Y}, \code{N}, \code{J}, \code{K}, \code{L}.
#' @noRd
ahmq_post_idx <- function(result, burn_in = NULL, chain_length = NULL)
{
  cs <- result[[1]]$chain_sample
  stored_length <- ncol(cs$s)
  data <- result[[1]]$data
  keep_burnin <- isTRUE(data$keep_burnin)

  if (is.null(burn_in)) {
    discard <- if (keep_burnin && !is.null(data$burn_in)) {
      as.integer(data$burn_in)
    } else {
      0L
    }
    chain_length <- stored_length
  } else if (!is.null(chain_length) && stored_length == chain_length - burn_in) {
    # AHMQ stores only post-burn-in draws; this branch supports callers who
    # pass the original MCMC chain_length/burn_in pair from the sampler.
    discard <- 0L
    chain_length <- stored_length
  } else {
    discard <- as.integer(burn_in)
    chain_length <- if (is.null(chain_length)) stored_length else min(chain_length, stored_length)
  }

  if (discard < 0L || discard >= chain_length) {
    stop("burn_in/discard leaves no posterior draws to summarize.", call. = FALSE)
  }

  list(
    post_idx = (discard + 1L):chain_length,
    sample_size = chain_length - discard,
    burn_in = discard,
    chain_length = chain_length,
    Y = result[[1]]$data$Y,
    N = nrow(result[[1]]$data$Y),
    J = ncol(result[[1]]$data$Y),
    K = nrow(cs$G),
    L = 2^nrow(cs$G)
  )
}

#' Select the MCMC chain with minimum DIC
#'
#' @param result Output of \code{AHMQ}.
#' @param post_idx Column indices of posterior draws to use.
#' @param sample_size Number of posterior iterations (for \code{compute_DIC}).
#' @return List with \code{best_chain}, vector \code{DIC}, scalar
#'   \code{DIC_value}, posterior means \code{s}, \code{g}, and arrays
#'   \code{Q_sample}, \code{G_sample}, \code{alpha_sample}, \code{pi_sample}.
#' @export
select_chain_by_DIC <- function(result, post_idx, sample_size)
{
  Y <- result[[1]]$data$Y
  N <- nrow(Y)
  J <- ncol(Y)
  K <- nrow(result[[1]]$chain_sample$G)
  n_chain <- length(result)
  DIC <- numeric(n_chain)
  known_Q <- isTRUE(result[[1]]$data$known_Q)
  Q_fixed <- result[[1]]$data$Q_fixed

  for (cl in seq_len(n_chain)) {
    cs <- result[[cl]]$chain_sample
    if (known_Q) {
      DIC[cl] <- compute_DIC_fixedQ(
        Y, N, J, K, sample_size,
        cs$s[, post_idx],
        cs$g[, post_idx],
        cs$alpha[, , post_idx],
        Q_fixed
      )
    } else {
      DIC[cl] <- compute_DIC(
        Y, N, J, K, sample_size,
        cs$s[, post_idx],
        cs$g[, post_idx],
        cs$alpha[, , post_idx],
        cs$Q[, , post_idx]
      )
    }
  }
  best_chain <- which.min(DIC)[1]
  cs <- result[[best_chain]]$chain_sample

  out <- list(
    best_chain = best_chain,
    DIC = DIC,
    DIC_value = DIC[best_chain],
    s = apply(cs$s[, post_idx], 1, mean),
    g = apply(cs$g[, post_idx], 1, mean),
    s_sd = apply(cs$s[, post_idx, drop = FALSE], 1, stats::sd),
    g_sd = apply(cs$g[, post_idx, drop = FALSE], 1, stats::sd),
    G_sample = cs$G[, , post_idx],
    alpha_sample = cs$alpha[, , post_idx],
    pi_sample = cs$pi[, post_idx]
  )
  if (known_Q) {
    out$Q_fixed <- Q_fixed
  } else {
    out$Q_sample <- cs$Q[, , post_idx]
  }
  out
}

#' Gelman-Rubin Rhat for AHMQ chains
#'
#' Computes potential scale reduction factors across MCMC chains for selected
#' numeric sample blocks. Matrices and arrays are flattened so each scalar
#' element receives its own univariate Rhat. A multivariate PSRF (\code{mpsrf})
#' is also returned for each parameter block, matching the role of
#' \code{coda::gelman.diag(..., multivariate = TRUE)$mpsrf} used in the
#' original analysis scripts.
#'
#' @param result Output of \code{\link{AHMQ}}.
#' @param post_idx Posterior column indices to use. If \code{NULL}, the
#'   recorded burn-in is discarded when the result was created with
#'   \code{keep_burnin = TRUE}; otherwise all stored draws are used.
#' @param parameters Character vector of chain sample components to diagnose.
#'   Only \code{"s"} and \code{"g"} are used. Other sampled blocks are omitted
#'   because label switching and discreteness make their scalar Rhat diagnostics
#'   unstable for this model.
#' @return A list with per-parameter univariate Rhat vectors, maximum
#'   univariate Rhat values, and multivariate \code{mpsrf} values.
#' @export
compute_Rhat <- function(result,
                         post_idx = NULL,
                         parameters = c("s", "g"))
{
  parameters <- intersect(parameters, c("s", "g"))
  n_chain <- length(result)
  if (n_chain < 2L) {
    empty <- stats::setNames(vector("list", length(parameters)), parameters)
    return(list(values = empty,
                max = stats::setNames(rep(NA_real_, length(parameters)), parameters),
                mpsrf = stats::setNames(rep(NA_real_, length(parameters)), parameters)))
  }
  if (is.null(post_idx)) {
    post_idx <- ahmq_post_idx(result)$post_idx
  }
  available <- vapply(parameters, function(param) {
    !is.null(result[[1]]$chain_sample[[param]])
  }, logical(1L))
  if (any(!available)) {
    warning("Skipping unavailable chain sample component(s): ",
            paste(parameters[!available], collapse = ", "), call. = FALSE)
    parameters <- parameters[available]
  }
  if (length(parameters) == 0L) {
    return(list(values = list(), max = numeric(0L), mpsrf = numeric(0L)))
  }

  rhat_one <- function(draw_list) {
    m <- length(draw_list)
    n <- ncol(draw_list[[1]])
    if (n < 2L) {
      return(rep(NA_real_, nrow(draw_list[[1]])))
    }

    chain_means <- vapply(draw_list, rowMeans, numeric(nrow(draw_list[[1]])))
    chain_vars <- vapply(draw_list, function(x) {
      apply(x, 1L, stats::var)
    }, numeric(nrow(draw_list[[1]])))
    W <- rowMeans(chain_vars)
    B <- n * apply(chain_means, 1L, stats::var)
    var_hat <- ((n - 1) / n) * W + B / n
    rhat <- sqrt(var_hat / W)
    rhat[!is.finite(rhat)] <- NA_real_
    rhat
  }

  mpsrf_one <- function(draw_list) {
    m <- length(draw_list)
    n <- ncol(draw_list[[1]])
    p <- nrow(draw_list[[1]])
    if (n < 2L || p < 2L) {
      return(NA_real_)
    }

    chain_means <- vapply(draw_list, rowMeans, numeric(p))
    grand_mean <- rowMeans(chain_means)
    centered_means <- sweep(t(chain_means), 2L, grand_mean)
    B_over_n <- stats::cov(centered_means)
    W <- Reduce(`+`, lapply(draw_list, function(x) stats::cov(t(x)))) / m

    eig <- tryCatch({
      # A tiny ridge keeps the diagnostic defined when binary/near-constant
      # chains make the within-chain covariance singular.
      W_ridge <- W + diag(1e-10, p)
      eigen(solve(W_ridge, B_over_n), only.values = TRUE)$values
    }, error = function(e) NA_real_)
    lambda_max <- max(Re(eig), na.rm = TRUE)
    if (!is.finite(lambda_max)) {
      return(NA_real_)
    }
    sqrt((n - 1) / n + (1 + 1 / m) * lambda_max)
  }

  values <- stats::setNames(vector("list", length(parameters)), parameters)
  mpsrf <- stats::setNames(rep(NA_real_, length(parameters)), parameters)
  for (param in parameters) {
    draw_list <- lapply(result, function(chain) {
      x <- chain$chain_sample[[param]]
      if (is.null(x)) {
        stop("Unknown chain sample component: ", param, call. = FALSE)
      }
      if (length(dim(x)) == 2L) {
        return(as.matrix(x[, post_idx, drop = FALSE]))
      }
      if (length(dim(x)) == 3L) {
        d <- dim(x)
        matrix(x[, , post_idx, drop = FALSE], nrow = d[1L] * d[2L])
      } else {
        stop("Unsupported sample shape for parameter: ", param, call. = FALSE)
      }
    })
    values[[param]] <- rhat_one(draw_list)
    mpsrf[[param]] <- mpsrf_one(draw_list)
  }

  list(values = values,
       max = vapply(values, function(x) {
         if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
       }, numeric(1L)),
       mpsrf = mpsrf)
}

#' Binarize aligned posterior samples and build restricted Q
#'
#' @param aligned List from \code{resolve_label_switch}.
#' @param s,g Slip and guess posterior means.
#' @param cut_value Threshold for binarizing \code{G}; \code{NULL} uses the
#'   package default \code{0.2}.
#' @param J,N,K Dimensions.
#' @return List with binarized \code{Q}, \code{G}, \code{alpha}, \code{pi},
#'   continuous \code{GG}, and restricted \code{Est_Q}.
#' @noRd
binarize_aligned_samples <- function(aligned, s, g, cut_value, J, N, K)
{
  Q_re_est <- aligned$Q
  G_re_est <- aligned$G
  alpha_re_est <- aligned$alpha
  pi_re_est <- aligned$pi

  QQ <- apply(Q_re_est, c(1, 2), mean)
  Q <- matrix(0, J, K)
  Q[which(QQ > 0.5)] <- 1

  GG <- apply(G_re_est, c(1, 2), mean)
  if (is.null(cut_value)) {
    cut_value <- 0.2
  }
  G <- matrix(0, K, K)
  G[which(GG > cut_value)] <- 1

  AA <- apply(alpha_re_est, c(1, 2), mean)
  alpha <- matrix(0, N, K)
  alpha[which(AA > 0.5)] <- 1

  pi <- apply(pi_re_est, 1, mean)
  pi <- pi / sum(pi)
  pi_sd <- apply(pi_re_est, 1, stats::sd)

  G <- Transitive(G, K)
  R <- Reachability(G, K)
  Est_Q <- Restricted_Q(Q, R, J, K)$Q_restrict

  list(Q = Q, G = G, alpha = alpha, pi = pi, pi_sd = pi_sd,
       GG = GG, Est_Q = Est_Q)
}

binarize_fixedQ_samples <- function(G_sample,
                                    alpha_sample,
                                    pi_sample,
                                    Q_fixed,
                                    s,
                                    g,
                                    cut_value,
                                    N,
                                    K)
{
  if (is.null(cut_value)) {
    cut_value <- 0.2
  }

  GG <- apply(G_sample, c(1, 2), mean)
  G <- matrix(0, K, K)
  G[which(GG > cut_value)] <- 1
  G <- Transitive(G, K)

  AA <- apply(alpha_sample, c(1, 2), mean)
  alpha <- matrix(0, N, K)
  alpha[which(AA > 0.5)] <- 1

  pi <- apply(pi_sample, 1, mean)
  pi <- pi / sum(pi)
  pi_sd <- apply(pi_sample, 1, stats::sd)

  list(s = s, g = g, pi = pi, pi_sd = pi_sd, G = G, GG = GG,
       alpha = alpha, Q_fixed = Q_fixed)
}

#' Point estimates from multi-chain \code{AHMQ} output
#'
#' Post-processes the list returned by \code{AHMQ}: selects the chain with
#' minimum DIC, resolves attribute label switching, binarizes parameters, and
#' returns hierarchy-consistent Q-matrix estimates.
#'
#' @param result Output of \code{AHMQ} (list of chains, each with \code{data}
#'   and \code{chain_sample}).
#' @param burn_in Burn-in columns to drop; \code{NULL} uses all stored draws.
#'   If the \code{AHMQ} result was created with \code{keep_burnin = TRUE},
#'   \code{NULL} automatically uses the burn-in value recorded in the result.
#' @param chain_length Total stored iterations; \code{NULL} inferred from data.
#' @param cut_value Threshold on posterior mean of \code{G}; \code{NULL} uses
#'   the package default \code{0.2}.
#' @param compute_rhat If \code{TRUE}, compute Rhat diagnostics across all
#'   chains before DIC chain selection.
#' @param rhat_parameters Chain sample components for \code{\link{compute_Rhat}}.
#'   Only \code{"s"} and \code{"g"} are used.
#' @param return_samples If \code{TRUE}, include the DIC-selected, relabeled
#'   posterior samples in the returned object.
#' @param verbose If \code{TRUE}, print label-switch progress.
#' @return A list with posterior means \code{Est_s}, \code{Est_g},
#'   \code{Est_pi}, posterior standard deviations \code{Est_s_sd},
#'   \code{Est_g_sd}, \code{Est_pi_sd}, estimated \code{Est_G},
#'   \code{Est_Q}, \code{Est_alpha}, \code{Est_GG}, model-fit values
#'   \code{DIC}, \code{DIC_all}, \code{best_chain}, \code{runtime}, and
#'   optionally \code{Rhat} and \code{posterior_samples}.
#' @seealso \code{\link{select_chain_by_DIC}}, \code{\link{resolve_label_switch}}
#' @export
Est_fun <- function(result,
                    burn_in = NULL,
                    chain_length = NULL,
                    cut_value = 0.2,
                    compute_rhat = TRUE,
                    rhat_parameters = c("s", "g"),
                    return_samples = FALSE,
                    verbose = TRUE)
{
  idx <- ahmq_post_idx(result, burn_in, chain_length)
  known_Q <- isTRUE(result[[1]]$data$known_Q)
  rhat_parameters <- intersect(rhat_parameters, c("s", "g"))
  rhat <- NULL
  if (compute_rhat) {
    rhat <- compute_Rhat(result, idx$post_idx, rhat_parameters)
  }
  sel <- select_chain_by_DIC(result, idx$post_idx, idx$sample_size)
  if (known_Q) {
    est <- binarize_fixedQ_samples(
      sel$G_sample, sel$alpha_sample, sel$pi_sample, sel$Q_fixed,
      sel$s, sel$g, cut_value, idx$N, idx$K
    )
  } else {
    aligned <- resolve_label_switch(
      sel$Q_sample, sel$G_sample, sel$alpha_sample, sel$pi_sample,
      K = idx$K,
      verbose = verbose
    )
    est <- binarize_aligned_samples(
      aligned, sel$s, sel$g, cut_value, idx$J, idx$N, idx$K
    )
  }

  runtime <- attr(result, 'runtime', exact = TRUE) %||% result[[1]]$data$runtime

  out <- list(Est_s = sel$s,
              Est_g = sel$g,
              Est_pi = est$pi,
              Est_s_sd = sel$s_sd,
              Est_g_sd = sel$g_sd,
              Est_pi_sd = est$pi_sd,
              Est_G = est$G,
              Est_Q = if (known_Q) NULL else est$Est_Q,
              Q_fixed = if (known_Q) est$Q_fixed else NULL,
              Est_alpha = est$alpha,
              Est_GG = est$GG,
              cut_value = cut_value,
              known_Q = known_Q,
              model = if (known_Q) "AHM" else "AHMQ",
              DIC = sel$DIC_value,
              DIC_all = sel$DIC,
              best_chain = sel$best_chain,
              Rhat = rhat,
              runtime = runtime)
  if (return_samples) {
    out$posterior_samples <- if (known_Q) {
      list(G = sel$G_sample, alpha = sel$alpha_sample, pi = sel$pi_sample)
    } else {
      aligned
    }
  }
  out
}







