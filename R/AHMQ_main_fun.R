default_minibatch_size <- function(N1, N)
{
  if (is.null(N1)) {
    return(max(1L, floor(as.integer(N) / 2L)))
  }
  N1 <- as.integer(N1)
  if (length(N1) != 1L || is.na(N1) || N1 < 1L) {
    stop("N1 must be NULL or a positive integer.", call. = FALSE)
  }
  min(N1, as.integer(N))
}

validate_p_add <- function(p_add)
{
  if (length(p_add) != 1L || is.na(p_add) || p_add <= 0 || p_add >= 1) {
    stop("p_add must be a single value in (0, 1).", call. = FALSE)
  }
  as.numeric(p_add)
}

#' Single-chain DINA estimation with AHM Q-matrix prior
#'
#' Runs one MCMC chain via the Rcpp sampler. Initial values for item
#' parameters, class probabilities, Q, G, and individual attributes are drawn
#' inside the C++ sampler.
#'
#' @param Y An \eqn{N \times J} binary response matrix.
#' @param K Number of attributes.
#' @param N1 Mini-batch sample size for updating item parameters, Q, and G.
#'   If \code{NULL}, the default is \code{floor(N / 2)}, where \code{N} is the
#'   number of examinees.
#' @param chain_length Total MCMC iterations (default 20000).
#' @param burn_in Burn-in iterations (default 10000).
#' @param keep_burnin If \code{TRUE}, keep all iterations in the returned
#'   chain samples. If \code{FALSE}, only post-burn-in draws are stored.
#' @param a_s0,b_s0,a_g0,b_g0 Beta prior hyperparameters for slip and guess.
#' @param p_add Proposal probability for adding an edge in the hierarchy
#'   \code{G} update; removing an edge is proposed with probability
#'   \code{1 - p_add} (default 0.5).
#' @param progress Whether to print per-chain iteration progress messages.
#' @param print_every Print progress every this many MCMC iterations.
#' @return A list with posterior draws: \code{alpha}, \code{s}, \code{g},
#'   \code{pi}, \code{Q}, and \code{G}.
#' @noRd
AHMQ_single <- function(Y, K,
                        N1 = NULL,
                        chain_length = 20000L,
                        burn_in = 10000L,
                        keep_burnin = FALSE,
                        a_s0 = 1.0, a_g0 = 1.0,
                        b_s0 = 1.0, b_g0 = 1.0,
                        p_add = 0.5,
                        chain_id = 1L,
                        progress = TRUE,
                        print_every = 1000L)
{
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be a binary response matrix containing only 0 and 1.", call. = FALSE)
  }
  if (K < 1L) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (burn_in >= chain_length) {
    stop("burn_in must be smaller than chain_length.", call. = FALSE)
  }

  N <- dim(Y)[1]
  N1 <- default_minibatch_size(N1, N)
  print_every <- max(1L, as.integer(print_every))
  p_add <- validate_p_add(p_add)
  stored_burn_in <- if (isTRUE(keep_burnin)) 0L else burn_in

  AHM_Q(Y, K, N1, chain_length, stored_burn_in,
        a_s0, a_g0, b_s0, b_g0, p_add,
        as.integer(chain_id), isTRUE(progress), print_every)
}

#' Multi-chain DINA estimation with AHM Q-matrix prior
#'
#' Runs \code{chain_num} independent MCMC chains. Each chain draws its initial
#' values inside the C++ sampler and returns posterior draws together with the
#' analysis settings used.
#'
#' @param Y An \eqn{N \times J} binary response matrix.
#' @param K Number of attributes.
#' @param N1 Mini-batch sample size for updating item parameters, Q, and G.
#'   If \code{NULL}, the default is \code{floor(N / 2)}, where \code{N} is the
#'   number of examinees.
#' @param chain_length MCMC length per chain (default 20000).
#' @param burn_in Burn-in per chain (default 10000).
#' @param keep_burnin If \code{TRUE}, retain burn-in draws in
#'   \code{chain_sample}. This is useful for diagnostic plots over the full
#'   chain. Point-estimation functions still discard \code{burn_in} draws by
#'   default when this flag is present.
#' @param chain_num Number of independent chains (default 4).
#' @param a_s0,b_s0,a_g0,b_g0 Beta prior hyperparameters for slip and guess.
#' @param p_add Proposal probability for adding an edge in the hierarchy
#'   \code{G} update; removing an edge is proposed with probability
#'   \code{1 - p_add} (default 0.5).
#' @param progress Whether to print per-chain iteration progress messages.
#' @param print_every Print progress every this many MCMC iterations.
#' @param parallel If \code{TRUE}, run independent chains in parallel with
#'   \code{parallel::mclapply} on Unix-like systems or a PSOCK cluster on
#'   Windows.
#' @param n_cores Number of worker processes used when \code{parallel = TRUE}.
#' @return A list of length \code{chain_num}. Each element is a list with
#'   \code{data} (settings) and \code{chain_sample} (posterior draws:
#'   \code{G}, \code{s}, \code{g}, \code{Q}, \code{pi}, \code{alpha}).
#' @export
AHMQ <- function(Y, K,
                 N1 = NULL,
                 chain_length = 20000L,
                 burn_in = 10000L,
                 keep_burnin = FALSE,
                 chain_num = 4L,
                 a_s0 = 1.0, a_g0 = 1.0,
                 b_s0 = 1.0, b_g0 = 1.0,
                 p_add = 0.5,
                 parallel = TRUE,
                 n_cores = min(chain_num, parallel::detectCores(logical = FALSE)),
                 progress = TRUE,
                 print_every = 1000L)
{
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be a binary response matrix containing only 0 and 1.", call. = FALSE)
  }
  if (burn_in >= chain_length) {
    stop("burn_in must be smaller than chain_length.", call. = FALSE)
  }
  if (chain_num < 1L) {
    stop("chain_num must be at least 1.", call. = FALSE)
  }
  N1_default <- is.null(N1)
  N1 <- default_minibatch_size(N1, nrow(Y))
  print_every <- max(1L, as.integer(print_every))
  p_add <- validate_p_add(p_add)

  # Shared settings stored with each chain (for reproducibility / labeling)
  data <- list(Y = Y, K = K, N1 = N1,
               chain_length = chain_length, burn_in = burn_in,
               known_Q = FALSE,
               keep_burnin = keep_burnin,
               N1_default = N1_default,
               stored_burn_in = if (isTRUE(keep_burnin)) 0L else burn_in,
               stored_draws = if (isTRUE(keep_burnin)) chain_length else chain_length - burn_in,
               chain_num = chain_num,
               a_s0 = a_s0, a_g0 = a_g0,
               b_s0 = b_s0, b_g0 = b_g0,
               p_add = p_add,
               parallel = isTRUE(parallel),
               n_cores = if (isTRUE(parallel) && chain_num > 1L) max(1L, min(as.integer(n_cores), chain_num)) else 1L,
               progress = isTRUE(progress),
               print_every = print_every)

  run_chain <- function(i) {
    if (isTRUE(progress)) {
      cat(sprintf("AHMQ chain %d/%d started (N1 = %d, iterations = %d).\n",
                  i, chain_num, N1, chain_length))
      flush.console()
    }
    r <- AHMQ_single(Y = Y, K = K,
                     N1 = N1,
                     chain_length = chain_length,
                     burn_in = burn_in,
                     keep_burnin = keep_burnin,
                     a_s0 = a_s0, a_g0 = a_g0,
                     b_s0 = b_s0, b_g0 = b_g0,
                     p_add = p_add,
                     chain_id = i,
                     progress = progress,
                     print_every = print_every)
    if (isTRUE(progress)) {
      cat(sprintf("AHMQ chain %d/%d finished.\n", i, chain_num))
      flush.console()
    }
    chain_sample <- list(G = r$G, s = r$s, g = r$g,
                         Q = r$Q, pi = r$pi, alpha = r$alpha)
    list(data = data, chain_sample = chain_sample)
  }

  if (isTRUE(parallel) && chain_num > 1L) {
    n_cores <- max(1L, min(as.integer(n_cores), chain_num))
    if (.Platform$OS.type == "windows") {
      cl <- if (isTRUE(progress)) {
        parallel::makeCluster(n_cores, outfile = "")
      } else {
        parallel::makeCluster(n_cores)
      }
      on.exit(parallel::stopCluster(cl), add = TRUE)
      lib_paths <- .libPaths()
      parallel::clusterExport(cl, "lib_paths", envir = environment())
      parallel::clusterEvalQ(cl, .libPaths(lib_paths))
      dev_dir <- Sys.getenv("AHM_DEV_DIR", "")
      if (nzchar(dev_dir) && dir.exists(dev_dir)) {
        parallel::clusterExport(cl, "dev_dir", envir = environment())
        parallel::clusterEvalQ(cl, {
          suppressPackageStartupMessages(library(Rcpp))
          dll_path <- file.path(dev_dir, "src", "AHM.dll")
          if (!file.exists(dll_path)) {
            stop("Cannot find development DLL at ", dll_path)
          }
          dyn.load(normalizePath(dll_path))
          regs <- getDLLRegisteredRoutines("AHM")$.Call
          for (nm in names(regs)) {
            assign(nm, getNativeSymbolInfo(nm, PACKAGE = "AHM"), envir = .GlobalEnv)
          }
          for (f in list.files(file.path(dev_dir, "R"), full.names = TRUE, pattern = "[.]R$")) {
            source(f, local = .GlobalEnv)
          }
          NULL
        })
      } else {
        parallel::clusterEvalQ(cl, library(AHM))
      }
      parallel::clusterExport(
        cl,
        varlist = c("AHMQ_single", "Y", "K", "N1", "chain_length", "burn_in",
                    "keep_burnin", "a_s0", "a_g0", "b_s0", "b_g0",
                    "p_add", "data", "chain_num", "progress", "print_every",
                    "default_minibatch_size"),
        envir = environment()
      )
      parallel::clusterSetRNGStream(cl)
      result <- parallel::parLapply(cl, seq_len(chain_num), run_chain)
    } else {
      result <- parallel::mclapply(seq_len(chain_num), run_chain,
                                   mc.cores = n_cores)
    }
  } else {
    result <- lapply(seq_len(chain_num), run_chain)
  }

  result
}

#' Single-chain DINA estimation with a fixed Q-matrix
#'
#' Runs one MCMC chain for the model where the Q-matrix is known. The sampler
#' updates item slip/guess parameters, individual attributes, class
#' probabilities, and the attribute hierarchy, while keeping Q fixed.
#'
#' @param Y An \eqn{N \times J} binary response matrix.
#' @param Q Known binary \eqn{J \times K} Q-matrix.
#' @param N1 Mini-batch sample size for updating item parameters and G. Q is
#'   fixed in this model. If \code{NULL}, the default is \code{floor(N / 2)},
#'   where \code{N} is the number of examinees.
#' @param chain_length Total MCMC iterations.
#' @param burn_in Burn-in iterations.
#' @param keep_burnin If \code{TRUE}, keep all iterations in returned chain
#'   samples for diagnostics. Point estimation still discards burn-in by
#'   default.
#' @param a_s0,b_s0,a_g0,b_g0 Beta prior hyperparameters for slip and guess.
#' @param p_add Proposal probability for adding an edge in the hierarchy
#'   \code{G} update; removing an edge is proposed with probability
#'   \code{1 - p_add} (default 0.5).
#' @param progress Whether to print per-chain iteration progress messages.
#' @param print_every Print progress every this many MCMC iterations.
#' @return A list with posterior draws \code{alpha}, \code{s}, \code{g},
#'   \code{pi}, and \code{G}. Q is fixed and is not sampled.
#' @noRd
AHM_single <- function(Y, Q,
                       N1 = NULL,
                       chain_length = 20000L,
                       burn_in = 10000L,
                       keep_burnin = FALSE,
                       a_s0 = 1.0, a_g0 = 1.0,
                       b_s0 = 1.0, b_g0 = 1.0,
                       p_add = 0.5,
                       chain_id = 1L,
                       progress = TRUE,
                       print_every = 1000L)
{
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  storage.mode(Y) <- "double"
  storage.mode(Q) <- "double"
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be a binary response matrix containing only 0 and 1.", call. = FALSE)
  }
  if (!all(Q %in% c(0, 1))) {
    stop("Q must be a binary matrix containing only 0 and 1.", call. = FALSE)
  }
  if (ncol(Y) != nrow(Q)) {
    stop("ncol(Y) must equal nrow(Q).", call. = FALSE)
  }
  if (burn_in >= chain_length) {
    stop("burn_in must be smaller than chain_length.", call. = FALSE)
  }

  stored_burn_in <- if (isTRUE(keep_burnin)) 0L else burn_in
  N1 <- default_minibatch_size(N1, nrow(Y))
  print_every <- max(1L, as.integer(print_every))
  p_add <- validate_p_add(p_add)
  AHM_fixedQ(Y, Q, N1, chain_length, stored_burn_in,
                 a_s0, a_g0, b_s0, b_g0, p_add,
                 as.integer(chain_id), isTRUE(progress), print_every)
}

#' Multi-chain DINA estimation with a known Q-matrix
#'
#' \code{AHM} is the fixed-Q counterpart of \code{\link{AHMQ}}. It estimates
#' item parameters, individual attributes, class probabilities, and the
#' attribute hierarchy while treating the supplied Q-matrix as known. Because Q
#' anchors the attribute labels, the default post-processing for \code{AHM}
#' output does not perform Q-based label-switch correction and does not report
#' Q recovery metrics.
#'
#' @param Y An \eqn{N \times J} binary response matrix.
#' @param Q Known binary \eqn{J \times K} Q-matrix.
#' @param N1 Mini-batch sample size for updating item parameters and G. Q is
#'   fixed in this model. If \code{NULL}, the default is \code{floor(N / 2)},
#'   where \code{N} is the number of examinees.
#' @param chain_length MCMC length per chain.
#' @param burn_in Burn-in per chain.
#' @param keep_burnin If \code{TRUE}, retain burn-in draws for diagnostic
#'   plots over the full chain. Point-estimation functions still discard
#'   \code{burn_in} draws by default.
#' @param chain_num Number of independent chains.
#' @param a_s0,b_s0,a_g0,b_g0 Beta prior hyperparameters for slip and guess.
#' @param p_add Proposal probability for adding an edge in the hierarchy
#'   \code{G} update; removing an edge is proposed with probability
#'   \code{1 - p_add} (default 0.5).
#' @param progress Whether to print per-chain iteration progress messages.
#' @param print_every Print progress every this many MCMC iterations.
#' @param parallel If \code{TRUE}, run chains in parallel.
#' @param n_cores Number of worker processes used when \code{parallel = TRUE}.
#' @return A list of chain outputs with model settings and posterior draws.
#' @export
AHM <- function(Y, Q,
                N1 = NULL,
                chain_length = 20000L,
                burn_in = 10000L,
                keep_burnin = FALSE,
                chain_num = 4L,
                a_s0 = 1.0, a_g0 = 1.0,
                b_s0 = 1.0, b_g0 = 1.0,
                       p_add = 0.5,
                parallel = TRUE,
                n_cores = min(chain_num, parallel::detectCores(logical = FALSE)),
                progress = TRUE,
                print_every = 1000L)
{
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  storage.mode(Y) <- "double"
  storage.mode(Q) <- "double"
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be a binary response matrix containing only 0 and 1.", call. = FALSE)
  }
  if (!all(Q %in% c(0, 1))) {
    stop("Q must be a binary matrix containing only 0 and 1.", call. = FALSE)
  }
  if (ncol(Y) != nrow(Q)) {
    stop("ncol(Y) must equal nrow(Q).", call. = FALSE)
  }
  if (burn_in >= chain_length) {
    stop("burn_in must be smaller than chain_length.", call. = FALSE)
  }
  if (chain_num < 1L) {
    stop("chain_num must be at least 1.", call. = FALSE)
  }
  K <- ncol(Q)
  N1_default <- is.null(N1)
  N1 <- default_minibatch_size(N1, nrow(Y))
  print_every <- max(1L, as.integer(print_every))
  p_add <- validate_p_add(p_add)

  data <- list(Y = Y, K = K, Q_fixed = Q, N1 = N1,
               chain_length = chain_length, burn_in = burn_in,
               known_Q = TRUE,
               keep_burnin = keep_burnin,
               N1_default = N1_default,
               stored_burn_in = if (isTRUE(keep_burnin)) 0L else burn_in,
               stored_draws = if (isTRUE(keep_burnin)) chain_length else chain_length - burn_in,
               chain_num = chain_num,
               a_s0 = a_s0, a_g0 = a_g0,
               b_s0 = b_s0, b_g0 = b_g0,
               p_add = p_add,
               parallel = isTRUE(parallel),
               n_cores = if (isTRUE(parallel) && chain_num > 1L) max(1L, min(as.integer(n_cores), chain_num)) else 1L,
               progress = isTRUE(progress),
               print_every = print_every)

  run_chain <- function(i) {
    if (isTRUE(progress)) {
      cat(sprintf("AHM chain %d/%d started (N1 = %d, iterations = %d).\n",
                  i, chain_num, N1, chain_length))
      flush.console()
    }
    r <- AHM_single(Y = Y, Q = Q,
                    N1 = N1,
                    chain_length = chain_length,
                    burn_in = burn_in,
                    keep_burnin = keep_burnin,
                    a_s0 = a_s0, a_g0 = a_g0,
                    b_s0 = b_s0, b_g0 = b_g0,
                    p_add = p_add,
                    chain_id = i,
                    progress = progress,
                    print_every = print_every)
    if (isTRUE(progress)) {
      cat(sprintf("AHM chain %d/%d finished.\n", i, chain_num))
      flush.console()
    }
    chain_sample <- list(G = r$G, s = r$s, g = r$g,
                         pi = r$pi, alpha = r$alpha)
    list(data = data, chain_sample = chain_sample)
  }

  if (isTRUE(parallel) && chain_num > 1L) {
    n_cores <- max(1L, min(as.integer(n_cores), chain_num))
    if (.Platform$OS.type == "windows") {
      cl <- if (isTRUE(progress)) {
        parallel::makeCluster(n_cores, outfile = "")
      } else {
        parallel::makeCluster(n_cores)
      }
      on.exit(parallel::stopCluster(cl), add = TRUE)
      lib_paths <- .libPaths()
      parallel::clusterExport(cl, "lib_paths", envir = environment())
      parallel::clusterEvalQ(cl, .libPaths(lib_paths))
      dev_dir <- Sys.getenv("AHM_DEV_DIR", "")
      if (nzchar(dev_dir) && dir.exists(dev_dir)) {
        parallel::clusterExport(cl, "dev_dir", envir = environment())
        parallel::clusterEvalQ(cl, {
          suppressPackageStartupMessages(library(Rcpp))
          dll_path <- file.path(dev_dir, "src", "AHM.dll")
          if (!file.exists(dll_path)) {
            stop("Cannot find development DLL at ", dll_path)
          }
          dyn.load(normalizePath(dll_path))
          regs <- getDLLRegisteredRoutines("AHM")$.Call
          for (nm in names(regs)) {
            assign(nm, getNativeSymbolInfo(nm, PACKAGE = "AHM"), envir = .GlobalEnv)
          }
          for (f in list.files(file.path(dev_dir, "R"), full.names = TRUE, pattern = "[.]R$")) {
            source(f, local = .GlobalEnv)
          }
          NULL
        })
      } else {
        parallel::clusterEvalQ(cl, library(AHM))
      }
      parallel::clusterExport(
        cl,
        varlist = c("AHM_single", "Y", "Q", "N1", "chain_length", "burn_in",
                    "keep_burnin", "a_s0", "a_g0", "b_s0", "b_g0",
                    "p_add", "data", "chain_num", "progress", "print_every",
                    "default_minibatch_size"),
        envir = environment()
      )
      parallel::clusterSetRNGStream(cl)
      result <- parallel::parLapply(cl, seq_len(chain_num), run_chain)
    } else {
      result <- parallel::mclapply(seq_len(chain_num), run_chain,
                                   mc.cores = n_cores)
    }
  } else {
    result <- lapply(seq_len(chain_num), run_chain)
  }

  result
}

