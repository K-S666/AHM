#' Compute DINA ideal responses for one examinee
#'
#' For one attribute profile \eqn{\alpha_i} and a Q-matrix, this helper returns
#' the deterministic DINA ideal response vector \eqn{\eta_i}. An item has
#' \eqn{\eta_{ij}=1} only when the examinee masters every attribute required by
#' that item's Q-vector.
#'
#' @param alpha Binary attribute vector of length \code{K}.
#' @param Q Binary Q-matrix with dimension \eqn{J \times K}.
#' @return A binary vector of length \code{J}.
#' @noRd
eta_onesubject_dina <- function(alpha, Q)
{
  J <- nrow(Q)
  K <- ncol(Q)
  alphamatrix <- matrix(rep(alpha, J), J, K, byrow = TRUE)^Q
  apply(alphamatrix, 1L, prod)
}

#' Simulate response data from the DINA model
#'
#' Generates a binary response matrix from fixed DINA item parameters,
#' Q-matrix, and examinee attribute profiles. For item \eqn{j}, a capable
#' examinee answers correctly with probability \eqn{1-s_j}; an incapable
#' examinee answers correctly with probability \eqn{g_j}.
#'
#' @param N Number of examinees.
#' @param J Number of items.
#' @param K Number of attributes.
#' @param s Numeric vector of item slipping parameters, length \code{J}.
#' @param g Numeric vector of item guessing parameters, length \code{J}.
#' @param Q Binary Q-matrix with dimension \eqn{J \times K}.
#' @param alpha Binary attribute profile matrix with dimension \eqn{N \times K}.
#' @return Binary response matrix with dimension \eqn{N \times J}.
#' @export
DINA_data <- function(N, J, K, s, g, Q, alpha)
{
  Q <- as.matrix(Q)
  alpha <- as.matrix(alpha)
  storage.mode(Q) <- "numeric"
  storage.mode(alpha) <- "numeric"

  if (!identical(dim(Q), c(as.integer(J), as.integer(K)))) {
    stop("Q must have dimension J x K.", call. = FALSE)
  }
  if (!identical(dim(alpha), c(as.integer(N), as.integer(K)))) {
    stop("alpha must have dimension N x K.", call. = FALSE)
  }
  if (length(s) != J || length(g) != J) {
    stop("s and g must both have length J.", call. = FALSE)
  }
  if (!all(Q %in% c(0, 1)) || !all(alpha %in% c(0, 1))) {
    stop("Q and alpha must be binary.", call. = FALSE)
  }

  r <- matrix(NA_real_, N, J)
  for (i in seq_len(N)) {
    eta <- eta_onesubject_dina(alpha[i, ], Q)
    p <- (1 - s)^eta * g^(1 - eta)
    r[i, ] <- stats::rbinom(J, 1L, p)
  }
  r
}

#' Restrict Q-matrix rows under an attribute hierarchy
#'
#' If attribute \eqn{k} requires prerequisite attributes, every item requiring
#' \eqn{k} must also require those prerequisites. This function closes each row
#' of a candidate Q-matrix under a reachability matrix.
#'
#' @param Q_all Candidate binary Q-matrix with dimension \eqn{J \times K}.
#' @param R Reachability matrix of an attribute hierarchy. Diagonal entries
#'   should be one; \code{R[a, b] = 1} means attribute \code{a} is a prerequisite
#'   of attribute \code{b}.
#' @param J Number of items.
#' @param K Number of attributes.
#' @return A list with \code{Q_restrict}, row decimal codes
#'   \code{Q_restrict_code}, and unique sorted row codes
#'   \code{Q_restrict_index}.
#' @export
Restricted_Q <- function(Q_all, R, J, K)
{
  Q_restrict <- as.matrix(Q_all)
  R <- as.matrix(R)
  if (!identical(dim(Q_restrict), c(as.integer(J), as.integer(K)))) {
    stop("Q_all must have dimension J x K.", call. = FALSE)
  }
  if (!identical(dim(R), c(as.integer(K), as.integer(K)))) {
    stop("R must have dimension K x K.", call. = FALSE)
  }

  for (k in seq_len(K)) {
    pre_index <- which(R[, k] == 1)
    index <- which(Q_restrict[, k] == 1)
    Q_restrict[index, pre_index] <- 1
  }
  Q_restrict_code <- apply(Q_restrict, 1L, twoToten)
  list(
    Q_restrict = Q_restrict,
    Q_restrict_code = Q_restrict_code,
    Q_restrict_index = sort(unique(Q_restrict_code))
  )
}

#' Count rows equal to a Q-vector
#'
#' Counts how many rows of a Q-matrix match a supplied item requirement vector.
#'
#' @param Q Binary Q-matrix.
#' @param q Binary vector with length \code{ncol(Q)}.
#' @return Integer count of matching rows.
#' @export
samevec <- function(Q, q)
{
  Q <- as.matrix(Q)
  sum(apply(Q, 1L, function(row) all(row == q)))
}

ahmq_alpha_space <- function(G)
{
  G <- as.matrix(G)
  K <- ncol(G)
  R <- Reachability(G, K)
  alpha_all <- Trans_10to2_mat(K, 0:(2^K - 1L))
  reduced <- Reduced_alpha(alpha_all, R, K)
  list(
    R = R,
    alpha_all = alpha_all,
    admissible_codes = as.integer(reduced$index),
    reduced_alpha = reduced$alpha_current
  )
}

ahmq_reduce_alpha <- function(alpha, G)
{
  alpha <- as.matrix(alpha)
  K <- ncol(alpha)
  space <- ahmq_alpha_space(G)
  codes <- apply(alpha, 1L, Trans_2to10, K = K)
  space$reduced_alpha[codes + 1L, , drop = FALSE]
}

#' Generate common attribute hierarchy structures
#'
#' Returns binary adjacency matrices for the four hierarchy structures used in
#' the simulation design. Matrix entries follow the convention
#' \code{G[prerequisite, target] = 1}; for example, \code{G[1, 2] = 1}
#' represents \code{A1 -> A2}. The named G1--G4 presets are defined for
#' \code{K = 4} attributes.
#'
#' @param type Hierarchy type. \code{"linear"} is G1
#'   (A1 -> A2 -> A3 -> A4); \code{"convergent"} is G2
#'   (A1 -> A2, A1 -> A3, A2 -> A4, A3 -> A4); \code{"divergent"}
#'   is G3 (A1 -> A2, A1 -> A3, A1 -> A4); \code{"partially_structured"} is G4
#'   partially structured (A1 -> A2, A1 -> A3, A3 -> A4); and \code{"none"} has no edges.
#'   The legacy value \code{"unstructured"} is accepted as an alias for
#'   \code{"partially_structured"}.
#' @param K Number of attributes. The G1--G4 presets require \code{K = 4};
#'   \code{"none"} is available for any positive \code{K}.
#' @return A binary \eqn{K \times K} adjacency matrix.
#' @export
simu_G <- function(type = c("linear", "convergent", "divergent",
                         "partially_structured", "none"),
                   K = 4L)
{
  if (missing(type)) {
    type <- "linear"
  } else {
    type <- match.arg(type, c("linear", "convergent", "divergent",
                             "partially_structured", "none",
                             "unstructured"))
  }
  if (identical(type, "unstructured")) {
    warning(
      "simu_G(\"unstructured\") is deprecated; use ",
      "simu_G(\"partially_structured\") instead.",
      call. = FALSE
    )
    type <- "partially_structured"
  }
  K <- as.integer(K)
  if (length(K) != 1L || is.na(K) || K < 1L) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (type != "none" && K != 4L) {
    stop("The G1--G4 hierarchy presets are defined for K = 4.", call. = FALSE)
  }

  G <- matrix(0, K, K)
  if (type == "linear") {
    G[1L, 2L] <- 1
    G[2L, 3L] <- 1
    G[3L, 4L] <- 1
  } else if (type == "convergent") {
    G[1L, 2L] <- 1
    G[1L, 3L] <- 1
    G[2L, 4L] <- 1
    G[3L, 4L] <- 1
  } else if (type == "divergent") {
    G[1L, 2L] <- 1
    G[1L, 3L] <- 1
    G[1L, 4L] <- 1
  } else if (type == "partially_structured") {
    G[1L, 2L] <- 1
    G[1L, 3L] <- 1
    G[3L, 4L] <- 1
  }
  rownames(G) <- colnames(G) <- paste0("A", seq_len(K))
  G
}
#' Simulate attribute profiles under an attribute hierarchy
#'
#' Generates an examinee-by-attribute matrix under a fixed hierarchy. The
#' \code{"balanced"} condition samples uniformly from the hierarchy-permissible
#' attribute patterns. The \code{"unbalanced"} condition follows the simulation
#' design in the manuscript: draw a latent multivariate normal vector with
#' common off-diagonal correlation \code{sigma}, threshold it at zero, and then
#' map impermissible patterns to their hierarchy-permissible reduced form.
#'
#' @param N Number of examinees.
#' @param G Attribute hierarchy adjacency matrix.
#' @param distribution Either \code{"balanced"} or \code{"unbalanced"}.
#' @param sigma Common off-diagonal covariance/correlation used for the
#'   unbalanced multivariate normal generator.
#' @param seed Optional random seed.
#' @return A list with binary \code{alpha}, class probabilities \code{pi} over
#'   all \eqn{2^K} patterns, decimal pattern \code{codes}, reachability matrix
#'   \code{R}, and the hierarchy-admissible patterns.
#' @export
simu_alpha <- function(N,
                       G,
                       distribution = c("balanced", "unbalanced"),
                       sigma = 0.5,
                       seed = NULL)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  distribution <- match.arg(distribution)
  G <- as.matrix(G)
  K <- ncol(G)
  if (!identical(dim(G), c(K, K))) {
    stop("G must be a square K x K adjacency matrix.", call. = FALSE)
  }

  space <- ahmq_alpha_space(G)
  admissible <- space$alpha_all[space$admissible_codes + 1L, , drop = FALSE]

  if (distribution == "balanced") {
    sampled <- sample(seq_len(nrow(admissible)), N, replace = TRUE)
    alpha <- admissible[sampled, , drop = FALSE]
  } else {
    Sigma <- matrix(sigma, K, K)
    diag(Sigma) <- 1
    z <- matrix(stats::rnorm(N * K), N, K) %*% chol(Sigma)
    alpha_raw <- (z > 0) * 1
    alpha <- ahmq_reduce_alpha(alpha_raw, G)
  }

  codes <- apply(alpha, 1L, Trans_2to10, K = K)
  pi <- tabulate(codes + 1L, nbins = 2^K) / N
  list(
    alpha = alpha,
    pi = pi,
    codes = codes,
    R = space$R,
    admissible_codes = space$admissible_codes,
    admissible_alpha = admissible,
    distribution = distribution
  )
}

densify_Q <- function(Q, G)
{
  Q <- as.matrix(Q)
  K <- ncol(Q)
  R <- Reachability(as.matrix(G), K)
  Restricted_Q(Q, R, nrow(Q), K)$Q_restrict
}

sparsify_Q <- function(Q, G)
{
  Q <- as.matrix(Q)
  G <- as.matrix(G)
  K <- ncol(Q)
  R <- Reachability(G, K)
  out <- Q
  for (j in seq_len(nrow(out))) {
    for (k in seq_len(K)) {
      descendants <- which(R[k, ] == 1)
      descendants <- setdiff(descendants, k)
      if (out[j, k] == 1 && any(out[j, descendants] == 1)) {
        out[j, k] <- 0
      }
    }
  }
  out
}

check_Q_identifiable <- function(Q, G)
{
  Q <- as.matrix(Q)
  K <- ncol(Q)
  if (nrow(Q) < K) {
    return(FALSE)
  }
  Q0 <- Q[seq_len(K), , drop = FALSE]
  Q_star <- Q[-seq_len(K), , drop = FALSE]
  Q_sparse <- sparsify_Q(Q, G)
  Q_star_dense <- densify_Q(Q_star, G)

  cond_A <- all(sparsify_Q(Q0, G) == diag(K))
  cond_B <- all(colSums(Q_sparse) >= 3)
  cond_C <- nrow(unique(t(Q_star_dense), MARGIN = 1L)) == K
  cond_A && cond_B && cond_C
}

#' Simulate an identifiable Q-matrix under an attribute hierarchy
#'
#' Generates a structured Q-matrix satisfying the identifiability conditions
#' described in the manuscript: \code{Q0} is equivalent to the identity matrix
#' under \code{G}; the sparsified Q-matrix has at least three ones in every
#' column; and the densified \code{Q_star} part contains \code{K} distinct
#' column vectors. The returned Q-matrix is the densified, hierarchy-structured
#' version used for data generation.
#'
#' @param J Number of items.
#' @param G Attribute hierarchy adjacency matrix.
#' @param seed Optional random seed.
#' @param max_tries Maximum number of random completion attempts.
#' @return A list with structured \code{Q}, sparsified \code{Q_sparse},
#'   reachability matrix \code{R}, and a logical \code{identifiable}.
#' @export
simu_Q <- function(J, G, seed = NULL, max_tries = 1000L)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  G <- as.matrix(G)
  K <- ncol(G)
  if (!identical(dim(G), c(K, K))) {
    stop("G must be a square K x K adjacency matrix.", call. = FALSE)
  }
  if (J < 3L * K) {
    stop("J must be at least 3*K to guarantee condition (B).",
         call. = FALSE)
  }

  R <- Reachability(G, K)
  candidate_sparse <- Trans_10to2_mat(K, 1:(2^K - 1L))
  Q0_sparse <- diag(K)
  fixed_star_sparse <- rbind(diag(K), diag(K))
  remaining <- J - K - nrow(fixed_star_sparse)

  for (try in seq_len(max_tries)) {
    random_star <- if (remaining > 0L) {
      candidate_sparse[
        sample(seq_len(nrow(candidate_sparse)), remaining, replace = TRUE),
        ,
        drop = FALSE
      ]
    } else {
      matrix(numeric(0), 0L, K)
    }
    Q_sparse <- rbind(Q0_sparse, fixed_star_sparse, random_star)
    Q <- densify_Q(Q_sparse, G)
    if (check_Q_identifiable(Q, G)) {
      return(list(Q = Q,
                  Q_sparse = sparsify_Q(Q, G),
                  R = R,
                  identifiable = TRUE,
                  tries = try))
    }
  }

  stop("Failed to generate an identifiable Q-matrix within max_tries.",
       call. = FALSE)
}

#' Simulate a complete AHMQ/DINA data set
#'
#' Creates a self-contained synthetic data set for testing the Bayesian DINA
#' estimator with an unknown Q-matrix and unknown attribute hierarchy. The
#' function can either use supplied true parameters or generate reasonable
#' defaults: a random Q-matrix, a zero hierarchy, uniform class probabilities
#' over hierarchy-admissible attribute profiles, and moderate slip/guess
#' parameters.
#'
#' @param N Number of examinees.
#' @param J Number of items.
#' @param K Number of attributes.
#' @param Q Optional true Q-matrix. If \code{NULL}, \code{random_Q(J, K)} is
#'   used and then restricted by \code{G}.
#' @param G Optional true attribute hierarchy adjacency matrix. If \code{NULL},
#'   the no-prerequisite hierarchy is used.
#' @param s Optional slipping vector. If scalar, recycled to length \code{J};
#'   if \code{NULL}, sampled uniformly from \code{s_range}.
#' @param g Optional guessing vector. If scalar, recycled to length \code{J};
#'   if \code{NULL}, sampled uniformly from \code{g_range}.
#' @param pi Optional latent class probabilities over all \eqn{2^K} classes.
#'   When \code{NULL}, admissible profiles receive equal probability and
#'   inadmissible profiles receive zero probability.
#' @param alpha Optional true attribute profile matrix. If supplied, \code{pi}
#'   is ignored for profile generation.
#' @param alpha_distribution Attribute profile distribution used when
#'   \code{alpha} is not supplied.
#' @param sigma Common off-diagonal covariance/correlation for unbalanced
#'   \code{\link{simu_alpha}} generation.
#' @param s_range,g_range Length-two numeric ranges used when \code{s} or
#'   \code{g} are generated.
#' @param seed Optional random seed.
#' @return A list containing \code{Y}, \code{Q}, \code{G}, \code{R},
#'   \code{alpha}, \code{pi}, \code{s}, \code{g}, \code{N}, \code{J}, and
#'   \code{K}. \code{Y} can be passed directly to \code{\link{AHMQ}}.
#' @examples
#' set.seed(1)
#' dat <- simulate_ahmq_data(N = 50, J = 12, K = 3)
#' dim(dat$Y)
#' dat$Q
#' @export
simulate_ahmq_data <- function(N, J, K,
                               Q = NULL,
                               G = NULL,
                               s = NULL,
                               g = NULL,
                               pi = NULL,
                               alpha = NULL,
                               alpha_distribution = c("balanced", "unbalanced"),
                               sigma = 0.5,
                               s_range = c(0.05, 0.25),
                               g_range = c(0.05, 0.25),
                               seed = NULL)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  N <- as.integer(N)
  J <- as.integer(J)
  K <- as.integer(K)
  alpha_distribution <- match.arg(alpha_distribution)
  if (N < 1L || J < 1L || K < 1L) {
    stop("N, J, and K must be positive integers.", call. = FALSE)
  }

  if (is.null(G)) {
    G <- matrix(0, K, K)
  } else {
    G <- as.matrix(G)
    if (!identical(dim(G), c(K, K))) {
      stop("G must have dimension K x K.", call. = FALSE)
    }
    if (!all(G %in% c(0, 1))) {
      stop("G must be binary.", call. = FALSE)
    }
    diag(G) <- 0
    G <- Transitive(G, K)
  }
  R <- Reachability(G, K)

  if (is.null(Q)) {
    Q_info <- simu_Q(J, G)
    Q <- Q_info$Q
  } else {
    Q <- as.matrix(Q)
    if (!identical(dim(Q), c(J, K))) {
      stop("Q must have dimension J x K.", call. = FALSE)
    }
    if (!all(Q %in% c(0, 1))) {
      stop("Q must be binary.", call. = FALSE)
    }
  }
  Q <- Restricted_Q(Q, R, J, K)$Q_restrict

  make_item_param <- function(x, range, name) {
    if (is.null(x)) {
      return(stats::runif(J, range[1L], range[2L]))
    }
    if (length(x) == 1L) {
      x <- rep(x, J)
    }
    if (length(x) != J) {
      stop(name, " must be NULL, scalar, or length J.", call. = FALSE)
    }
    as.numeric(x)
  }
  s <- make_item_param(s, s_range, "s")
  g <- make_item_param(g, g_range, "g")
  if (any(s < 0 | s >= 1) || any(g < 0 | g >= 1)) {
    stop("s and g must lie in [0, 1).", call. = FALSE)
  }

  alpha_info <- simu_alpha(N, G, distribution = alpha_distribution,
                           sigma = sigma)
  admissible <- alpha_info$admissible_alpha

  if (is.null(alpha)) {
    if (is.null(pi)) {
      alpha <- alpha_info$alpha
      pi <- alpha_info$pi
    } else {
      pi <- as.numeric(pi)
      if (length(pi) != 2^K) {
        stop("pi must have length 2^K.", call. = FALSE)
      }
      inadmissible <- setdiff(seq_len(2^K), alpha_info$admissible_codes + 1L)
      pi[inadmissible] <- 0
      if (sum(pi) <= 0) {
        stop("pi must assign positive mass to admissible profiles.", call. = FALSE)
      }
      pi <- pi / sum(pi)
      alpha_all <- Trans_10to2_mat(K, 0:(2^K - 1L))
      class_id <- sample(seq_len(2^K), N, replace = TRUE, prob = pi)
      alpha <- alpha_all[class_id, , drop = FALSE]
    }
  } else {
    alpha <- as.matrix(alpha)
    if (!identical(dim(alpha), c(N, K))) {
      stop("alpha must have dimension N x K.", call. = FALSE)
    }
    if (!all(alpha %in% c(0, 1))) {
      stop("alpha must be binary.", call. = FALSE)
    }
    if (is.null(pi)) {
      codes <- apply(alpha, 1L, Trans_2to10, K = K)
      pi <- tabulate(codes + 1L, nbins = 2^K) / N
    }
  }

  Y <- DINA_data(N, J, K, s, g, Q, alpha)
  list(Y = Y, Q = Q, G = G, R = R, alpha = alpha, pi = pi,
       s = s, g = g, N = N, J = J, K = K,
       admissible_alpha = admissible,
       alpha_distribution = alpha_distribution)
}




