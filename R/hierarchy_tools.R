#' Transitive reduction of an attribute hierarchy
#'
#' Returns the minimal directed acyclic graph (DAG) that preserves the same
#' reachability relation among attributes. In an attribute hierarchy matrix,
#' \code{G[a, b] = 1} means attribute \code{a} is a prerequisite of
#' attribute \code{b}. If an edge can be removed without changing which
#' attributes are reachable from which other attributes, that edge is removed.
#'
#' @param G Binary square adjacency matrix for an attribute hierarchy.
#' @param check_acyclic If \code{TRUE}, stop when \code{G} contains a directed
#'   cycle. Attribute hierarchies should be DAGs.
#' @return A binary square matrix with redundant transitive edges removed.
#' @examples
#' G <- matrix(0, 4, 4)
#' G[1, 2] <- G[2, 3] <- G[1, 3] <- G[3, 4] <- 1
#' transitive_reduction_G(G)
#' @export
transitive_reduction_G <- function(G, check_acyclic = TRUE)
{
  G <- as.matrix(G)
  if (nrow(G) != ncol(G)) {
    stop("G must be a square adjacency matrix.", call. = FALSE)
  }
  if (!all(G %in% c(0, 1))) {
    stop("G must be binary.", call. = FALSE)
  }
  K <- nrow(G)
  diag(G) <- 0
  R <- Reachability(G, K)
  off_diag <- row(R) != col(R)
  if (isTRUE(check_acyclic) && any(R[off_diag] == 1 & t(R)[off_diag] == 1)) {
    stop("G must be a directed acyclic graph.", call. = FALSE)
  }
  target_R <- R
  out <- G
  edges <- which(out == 1, arr.ind = TRUE)
  if (nrow(edges) > 0L) {
    for (edge_id in seq_len(nrow(edges))) {
      candidate <- out
      candidate[edges[edge_id, 1L], edges[edge_id, 2L]] <- 0
      candidate_R <- Reachability(candidate, K)
      if (all(candidate_R == target_R)) {
        out <- candidate
      }
    }
  }
  rownames(out) <- rownames(G)
  colnames(out) <- colnames(G)
  storage.mode(out) <- "numeric"
  out
}
