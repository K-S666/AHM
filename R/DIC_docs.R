#' Compute DIC for fitted DINA MCMC chains
#'
#' These low-level helpers compute deviance information criterion (DIC) values
#' from posterior samples. They are primarily called by \code{Est_fun()} and
#' \code{select_chain_by_DIC()}, but are exported for users who need direct
#' access to the model comparison calculation.
#'
#' @param y Binary response matrix with dimension \eqn{N \times J}.
#' @param N Number of examinees.
#' @param J Number of items.
#' @param K Number of attributes.
#' @param chain Number of posterior draws used in the calculation.
#' @param s Posterior draws of slip parameters, stored as a \eqn{J \times T}
#'   matrix.
#' @param g Posterior draws of guess parameters, stored as a \eqn{J \times T}
#'   matrix.
#' @param alpha Posterior draws of attribute profiles, stored as an
#'   \eqn{N \times K \times T} array.
#' @param Q For \code{compute_DIC()}, posterior draws of Q, stored as a
#'   \eqn{J \times K \times T} array. For \code{compute_DIC_fixedQ()}, the
#'   fixed \eqn{J \times K} Q-matrix.
#' @return A scalar DIC value.
#' @name compute_DIC
#' @aliases compute_DIC_fixedQ
#' @usage
#' compute_DIC(y, N, J, K, chain, s, g, alpha, Q)
#' compute_DIC_fixedQ(y, N, J, K, chain, s, g, alpha, Q)
NULL



