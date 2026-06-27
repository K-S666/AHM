#' AHM: Bayesian Estimation via MCMC for the DINA Model
#'
#' The package implements Bayesian estimation via MCMC for the DINA model with
#' an unknown attribute hierarchy, allowing the Q-matrix to be either unknown
#' and estimated jointly or supplied as known. It provides multi-chain samplers,
#' DIC-based chain selection, Rhat diagnostics for item parameters,
#' label-switching post-processing, simulation utilities, and an optional Shiny
#' interface.
#'
#' @keywords internal
#' @useDynLib AHM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils flush.console
#' @rawNamespace export(compute_DIC)
#' @rawNamespace export(compute_DIC_fixedQ)
"_PACKAGE"



