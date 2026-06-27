#' Compute Rhat trajectories over chain length
#'
#' Computes Gelman-Rubin Rhat repeatedly on prefixes of MCMC draws. By default
#' the recorded burn-in is excluded, matching the final Rhat diagnostics used
#' in summaries. When \code{include_burnin = TRUE} and \code{\link{AHMQ}} was
#' run with \code{keep_burnin = TRUE}, prefixes start from iteration 1 of the
#' full chain and the returned object records the burn-in boundary for plotting.
#'
#' @param result Output from \code{\link{AHMQ}}.
#' @param parameters Chain components passed to \code{\link{compute_Rhat}}.
#'   Defaults to \code{c("s", "g")}.
#' @param step Prefix spacing in stored draws.
#' @param min_draws First prefix length to evaluate.
#' @param max_draws Last prefix length; \code{NULL} uses all stored draws.
#' @param include_burnin Whether to include burn-in draws when full chains were
#'   saved with \code{keep_burnin = TRUE}. The default \code{FALSE} computes the
#'   Rhat trajectory from post-burn-in draws.
#' @return A list with \code{max}, a data frame of maximum Rhat by prefix, and
#'   \code{item}, a data frame of scalar Rhat values for each flattened
#'   component.
#' @export
compute_Rhat_curve <- function(result,
                               parameters = c("s", "g"),
                               step = 1000L,
                               min_draws = step,
                               max_draws = NULL,
                               include_burnin = FALSE)
{
  stored_draws <- ncol(result[[1]]$chain_sample$s)
  data <- result[[1]]$data
  stored_is_full <- isTRUE(data$keep_burnin)
  burn_in <- if (!is.null(data$burn_in)) as.integer(data$burn_in) else 0L
  start <- 1L
  if (!isTRUE(include_burnin) && stored_is_full) {
    start <- burn_in + 1L
  }
  if (start > stored_draws) {
    stop("burn-in leaves no stored draws for Rhat trajectory.", call. = FALSE)
  }
  available_draws <- stored_draws - start + 1L
  if (is.null(max_draws)) {
    max_draws <- available_draws
  }
  max_draws <- min(as.integer(max_draws), available_draws)
  min_draws <- min(as.integer(min_draws), max_draws)
  step <- max(1L, as.integer(step))
  draws <- unique(c(seq(min_draws, max_draws, by = step), max_draws))
  draws <- draws[draws >= 2L & draws <= available_draws]
  end_idx <- start + draws - 1L

  max_table <- data.frame(
    stored_draws = draws,
    start_iteration = start,
    end_iteration = end_idx,
    total_iterations = if (stored_is_full) end_idx else burn_in + draws
  )
  item_tables <- vector("list", length(draws))

  for (i in seq_along(draws)) {
    n_draw <- draws[i]
    rhat <- compute_Rhat(result, post_idx = start:end_idx[i],
                         parameters = parameters)
    for (param in parameters) {
      max_table[[paste0(param, "_max_Rhat")]][i] <- rhat$max[[param]]
      max_table[[paste0(param, "_mpsrf")]][i] <- rhat$mpsrf[[param]]
    }

    scalar_rows <- lapply(parameters, function(param) {
      values <- as.numeric(rhat$values[[param]])
      data.frame(
        stored_draws = n_draw,
        start_iteration = start,
        end_iteration = end_idx[i],
        total_iterations = max_table$total_iterations[i],
        parameter = param,
        element = seq_along(values),
        Rhat = values
      )
    })
    item_tables[[i]] <- do.call(rbind, scalar_rows)
  }

  out <- list(max = max_table,
              item = do.call(rbind, item_tables),
              parameters = parameters,
              include_burnin = isTRUE(include_burnin) && stored_is_full,
              burn_in = if (stored_is_full) burn_in else NULL)
  class(out) <- "AHMQ_Rhat_curve"
  out
}

#' Plot Rhat trajectories
#'
#' @param curve Output from \code{\link{compute_Rhat_curve}}.
#' @param parameters Parameters to plot. \code{NULL} uses
#'   \code{curve$parameters}.
#' @param threshold Horizontal convergence reference line. The default 1.1 is
#'   the usual Rhat convergence threshold.
#' @param cutoff Deprecated alias for \code{threshold}.
#' @param file Optional file path. If supplied, a PNG file is written.
#' @param width,height,res PNG device settings used when \code{file} is not
#'   \code{NULL}.
#' @return Invisibly returns \code{curve}.
#' @export
plot_Rhat_curve <- function(curve,
                            parameters = NULL,
                            threshold = 1.1,
                            cutoff = NULL,
                            file = NULL,
                            width = 1200,
                            height = 800,
                            res = 130)
{
  if (!is.null(cutoff)) {
    threshold <- cutoff
  }
  if (is.null(parameters)) {
    parameters <- curve$parameters
  }
  if (!is.null(file)) {
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  cols <- grDevices::rainbow(length(parameters))
  y_values <- unlist(lapply(parameters, function(p) {
    curve$max[[paste0(p, "_max_Rhat")]]
  }), use.names = FALSE)
  y_max <- max(y_values, threshold, na.rm = TRUE)
  graphics::plot(
    NA,
    xlim = range(curve$max$total_iterations),
    ylim = c(1, y_max * 1.03),
    xlab = "Total MCMC iterations per chain",
    ylab = "Maximum scalar Rhat",
    main = "Rhat over chain length"
  )
  for (i in seq_along(parameters)) {
    p <- parameters[i]
    graphics::lines(curve$max$total_iterations,
                    curve$max[[paste0(p, "_max_Rhat")]],
                    col = cols[i], lwd = 2)
  }
  if (isTRUE(curve$include_burnin) && !is.null(curve$burn_in) &&
      is.finite(curve$burn_in) && curve$burn_in > min(curve$max$total_iterations) &&
      curve$burn_in < max(curve$max$total_iterations)) {
    graphics::abline(v = curve$burn_in, lty = 3, col = "gray50")
  }
  graphics::abline(h = threshold, lty = 2, col = "gray40")
  legend_text <- c(paste0(parameters, " max Rhat"), paste0("threshold = ", threshold))
  legend_col <- c(cols, "gray40")
  legend_lty <- c(rep(1, length(parameters)), 2)
  legend_lwd <- c(rep(2, length(parameters)), 1)
  if (isTRUE(curve$include_burnin) && !is.null(curve$burn_in) &&
      is.finite(curve$burn_in)) {
    legend_text <- c(legend_text, paste0("burn-in = ", curve$burn_in))
    legend_col <- c(legend_col, "gray50")
    legend_lty <- c(legend_lty, 3)
    legend_lwd <- c(legend_lwd, 1)
  }
  graphics::legend("topright",
                   legend = legend_text,
                   col = legend_col,
                   lty = legend_lty,
                   lwd = legend_lwd,
                   bty = "n")
  invisible(curve)
}

#' Plot an estimated hierarchy
#'
#' Draws the estimated attribute hierarchy using base R graphics. The function
#' can draw a directed graph, a posterior edge-probability heatmap, or both in
#' one figure. The hierarchy matrix is interpreted as
#' \code{G[prerequisite, target] = 1}. If an \code{\link{Est_fun}} or
#' \code{\link{extract_estimates}} result is supplied, \code{G} and
#' \code{G_posterior} are read from the object automatically.
#'
#' @param G Binary hierarchy adjacency matrix, posterior probability matrix,
#'   or an object returned by \code{\link{Est_fun}} or
#'   \code{\link{extract_estimates}}.
#' @param G_posterior Optional posterior edge-probability matrix. If omitted
#'   and \code{G} is an estimate object, the stored posterior mean of G is used.
#' @param type Plot type: \code{"graph"}, \code{"heatmap"}, or \code{"both"}.
#' @param cut_value Threshold used to binarize \code{G_posterior} for the graph
#'   when a binary \code{G} is not supplied.
#' @param labels Optional attribute labels.
#' @param file Optional PNG file path. If \code{NULL}, the active graphics
#'   device is used.
#' @param width,height,res PNG device settings used when \code{file} is not
#'   \code{NULL}.
#' @param main Plot title.
#' @param node_col Node fill color.
#' @param edge_col Arrow color.
#' @return Invisibly returns a data frame with node coordinates.
#' @export
plot_G_graph <- function(G,
                         G_posterior = NULL,
                         type = c("graph", "heatmap", "both"),
                         cut_value = 0.2,
                         labels = NULL,
                         file = NULL,
                         width = 1000,
                         height = 700,
                         res = 130,
                         main = "Estimated hierarchy G",
                         node_col = "white",
                         edge_col = "gray30")
{
  type <- match.arg(type)

  if (inherits(G, "simulation_summary")) {
    if (is.null(G_posterior)) {
      G_posterior <- G$aligned_estimates$Est_GG
    }
    if (is.null(cut_value) && !is.null(G$cut_value)) {
      cut_value <- G$cut_value
    }
    G <- G$aligned_estimates$Est_G
  } else if (is.list(G) && all(c("Est_G", "Est_GG") %in% names(G))) {
    if (is.null(G_posterior)) {
      G_posterior <- G$Est_GG
    }
    if (is.null(cut_value) && !is.null(G$cut_value)) {
      cut_value <- G$cut_value
    }
    G <- G$Est_G
  } else if (is.list(G) && all(c("est_G", "est_G_posterior") %in% names(G))) {
    if (is.null(G_posterior)) {
      G_posterior <- G$est_G_posterior
    }
    if (is.null(cut_value) && !is.null(G$cut_value)) {
      cut_value <- G$cut_value
    }
    G <- G$est_G
  }

  G_input <- as.matrix(G)
  K <- nrow(G_input)
  if (!identical(dim(G_input), c(K, K))) {
    stop("G must be a square matrix.", call. = FALSE)
  }
  if (is.null(G_posterior) && any(G_input > 0 & G_input < 1, na.rm = TRUE)) {
    G_posterior <- G_input
  }
  if (!is.null(G_posterior)) {
    G_posterior <- as.matrix(G_posterior)
    if (!identical(dim(G_posterior), c(K, K))) {
      stop("G_posterior must have the same dimensions as G.", call. = FALSE)
    }
  }

  G_binary <- G_input
  if (any(G_binary > 0 & G_binary < 1, na.rm = TRUE)) {
    G_binary <- (G_binary > cut_value) * 1
  }
  G_binary[G_binary != 0] <- 1

  if (is.null(labels)) {
    labels <- paste0("A", seq_len(K))
  }

  dag_layout <- function(G_binary) {
    if (K == 4L && identical(unname(G_binary), unname(simu_G("linear")))) {
      return(data.frame(x = c(0, 0, 0, 0), y = c(3, 2, 1, 0), depth = c(0, 1, 2, 3)))
    }
    if (K == 4L && identical(unname(G_binary), unname(simu_G("convergent")))) {
      return(data.frame(x = c(0, -1, 1, 0), y = c(2, 1, 1, 0), depth = c(0, 1, 1, 2)))
    }
    if (K == 4L && identical(unname(G_binary), unname(simu_G("divergent")))) {
      return(data.frame(x = c(0, -1, 0, 1), y = c(1, 0, 0, 0), depth = c(0, 1, 1, 1)))
    }
    if (K == 4L && identical(unname(G_binary), unname(simu_G("unstructured")))) {
      return(data.frame(x = c(0, -1, 1, 1), y = c(2, 1, 1, 0), depth = c(0, 1, 1, 2)))
    }

    roots <- which(colSums(G_binary) == 0)
    if (length(roots) == 0L) {
      roots <- seq_len(K)
    }
    depth <- rep(NA_integer_, K)
    depth[roots] <- 0L
    changed <- TRUE
    while (changed) {
      changed <- FALSE
      for (parent in seq_len(K)) {
        if (is.na(depth[parent])) next
        for (child in which(G_binary[parent, ] == 1)) {
          new_depth <- depth[parent] + 1L
          if (is.na(depth[child]) || new_depth > depth[child]) {
            depth[child] <- new_depth
            changed <- TRUE
          }
        }
      }
    }
    depth[is.na(depth)] <- 0L

    x <- numeric(K)
    y <- -depth
    for (d in sort(unique(depth))) {
      idx <- which(depth == d)
      x[idx] <- seq(-(length(idx) - 1) / 2, (length(idx) - 1) / 2, length.out = length(idx))
    }
    data.frame(x = x, y = y, depth = depth)
  }

  layout <- dag_layout(G_binary)
  x <- layout$x
  y <- layout$y
  depth <- layout$depth

  if (!is.null(file)) {
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  draw_graph <- function(plot_main) {
    graphics::par(mar = c(1, 1, 3, 1))
    graphics::plot(NA,
                   xlim = range(x) + c(-0.8, 0.8),
                   ylim = range(y) + c(-0.8, 0.8),
                   axes = FALSE,
                   xlab = "",
                   ylab = "",
                   main = plot_main,
                   asp = 1)

    edge_index <- which(G_binary == 1, arr.ind = TRUE)
    if (nrow(edge_index) > 0L) {
      for (i in seq_len(nrow(edge_index))) {
        parent <- edge_index[i, 1L]
        child <- edge_index[i, 2L]
        dx <- x[child] - x[parent]
        dy <- y[child] - y[parent]
        edge_len <- sqrt(dx^2 + dy^2)
        radius <- 0.18
        if (edge_len > 0) {
          graphics::arrows(x[parent] + radius * dx / edge_len,
                           y[parent] + radius * dy / edge_len,
                           x[child] - radius * dx / edge_len,
                           y[child] - radius * dy / edge_len,
                           length = 0.12,
                           lwd = 1.5,
                           col = edge_col)
        }
      }
    }

    radius <- 0.18
    graphics::symbols(x, y,
                      circles = rep(radius, K),
                      inches = FALSE,
                      add = TRUE,
                      bg = node_col,
                      fg = "black")
    graphics::text(x, y, labels = labels, cex = 0.9)
  }

  draw_heatmap <- function(plot_main) {
    if (is.null(G_posterior)) {
      heat <- G_binary
    } else {
      heat <- G_posterior
    }
    heat <- matrix(pmax(0, pmin(1, as.numeric(heat))), nrow = K, ncol = K)
    graphics::par(mar = c(4, 4, 3, 5))
    image_cols <- grDevices::colorRampPalette(c("white", "#2166ac"))(101)
    graphics::plot(NA,
                   xlim = c(0.5, K + 1.2),
                   ylim = c(0.5, K + 0.5),
                   axes = FALSE,
                   xlab = "Target attribute",
                   ylab = "Prerequisite attribute",
                   main = plot_main,
                   xaxs = "i",
                   yaxs = "i")
    graphics::axis(1, at = seq_len(K), labels = labels)
    graphics::axis(2, at = seq_len(K), labels = rev(labels), las = 1)
    for (parent in seq_len(K)) {
      y_pos <- K - parent + 1L
      for (child in seq_len(K)) {
        prob <- heat[parent, child]
        col_id <- max(1L, min(101L, floor(prob * 100) + 1L))
        graphics::rect(child - 0.5, y_pos - 0.5,
                       child + 0.5, y_pos + 0.5,
                       col = image_cols[col_id],
                       border = "gray80")
        graphics::text(child, y_pos, labels = sprintf("%.2f", prob),
                       cex = 0.75)
      }
    }
    if (!is.null(cut_value)) {
      graphics::mtext(paste0("cut_value = ", cut_value), side = 3,
                      line = 0.2, cex = 0.8)
    }
    graphics::box()
    legend_x0 <- K + 0.75
    legend_x1 <- K + 0.95
    legend_y <- seq(0.5, K + 0.5, length.out = 102)
    for (i in seq_len(101)) {
      graphics::rect(legend_x0, legend_y[i],
                     legend_x1, legend_y[i + 1L],
                     col = image_cols[i], border = NA)
    }
    graphics::rect(legend_x0, 0.5, legend_x1, K + 0.5, border = "gray50")
    graphics::axis(4, at = c(0.5, K + 0.5), labels = c("0", "1"), las = 1)
    graphics::mtext("Posterior probability", side = 4, line = 2.2, cex = 0.75)
  }

  if (type == "both") {
    graphics::par(mfrow = c(1, 2))
    draw_graph("Estimated G")
    draw_heatmap("Posterior mean of G")
  } else if (type == "heatmap") {
    draw_heatmap(main)
  } else {
    draw_graph(main)
  }

  invisible(data.frame(attribute = seq_len(K),
                       label = labels,
                       x = x,
                       y = y,
                       depth = depth))
}

#' Trace plots for AHMQ MCMC chains
#'
#' Draws trace plots for selected scalar parameters. By default, all item-level
#' \code{s} and \code{g} chains are plotted, with \code{s} and \code{g} kept
#' on separate pages. When the number of panels exceeds \code{panels_per_page},
#' multiple pages are produced automatically.
#'
#' @param result Output from \code{\link{AHMQ}}.
#' @param parameter Parameter block to plot. Supported values are \code{"s"}
#'   and \code{"g"}. The default \code{c("s", "g")} plots both, on separate
#'   pages.
#' @param items Item indices to plot. \code{NULL} plots all items.
#' @param chains Chain indices to include. \code{NULL} includes all chains.
#' @param post_idx Optional stored draw indices. \code{NULL} uses all stored
#'   draws selected by \code{include_burnin}. If supplied, \code{post_idx}
#'   overrides \code{include_burnin}.
#' @param include_burnin Whether to include burn-in draws when the result was
#'   created with \code{keep_burnin = TRUE}. If \code{TRUE}, a vertical dashed
#'   line marks the burn-in boundary.
#' @param discard_burnin Deprecated alias. Use \code{include_burnin = FALSE}.
#' @param panels_per_page Maximum number of item panels per page.
#' @param file Optional PNG file prefix. If supplied, files are written as
#'   \code{<prefix>_<parameter>_page<page>.png}. If \code{NULL}, plots are
#'   drawn to the active graphics device.
#' @param width,height,res PNG device settings used when \code{file} is not
#'   \code{NULL}.
#' @return Invisibly returns a data frame listing produced pages.
#' @export
traceplot_AHMQ <- function(result,
                           parameter = c("s", "g"),
                           items = NULL,
                           chains = NULL,
                           post_idx = NULL,
                           include_burnin = TRUE,
                           discard_burnin = FALSE,
                           panels_per_page = 12L,
                           file = NULL,
                           width = 1400,
                           height = 900,
                           res = 130)
{
  parameter <- match.arg(parameter, choices = c("s", "g"), several.ok = TRUE)
  n_chain <- length(result)
  if (is.null(chains)) {
    chains <- seq_len(n_chain)
  }
  chains <- as.integer(chains)
  if (any(chains < 1L | chains > n_chain)) {
    stop("chains contains invalid chain indices.", call. = FALSE)
  }

  stored_draws <- ncol(result[[1L]]$chain_sample$s)
  saved_full_chain <- isTRUE(result[[1L]]$data$keep_burnin)
  recorded_burn_in <- if (!is.null(result[[1L]]$data$burn_in)) {
    as.integer(result[[1L]]$data$burn_in)
  } else {
    0L
  }
  if (isTRUE(discard_burnin)) {
    include_burnin <- FALSE
  }
  if (is.null(post_idx)) {
    start <- 1L
    if (!isTRUE(include_burnin) && saved_full_chain) {
      start <- recorded_burn_in + 1L
    }
    post_idx <- start:stored_draws
  }
  post_idx <- as.integer(post_idx)
  if (any(post_idx < 1L | post_idx > stored_draws)) {
    stop("post_idx contains invalid stored draw indices.", call. = FALSE)
  }

  J <- nrow(result[[1L]]$chain_sample$s)
  if (is.null(items)) {
    items <- seq_len(J)
  }
  items <- as.integer(items)
  if (any(items < 1L | items > J)) {
    stop("items contains invalid item indices.", call. = FALSE)
  }

  panels_per_page <- max(1L, as.integer(panels_per_page))
  chain_cols <- grDevices::rainbow(length(chains))
  pages <- data.frame(parameter = character(),
                      page = integer(),
                      file = character(),
                      stringsAsFactors = FALSE)

  for (param in parameter) {
    item_pages <- split(items, ceiling(seq_along(items) / panels_per_page))
    for (page_id in seq_along(item_pages)) {
      page_items <- item_pages[[page_id]]
      out_file <- NA_character_
      if (!is.null(file)) {
        out_file <- paste0(file, "_", param, "_page", page_id, ".png")
        grDevices::png(out_file, width = width, height = height, res = res)
      }

      nr <- ceiling(sqrt(length(page_items)))
      nc <- ceiling(length(page_items) / nr)
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))

      for (item in page_items) {
        y_range <- range(unlist(lapply(chains, function(cl) {
          result[[cl]]$chain_sample[[param]][item, post_idx]
        })), na.rm = TRUE)
        graphics::plot(post_idx,
                       result[[chains[1L]]]$chain_sample[[param]][item, post_idx],
                       type = "l",
                       col = chain_cols[1L],
                       ylim = y_range,
                       xlab = "Stored draw",
                       ylab = param,
                       main = paste0(param, "[", item, "]"))
        if (saved_full_chain && isTRUE(include_burnin) &&
            recorded_burn_in > min(post_idx) &&
            recorded_burn_in < max(post_idx)) {
          graphics::abline(v = recorded_burn_in, lty = 2, col = "gray50")
        }
        if (length(chains) > 1L) {
          for (i in 2:length(chains)) {
            graphics::lines(post_idx,
                            result[[chains[i]]]$chain_sample[[param]][item, post_idx],
                            col = chain_cols[i])
          }
        }
      }
      graphics::mtext(paste0("Trace plots for ", param,
                             " (page ", page_id, "/", length(item_pages), ")"),
                      outer = TRUE, cex = 1.1)
      if (length(chains) <= 6L) {
        graphics::legend("topright",
                         legend = paste0("chain ", chains),
                         col = chain_cols,
                         lty = 1,
                         bty = "n",
                         cex = 0.7)
      }
      if (!is.null(file)) {
        grDevices::dev.off()
      }
      pages <- rbind(pages,
                     data.frame(parameter = param,
                                page = page_id,
                                file = out_file,
                                stringsAsFactors = FALSE))
    }
  }
  invisible(pages)
}




