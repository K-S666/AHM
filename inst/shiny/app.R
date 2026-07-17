if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The cdmArch Shiny app requires the 'shiny' package.")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_matrix_text <- function(text, nrow = NULL, ncol = NULL, name = "matrix")
{
  text <- trimws(text %||% "")
  if (!nzchar(text)) {
    return(NULL)
  }
  rows <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  rows <- rows[nzchar(trimws(rows))]
  vals <- lapply(rows, function(row) {
    as.numeric(strsplit(trimws(row), "[,;[:space:]]+")[[1L]])
  })
  lens <- vapply(vals, length, integer(1L))
  if (length(unique(lens)) != 1L || any(!is.finite(unlist(vals)))) {
    stop(name, " must be a rectangular numeric matrix.", call. = FALSE)
  }
  out <- do.call(rbind, vals)
  storage.mode(out) <- "double"
  if (!is.null(nrow) && nrow(out) != nrow) {
    stop(name, " must have ", nrow, " rows.", call. = FALSE)
  }
  if (!is.null(ncol) && ncol(out) != ncol) {
    stop(name, " must have ", ncol, " columns.", call. = FALSE)
  }
  out
}

read_uploaded_matrix <- function(file, header = TRUE, drop_first_col = FALSE,
                                 name = "uploaded matrix")
{
  if (is.null(file)) {
    return(NULL)
  }
  dat <- utils::read.csv(file$datapath, header = header, check.names = FALSE)
  if (isTRUE(drop_first_col)) {
    dat <- dat[, -1L, drop = FALSE]
  }
  out <- as.matrix(dat)
  storage.mode(out) <- "double"
  if (any(!is.finite(out))) {
    stop(name, " must contain only finite numeric values.", call. = FALSE)
  }
  out
}

make_G_preset <- function(type, K)
{
  cdmArch::simu_G(type = type, K = K)
}

binary_check <- function(x, name)
{
  if (!all(x %in% c(0, 1))) {
    stop(name, " must contain only 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

matrix_to_df <- function(x, row_prefix = "row")
{
  x <- as.matrix(x)
  if (all(is.na(x) | x %in% c(0, 1))) {
    storage.mode(x) <- "integer"
  }
  out <- as.data.frame(x)
  names(out) <- paste0("A", seq_len(ncol(out)))
  out <- cbind(id = paste0(row_prefix, seq_len(nrow(out))), out)
  row.names(out) <- NULL
  out
}

ui <- shiny::navbarPage(
  title = "cdmArch Workbench",
  id = "main_nav",
  shiny::tabPanel(
    "Data",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::radioButtons(
          "data_mode", "Data source",
          choices = c("Simulated data" = "simulation",
                      "Real data" = "real"),
          selected = "simulation"
        ),
        shiny::conditionalPanel(
          "input.data_mode == 'simulation'",
          shiny::h4("Simulation conditions"),
          shiny::numericInput("sim_N", "Sample size N", 500, min = 1, step = 10),
          shiny::numericInput("sim_J", "Number of items J", 20, min = 2, step = 1),
          shiny::helpText("Number of attributes K is fixed at 4 for simulation examples."),
          shiny::selectInput(
            "sim_G_type", "True hierarchy G",
            choices = c("No hierarchy" = "none",
                        "G1 Linear" = "linear",
                        "G2 Convergent" = "convergent",
                        "G3 Divergent" = "divergent",
                        "G4 Partially Structured" = "partially_structured",
                        "Custom matrix" = "custom"),
            selected = "linear"
          ),
          shiny::conditionalPanel(
            "input.sim_G_type == 'custom'",
            shiny::textAreaInput(
              "sim_G_text", "Custom G adjacency matrix",
              value = "0 1 0 0\n0 0 1 0\n0 0 0 1\n0 0 0 0",
              rows = 5
            )
          ),
          shiny::selectInput(
            "alpha_distribution", "Attribute profile distribution",
            choices = c("Balanced" = "balanced",
                        "Unbalanced" = "unbalanced"),
            selected = "balanced"
          ),
          shiny::conditionalPanel(
            "input.alpha_distribution == 'unbalanced'",
            shiny::numericInput("sigma", "Sigma for unbalanced profiles", 0.5,
                                min = -0.9, max = 0.99, step = 0.05)
          ),
          shiny::radioButtons(
            "sg_mode", "Slip/guess generation",
            choices = c("Fixed s_j = g_j" = "fixed",
                        "Uniform ranges" = "range"),
            selected = "fixed"
          ),
          shiny::conditionalPanel(
            "input.sg_mode == 'fixed'",
            shiny::numericInput("sim_s_fixed", "Fixed slip s_j", 0.2,
                                min = 0, max = 0.99, step = 0.01),
            shiny::numericInput("sim_g_fixed", "Fixed guess g_j", 0.2,
                                min = 0, max = 0.99, step = 0.01)
          ),
          shiny::conditionalPanel(
            "input.sg_mode == 'range'",
            shiny::numericInput("sim_s_min", "Slip lower", 0.05,
                                min = 0, max = 0.99, step = 0.01),
            shiny::numericInput("sim_s_max", "Slip upper", 0.25,
                                min = 0, max = 0.99, step = 0.01),
            shiny::numericInput("sim_g_min", "Guess lower", 0.05,
                                min = 0, max = 0.99, step = 0.01),
            shiny::numericInput("sim_g_max", "Guess upper", 0.25,
                                min = 0, max = 0.99, step = 0.01)
          ),
          shiny::checkboxInput("sim_q_known", "Treat true Q as known and run AHM()", FALSE),
          shiny::numericInput("sim_seed", "Simulation seed", 1, min = 1, step = 1),
          shiny::actionButton("prepare_sim", "Generate simulated data")
        ),
        shiny::conditionalPanel(
          "input.data_mode == 'real'",
          shiny::h4("Real data"),
          shiny::radioButtons(
            "real_source", "Real-data source",
            choices = c("Bundled ECPE data" = "ecpe",
                        "Upload CSV files" = "upload"),
            selected = "ecpe"
          ),
          shiny::conditionalPanel(
            "input.real_source == 'ecpe'",
            shiny::checkboxInput("ecpe_use_q", "Use bundled ECPE Q as known", TRUE)
          ),
          shiny::conditionalPanel(
            "input.real_source == 'upload'",
            shiny::fileInput("upload_Y", "Upload response matrix Y (.csv)"),
            shiny::checkboxInput("upload_Y_header", "Y file has header", TRUE),
            shiny::checkboxInput("upload_Y_drop_first", "Drop first Y column as ID", FALSE),
            shiny::fileInput("upload_Q", "Optional known Q-matrix (.csv)"),
            shiny::checkboxInput("upload_Q_header", "Q file has header", TRUE),
            shiny::numericInput("upload_K", "K if Q is not uploaded", 4,
                                min = 1, step = 1)
          ),
          shiny::actionButton("prepare_real", "Load real data")
        )
      ),
      shiny::mainPanel(
        shiny::h4("Prepared data"),
        shiny::verbatimTextOutput("data_status"),
        shiny::plotOutput("item_prop_plot", height = 260),
        shiny::h4("Y preview"),
        shiny::tableOutput("Y_preview"),
        shiny::h4("Q preview"),
        shiny::tableOutput("Q_preview"),
        shiny::h4("True G DAG"),
        shiny::plotOutput("true_G_plot", height = 320)
      )
    )
  ),
  shiny::tabPanel(
    "MCMC",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::h4("Model settings"),
        shiny::checkboxInput("use_default_N1", "Use default minibatch: floor(N/2)", TRUE),
        shiny::conditionalPanel(
          "!input.use_default_N1",
          shiny::numericInput("N1", "minibatch N1", 256, min = 1, step = 1)
        ),
        shiny::checkboxInput("use_default_mcmc", "Use default MCMC length", TRUE),
        shiny::conditionalPanel(
          "input.use_default_mcmc",
          shiny::helpText("Default: chain_length = 20000, burn_in = 10000")
        ),
        shiny::conditionalPanel(
          "!input.use_default_mcmc",
          shiny::numericInput("chain_length", "chain_length", 20000, min = 2, step = 1000),
          shiny::numericInput("burn_in", "burn_in", 10000, min = 0, step = 1000)
        ),
        shiny::numericInput("chain_num", "Number of chains", 4, min = 1, step = 1),
        shiny::checkboxInput("parallel", "Run chains in parallel", TRUE),
        shiny::numericInput("n_cores", "Number of cores", 4, min = 1, step = 1),
        shiny::checkboxInput("keep_burnin", "Keep burn-in draws for diagnostics", FALSE),
        shiny::numericInput("cut_value", "G posterior cut_value", 0.2,
                            min = 0, max = 1, step = 0.05),
        shiny::checkboxInput("compute_rhat", "Compute Rhat diagnostics", TRUE),
        shiny::checkboxInput("return_samples", "Keep posterior samples for simulation label alignment", TRUE),
        shiny::h4("Priors and proposals"),
        shiny::numericInput("a_s0", "a_s0", 1, min = 0.001, step = 0.1),
        shiny::numericInput("b_s0", "b_s0", 1, min = 0.001, step = 0.1),
        shiny::numericInput("a_g0", "a_g0", 1, min = 0.001, step = 0.1),
        shiny::numericInput("b_g0", "b_g0", 1, min = 0.001, step = 0.1),
        shiny::numericInput("p_add", "p_add (add-edge proposal prob.)", 0.5,
                            min = 0.01, max = 0.99, step = 0.05),
        shiny::checkboxInput("show_progress", "Show per-chain iteration progress", TRUE),
        shiny::numericInput("print_every", "Progress print_every", 1000, min = 1, step = 100),
        shiny::actionButton("run_mcmc", "Run MCMC", class = "btn-primary")
      ),
      shiny::mainPanel(
        shiny::h4("Run status"),
        shiny::verbatimTextOutput("run_status"),
        shiny::h4("DIC by chain"),
        shiny::tableOutput("DIC_table"),
        shiny::h4("Rhat summary"),
        shiny::tableOutput("Rhat_table")
      )
    )
  ),
  shiny::tabPanel(
    "Estimates",
    shiny::h4("Posterior summary"),
    shiny::verbatimTextOutput("summary_text"),
    shiny::h4("Item parameters"),
    shiny::tableOutput("item_param_table"),
    shiny::h4("G posterior mean"),
    shiny::tableOutput("G_post_table"),
    shiny::h4("Estimated G"),
    shiny::tableOutput("G_table"),
    shiny::h4("Q matrix"),
    shiny::tableOutput("Q_table"),
    shiny::h4("First estimated alpha profiles"),
    shiny::tableOutput("alpha_table")
  ),
  shiny::tabPanel(
    "Diagnostics",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("G_plot_type", "G plot type",
                           choices = c("Graph" = "graph",
                                       "Posterior heatmap" = "heatmap",
                                       "Both" = "both"),
                           selected = "both"),
        shiny::selectInput("trace_param", "Trace parameter",
                           choices = c("s", "g"), selected = "s"),
        shiny::textInput("trace_items", "Trace item indices", "1:6"),
        shiny::numericInput("trace_panels", "Panels per page", 6, min = 1, step = 1),
        shiny::numericInput("rhat_step", "Rhat curve step", 1000, min = 2, step = 100),
        shiny::checkboxInput("rhat_include_burnin", "Rhat curve includes burn-in if stored", FALSE)
      ),
      shiny::mainPanel(
        shiny::h4("G visualization"),
        shiny::plotOutput("G_plot", height = 420),
        shiny::h4("Trace plot"),
        shiny::plotOutput("trace_plot", height = 620),
        shiny::h4("Rhat over chain length"),
        shiny::plotOutput("rhat_curve_plot", height = 420)
      )
    )
  ),
  shiny::tabPanel(
    "Simulation Evaluation",
    shiny::h4("Simulation-only recovery metrics"),
    shiny::verbatimTextOutput("sim_summary_text"),
    shiny::h4("Bias for s and g"),
    shiny::tableOutput("sim_bias_table"),
    shiny::h4("G recovery"),
    shiny::tableOutput("sim_G_table"),
    shiny::h4("Q recovery"),
    shiny::tableOutput("sim_Q_table"),
    shiny::h4("Truth-aligned estimated G DAG"),
    shiny::plotOutput("sim_G_plot", height = 420)
  ),
  shiny::tabPanel(
    "Export",
    shiny::p("Download fitted objects and key tables."),
    shiny::downloadButton("download_fit", "Download fit.rds"),
    shiny::downloadButton("download_est", "Download estimates.rds"),
    shiny::downloadButton("download_item", "Download item parameters CSV"),
    shiny::downloadButton("download_alpha", "Download alpha estimates CSV"),
    shiny::downloadButton("download_G", "Download G estimate CSV"),
    shiny::downloadButton("download_Q", "Download Q matrix CSV")
  )
)

server <- function(input, output, session)
{
  rv <- shiny::reactiveValues(
    data = NULL,
    truth = NULL,
    fit = NULL,
    est = NULL,
    summary = NULL,
    sim_summary = NULL,
    runtime = NULL,
    status = "No data loaded."
  )

  prepared_model_name <- shiny::reactive({
    shiny::req(rv$data)
    if (isTRUE(rv$data$q_known)) "AHM (known Q)" else "AHMQ (unknown Q)"
  })

  shiny::observeEvent(input$prepare_sim, {
    tryCatch({
      K <- 4L
      G <- if (input$sim_G_type == "custom") {
        parse_matrix_text(input$sim_G_text, nrow = K, ncol = K, name = "Custom G")
      } else {
        make_G_preset(input$sim_G_type, K)
      }
      binary_check(G, "G")
      if (input$sg_mode == "fixed") {
        s <- input$sim_s_fixed
        g <- input$sim_g_fixed
        s_range <- c(0.05, 0.25)
        g_range <- c(0.05, 0.25)
      } else {
        s <- NULL
        g <- NULL
        s_range <- c(input$sim_s_min, input$sim_s_max)
        g_range <- c(input$sim_g_min, input$sim_g_max)
      }
      dat <- cdmArch::simulate_ahmq_data(
        N = input$sim_N,
        J = input$sim_J,
        K = K,
        G = G,
        s = s,
        g = g,
        alpha_distribution = input$alpha_distribution,
        sigma = input$sigma,
        s_range = s_range,
        g_range = g_range,
        seed = input$sim_seed
      )
      rv$data <- list(Y = dat$Y, Q = dat$Q, K = dat$K,
                      q_known = isTRUE(input$sim_q_known),
                      source = "simulation")
      rv$truth <- dat
      rv$fit <- rv$est <- rv$summary <- rv$sim_summary <- NULL
      rv$status <- "Simulated data generated."
    }, error = function(e) {
      shiny::showNotification(conditionMessage(e), type = "error", duration = 8)
    })
  })

  shiny::observeEvent(input$prepare_real, {
    tryCatch({
      if (input$real_source == "ecpe") {
        dat <- cdmArch::ecpe_ahmq_data()
        rv$data <- list(Y = dat$Y, Q = dat$Q, K = dat$K,
                        q_known = isTRUE(input$ecpe_use_q),
                        source = "real_ecpe")
      } else {
        Y <- read_uploaded_matrix(input$upload_Y,
                                  header = input$upload_Y_header,
                                  drop_first_col = input$upload_Y_drop_first,
                                  name = "Y")
        if (is.null(Y)) {
          stop("Please upload a response matrix Y.", call. = FALSE)
        }
        binary_check(Y, "Y")
        Q <- read_uploaded_matrix(input$upload_Q,
                                  header = input$upload_Q_header,
                                  drop_first_col = FALSE,
                                  name = "Q")
        if (!is.null(Q)) {
          binary_check(Q, "Q")
          if (nrow(Q) != ncol(Y)) {
            stop("nrow(Q) must equal ncol(Y).", call. = FALSE)
          }
          K <- ncol(Q)
          q_known <- TRUE
        } else {
          K <- as.integer(input$upload_K)
          q_known <- FALSE
        }
        rv$data <- list(Y = Y, Q = Q, K = K,
                        q_known = q_known,
                        source = "real_upload")
      }
      rv$truth <- NULL
      rv$fit <- rv$est <- rv$summary <- rv$sim_summary <- NULL
      rv$status <- "Real data loaded."
    }, error = function(e) {
      shiny::showNotification(conditionMessage(e), type = "error", duration = 8)
    })
  })

  output$data_status <- shiny::renderPrint({
    if (is.null(rv$data)) {
      cat(rv$status, "\n")
      return(invisible(NULL))
    }
    Y <- rv$data$Y
    cat(rv$status, "\n")
    cat("Source:", rv$data$source, "\n")
    cat("Model:", prepared_model_name(), "\n")
    cat("N:", nrow(Y), " J:", ncol(Y), " K:", rv$data$K, "\n")
    cat("Missing values:", sum(is.na(Y)), "\n")
    cat("Q available:", !is.null(rv$data$Q), "\n")
  })

  output$Y_preview <- shiny::renderTable({
    shiny::req(rv$data)
    Y <- rv$data$Y
    storage.mode(Y) <- "integer"
    utils::head(as.data.frame(Y), 8)
  }, rownames = TRUE)

  output$Q_preview <- shiny::renderTable({
    shiny::req(rv$data)
    if (is.null(rv$data$Q)) return(data.frame(message = "No Q supplied. AHMQ will estimate Q."))
    utils::head(matrix_to_df(rv$data$Q, "item"), 12)
  })

  output$item_prop_plot <- shiny::renderPlot({
    shiny::req(rv$data)
    p <- colMeans(rv$data$Y, na.rm = TRUE)
    graphics::barplot(p, ylim = c(0, 1), col = "gray70",
                      xlab = "Item", ylab = "Proportion correct",
                      main = "Item correct proportions")
  })
  output$true_G_plot <- shiny::renderPlot({
    shiny::req(rv$truth)
    cdmArch::plot_G_graph(rv$truth$G, type = "graph", labels = seq_len(rv$truth$K),
                      main = "True hierarchy G")
  })

  shiny::observeEvent(input$run_mcmc, {
    shiny::req(rv$data)
    tryCatch({
      Y <- rv$data$Y
      Q <- rv$data$Q
      K <- rv$data$K
      N1 <- if (isTRUE(input$use_default_N1)) NULL else as.integer(input$N1)
      chain_length <- if (isTRUE(input$use_default_mcmc)) 20000L else as.integer(input$chain_length)
      burn_in <- if (isTRUE(input$use_default_mcmc)) 10000L else as.integer(input$burn_in)
      if (burn_in >= chain_length) {
        stop("burn_in must be smaller than chain_length.", call. = FALSE)
      }
      keep_samples <- isTRUE(input$return_samples) && rv$data$source == "simulation"
      with_shiny_progress <- shiny::Progress$new(session, min = 0, max = 3)
      on.exit(with_shiny_progress$close(), add = TRUE)
      with_shiny_progress$set(message = "Running MCMC", value = 0)
      shiny_progress_callback <- NULL
      if (isTRUE(input$show_progress)) {
        shiny_progress_callback <- function(current, total, model, final = FALSE) {
          pct <- pmin(100, 100 * current / total)
          detail <- paste(
            sprintf("chain %d: %d/%d, %.1f%%", seq_along(current), current, total, pct),
            collapse = "\n"
          )
          with_shiny_progress$set(
            message = paste(model, "chain progress", if (final) "(final)" else ""),
            detail = detail,
            value = mean(pct) / 100
          )
          try(session$flushReact(), silent = TRUE)
        }
      }
      t0 <- proc.time()
      if (isTRUE(rv$data$q_known)) {
        if (is.null(Q)) {
          stop("Known-Q model selected but Q is not available.", call. = FALSE)
        }
        fit <- cdmArch::AHM(
          Y, Q,
          N1 = N1,
          chain_length = chain_length,
          burn_in = burn_in,
          keep_burnin = isTRUE(input$keep_burnin),
          chain_num = as.integer(input$chain_num),
          a_s0 = input$a_s0, b_s0 = input$b_s0,
          a_g0 = input$a_g0, b_g0 = input$b_g0,
          p_add = input$p_add,
          parallel = isTRUE(input$parallel),
          n_cores = as.integer(input$n_cores),
          progress = isTRUE(input$show_progress),
          print_every = as.integer(input$print_every),
          progress_callback = shiny_progress_callback
        )
      } else {
        fit <- cdmArch::AHMQ(
          Y, K,
          N1 = N1,
          chain_length = chain_length,
          burn_in = burn_in,
          keep_burnin = isTRUE(input$keep_burnin),
          chain_num = as.integer(input$chain_num),
          a_s0 = input$a_s0, b_s0 = input$b_s0,
          a_g0 = input$a_g0, b_g0 = input$b_g0,
          p_add = input$p_add,
          parallel = isTRUE(input$parallel),
          n_cores = as.integer(input$n_cores),
          progress = isTRUE(input$show_progress),
          print_every = as.integer(input$print_every),
          progress_callback = shiny_progress_callback
        )
      }
      runtime <- proc.time() - t0
      with_shiny_progress$set(message = "Post-processing estimates", value = 1)
      est <- cdmArch::Est_fun(
        fit,
        cut_value = input$cut_value,
        compute_rhat = isTRUE(input$compute_rhat),
        return_samples = keep_samples,
        verbose = FALSE
      )
      est$runtime <- as.difftime(as.numeric(runtime["elapsed"]), units = "secs")
      summ <- cdmArch::summary_est(est)
      sim_summ <- NULL
      if (rv$data$source == "simulation") {
        with_shiny_progress$set(message = "Computing simulation recovery metrics", value = 2)
        sim_summ <- cdmArch::simu_result_summary(
          est,
          truth = rv$truth,
          align_labels = !isTRUE(est$known_Q) && keep_samples,
          verbose = FALSE
        )
      }
      rv$fit <- fit
      rv$est <- est
      rv$summary <- summ
      rv$sim_summary <- sim_summ
      rv$runtime <- runtime
      rv$status <- "MCMC finished."
      with_shiny_progress$set(message = "Done", value = 3)
    }, error = function(e) {
      shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })
  output$run_status <- shiny::renderPrint({
    cat(rv$status, "\n")
    if (!is.null(rv$fit)) {
      dat <- rv$fit[[1L]]$data
      cat("Model:", if (isTRUE(dat$known_Q)) "AHM" else "AHMQ", "\n")
      cat("Chains:", dat$chain_num, "\n")
      cat("N1:", dat$N1, if (isTRUE(dat$N1_default)) "(default)" else "", "\n")
      cat("chain_length:", dat$chain_length, " burn_in:", dat$burn_in, "\n")
      cat("parallel:", dat$parallel, " n_cores:", dat$n_cores, "\n")
      if (!is.null(rv$runtime)) {
        cat("Runtime elapsed:", round(as.numeric(rv$runtime["elapsed"]), 2), "seconds\n")
      }
    }
  })

  output$DIC_table <- shiny::renderTable({
    shiny::req(rv$est)
    data.frame(chain = seq_along(rv$est$DIC_all),
               DIC = round(rv$est$DIC_all, 3),
               best = seq_along(rv$est$DIC_all) == rv$est$best_chain)
  })

  output$Rhat_table <- shiny::renderTable({
    shiny::req(rv$est)
    if (is.null(rv$est$Rhat)) return(data.frame(message = "Rhat was not computed."))
    data.frame(parameter = names(rv$est$Rhat$max),
               max_Rhat = round(as.numeric(rv$est$Rhat$max), 3),
               mpsrf = round(as.numeric(rv$est$Rhat$mpsrf), 3))
  })

  output$summary_text <- shiny::renderPrint({
    shiny::req(rv$summary)
    print(rv$summary)
  })

  output$item_param_table <- shiny::renderTable({
    shiny::req(rv$est)
    data.frame(item = seq_along(rv$est$Est_s),
               slip = round(rv$est$Est_s, 4),
               guess = round(rv$est$Est_g, 4))
  })

  output$G_post_table <- shiny::renderTable({
    shiny::req(rv$est)
    round(rv$est$Est_GG, 4)
  }, rownames = TRUE)

  output$G_table <- shiny::renderTable({
    shiny::req(rv$est)
    G <- rv$est$Est_G
    storage.mode(G) <- "integer"
    G
  }, rownames = TRUE)

  output$Q_table <- shiny::renderTable({
    shiny::req(rv$est)
    Q <- if (isTRUE(rv$est$known_Q)) rv$est$Q_fixed else rv$est$Est_Q
    matrix_to_df(Q, "item")
  })

  output$alpha_table <- shiny::renderTable({
    shiny::req(rv$est)
    alpha <- rv$est$Est_alpha
    storage.mode(alpha) <- "integer"
    utils::head(matrix_to_df(alpha, "person"), 20)
  })

  parse_indices <- shiny::reactive({
    txt <- trimws(input$trace_items %||% "")
    if (!nzchar(txt)) return(NULL)
    eval(parse(text = txt), envir = baseenv())
  })

  output$G_plot <- shiny::renderPlot({
    shiny::req(rv$est)
    cdmArch::plot_G_graph(rv$est, type = input$G_plot_type,
                      cut_value = input$cut_value)
  })

  output$trace_plot <- shiny::renderPlot({
    shiny::req(rv$fit)
    cdmArch::traceplot_AHMQ(
      rv$fit,
      parameter = input$trace_param,
      items = parse_indices(),
      panels_per_page = input$trace_panels,
      file = NULL
    )
  })

  output$rhat_curve_plot <- shiny::renderPlot({
    shiny::req(rv$fit)
    if (length(rv$fit) < 2L) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Rhat curve requires at least two chains.")
      return(invisible(NULL))
    }
    params <- c("s", "g")
    curve <- cdmArch::compute_Rhat_curve(
      rv$fit,
      parameters = params,
      step = as.integer(input$rhat_step),
      include_burnin = isTRUE(input$rhat_include_burnin)
    )
    cdmArch::plot_Rhat_curve(curve, threshold = 1.1)
  })

  output$sim_summary_text <- shiny::renderPrint({
    if (is.null(rv$sim_summary)) {
      cat("Simulation summary is available only for simulated data after MCMC.\n")
    } else {
      print(rv$sim_summary)
    }
  })

  output$sim_bias_table <- shiny::renderTable({
    shiny::req(rv$sim_summary)
    out <- rv$sim_summary$s_g_metrics[, c("parameter", "mean_bias", "mean_abs_bias")]
    num_cols <- vapply(out, is.numeric, logical(1L))
    out[num_cols] <- lapply(out[num_cols], round, digits = 3)
    out
  })

  output$sim_G_table <- shiny::renderTable({
    shiny::req(rv$sim_summary)
    out <- as.data.frame(rv$sim_summary$G_metrics)
    num_cols <- vapply(out, is.numeric, logical(1L))
    out[num_cols] <- lapply(out[num_cols], round, digits = 3)
    out
  })

  output$sim_Q_table <- shiny::renderTable({
    shiny::req(rv$sim_summary)
    if (is.null(rv$sim_summary$Q_metrics)) {
      return(data.frame(message = "Q is fixed/known; Q recovery is not evaluated."))
    }
    out <- data.frame(t(rv$sim_summary$Q_metrics$QRR),
               AQRR = rv$sim_summary$Q_metrics$AQRR,
               Q_exact = rv$sim_summary$Q_metrics$Q_exact,
               check.names = FALSE)
    num_cols <- vapply(out, is.numeric, logical(1L))
    out[num_cols] <- lapply(out[num_cols], round, digits = 3)
    out
  })
  output$sim_G_plot <- shiny::renderPlot({
    shiny::req(rv$sim_summary)
    cdmArch::plot_G_graph(rv$sim_summary, type = "graph",
                      cut_value = input$cut_value,
                      labels = seq_len(nrow(rv$sim_summary$aligned_estimates$Est_G)),
                      main = "Truth-aligned estimated hierarchy G")
  })

  output$download_fit <- shiny::downloadHandler(
    filename = function() "cdmArch_fit.rds",
    content = function(file) saveRDS(rv$fit, file)
  )
  output$download_est <- shiny::downloadHandler(
    filename = function() "cdmArch_estimates.rds",
    content = function(file) saveRDS(rv$est, file)
  )
  output$download_item <- shiny::downloadHandler(
    filename = function() "cdmArch_item_parameters.csv",
    content = function(file) {
      shiny::req(rv$est)
      utils::write.csv(data.frame(item = seq_along(rv$est$Est_s),
                                  slip = rv$est$Est_s,
                                  guess = rv$est$Est_g),
                       file, row.names = FALSE)
    }
  )
  output$download_alpha <- shiny::downloadHandler(
    filename = function() "cdmArch_alpha_estimates.csv",
    content = function(file) {
      shiny::req(rv$est)
      alpha <- rv$est$Est_alpha
      storage.mode(alpha) <- "integer"
      utils::write.csv(alpha, file, row.names = FALSE)
    }
  )
  output$download_G <- shiny::downloadHandler(
    filename = function() "cdmArch_G_estimate.csv",
    content = function(file) {
      shiny::req(rv$est)
      G <- rv$est$Est_G
      storage.mode(G) <- "integer"
      utils::write.csv(G, file, row.names = FALSE)
    }
  )
  output$download_Q <- shiny::downloadHandler(
    filename = function() "cdmArch_Q_matrix.csv",
    content = function(file) {
      shiny::req(rv$est)
      Q <- if (isTRUE(rv$est$known_Q)) rv$est$Q_fixed else rv$est$Est_Q
      storage.mode(Q) <- "integer"
      utils::write.csv(Q, file, row.names = FALSE)
    }
  )
}

shiny::shinyApp(ui, server)
