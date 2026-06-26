# AHM

**Bayesian DINA estimation with unknown Q-matrix and attribute hierarchy**

[![GitHub](https://img.shields.io/badge/GitHub-K--S666%2FAHM-blue)](https://github.com/K-S666/AHM)

**AHM** is an R package for Bayesian MCMC estimation of the [DINA](https://en.wikipedia.org/wiki/DINA_model) (Deterministic Input, Noisy “And” gate) cognitive diagnosis model when both the **Q-matrix** and the **attribute hierarchy** are unknown. The sampler jointly estimates item slipping and guessing parameters, examinee attribute profiles, latent class probabilities, the Q-matrix, and a transitive-reduced attribute hierarchy.

A fixed-Q variant (`AHM()`) is also provided when the Q-matrix is known. The package includes post-processing for attribute label switching, Gelman–Rubin R̂ diagnostics, DIC-based chain selection, posterior summaries, simulation utilities, a bundled ECPE real-data example, and an optional Shiny analysis app.

## Installation

### Prerequisites

- R (≥ 2.10)
- A C++ compiler supported by **Rcpp** (on Windows, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required to build from source)
- System dependencies for compiling **RcppArmadillo**

### Install from GitHub

The recommended way to install the development version is with [**remotes**](https://cran.r-project.org/package=remotes) or [**devtools**](https://cran.r-project.org/package=devtools):

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("K-S666/AHM")
```

Alternatively:

```r
devtools::install_github("K-S666/AHM")
```

**Rcpp**, **RcppArmadillo**, and **rgen** (compile-time dependency) are installed automatically when missing.

### Optional: Shiny app

The interactive workbench requires the **shiny** package:

```r
install.packages("shiny")
```

## Quick start

### Unknown Q-matrix (`AHMQ`)

```r
library(AHM)

# Bundled ECPE real-data example
ecpe <- ecpe_ahmq_data()
fit <- AHMQ(
  Y = ecpe$Y,
  K = ecpe$K,
  N1 = 128,
  chain_length = 20000,
  burn_in = 10000,
  chain_num = 4
)

# Post-process: Rhat diagnostics, DIC chain selection, label switching
est <- Est_fun(fit, cut_value = 0.2)
summary_AHMQ(est)
```

### Known Q-matrix (`AHM`)

When the Q-matrix is supplied and fixed, use `AHM()` instead:

```r
dat <- simulate_ahmq_data(N = 100, J = 12, K = 3, seed = 1)

fit <- AHM(
  Y = dat$Y,
  Q = dat$Q,
  chain_length = 20000,
  burn_in = 10000,
  chain_num = 4
)

est <- Est_fun(fit)
summary_AHMQ(est)
```

### Simulate data

```r
dat <- simulate_ahmq_data(N = 200, J = 15, K = 4, seed = 42)
head(dat$Y)
dat$Q   # true Q-matrix
dat$G   # true attribute hierarchy
```

## Interactive Shiny app

Launch the bundled analysis workbench for simulated data, the ECPE example, or user-uploaded CSV files:

```r
library(AHM)
run_AHM_app()
```

## Main functions

| Function | Description |
|----------|-------------|
| `AHMQ()` | Multi-chain MCMC when Q-matrix and hierarchy are unknown |
| `AHM()` | Multi-chain MCMC with a known, fixed Q-matrix |
| `Est_fun()` | Post-process chains: R̂, DIC selection, label switching, point estimates |
| `summary_AHMQ()` | Print a concise summary of `Est_fun()` output |
| `compute_DIC()` / `select_chain_by_DIC()` | DIC computation and best-chain selection |
| `compute_Rhat()` / `plot_Rhat_curve()` | Gelman–Rubin convergence diagnostics |
| `resolve_label_switch()` / `align_estimates_to_truth()` | Attribute label-switching utilities |
| `simulate_ahmq_data()` | Simulate DINA data under AHM priors |
| `ecpe_ahmq_data()` | Prepare bundled ECPE data for `AHMQ()` |
| `run_AHM_app()` | Launch the Shiny analysis app |
| `traceplot_AHMQ()` / `plot_G_graph()` | Diagnostic and hierarchy visualization |

See `?AHMQ`, `?AHM`, and `?Est_fun` for full argument lists and details.

## Package structure

- **`R/`** — R wrappers, post-processing, simulation, and the Shiny launcher
- **`src/`** — C++ MCMC samplers (Rcpp / RcppArmadillo)
- **`inst/shiny/`** — Interactive Shiny application
- **`man/`** — Function documentation

## License

GPL (≥ 2)

## Author

**Xue Wang** — [wangx625@nenu.edu.cn](mailto:wangx625@nenu.edu.cn)

## Repository

https://github.com/K-S666/AHM
