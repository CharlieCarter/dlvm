#!/usr/bin/env Rscript
# =============================================================================
# sim_compare_obs_models.R
#
# Simulation comparison of four observation noise models for DLVM means mode:
#   Baseline : dlvm.stan         — global sigma, SE=0 for n=1 (original bug)
#   Option A : dlvm_obs_n.stan   — sigma/sqrt(n), zero new params
#   Option B : dlvm_hier_sigma   — hierarchical per-unit sigma
#   Option C : dlvm_hier_sigma_n — hierarchical sigma + sqrt(n)
#
# Usage:
#   Rscript tests/sim_compare_obs_models.R           # run all scenarios
#   Rscript tests/sim_compare_obs_models.R quick      # 1 rep, 2 scenarios
#
# Designed to run from the dlvm/ directory.
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(data.table)
})

# Prevent TBB conflict: CmdStan models compiled with TBB can corrupt
# data.table's internal threading on macOS. Stan sampling runs in a
# separate process so this doesn't affect Stan's performance.
setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
QUICK <- "quick" %in% args

cat("=== DLVM Observation Noise Model Comparison ===\n")
cat(sprintf("Mode: %s\n\n", if (QUICK) "QUICK (reduced)" else "FULL"))

# =============================================================================
# 1. Model paths & compilation
# =============================================================================

stan_dir <- file.path("inst", "stan")
if (!dir.exists(stan_dir)) stan_dir <- file.path("stan")
if (!dir.exists(stan_dir)) stop("Cannot find Stan model directory")

model_specs <- list(
  baseline = list(stan = file.path(stan_dir, "dlvm.stan"),          label = "Baseline (SE=0)"),
  obs_n    = list(stan = file.path(stan_dir, "dlvm_obs_n.stan"),    label = "A: sigma/sqrt(n)"),
  hier     = list(stan = file.path(stan_dir, "dlvm_hier_sigma.stan"),     label = "B: Hier sigma"),
  hier_n   = list(stan = file.path(stan_dir, "dlvm_hier_sigma_n.stan"),   label = "C: Hier + sqrt(n)")
)

cat("Compiling models...\n")

# Resolve TBB path for runtime use (set only during model$sample() to avoid
# conflicting with data.table's own threading)
TBB_PATH <- file.path(cmdstanr::cmdstan_path(), "stan", "lib", "stan_math", "lib", "tbb")
if (!dir.exists(TBB_PATH)) TBB_PATH <- NULL
if (!is.null(TBB_PATH)) cat(sprintf("  TBB path found: %s\n", TBB_PATH))

models <- list()
for (nm in names(model_specs)) {
  cat(sprintf("  %s: %s\n", nm, model_specs[[nm]]$stan))
  models[[nm]] <- cmdstan_model(
    model_specs[[nm]]$stan,
    cpp_options = list(stan_threads = TRUE),
    quiet = TRUE
  )
}
cat("All models compiled.\n\n")

# =============================================================================
# 2. Data generating process
# =============================================================================

simulate_panel <- function(n_units, n_times, innov_true, sigma_true,
                           obs_counts, sigma_by_unit = NULL,
                           missing_frac = 0, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  nu <- 4       # degrees of freedom matching the model default
  scale <- 4    # scale_state

  # Per-unit sigma
  if (is.null(sigma_by_unit)) {
    sigma_by_unit <- rep(sigma_true, n_units)
  }

  # Generate true latent trajectories
  thetas <- matrix(NA_real_, nrow = n_units, ncol = n_times)
  for (u in 1:n_units) {
    thetas[u, 1] <- rt(1, df = nu) * scale
    for (t in 2:n_times) {
      z <- rt(1, df = nu) * scale
      thetas[u, t] <- thetas[u, t - 1] + innov_true * z
    }
  }

  # Resolve obs_counts to a matrix
  if (is.matrix(obs_counts)) {
    n_mat <- obs_counts
  } else {
    n_mat <- matrix(as.integer(obs_counts), nrow = n_units, ncol = n_times)
  }

  # Apply missingness (set some cells to 0)
  if (missing_frac > 0) {
    n_cells <- n_units * n_times
    n_missing <- round(n_cells * missing_frac)
    eligible <- which(n_mat > 0)
    if (n_missing > 0 && length(eligible) > n_missing) {
      drop_idx <- sample(eligible, n_missing)
      n_mat[drop_idx] <- 0L
      for (u in 1:n_units) {
        if (all(n_mat[u, ] == 0)) {
          t_restore <- sample(n_times, 1)
          n_mat[u, t_restore] <- max(1L, obs_counts[1])
        }
      }
    }
  }

  # Generate observations
  rows <- list()
  for (u in 1:n_units) {
    for (t in 1:n_times) {
      n_obs_cell <- n_mat[u, t]
      if (n_obs_cell == 0) next
      for (k in 1:n_obs_cell) {
        y <- thetas[u, t] + sigma_by_unit[u] * rt(1, df = nu)
        rows[[length(rows) + 1]] <- data.frame(
          unit = paste0("U", u), time = t, value = y,
          unit_idx = u, stringsAsFactors = FALSE
        )
      }
    }
  }

  df <- as.data.table(do.call(rbind, rows))

  # Aggregate to means
  agg <- df[, .(mean_val = mean(value),
                n_obs = .N,
                unit_idx = unit_idx[1]),
            by = .(unit, time)]

  list(
    raw_data = df,
    agg_data = agg,
    thetas_true = thetas,
    innov_true = innov_true,
    sigma_true = sigma_true,
    sigma_by_unit = sigma_by_unit,
    obs_counts = n_mat,
    n_units = n_units,
    n_times = n_times
  )
}

# =============================================================================
# 3. Stan data builders
# =============================================================================

build_state_grid <- function(agg_dt, n_units, n_times) {
  #' Build state grid and mappings (mimics dlvm_prepare internals)
  units <- sort(unique(agg_dt$unit))
  times <- sort(unique(agg_dt$time))

  # Full grid: all units x all observed times
  all_times <- 1:n_times
  sg <- CJ(unit = units, time = all_times)
  sg[, state_id := .I]
  sg[, prev_state_id := 0L]
  sg[, delta_t := 0.0]

  for (u in units) {
    idx <- which(sg$unit == u)
    if (length(idx) > 1) {
      sg$prev_state_id[idx[-1]] <- sg$state_id[idx[-length(idx)]]
      sg$delta_t[idx[-1]] <- diff(sg$time[idx])
    }
  }

  # Map observations to states
  agg_dt[sg, state_id := i.state_id, on = .(unit, time)]

  # Unit index mapping for hierarchical models
  unit_map <- data.table(unit = units, unit_num = seq_along(units))
  agg_dt[unit_map, unit_num := i.unit_num, on = "unit"]

  list(
    state_grid = sg,
    n_states = nrow(sg),
    state_prev = as.integer(sg$prev_state_id),
    delta_t = sg$delta_t
  )
}

make_stan_data <- function(sim, variant = "baseline") {
  #' Build Stan data list for a given model variant.
  agg <- copy(sim$agg_data)
  grid <- build_state_grid(agg, sim$n_units, sim$n_times)

  base <- list(
    n_states   = grid$n_states,
    n_obs      = nrow(agg),
    n_units    = sim$n_units,
    y          = agg$mean_val,
    obs_to_state = as.integer(agg$state_id),
    has_se     = 0L,
    se         = rep(0, nrow(agg)),
    obs_n      = rep(1L, nrow(agg)),  # default: no sqrt(n) scaling
    state_prev = grid$state_prev,
    delta_t    = grid$delta_t,
    nu_obs     = 4,
    nu_state   = 4,
    scale_state = 4,
    grainsize  = 1L,
    compute_gq = 0L
  )

  # Models that use observation-count scaling get the real counts
  if (variant == "obs_n" || variant == "hier_n") {
    base$obs_n <- as.integer(agg$n_obs)
  }
  if (variant == "hier" || variant == "hier_n") {
    base$obs_to_unit <- as.integer(agg$unit_num)
  }

  base
}

# =============================================================================
# 4. Fitting & evaluation
# =============================================================================

fit_model <- function(model, stan_data, warmup = 500, sampling = 500,
                      chains = 2, seed = 42) {
  # Set TBB path only for sampling (avoid data.table conflict)
  old_dyld <- Sys.getenv("DYLD_LIBRARY_PATH", unset = NA)
  if (!is.null(TBB_PATH)) {
    if (is.na(old_dyld) || old_dyld == "") {
      Sys.setenv(DYLD_LIBRARY_PATH = TBB_PATH)
    } else {
      Sys.setenv(DYLD_LIBRARY_PATH = paste(TBB_PATH, old_dyld, sep = ":"))
    }
  }

  t0 <- proc.time()["elapsed"]
  fit <- tryCatch(
    model$sample(
      data = stan_data,
      seed = seed,
      chains = chains,
      parallel_chains = chains,
      threads_per_chain = 1L,
      iter_warmup = warmup,
      iter_sampling = sampling,
      adapt_delta = 0.9,
      refresh = 0,
      show_messages = FALSE,
      show_exceptions = FALSE
    ),
    finally = {
      # Restore DYLD_LIBRARY_PATH
      if (is.na(old_dyld)) Sys.unsetenv("DYLD_LIBRARY_PATH")
      else Sys.setenv(DYLD_LIBRARY_PATH = old_dyld)
    }
  )
  elapsed <- proc.time()["elapsed"] - t0
  list(fit = fit, elapsed = elapsed)
}

evaluate_fit <- function(fit_result, sim, variant_label) {
  fit <- fit_result$fit

  # Scalar parameter recovery
  innov_draws <- as.vector(fit$draws("innov", format = "matrix"))
  innov_med <- median(innov_draws)
  innov_ci <- quantile(innov_draws, c(0.05, 0.95))
  innov_covered <- sim$innov_true >= innov_ci[1] && sim$innov_true <= innov_ci[2]

  # Sigma recovery — for hier models, use exp(log_sigma_mu) as the global equiv
  if (variant_label %in% c("B: Hier sigma", "C: Hier + sqrt(n)")) {
    sigma_draws <- as.vector(fit$draws("sigma", format = "matrix"))
  } else {
    sigma_draws <- as.vector(fit$draws("sigma", format = "matrix"))
  }
  sigma_med <- median(sigma_draws)
  sigma_ci <- quantile(sigma_draws, c(0.05, 0.95))
  sigma_covered <- sim$sigma_true >= sigma_ci[1] && sim$sigma_true <= sigma_ci[2]

  # Theta recovery
  theta_draws <- fit$draws("theta", format = "matrix")
  theta_med <- apply(theta_draws, 2, median)
  theta_lo <- apply(theta_draws, 2, quantile, 0.05)
  theta_hi <- apply(theta_draws, 2, quantile, 0.95)

  true_flat <- as.vector(t(sim$thetas_true))  # flatten row-major (unit × time)
  n_states <- length(true_flat)

  # Overall theta metrics
  theta_bias <- mean(theta_med - true_flat)
  theta_rmse <- sqrt(mean((theta_med - true_flat)^2))
  theta_coverage <- mean(true_flat >= theta_lo & true_flat <= theta_hi)

  # Identify which states have n=1 vs n>1 observations.
  # Use sim$obs_counts matrix directly: it's [n_units x n_times], same order
  # as true_flat = as.vector(t(thetas_true)).
  obs_n_by_state <- as.vector(t(sim$obs_counts))

  n1_idx <- which(obs_n_by_state == 1)
  ngt1_idx <- which(obs_n_by_state > 1)
  n0_idx <- which(obs_n_by_state == 0)

  # n=1 specific metrics
  if (length(n1_idx) > 0) {
    theta_rmse_n1 <- sqrt(mean((theta_med[n1_idx] - true_flat[n1_idx])^2))
    theta_cov_n1 <- mean(true_flat[n1_idx] >= theta_lo[n1_idx] &
                         true_flat[n1_idx] <= theta_hi[n1_idx])
    theta_width_n1 <- mean(theta_hi[n1_idx] - theta_lo[n1_idx])
  } else {
    theta_rmse_n1 <- NA; theta_cov_n1 <- NA; theta_width_n1 <- NA
  }

  # n>1 specific metrics
  if (length(ngt1_idx) > 0) {
    theta_rmse_ngt1 <- sqrt(mean((theta_med[ngt1_idx] - true_flat[ngt1_idx])^2))
    theta_cov_ngt1 <- mean(true_flat[ngt1_idx] >= theta_lo[ngt1_idx] &
                           true_flat[ngt1_idx] <= theta_hi[ngt1_idx])
    theta_width_ngt1 <- mean(theta_hi[ngt1_idx] - theta_lo[ngt1_idx])
  } else {
    theta_rmse_ngt1 <- NA; theta_cov_ngt1 <- NA; theta_width_ngt1 <- NA
  }

  # n=0 (unobserved) metrics
  if (length(n0_idx) > 0) {
    theta_rmse_n0 <- sqrt(mean((theta_med[n0_idx] - true_flat[n0_idx])^2))
    theta_cov_n0 <- mean(true_flat[n0_idx] >= theta_lo[n0_idx] &
                         true_flat[n0_idx] <= theta_hi[n0_idx])
  } else {
    theta_rmse_n0 <- NA; theta_cov_n0 <- NA
  }

  # Diagnostics
  diag <- fit$diagnostic_summary(quiet = TRUE)
  n_div <- sum(diag$num_divergent)

  data.table(
    model = variant_label,
    innov_med = innov_med,
    innov_true = sim$innov_true,
    innov_bias = innov_med - sim$innov_true,
    innov_covered = innov_covered,
    sigma_med = sigma_med,
    sigma_true = sim$sigma_true,
    sigma_bias = sigma_med - sim$sigma_true,
    sigma_covered = sigma_covered,
    theta_rmse = theta_rmse,
    theta_coverage = theta_coverage,
    theta_rmse_n1 = theta_rmse_n1,
    theta_cov_n1 = theta_cov_n1,
    theta_width_n1 = theta_width_n1,
    theta_rmse_ngt1 = theta_rmse_ngt1,
    theta_cov_ngt1 = theta_cov_ngt1,
    theta_width_ngt1 = theta_width_ngt1,
    theta_rmse_n0 = theta_rmse_n0,
    theta_cov_n0 = theta_cov_n0,
    n_divergent = n_div,
    elapsed_s = fit_result$elapsed,
    n_n1_cells = length(n1_idx),
    n_ngt1_cells = length(ngt1_idx),
    n_n0_cells = length(n0_idx)
  )
}

# =============================================================================
# 5. Scenario definitions
# =============================================================================

make_mixed_obs_counts <- function(n_units, n_times, frac_n1 = 0.5,
                                  n_sparse = 2, n_dense = 10, seed = NULL) {
  #' Create a matrix of observation counts with a mix of sparse and dense cells.
  if (!is.null(seed)) set.seed(seed)
  n_mat <- matrix(n_dense, nrow = n_units, ncol = n_times)
  n_cells <- n_units * n_times
  n_sparse_cells <- round(n_cells * frac_n1)
  sparse_idx <- sample(n_cells, n_sparse_cells)
  n_mat[sparse_idx] <- n_sparse
  n_mat
}

scenarios <- list(
  # Scenario 1: Balanced, homogeneous — all cells have n=5
  balanced_homo = list(
    n_units = 8, n_times = 15, innov = 0.3, sigma = 0.25,
    obs_counts = 5L,
    sigma_by_unit = NULL,
    missing_frac = 0,
    desc = "Balanced (n=5), homogeneous sigma"
  ),

  # Scenario 2: Mixed n (many n=1), homogeneous — tests n=1 handling
  mixed_homo = list(
    n_units = 8, n_times = 15, innov = 0.3, sigma = 0.25,
    obs_counts = "mixed_n1",
    sigma_by_unit = NULL,
    missing_frac = 0,
    desc = "Mixed n (60% n=1, 40% n=8), homogeneous sigma"
  ),

  # Scenario 3: Heavy sparsity with missingness — near real-world conditions
  sparse_missing = list(
    n_units = 10, n_times = 15, innov = 0.3, sigma = 0.25,
    obs_counts = "mostly_n1",
    sigma_by_unit = NULL,
    missing_frac = 0.3,
    desc = "Sparse (80% n=1) + 30% missing cells"
  ),

  # Scenario 4: Balanced, heterogeneous sigma — tests hierarchical need
  balanced_hetero = list(
    n_units = 8, n_times = 15, innov = 0.3, sigma = 0.25,
    obs_counts = 5L,
    sigma_by_unit = c(0.1, 0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.5),
    missing_frac = 0,
    desc = "Balanced (n=5), heterogeneous sigma (0.1-0.5)"
  ),

  # Scenario 5: Mixed n + heterogeneous sigma — full complexity
  mixed_hetero = list(
    n_units = 8, n_times = 15, innov = 0.3, sigma = 0.25,
    obs_counts = "mixed_n1",
    sigma_by_unit = c(0.1, 0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.5),
    missing_frac = 0,
    desc = "Mixed n + heterogeneous sigma — full complexity"
  ),

  # Scenario 6: Extreme sparsity — some units with only 1 obs total
  extreme_sparse = list(
    n_units = 10, n_times = 15, innov = 0.3, sigma = 0.25,
    obs_counts = "extreme",
    sigma_by_unit = NULL,
    missing_frac = 0,
    desc = "Extreme: 3 units have only 1 obs total"
  )
)

resolve_obs_counts <- function(spec, seed) {
  nu <- spec$n_units; nt <- spec$n_times
  oc <- spec$obs_counts
  if (is.numeric(oc) && length(oc) == 1) {
    return(matrix(as.integer(oc), nrow = nu, ncol = nt))
  }
  set.seed(seed + 999)
  if (oc == "mixed_n1") {
    # 60% of cells get n=1, rest get n=8
    n_mat <- matrix(8L, nrow = nu, ncol = nt)
    sparse_count <- round(nu * nt * 0.6)
    sparse_idx <- sample(nu * nt, sparse_count)
    n_mat[sparse_idx] <- 1L
    return(n_mat)
  }
  if (oc == "mostly_n1") {
    # 80% n=1, 15% n=2-3, 5% n=5
    n_mat <- matrix(1L, nrow = nu, ncol = nt)
    cells <- nu * nt
    n2_count <- round(cells * 0.15)
    n5_count <- round(cells * 0.05)
    idx <- sample(cells, n2_count + n5_count)
    n_mat[idx[1:n2_count]] <- sample(2:3, n2_count, replace = TRUE)
    n_mat[idx[(n2_count + 1):length(idx)]] <- 5L
    return(n_mat)
  }
  if (oc == "extreme") {
    # Units 1-5: good coverage (n=5-10); units 6-7: sparse (n=1-2)
    # Units 8-10: only 1 observation in entire time series
    n_mat <- matrix(0L, nrow = nu, ncol = nt)
    for (u in 1:5) n_mat[u, ] <- sample(5:10, nt, replace = TRUE)
    for (u in 6:7) {
      obs_times <- sample(nt, round(nt * 0.4))
      n_mat[u, obs_times] <- sample(1:2, length(obs_times), replace = TRUE)
    }
    for (u in 8:10) {
      single_t <- sample(nt, 1)
      n_mat[u, single_t] <- 1L
    }
    return(n_mat)
  }
  stop("Unknown obs_counts spec: ", oc)
}

# =============================================================================
# 6. Main runner
# =============================================================================

N_REPS <- if (QUICK) 1 else 3
scenarios_to_run <- if (QUICK) c("mixed_homo", "extreme_sparse") else names(scenarios)

cat(sprintf("Running %d scenario(s) x %d rep(s) x 4 models = %d fits\n\n",
            length(scenarios_to_run), N_REPS,
            length(scenarios_to_run) * N_REPS * 4))

all_results <- list()

for (sc_name in scenarios_to_run) {
  sc <- scenarios[[sc_name]]
  cat(sprintf("--- Scenario: %s ---\n", sc$desc)); flush(stdout())

  for (rep_i in 1:N_REPS) {
    seed <- 1000 * which(names(scenarios) == sc_name) + rep_i
    cat(sprintf("  Resolving obs counts (seed=%d)...\n", seed)); flush(stdout())
    obs_mat <- resolve_obs_counts(sc, seed)

    cat("  Generating simulated data...\n"); flush(stdout())
    sim <- simulate_panel(
      n_units = sc$n_units,
      n_times = sc$n_times,
      innov_true = sc$innov,
      sigma_true = sc$sigma,
      obs_counts = obs_mat,
      sigma_by_unit = sc$sigma_by_unit,
      missing_frac = sc$missing_frac,
      seed = seed
    )
    cat(sprintf("  Data generated: %d raw obs, %d aggregated cells\n",
                nrow(sim$raw_data), nrow(sim$agg_data))); flush(stdout())

    for (variant in names(models)) {
      label <- model_specs[[variant]]$label
      cat(sprintf("  [rep %d] Building data for %s ... ", rep_i, label)); flush(stdout())

      stan_data <- tryCatch(
        make_stan_data(sim, variant),
        error = function(e) { cat(sprintf("DATA ERROR: %s\n", e$message)); flush(stdout()); NULL }
      )
      if (is.null(stan_data)) next
      cat(sprintf("OK (n_states=%d, n_obs=%d)\n", stan_data$n_states, stan_data$n_obs)); flush(stdout())

      cat(sprintf("  [rep %d] Fitting %s ... ", rep_i, label)); flush(stdout())
      fit_result <- tryCatch(
        fit_model(models[[variant]], stan_data, seed = seed),
        error = function(e) { cat(sprintf("FIT ERROR: %s\n", e$message)); flush(stdout()); NULL }
      )
      if (is.null(fit_result)) next

      cat("evaluating... "); flush(stdout())
      metrics <- tryCatch(
        evaluate_fit(fit_result, sim, label),
        error = function(e) {
          cat(sprintf("EVAL ERROR: %s\n", e$message)); flush(stdout())
          NULL
        }
      )
      if (is.null(metrics)) next

      metrics[, `:=`(scenario = sc_name, scenario_desc = sc$desc, rep = rep_i)]
      all_results[[length(all_results) + 1]] <- metrics

      cat(sprintf("done (%.0fs, %d div)\n", fit_result$elapsed, metrics$n_divergent)); flush(stdout())

      # Cleanup to free memory
      rm(fit_result, stan_data, metrics)
      gc(verbose = FALSE)
    }
  }
  cat("\n")
}

# =============================================================================
# 7. Summary tables
# =============================================================================

results <- rbindlist(all_results)

cat("\n=== PARAMETER RECOVERY (averaged over replications) ===\n\n")

param_summary <- results[, .(
  innov_bias = mean(innov_bias),
  innov_cov  = mean(innov_covered),
  sigma_bias = mean(sigma_bias),
  sigma_cov  = mean(sigma_covered),
  n_div      = mean(n_divergent),
  time_s     = mean(elapsed_s)
), by = .(scenario_desc, model)]

# Print wide format
for (sc_desc in unique(param_summary$scenario_desc)) {
  cat(sprintf("\n  %s\n", sc_desc))
  sub <- param_summary[scenario_desc == sc_desc]
  cat(sprintf("  %-22s %8s %8s %8s %8s %6s %6s\n",
              "Model", "fInnBias", "InnCov", "fSigBias", "SigCov", "Divs", "Time"))
  cat(sprintf("  %s\n", strrep("-", 80)))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-22s %+8.4f %7.0f%% %+8.4f %7.0f%% %6.1f %5.0fs\n",
                sub$model[i],
                sub$innov_bias[i], sub$innov_cov[i] * 100,
                sub$sigma_bias[i], sub$sigma_cov[i] * 100,
                sub$n_div[i], sub$time_s[i]))
  }
}

cat("\n\n=== THETA RECOVERY BY OBSERVATION COUNT ===\n\n")

theta_summary <- results[, .(
  rmse_overall = mean(theta_rmse),
  cov_overall  = mean(theta_coverage),
  rmse_n1      = mean(theta_rmse_n1, na.rm = TRUE),
  cov_n1       = mean(theta_cov_n1, na.rm = TRUE),
  width_n1     = mean(theta_width_n1, na.rm = TRUE),
  rmse_ngt1    = mean(theta_rmse_ngt1, na.rm = TRUE),
  cov_ngt1     = mean(theta_cov_ngt1, na.rm = TRUE),
  width_ngt1   = mean(theta_width_ngt1, na.rm = TRUE),
  rmse_n0      = mean(theta_rmse_n0, na.rm = TRUE),
  cov_n0       = mean(theta_cov_n0, na.rm = TRUE)
), by = .(scenario_desc, model)]

for (sc_desc in unique(theta_summary$scenario_desc)) {
  cat(sprintf("\n  %s\n", sc_desc))
  sub <- theta_summary[scenario_desc == sc_desc]

  cat(sprintf("  %-22s %8s %8s | %8s %8s %8s | %8s %8s %8s\n",
              "Model", "RMSE_all", "Cov_all",
              "RMSE_n1", "Cov_n1", "Width_n1",
              "RMSE_n>1", "Cov_n>1", "Wid_n>1"))
  cat(sprintf("  %s\n", strrep("-", 105)))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-22s %8.4f %7.0f%% | %8s %8s %8s | %8s %8s %8s\n",
                sub$model[i],
                sub$rmse_overall[i], sub$cov_overall[i] * 100,
                ifelse(is.na(sub$rmse_n1[i]), "  N/A   ", sprintf("%8.4f", sub$rmse_n1[i])),
                ifelse(is.na(sub$cov_n1[i]),  "  N/A   ", sprintf("%6.0f%%", sub$cov_n1[i] * 100)),
                ifelse(is.na(sub$width_n1[i]),"  N/A   ", sprintf("%8.4f", sub$width_n1[i])),
                ifelse(is.na(sub$rmse_ngt1[i]), " N/A   ", sprintf("%8.4f", sub$rmse_ngt1[i])),
                ifelse(is.na(sub$cov_ngt1[i]),  " N/A   ", sprintf("%6.0f%%", sub$cov_ngt1[i] * 100)),
                ifelse(is.na(sub$width_ngt1[i])," N/A   ", sprintf("%8.4f", sub$width_ngt1[i]))))
  }
}

# Save results
out_path <- file.path("tests", "sim_results_obs_models.rds")
saveRDS(results, out_path)
cat(sprintf("\n\nResults saved to: %s\n", out_path))
cat("=== Simulation complete ===\n")
