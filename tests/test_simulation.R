# dlvm/tests/test_simulation.R
# Simulation-based calibration: generate known data, fit, check recovery
# Requires CmdStan — will be skipped if not available.

cat("  --- Simulation-Based Calibration Tests ---\n")

# Check CmdStan availability
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("  SKIPPED: cmdstanr not available\n")
  return(invisible(NULL))
}
tryCatch(cmdstanr::cmdstan_path(), error = function(e) {
  cat("  SKIPPED: CmdStan backend not installed\n")
  return(invisible(NULL))
})

# --- Simulate data from the generative model ---
simulate_dlvm_data <- function(n_units = 5,
                               n_times = 20,
                               innov_true = 0.3,
                               sigma_true = 0.2,
                               nu = 4,
                               obs_per_time = 1,
                               seed = 12345) {
  set.seed(seed)
  units <- paste0("unit_", seq_len(n_units))
  times <- seq_len(n_times)

  # Generate latent trajectories
  thetas <- matrix(NA_real_, nrow = n_units, ncol = n_times)
  for (u in seq_len(n_units)) {
    thetas[u, 1] <- rt(1, df = nu) * 4  # Draw from init prior: student_t(4, 0, 4)
    for (t in 2:n_times) {
      # Innovation: theta[t] = theta[t-1] + innov * z, z ~ student_t(nu, 0, 4)
      # For delta_t = 1
      z <- rt(1, df = nu) * 4
      thetas[u, t] <- thetas[u, t - 1] + innov_true * z
    }
  }

  # Generate observations from latent states
  rows <- list()
  for (u in seq_len(n_units)) {
    for (t in times) {
      for (k in seq_len(obs_per_time)) {
        # y ~ student_t(nu, theta, sigma)
        y <- thetas[u, t] + sigma_true * rt(1, df = nu)
        rows[[length(rows) + 1]] <- data.frame(
          unit = units[u],
          time = t,
          value = y,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  df <- do.call(rbind, rows)

  return(list(
    data = df,
    thetas_true = thetas,
    innov_true = innov_true,
    sigma_true = sigma_true,
    units = units,
    times = times
  ))
}

# --- Test 1: Single simulation recovery ---
test_that("parameter recovery from simulated data (single run)", {
  sim <- simulate_dlvm_data(n_units = 3, n_times = 15, innov_true = 0.3, sigma_true = 0.2)

  prep <- dlvm_prepare(sim$data, "unit", "time", "value", mode = "means")
  fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                  iter_warmup = 500, iter_sampling = 1000, show_messages = FALSE)

  # Extract scalar parameter posteriors
  innov_draws <- fit$fit$draws("innov", format = "matrix")
  sigma_draws <- fit$fit$draws("sigma", format = "matrix")

  innov_ci <- quantile(innov_draws, c(0.05, 0.95))
  sigma_ci <- quantile(sigma_draws, c(0.05, 0.95))

  cat(sprintf("    innov: true=%.2f, 90%% CI=[%.3f, %.3f]\n",
              sim$innov_true, innov_ci[1], innov_ci[2]))
  cat(sprintf("    sigma: true=%.2f, 90%% CI=[%.3f, %.3f]\n",
              sim$sigma_true, sigma_ci[1], sigma_ci[2]))

  # Check that true values fall within 90% CI
  expect_true(innov_ci[1] <= sim$innov_true && sim$innov_true <= innov_ci[2])
  expect_true(sigma_ci[1] <= sim$sigma_true && sim$sigma_true <= sigma_ci[2])

  # Check no divergences
  diag <- fit$fit$diagnostic_summary(quiet = TRUE)
  expect_equal(sum(diag$num_divergent), 0)
})

# --- Test 2: Multiple-obs-per-time recovery ---
test_that("individual-obs mode recovers parameters", {
  sim <- simulate_dlvm_data(n_units = 3, n_times = 10,
                            innov_true = 0.25, sigma_true = 0.3,
                            obs_per_time = 3, seed = 54321)

  prep <- dlvm_prepare(sim$data, "unit", "time", "value", mode = "individual")
  expect_equal(prep$metadata$mode, "individual")

  fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                  iter_warmup = 500, iter_sampling = 1000, show_messages = FALSE)

  # Parameter recovery
  innov_draws <- fit$fit$draws("innov", format = "matrix")
  sigma_draws <- fit$fit$draws("sigma", format = "matrix")

  innov_ci <- quantile(innov_draws, c(0.05, 0.95))
  sigma_ci <- quantile(sigma_draws, c(0.05, 0.95))

  cat(sprintf("    innov: true=%.2f, 90%% CI=[%.3f, %.3f]\n",
              sim$innov_true, innov_ci[1], innov_ci[2]))
  cat(sprintf("    sigma: true=%.2f, 90%% CI=[%.3f, %.3f]\n",
              sim$sigma_true, sigma_ci[1], sigma_ci[2]))

  expect_true(innov_ci[1] <= sim$innov_true && sim$innov_true <= innov_ci[2])
  expect_true(sigma_ci[1] <= sim$sigma_true && sim$sigma_true <= sigma_ci[2])
})

# --- Test 3: Latent trajectory corridor coverage ---
test_that("latent trajectory recovery (corridor coverage)", {
  sim <- simulate_dlvm_data(n_units = 3, n_times = 15, innov_true = 0.3, sigma_true = 0.2)
  prep <- dlvm_prepare(sim$data, "unit", "time", "value", mode = "means")
  fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                  iter_warmup = 500, iter_sampling = 1000, show_messages = FALSE)

  summ <- dlvm_summary(fit, prob = 0.90)

  # Check what proportion of true thetas fall within 90% CI
  true_flat <- as.vector(t(sim$thetas_true))
  n_states <- length(true_flat)

  covered <- sum(true_flat >= summ$lower & true_flat <= summ$upper, na.rm = TRUE)
  coverage <- covered / n_states

  cat(sprintf("    Corridor coverage: %.1f%% (target: >= 70%%)\n", coverage * 100))

  # 90% CI should cover at least 70% of true values
  # (relaxed from 90% because small data + heavy tails)
  expect_gte(coverage, 0.70)
})

# --- Test 4: LOO-CV computes without error ---
test_that("LOO-CV runs on simulated fit", {
  sim <- simulate_dlvm_data(n_units = 2, n_times = 10, innov_true = 0.3, sigma_true = 0.2)
  prep <- dlvm_prepare(sim$data, "unit", "time", "value")
  fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                  iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)

  loo_result <- dlvm_loo(fit)
  expect_true(!is.null(loo_result$estimates))
  cat(sprintf("    LOO ELPD: %.2f (SE: %.2f)\n",
              loo_result$estimates["elpd_loo", "Estimate"],
              loo_result$estimates["elpd_loo", "SE"]))
})

cat("  --- Simulation Tests Complete ---\n")
