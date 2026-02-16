#!/usr/bin/env Rscript
# tests/simulation_cureton_effect.R
# Compare Cureton-adjusted vs Standard SEM on small-N, high-variance data
# With Targeted "Anomaly" Injection + Posterior Width Analysis

suppressPackageStartupMessages({
  library(dlvm)
  library(cmdstanr)
  library(data.table)
  library(ggplot2)
})

set.seed(2023)

cat("=== Cureton Adjustment Simulation (Posterior Width Analysis) ===\n")
cat("Parameters:\n")
cat("  Units: 5\n  Time points: 30\n  True Innovation: 0.5\n  True Sigma: 1.0\n  N per time: 2 (Outlier), 3 (Normal)\n\n")

# --- 1. Data Generation ---
n_units <- 5
n_times <- 30
innov_true <- 0.5
sigma_true <- 1.0
units <- paste0("unit_", 1:n_units)

# Generate latent trajectories
thetas <- matrix(NA_real_, nrow = n_units, ncol = n_times)
for (u in 1:n_units) {
  thetas[u, 1] <- rnorm(1, 0, 2)
  for (t in 2:n_times) {
    z <- rnorm(1, 0, 1)
    thetas[u, t] <- thetas[u, t - 1] + innov_true * z
  }
}

# Generate observations
obs_list <- list()
anomaly_indices <- list() # Store indices of anomalous observations

for (u in 1:n_units) {
  for (t in 1:n_times) {
    # Is this an anomaly time point?
    # Inject anomaly at t=15 for unit 1, t=25 for unit 2
    is_anomaly <- (u == 1 && t == 15) || (u == 2 && t == 25)

    if (is_anomaly) {
      target_val <- thetas[u, t] + 5.0
      vals <- c(target_val, target_val + 0.1) # N=2, tight cluster, far from truth
      n_obs <- 2
      anomaly_indices[[length(anomaly_indices) + 1]] <- list(u = u, t = t, val = mean(vals))
    } else {
      n_obs <- 3
      vals <- rnorm(n_obs, mean = thetas[u, t], sd = sigma_true)
    }

    obs_list[[length(obs_list) + 1]] <- data.frame(
      unit = units[u],
      time = t,
      value = vals,
      stringsAsFactors = FALSE
    )
  }
}
sim_data <- do.call(rbind, obs_list)

cat(sprintf("Generated %d total observations.\n", nrow(sim_data)))

# --- 2. Fit Models ---

# Model A: With Cureton Correction
cat("\nFitting Model A: With Cureton Correction...\n")
prep_cureton <- dlvm_prepare(sim_data, "unit", "time", "value", mode = "means", se_correction = TRUE)
fit_cureton <- dlvm_fit(prep_cureton, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                        iter_warmup = 300, iter_sampling = 500,
                        show_messages = FALSE, refresh = 0)

# Model B: Standard SEM
cat("\nFitting Model B: Standard SEM (Uncorrected)...\n")
prep_std <- dlvm_prepare(sim_data, "unit", "time", "value", mode = "means", se_correction = FALSE)
fit_std <- dlvm_fit(prep_std, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500,
                    show_messages = FALSE, refresh = 0)


# --- 3. Compare Results ---
get_theta_stats <- function(fit) {
  tryCatch({
    theta_draws <- fit$fit$draws("theta", format = "matrix")
    list(
      mean = colMeans(theta_draws),
      sd = apply(theta_draws, 2, sd)
    )
  }, error = function(e) return(NULL))
}

stats_cureton <- get_theta_stats(fit_cureton)
stats_std     <- get_theta_stats(fit_std)

if (is.null(stats_cureton) || is.null(stats_std)) stop("Model fitting failed.")

# Align with Truth
true_vec <- as.vector(t(thetas))

cat("\n=== Posterior Stability Analysis ===\n")
cat(sprintf("Average Posterior SD (Uncertainty):\n"))
cat(sprintf("  Cureton:  %.4f\n", mean(stats_cureton$sd)))
cat(sprintf("  Standard: %.4f\n", mean(stats_std$sd)))
diff_sd <- mean(stats_cureton$sd) - mean(stats_std$sd)
cat(sprintf("  Difference: %.4f (%.1f%%)\n", diff_sd, (diff_sd / mean(stats_std$sd)) * 100))


cat("\n=== Anomaly Point Analysis ===\n")
cat(sprintf("%-6s %-6s %-10s %-10s %-10s %-10s %-10s\n",
            "Unit", "Time", "Type", "Est", "PostSD", "InputSE", "True"))

for (a in anomaly_indices) {
  idx <- (a$u - 1) * n_times + a$t

  # Get input SEs from prep objects
  # We need to find the row in prep$metadata$obs_data corresponding to this unit/time
  # But since we have balanced data, the row index in obs_data matches idx (roughly)
  # Actually, prep objects sort by unit then time, so idx is correct for obs_data too.
  se_cureton <- prep_cureton$metadata$obs_data$se[idx]
  se_std     <- prep_std$metadata$obs_data$se[idx]

  cat(sprintf("%-6d %-6d %-10s %-10.4f %-10.4f %-10.4f %-10.4f\n",
              a$u, a$t, "Cureton", stats_cureton$mean[idx], stats_cureton$sd[idx], se_cureton, thetas[a$u, a$t]))
  cat(sprintf("%-6d %-6d %-10s %-10.4f %-10.4f %-10.4f %-10.4f\n",
              a$u, a$t, "Standard", stats_std$mean[idx], stats_std$sd[idx], se_std, thetas[a$u, a$t]))
  cat("-------------------------------------------------------------\n")
}

cat("\nDone.\n")
