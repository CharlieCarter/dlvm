#!/usr/bin/env Rscript
# dlvm/tests/test_simulation_families.R
#
# Simulation-based parameter recovery for all four observation families.
# Generates data from known DGPs, fits the corresponding DLVM model,
# and checks whether the true parameters fall within the 90% credible intervals.
#
# Usage (from dlvm/ root):
#   Rscript tests/test_simulation_families.R
#
# Requires: CmdStan compiled models (run dlvm_compile() for each family first)

cat("=== Simulation-Based Parameter Recovery ===\n")
cat("Loading source files...\n")
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

set.seed(42)

# ============================================================================
# Helper: simulate latent theta via random walk
# ============================================================================
simulate_theta <- function(n_units, n_times, innov_true, nu_state = 4) {
  theta <- matrix(NA, n_units, n_times)
  for (u in 1:n_units) {
    # First time point: draw from student_t(4, 0, 4)
    theta[u, 1] <- rt(1, df = 4) * 4
    for (t in 2:n_times) {
      # Random walk with student_t innovations (delta_t = 1)
      theta[u, t] <- theta[u, t-1] + rt(1, df = nu_state) * innov_true
    }
  }
  theta
}

# ============================================================================
# Helper: check parameter recovery
# ============================================================================
check_recovery <- function(fit, param_name, true_value, prob = 0.90) {
  draws <- fit$fit$draws(param_name, format = "matrix")
  lower <- quantile(draws, (1 - prob) / 2)
  upper <- quantile(draws, 1 - (1 - prob) / 2)
  median_est <- median(draws)
  covered <- (true_value >= lower) && (true_value <= upper)
  cat(sprintf("    %s: true=%.3f, median=%.3f, 90%% CI=[%.3f, %.3f] %s\n",
              param_name, true_value, median_est, lower, upper,
              if (covered) "RECOVERED" else "MISSED"))
  covered
}

check_theta_recovery <- function(fit, theta_true, prob = 0.90) {
  # Check what fraction of true theta values fall within 90% CIs
  theta_draws <- fit$fit$draws("theta", format = "matrix")
  n_states <- ncol(theta_draws)
  covered <- 0
  for (s in 1:n_states) {
    lower <- quantile(theta_draws[, s], (1 - prob) / 2)
    upper <- quantile(theta_draws[, s], 1 - (1 - prob) / 2)
    if (theta_true[s] >= lower && theta_true[s] <= upper) covered <- covered + 1
  }
  coverage <- covered / n_states
  cat(sprintf("    theta: %d/%d states covered (%.0f%% vs nominal %.0f%%)\n",
              covered, n_states, coverage * 100, prob * 100))
  coverage
}

# Sampling settings (fast but sufficient for recovery checks)
CHAINS <- 4
ITER_WARMUP <- 500
ITER_SAMPLING <- 1000
THREADS <- 2

# ============================================================================
# 1. Student-t family
# ============================================================================
cat("\n--- 1. Student-t Family ---\n")

n_units <- 5; n_times <- 10
innov_true <- 0.5; sigma_true <- 1.0

theta <- simulate_theta(n_units, n_times, innov_true)

# Generate observations: y ~ student_t(4, theta, sigma)
df_list <- list()
for (u in 1:n_units) {
  for (t in 1:n_times) {
    n_obs_per <- sample(3:6, 1)
    for (j in 1:n_obs_per) {
      df_list[[length(df_list) + 1]] <- data.frame(
        unit = paste0("U", u), time = 2010 + t - 1,
        value = theta[u, t] + rt(1, df = 4) * sigma_true
      )
    }
  }
}
df_st <- do.call(rbind, df_list)

cat("  Preparing data...\n")
prep <- dlvm_prepare(df_st, "unit", "time", "value",
                     family = "student_t", mode = "individual",
                     compute_gq = TRUE)

cat("  Fitting model...\n")
fit_st <- dlvm_fit(prep, chains = CHAINS, threads_per_chain = THREADS,
                   iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                   show_messages = FALSE)

cat("  Checking recovery:\n")
innov_ok <- check_recovery(fit_st, "innov", innov_true)
sigma_ok <- check_recovery(fit_st, "sigma", sigma_true)
theta_vec <- as.vector(t(theta))  # flatten by unit then time
theta_cov <- check_theta_recovery(fit_st, theta_vec)

cat(sprintf("  Student-t result: innov=%s, sigma=%s, theta_coverage=%.0f%%\n",
            if(innov_ok) "OK" else "MISS", if(sigma_ok) "OK" else "MISS",
            theta_cov * 100))

# ============================================================================
# 2. Bernoulli family
# ============================================================================
cat("\n--- 2. Bernoulli Family ---\n")

n_units <- 5; n_times <- 10
innov_true <- 0.3

theta <- simulate_theta(n_units, n_times, innov_true)
# Keep theta in a reasonable range for logit
theta <- theta * 0.5

df_list <- list()
for (u in 1:n_units) {
  for (t in 1:n_times) {
    n_obs_per <- sample(10:20, 1)  # need more obs for binary
    prob <- plogis(theta[u, t])
    for (j in 1:n_obs_per) {
      df_list[[length(df_list) + 1]] <- data.frame(
        unit = paste0("U", u), time = 2010 + t - 1,
        value = rbinom(1, 1, prob)
      )
    }
  }
}
df_bern <- do.call(rbind, df_list)

cat("  Preparing data...\n")
prep <- dlvm_prepare(df_bern, "unit", "time", "value",
                     family = "bernoulli", mode = "individual",
                     compute_gq = TRUE)

cat("  Fitting model...\n")
fit_bern <- dlvm_fit(prep, chains = CHAINS, threads_per_chain = THREADS,
                     iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                     show_messages = FALSE)

cat("  Checking recovery:\n")
innov_ok <- check_recovery(fit_bern, "innov", innov_true * 0.5)  # scaled
theta_vec <- as.vector(t(theta))
theta_cov <- check_theta_recovery(fit_bern, theta_vec)

cat(sprintf("  Bernoulli result: innov=%s, theta_coverage=%.0f%%\n",
            if(innov_ok) "OK" else "MISS", theta_cov * 100))

# ============================================================================
# 3. Binomial family
# ============================================================================
cat("\n--- 3. Binomial Family ---\n")

n_units <- 5; n_times <- 10
innov_true <- 0.3

theta <- simulate_theta(n_units, n_times, innov_true)
theta <- theta * 0.5  # keep in reasonable logit range

df_list <- list()
for (u in 1:n_units) {
  for (t in 1:n_times) {
    n_trials <- sample(20:50, 1)
    prob <- plogis(theta[u, t])
    successes <- rbinom(1, n_trials, prob)
    df_list[[length(df_list) + 1]] <- data.frame(
      unit = paste0("U", u), time = 2010 + t - 1,
      value = successes, trials = n_trials
    )
  }
}
df_binom <- do.call(rbind, df_list)

cat("  Preparing data...\n")
prep <- dlvm_prepare(df_binom, "unit", "time", "value",
                     trials_col = "trials", family = "binomial",
                     compute_gq = TRUE)

cat("  Fitting model...\n")
fit_binom <- dlvm_fit(prep, chains = CHAINS, threads_per_chain = THREADS,
                      iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                      show_messages = FALSE)

cat("  Checking recovery:\n")
innov_ok <- check_recovery(fit_binom, "innov", innov_true * 0.5)  # scaled
theta_vec <- as.vector(t(theta))
theta_cov <- check_theta_recovery(fit_binom, theta_vec)

cat(sprintf("  Binomial result: innov=%s, theta_coverage=%.0f%%\n",
            if(innov_ok) "OK" else "MISS", theta_cov * 100))

# ============================================================================
# 4. Cumulative (ordinal) family
# ============================================================================
cat("\n--- 4. Cumulative (Ordinal) Family ---\n")

n_units <- 5; n_times <- 10
innov_true <- 0.4
cutpoints_true <- c(-1.5, -0.3, 0.5, 1.5)  # K=5 categories, 4 cutpoints

theta <- simulate_theta(n_units, n_times, innov_true)
theta <- theta * 0.3  # keep in moderate range relative to cutpoints

# Generate ordinal data
rordinal <- function(eta, cutpoints) {
  K <- length(cutpoints) + 1
  probs <- numeric(K)
  for (k in 1:K) {
    upper <- if (k < K) plogis(cutpoints[k] - eta) else 1
    lower <- if (k > 1) plogis(cutpoints[k-1] - eta) else 0
    probs[k] <- upper - lower
  }
  probs[probs < 0] <- 0  # numerical safety
  probs <- probs / sum(probs)
  sample(1:K, 1, prob = probs)
}

df_list <- list()
for (u in 1:n_units) {
  for (t in 1:n_times) {
    n_obs_per <- sample(8:15, 1)
    for (j in 1:n_obs_per) {
      y_val <- rordinal(theta[u, t], cutpoints_true)
      df_list[[length(df_list) + 1]] <- data.frame(
        unit = paste0("U", u), time = 2010 + t - 1,
        value = y_val
      )
    }
  }
}
df_ord <- do.call(rbind, df_list)

cat("  Preparing data...\n")
prep <- dlvm_prepare(df_ord, "unit", "time", "value",
                     family = "cumulative", mode = "individual",
                     compute_gq = TRUE)

cat("  Fitting model...\n")
fit_ord <- dlvm_fit(prep, chains = CHAINS, threads_per_chain = THREADS,
                    iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
                    show_messages = FALSE)

cat("  Checking recovery:\n")
innov_ok <- check_recovery(fit_ord, "innov", innov_true * 0.3)  # scaled

# Check cutpoints
cutpoint_covered <- 0
for (k in 1:(length(cutpoints_true))) {
  cp_name <- sprintf("cutpoints[%d]", k)
  if (check_recovery(fit_ord, cp_name, cutpoints_true[k])) {
    cutpoint_covered <- cutpoint_covered + 1
  }
}

theta_vec <- as.vector(t(theta))
theta_cov <- check_theta_recovery(fit_ord, theta_vec)

cat(sprintf("  Cumulative result: innov=%s, cutpoints=%d/%d, theta_coverage=%.0f%%\n",
            if(innov_ok) "OK" else "MISS",
            cutpoint_covered, length(cutpoints_true),
            theta_cov * 100))

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Simulation Summary ===\n")
cat("All four families tested with parameter recovery checks.\n")
cat("A 90% coverage rate >= 80% for theta is considered acceptable\n")
cat("(accounting for model approximations and finite samples).\n")
cat("=== Done ===\n")
