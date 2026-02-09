# dlvm/tests/test_edge_cases.R
# Edge case handling: unusual inputs that should succeed or fail gracefully.
# Requires CmdStan for fitting tests.

cat("  --- Edge Case Tests ---\n")

# Check CmdStan availability
has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
  !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))

# --- Single unit ---
test_that("single unit works", {
  df <- data.frame(unit = rep("A", 10), time = 1:10, value = cumsum(rnorm(10, 0, 0.5)))
  prep <- dlvm_prepare(df, "unit", "time", "value")
  expect_equal(prep$stan_data$n_units, 1L)

  if (has_cmdstan) {
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)
    diag <- fit$fit$diagnostic_summary(quiet = TRUE)
    expect_equal(sum(diag$num_divergent), 0)
    cat("    Single unit: fit OK\n")
  }
})

# --- Two time points ---
test_that("two time points works", {
  df <- data.frame(unit = "A", time = c(1, 2), value = c(0.5, 0.8))
  prep <- dlvm_prepare(df, "unit", "time", "value")
  expect_equal(prep$stan_data$n_states, 2L)

  if (has_cmdstan) {
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)
    cat("    Two time points: fit OK\n")
  }
})

# --- Extreme outliers ---
test_that("extreme outliers handled by student-t", {
  set.seed(77)
  df <- data.frame(
    unit = rep("A", 20),
    time = 1:20,
    value = c(rnorm(18, 0, 0.5), 50, -50)  # Two massive outliers
  )
  prep <- dlvm_prepare(df, "unit", "time", "value")

  if (has_cmdstan) {
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)
    diag <- fit$fit$diagnostic_summary(quiet = TRUE)
    # Student-t should handle this without many divergences
    expect_lt(sum(diag$num_divergent), 10)  # Some acceptable with extreme data
    cat(sprintf("    Extreme outliers: %d divergences\n", sum(diag$num_divergent)))
  }
})

# --- Many units, short series ---
test_that("many units with short time series", {
  set.seed(33)
  df <- data.frame(
    unit = rep(paste0("U", 1:50), each = 3),
    time = rep(1:3, 50),
    value = rnorm(150)
  )
  prep <- dlvm_prepare(df, "unit", "time", "value")
  expect_equal(prep$stan_data$n_units, 50L)
  expect_equal(prep$stan_data$n_states, 150L)

  if (has_cmdstan) {
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)
    cat(sprintf("    50 units × 3 times: %s\n", fmt_duration(fit$timing$total_secs)))
  }
})

# --- All-zero SEs ---
test_that("all-zero SEs treated as no-SE mode", {
  df <- data.frame(
    unit = rep("A", 5),
    time = 1:5,
    value = rnorm(5),
    se = rep(0, 5)
  )
  prep <- dlvm_prepare(df, "unit", "time", "value", se_col = "se")
  # Should detect that all SEs are zero and set has_se = 0
  expect_equal(prep$stan_data$has_se, 0L)
})

# --- Mixed SEs (some zero, some positive) ---
test_that("mixed SEs work correctly", {
  df <- data.frame(
    unit = rep("A", 5),
    time = 1:5,
    value = rnorm(5),
    se = c(0.1, 0, 0.3, 0, 0.2)
  )
  prep <- dlvm_prepare(df, "unit", "time", "value", se_col = "se")
  expect_equal(prep$stan_data$has_se, 1L)
  # Zero SEs should be kept as-is (handled in Stan via sqrt(sigma^2 + 0))
})

# --- Summary and plot functions ---
if (has_cmdstan) {
  test_that("dlvm_summary returns valid output", {
    df <- data.frame(unit = rep(c("A", "B"), each = 5), time = rep(1:5, 2), value = rnorm(10))
    prep <- dlvm_prepare(df, "unit", "time", "value")
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)
    summ <- dlvm_summary(fit)
    expect_true(nrow(summ) == 10)
    expect_true(all(c("estimate", "lower", "upper") %in% names(summ)))
    expect_true(all(summ$lower <= summ$estimate, na.rm = TRUE))
    expect_true(all(summ$estimate <= summ$upper, na.rm = TRUE))
  })

  test_that("dlvm_diagnostics runs without error", {
    df <- data.frame(unit = rep("A", 5), time = 1:5, value = rnorm(5))
    prep <- dlvm_prepare(df, "unit", "time", "value")
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)
    diag_out <- dlvm_diagnostics(fit)
    expect_true(!is.null(diag_out$theta_max_rhat))
  })

  test_that("dlvm_plot returns ggplot", {
    df <- data.frame(unit = rep(c("A", "B"), each = 5), time = rep(1:5, 2), value = rnorm(10))
    prep <- dlvm_prepare(df, "unit", "time", "value")
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)
    p <- dlvm_plot(fit)
    expect_s3_class(p, "ggplot")
  })
}

cat("  --- Edge Case Tests Complete ---\n")
