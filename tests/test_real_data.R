# dlvm/tests/test_real_data.R
# Integration test with real public longitudinal data
# Uses gapminder (built-in, no API calls) as primary test data.
# Requires CmdStan.

cat("  --- Real Data Integration Tests ---\n")

# Check CmdStan availability
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("  SKIPPED: cmdstanr not available\n")
  return(invisible(NULL))
}
tryCatch(cmdstanr::cmdstan_path(), error = function(e) {
  cat("  SKIPPED: CmdStan backend not installed\n")
  return(invisible(NULL))
})

# --- Try gapminder first (built-in data, no downloads) ---
use_gapminder <- requireNamespace("gapminder", quietly = TRUE)

if (use_gapminder) {
  cat("  Using gapminder data (life expectancy by country-year)\n")

  test_that("gapminder: fit and diagnostics pass", {
    data("gapminder", package = "gapminder")
    gm <- as.data.frame(gapminder::gapminder)

    # Select a manageable subset: 10 countries across all years
    top_countries <- c("United States", "United Kingdom", "France", "Japan",
                       "Brazil", "India", "South Africa", "Egypt", "Australia", "China")
    gm_sub <- gm[gm$country %in% top_countries, ]

    # Standardise life expectancy for numerical stability
    gm_sub$lifeExp_z <- (gm_sub$lifeExp - mean(gm_sub$lifeExp)) / sd(gm_sub$lifeExp)

    cat(sprintf("    Data: %d observations, %d countries, years %d-%d\n",
                nrow(gm_sub), length(unique(gm_sub$country)),
                min(gm_sub$year), max(gm_sub$year)))

    # Prepare and fit
    prep <- dlvm_prepare(gm_sub, "country", "year", "lifeExp_z")
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 500, iter_sampling = 1000, show_messages = FALSE)

    # Check diagnostics
    diag <- fit$fit$diagnostic_summary(quiet = TRUE)
    n_div <- sum(diag$num_divergent)
    cat(sprintf("    Divergences: %d\n", n_div))

    theta_summ <- fit$fit$summary(variables = "theta")
    max_rhat <- max(theta_summ$rhat, na.rm = TRUE)
    min_ess <- min(theta_summ$ess_bulk, na.rm = TRUE)
    cat(sprintf("    Max R-hat: %.4f, Min ESS: %.0f\n", max_rhat, min_ess))

    expect_equal(n_div, 0)
    expect_lt(max_rhat, 1.05)  # Relaxed for short chains
    expect_gt(min_ess, 100)    # Relaxed for short chains

    # Check that latent estimates are reasonable (bounded, not extreme)
    summ <- dlvm_summary(fit)
    expect_true(all(abs(summ$estimate) < 20, na.rm = TRUE))

    # Timing
    cat(sprintf("    Timing: %s\n", fmt_duration(fit$timing$total_secs)))
  })

} else {
  cat("  gapminder package not available. Install with: install.packages('gapminder')\n")

  # Fallback: use synthetic data that mimics real panel structure
  test_that("synthetic panel data: fit completes", {
    set.seed(42)
    n_countries <- 8
    years <- seq(1990, 2020, by = 1)
    df <- expand.grid(
      country = paste0("Country_", LETTERS[1:n_countries]),
      year = years
    )
    # Generate smooth trends with noise (simulating real panel data)
    df$value <- unlist(lapply(1:n_countries, function(i) {
      trend <- seq(from = rnorm(1, 0, 1), to = rnorm(1, 0, 1), length.out = length(years))
      trend + rnorm(length(years), 0, 0.3)
    }))

    prep <- dlvm_prepare(df, "country", "year", "value")
    fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)

    diag <- fit$fit$diagnostic_summary(quiet = TRUE)
    expect_equal(sum(diag$num_divergent), 0)
    cat(sprintf("    Timing: %s\n", fmt_duration(fit$timing$total_secs)))
  })
}

# --- Test irregular time spacing with real-like data ---
test_that("irregular time spacing works", {
  # Simulate data observed at irregular intervals (like survey waves)
  set.seed(99)
  df <- data.frame(
    country = rep(c("Alpha", "Beta"), each = 6),
    year = rep(c(1990, 1995, 2000, 2005, 2010, 2020), 2),
    score = rnorm(12)
  )
  prep <- dlvm_prepare(df, "country", "year", "score")

  # Verify delta_t captures the gaps
  dt <- prep$stan_data$delta_t
  expect_true(any(dt == 5))   # 5-year gaps
  expect_true(any(dt == 10))  # 10-year gap (2010 to 2020)

  fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                  iter_warmup = 300, iter_sampling = 500, show_messages = FALSE)

  # Wider intervals expected for larger time gaps
  summ <- dlvm_summary(fit)
  cat(sprintf("    Latent states recovered (%d states, irregular spacing)\n", nrow(summ)))
})

cat("  --- Real Data Tests Complete ---\n")
