# dlvm/tests/test_families.R
# Tests for multi-family support in dlvm_prepare()
# These tests do NOT require Stan compilation.

library(data.table)

cat("  --- Family Dispatch Tests ---\n")

# ============================================================================
# Family resolution
# ============================================================================

test_that("aliases resolve to canonical family names", {
  # continuous → student_t
  df <- data.frame(country = rep("USA", 5), year = 2010:2014, score = rnorm(5))
  prep <- dlvm_prepare(df, "country", "year", "score", family = "continuous")
  expect_equal(prep$metadata$family, "student_t")

  # ordinal → cumulative
  df_ord <- data.frame(country = rep("USA", 10), year = rep(2010:2014, 2), rating = sample(1:5, 10, replace = TRUE))
  prep <- dlvm_prepare(df_ord, "country", "year", "rating", family = "ordinal", mode = "individual")
  expect_equal(prep$metadata$family, "cumulative")

  # binary → bernoulli
  df_bin <- data.frame(country = rep("USA", 10), year = rep(2010:2014, 2), outcome = sample(0:1, 10, replace = TRUE))
  prep <- dlvm_prepare(df_bin, "country", "year", "outcome", family = "binary", mode = "individual")
  expect_equal(prep$metadata$family, "bernoulli")
})

test_that("gaussian family is rejected with helpful error", {
  df <- data.frame(country = "USA", year = 2010, score = 1.5)
  expect_error(dlvm_prepare(df, "country", "year", "score", family = "gaussian"))
})

test_that("unknown family is rejected", {
  df <- data.frame(country = "USA", year = 2010, score = 1.5)
  expect_error(dlvm_prepare(df, "country", "year", "score", family = "poisson"))
})

# ============================================================================
# student_t family (backward compatibility)
# ============================================================================

test_that("student_t family works identically to default", {
  set.seed(42)
  df <- data.frame(
    country = rep(c("USA", "GBR"), each = 5),
    year = rep(2010:2014, 2),
    score = rnorm(10)
  )
  prep_default <- dlvm_prepare(df, "country", "year", "score")
  prep_explicit <- dlvm_prepare(df, "country", "year", "score", family = "student_t")

  expect_equal(prep_default$stan_data$n_states, prep_explicit$stan_data$n_states)
  expect_equal(prep_default$stan_data$n_obs, prep_explicit$stan_data$n_obs)
  expect_equal(prep_default$stan_data$y, prep_explicit$stan_data$y)
  expect_equal(prep_default$metadata$family, "student_t")
})

test_that("student_t family includes expected stan_data fields", {
  df <- data.frame(country = "USA", year = 2010:2014, score = rnorm(5))
  prep <- dlvm_prepare(df, "country", "year", "score", family = "student_t")
  expect_true("sigma" %in% names(prep$stan_data) == FALSE)  # sigma is a parameter, not data
  expect_true("nu_obs" %in% names(prep$stan_data))
  expect_true("has_se" %in% names(prep$stan_data))
  expect_true("se" %in% names(prep$stan_data))
  expect_true("obs_n" %in% names(prep$stan_data))
})

# ============================================================================
# cumulative family (ordinal)
# ============================================================================

test_that("cumulative family accepts integer ordinal data", {
  df <- data.frame(
    country = rep("USA", 20),
    year = rep(2010:2014, each = 4),
    rating = sample(1:5, 20, replace = TRUE)
  )
  prep <- dlvm_prepare(df, "country", "year", "rating", family = "cumulative", mode = "individual")
  expect_equal(prep$metadata$family, "cumulative")
  expect_equal(prep$stan_data$K, 5L)
  expect_true(is.integer(prep$stan_data$y))
  expect_true(all(prep$stan_data$y >= 1 & prep$stan_data$y <= 5))
})

test_that("cumulative family accepts factor data", {
  df <- data.frame(
    country = rep("USA", 10),
    year = rep(2010:2014, 2),
    rating = factor(c("low", "medium", "high", "low", "medium",
                       "high", "low", "medium", "high", "low"),
                    levels = c("low", "medium", "high"), ordered = TRUE)
  )
  prep <- dlvm_prepare(df, "country", "year", "rating", family = "cumulative", mode = "individual")
  expect_equal(prep$stan_data$K, 3L)
  expect_true(all(prep$stan_data$y %in% 1:3))
})

test_that("cumulative family detects K from data", {
  df <- data.frame(
    country = rep("USA", 15),
    year = rep(2010:2014, 3),
    rating = rep(c(1, 2, 7), 5)  # K = 7 (max value)
  )
  prep <- dlvm_prepare(df, "country", "year", "rating", family = "cumulative", mode = "individual")
  expect_equal(prep$stan_data$K, 7L)
})

test_that("cumulative family rejects means mode", {
  df <- data.frame(
    country = rep("USA", 10),
    year = rep(2010:2014, 2),
    rating = sample(1:5, 10, replace = TRUE)
  )
  expect_error(dlvm_prepare(df, "country", "year", "rating", family = "cumulative", mode = "means"))
})

test_that("cumulative family rejects non-integer numeric data", {
  df <- data.frame(country = "USA", year = 2010, rating = 1.5)
  expect_error(dlvm_prepare(df, "country", "year", "rating", family = "cumulative"))
})

test_that("cumulative family rejects K < 2", {
  df <- data.frame(country = rep("USA", 5), year = 2010:2014, rating = rep(1L, 5))
  expect_error(dlvm_prepare(df, "country", "year", "rating", family = "cumulative", mode = "individual"))
})

test_that("cumulative family does not include sigma-related fields", {
  df <- data.frame(
    country = rep("USA", 10), year = rep(2010:2014, 2),
    rating = sample(1:3, 10, replace = TRUE)
  )
  prep <- dlvm_prepare(df, "country", "year", "rating", family = "cumulative", mode = "individual")
  expect_true(!"nu_obs" %in% names(prep$stan_data))
  expect_true(!"has_se" %in% names(prep$stan_data))
  expect_true(!"se" %in% names(prep$stan_data))
  expect_true("K" %in% names(prep$stan_data))
})

# ============================================================================
# bernoulli family (binary)
# ============================================================================

test_that("bernoulli family accepts 0/1 data", {
  df <- data.frame(
    country = rep("USA", 20),
    year = rep(2010:2014, each = 4),
    outcome = sample(0:1, 20, replace = TRUE)
  )
  prep <- dlvm_prepare(df, "country", "year", "outcome", family = "bernoulli", mode = "individual")
  expect_equal(prep$metadata$family, "bernoulli")
  expect_true(is.integer(prep$stan_data$y))
  expect_true(all(prep$stan_data$y %in% 0:1))
})

test_that("bernoulli family rejects non-binary data", {
  df <- data.frame(country = "USA", year = 2010, outcome = 2L)
  expect_error(dlvm_prepare(df, "country", "year", "outcome", family = "bernoulli"))
})

test_that("bernoulli means mode auto-converts to binomial", {
  set.seed(42)
  df <- data.frame(
    country = rep("USA", 20),
    year = rep(2010:2014, each = 4),
    outcome = sample(0:1, 20, replace = TRUE)
  )
  prep <- dlvm_prepare(df, "country", "year", "outcome", family = "bernoulli", mode = "means")
  # Should have been auto-converted to binomial
  expect_equal(prep$metadata$family, "binomial")
  expect_true("n_trials" %in% names(prep$stan_data))
  expect_equal(prep$stan_data$n_obs, 5L)  # 5 year cells
  expect_true(all(prep$stan_data$n_trials == 4L))  # 4 obs per cell
})

test_that("bernoulli family minimal fields", {
  df <- data.frame(
    country = rep("USA", 10), year = rep(2010:2014, 2),
    outcome = sample(0:1, 10, replace = TRUE)
  )
  prep <- dlvm_prepare(df, "country", "year", "outcome", family = "bernoulli", mode = "individual")
  expect_true(!"sigma" %in% names(prep$stan_data))
  expect_true(!"nu_obs" %in% names(prep$stan_data))
  expect_true(!"n_trials" %in% names(prep$stan_data))
})

# ============================================================================
# binomial family
# ============================================================================

test_that("binomial family accepts successes and trials", {
  df <- data.frame(
    country = rep("USA", 5),
    year = 2010:2014,
    successes = c(5, 10, 15, 8, 12),
    trials = c(20, 30, 40, 25, 35)
  )
  prep <- dlvm_prepare(df, "country", "year", "successes",
                       trials_col = "trials", family = "binomial")
  expect_equal(prep$metadata$family, "binomial")
  expect_true("n_trials" %in% names(prep$stan_data))
  expect_equal(prep$stan_data$n_trials, c(20L, 30L, 40L, 25L, 35L))
  expect_equal(prep$stan_data$y, c(5L, 10L, 15L, 8L, 12L))
})

test_that("binomial family requires trials_col", {
  df <- data.frame(country = "USA", year = 2010, successes = 5)
  expect_error(dlvm_prepare(df, "country", "year", "successes", family = "binomial"))
})

test_that("binomial family rejects successes > trials", {
  df <- data.frame(country = "USA", year = 2010, successes = 25, trials = 20)
  expect_error(dlvm_prepare(df, "country", "year", "successes",
                            trials_col = "trials", family = "binomial"))
})

test_that("binomial family rejects missing trials column", {
  df <- data.frame(country = "USA", year = 2010, successes = 5)
  expect_error(dlvm_prepare(df, "country", "year", "successes",
                            trials_col = "nonexistent", family = "binomial"))
})

test_that("binomial family rejects negative successes", {
  df <- data.frame(country = "USA", year = 2010, successes = -1, trials = 20)
  expect_error(dlvm_prepare(df, "country", "year", "successes",
                            trials_col = "trials", family = "binomial"))
})

# ============================================================================
# Metadata
# ============================================================================

test_that("metadata includes family field", {
  df <- data.frame(country = "USA", year = 2010:2012, score = rnorm(3))
  prep <- dlvm_prepare(df, "country", "year", "score", family = "student_t")
  expect_true("family" %in% names(prep$metadata))
  expect_equal(prep$metadata$family, "student_t")
})

test_that("print method shows family info", {
  df <- data.frame(country = "USA", year = 2010:2012, score = rnorm(3))
  prep <- dlvm_prepare(df, "country", "year", "score", family = "student_t")
  out <- capture.output(print(prep))
  expect_true(any(grepl("Family", out)))
  expect_true(any(grepl("student_t", out)))
})

cat("  --- Family Dispatch Tests Complete ---\n")
