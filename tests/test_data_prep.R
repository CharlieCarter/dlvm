# dlvm/tests/test_data_prep.R
# Tests for dlvm_prepare() — data validation and preparation
# These tests do NOT require Stan compilation.

library(data.table)

cat("  --- Data Preparation Tests ---\n")

# --- Basic valid data ---
test_that("valid data produces dlvm_prepared object", {
  df <- data.frame(
    country = rep(c("USA", "GBR", "FRA"), each = 10),
    year = rep(2005:2014, 3),
    score = rnorm(30)
  )
  prep <- dlvm_prepare(df, "country", "year", "score")
  expect_s3_class(prep, "dlvm_prepared")
  expect_equal(prep$stan_data$n_units, 3L)
  expect_equal(prep$stan_data$n_obs, 30L)
  expect_equal(prep$stan_data$n_states, 30L)  # 3 units × 10 times
})

# --- Missing required columns ---
test_that("missing column throws error", {
  df <- data.frame(country = "USA", year = 2010, score = 1)
  expect_error(dlvm_prepare(df, "country", "time", "score"))
  expect_error(dlvm_prepare(df, "nation", "year", "score"))
  expect_error(dlvm_prepare(df, "country", "year", "value"))
})

# --- Non-numeric value column ---
test_that("non-numeric value column throws error", {
  df <- data.frame(country = "USA", year = 2010, score = "high")
  expect_error(dlvm_prepare(df, "country", "year", "score"))
})

# --- Empty data ---
test_that("empty data throws error", {
  df <- data.frame(country = character(), year = numeric(), score = numeric())
  expect_error(dlvm_prepare(df, "country", "year", "score"))
})

# --- SE column handling ---
test_that("standard errors are included when provided", {
  df <- data.frame(
    country = rep("USA", 5),
    year = 2010:2014,
    score = rnorm(5),
    se = runif(5, 0.1, 0.5)
  )
  prep <- dlvm_prepare(df, "country", "year", "score", se_col = "se")
  expect_equal(prep$stan_data$has_se, 1L)
  expect_true(all(prep$stan_data$se > 0))
})

test_that("missing SE column name throws error", {
  df <- data.frame(country = "USA", year = 2010, score = 1)
  expect_error(dlvm_prepare(df, "country", "year", "score", se_col = "se"))
})

# --- Mode auto-detection ---
test_that("auto mode detects means for 1:1 data", {
  df <- data.frame(
    country = rep("USA", 5),
    year = 2010:2014,
    score = rnorm(5)
  )
  prep <- dlvm_prepare(df, "country", "year", "score", mode = "auto")
  expect_equal(prep$metadata$mode, "means")
})

test_that("auto mode detects individual for multi-obs", {
  df <- data.frame(
    country = rep("USA", 15),
    year = rep(2010:2014, each = 3),
    score = rnorm(15)
  )
  prep <- dlvm_prepare(df, "country", "year", "score", mode = "auto")
  expect_equal(prep$metadata$mode, "individual")
})

# --- Aggregation when mode=means with multi-obs ---
test_that("means mode aggregates multi-obs correctly", {
  set.seed(42)
  df <- data.frame(
    country = rep("USA", 20),
    year = rep(2010:2014, each = 4),
    score = rnorm(20)
  )
  prep <- dlvm_prepare(df, "country", "year", "score", mode = "means")
  expect_equal(prep$stan_data$n_obs, 5L)  # 5 year-means
  expect_equal(prep$stan_data$has_se, 1L) # SEM computed from aggregation
})

# --- delta_t computation ---
test_that("delta_t is 1.0 for annual data", {
  df <- data.frame(
    country = rep("USA", 5),
    year = 2010:2014,
    score = rnorm(5)
  )
  prep <- dlvm_prepare(df, "country", "year", "score")
  # First state has delta_t = 0 (no predecessor)
  dt_vals <- prep$stan_data$delta_t
  expect_equal(dt_vals[1], 0)
  expect_true(all(dt_vals[2:5] == 1))
})

test_that("delta_t handles irregular spacing", {
  df <- data.frame(
    country = rep("USA", 4),
    year = c(2000, 2002, 2005, 2010),
    score = rnorm(4)
  )
  prep <- dlvm_prepare(df, "country", "year", "score")
  dt_vals <- prep$stan_data$delta_t
  expect_equal(dt_vals[1], 0)
  expect_equal(dt_vals[2], 2)
  expect_equal(dt_vals[3], 3)
  expect_equal(dt_vals[4], 5)
})

# --- Predecessor array ---
test_that("state_prev is correctly structured", {
  df <- data.frame(
    country = rep(c("USA", "GBR"), each = 3),
    year = rep(2010:2012, 2),
    score = rnorm(6)
  )
  prep <- dlvm_prepare(df, "country", "year", "score")
  sp <- prep$stan_data$state_prev
  # Each unit's first state has prev = 0
  # States should be: GBR-2010(1), GBR-2011(2), GBR-2012(3), USA-2010(4), USA-2011(5), USA-2012(6)
  expect_equal(sp[1], 0L)  # GBR-2010: no predecessor
  expect_equal(sp[2], 1L)  # GBR-2011: predecessor is GBR-2010
  expect_equal(sp[3], 2L)  # GBR-2012: predecessor is GBR-2011
  expect_equal(sp[4], 0L)  # USA-2010: no predecessor (different unit)
})

# --- NA handling ---
test_that("NA values are dropped with message", {
  df <- data.frame(
    country = rep("USA", 5),
    year = 2010:2014,
    score = c(1, NA, 3, NA, 5)
  )
  prep <- dlvm_prepare(df, "country", "year", "score")
  expect_equal(prep$stan_data$n_obs, 3L)
})

# --- obs_to_state mapping ---
test_that("obs_to_state correctly maps observations to states", {
  df <- data.frame(
    country = rep("USA", 10),
    year = rep(2010:2014, each = 2),
    score = rnorm(10)
  )
  prep <- dlvm_prepare(df, "country", "year", "score", mode = "individual")
  # Each pair of observations should map to the same state
  ots <- prep$stan_data$obs_to_state
  expect_equal(ots[1], ots[2])
  expect_equal(ots[3], ots[4])
  expect_true(ots[1] != ots[3])  # Different years → different states
})

# --- Print method ---
test_that("print method works", {
  df <- data.frame(country = "USA", year = 2010:2012, score = rnorm(3))
  prep <- dlvm_prepare(df, "country", "year", "score")
  out <- capture.output(print(prep))
  expect_true(length(out) > 0)
})

cat("  --- Data Preparation Tests Complete ---\n")
