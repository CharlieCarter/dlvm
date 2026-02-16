# Tests for small-sample bias correction (Cureton/Holtzman)

# User-provided benchmarks:
# vals <- c(1:5) -> SD approx 1.682088 -> SEM = 1.682088/sqrt(5) = 0.752243
# vals <- c(2, 1.8) -> SD approx 0.1772454 -> SEM = 0.1772454/sqrt(2) = 0.12533

test_that("dlvm_sd_cureton calculates correct values", {
  vals1 <- c(1:5)
  sd_c1 <- dlvm_sd_cureton(vals1)
  expect_equal(sd_c1, 1.682088, tolerance = 1e-6)

  vals2 <- c(2, 1.8)
  sd_c2 <- dlvm_sd_cureton(vals2)
  expect_equal(sd_c2, 0.1772454, tolerance = 1e-6)
})

test_that("dlvm_se_cureton calculates correct values", {
  vals1 <- c(1:5)
  se_c1 <- dlvm_se_cureton(vals1)
  expect_equal(se_c1, 1.682088 / sqrt(5), tolerance = 1e-6)

  vals2 <- c(2, 1.8)
  se_c2 <- dlvm_se_cureton(vals2)
  expect_equal(se_c2, 0.1772454 / sqrt(2), tolerance = 1e-6)
})

test_that("dlvm_prepare uses correction by default", {
  # Create a small dataset that will be aggregated
  df <- data.frame(
    unit = c("A", "A"),
    year = c(2000, 2000),
    val = c(2, 1.8)
  )

  # Default behavior (se_correction = TRUE)
  prep_default <- dlvm_prepare(df, "unit", "year", "val", mode = "means")

  # The single observed row should have the corrected SE
  obs_se <- prep_default$metadata$obs_data$se[1]
  expected_se <- 0.1772454 / sqrt(2)

  expect_equal(obs_se, expected_se, tolerance = 1e-6)
})

test_that("dlvm_prepare can disable correction", {
  df <- data.frame(
    unit = c("A", "A"),
    year = c(2000, 2000),
    val = c(2, 1.8)
  )

  # Explicitly disable
  prep_uncorr <- dlvm_prepare(df, "unit", "year", "val", mode = "means", se_correction = FALSE)

  # Should match standard R sd()
  obs_se <- prep_uncorr$metadata$obs_data$se[1]
  expected_se <- sd(c(2, 1.8)) / sqrt(2) # 0.1 / 1.414... = 0.2 / 2 = 0.1
  # wait:
  # sd(c(2, 1.8)) -> sqrt(((2-1.9)^2 + (1.8-1.9)^2)/1) = sqrt(0.01+0.01) = sqrt(0.02) = 0.1414214
  # se = 0.1414214 / 1.414214 = 0.1

  expect_equal(obs_se, 0.1, tolerance = 1e-6)
})
