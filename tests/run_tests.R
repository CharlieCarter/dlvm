#!/usr/bin/env Rscript
# dlvm/tests/run_tests.R
# Test runner — execute from the dlvm/ directory
#
# Usage:
#   Rscript tests/run_tests.R all           # Run all tests
#   Rscript tests/run_tests.R data_prep     # Data prep only (no Stan needed)
#   Rscript tests/run_tests.R simulation    # Simulation-based calibration
#   Rscript tests/run_tests.R real_data     # Real data integration
#   Rscript tests/run_tests.R edge_cases    # Edge case handling

args <- commandArgs(trailingOnly = TRUE)
suite <- if (length(args) > 0) args[1] else "all"

# Resolve dlvm root (robust: works from dlvm/ or dlvm/tests/)
test_root <- tryCatch({
  # When source()'d, get the script's directory
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
}, error = function(e) {
  # When run via Rscript, check common locations
  for (cand in c(".", "..", "dlvm")) {
    if (file.exists(file.path(cand, "stan", "dlvm.stan"))) {
      return(normalizePath(cand))
    }
  }
  normalizePath(".")
})
cat("DLVM root:", test_root, "\n")

# Load DLVM — use library() if installed, source() for development
if (requireNamespace("dlvm", quietly = TRUE)) {
  library(dlvm)
  # Only use installed package location if we can't find tests locally
  if (!dir.exists(file.path(test_root, "tests"))) {
    test_root <- system.file(package = "dlvm")
  }
} else {
  for (f in list.files(file.path(test_root, "R"), pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }
}

# Simple test framework
.test_results <- list(passed = 0L, failed = 0L, errors = character())

test_that <- function(desc, expr) {
  cat(sprintf("  TEST: %s ... ", desc))
  result <- tryCatch({
    expr
    .test_results$passed <<- .test_results$passed + 1L
    cat("PASS\n")
    TRUE
  }, error = function(e) {
    .test_results$failed <<- .test_results$failed + 1L
    .test_results$errors <<- c(.test_results$errors, sprintf("%s: %s", desc, e$message))
    cat(sprintf("FAIL: %s\n", e$message))
    FALSE
  })
  invisible(result)
}

expect_true <- function(x) if (!isTRUE(x)) stop("Expected TRUE, got ", deparse(x))
expect_false <- function(x) if (!isFALSE(x)) stop("Expected FALSE, got ", deparse(x))
expect_equal <- function(a, b, tol = 1e-8, ...) {
  args <- list(...)
  if ("tolerance" %in% names(args)) tol <- args$tolerance
  if (is.numeric(a) && is.numeric(b)) {
    if (any(abs(a - b) > tol, na.rm = TRUE)) stop(sprintf("Expected equal (tol=%g), got diff=%g", tol, max(abs(a - b), na.rm = TRUE)))
  } else if (!identical(a, b)) {
    stop(sprintf("Expected %s, got %s", deparse(b), deparse(a)))
  }
}
expect_error <- function(expr) {
  caught <- tryCatch({ expr; FALSE }, error = function(e) TRUE)
  if (!caught) stop("Expected an error but none was thrown")
}
expect_s3_class <- function(x, cls) {
  if (!inherits(x, cls)) stop(sprintf("Expected class '%s', got '%s'", cls, paste(class(x), collapse = ", ")))
}
expect_gt <- function(a, b) if (!(a > b)) stop(sprintf("Expected %g > %g", a, b))
expect_lt <- function(a, b) if (!(a < b)) stop(sprintf("Expected %g < %g", a, b))
expect_gte <- function(a, b) if (!(a >= b)) stop(sprintf("Expected %g >= %g", a, b))

# Run selected test suites
test_files <- list(
  data_prep  = "test_data_prep.R",
  simulation = "test_simulation.R",
  real_data  = "test_real_data.R",
  edge_cases = "test_edge_cases.R",
  cureton    = "test_cureton.R"
)

if (suite == "all") {
  suites_to_run <- names(test_files)
} else if (suite %in% names(test_files)) {
  suites_to_run <- suite
} else {
  stop("Unknown test suite: '", suite, "'. Use one of: ", paste(c("all", names(test_files)), collapse = ", "))
}

for (s in suites_to_run) {
  cat(sprintf("\n=== Running: %s ===\n", s))
  test_file <- file.path(test_root, "tests", test_files[[s]])
  if (!file.exists(test_file)) {
    cat(sprintf("  SKIPPED: %s not found\n", test_file))
    next
  }
  source(test_file, local = TRUE)
}

# Summary
cat(sprintf("\n=== Results: %d passed, %d failed ===\n",
            .test_results$passed, .test_results$failed))
if (length(.test_results$errors) > 0) {
  cat("Failures:\n")
  for (e in .test_results$errors) cat(sprintf("  ✗ %s\n", e))
}

quit(status = if (.test_results$failed > 0) 1 else 0, save = "no")
