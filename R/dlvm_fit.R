# dlvm/R/dlvm_fit.R
# Model fitting function wrapping CmdStanR
#
# This file provides the main fitting function that takes prepared data
# from dlvm_prepare() and returns a fitted model object.

#' Fit the Dynamic Robust Latent Variable Model
#'
#' Compiles the Stan model (if not already cached), then runs CmdStan's
#' NUTS sampler with within-chain parallelisation via reduce_sum.
#'
#' @param prepared_data A dlvm_prepared object from dlvm_prepare().
#' @param chains Integer; number of Markov chains (default: 4).
#' @param parallel_chains Integer; number of chains to run in parallel (default: 4).
#' @param threads_per_chain Integer; threads per chain for reduce_sum (default: 2).
#'   Total CPU usage = parallel_chains × threads_per_chain.
#'   On an 8-core machine, try parallel_chains=4, threads_per_chain=2.
#'   On HPC with 32 cores, try parallel_chains=4, threads_per_chain=8.
#' @param iter_warmup Integer; warmup iterations per chain (default: 1000).
#' @param iter_sampling Integer; sampling iterations per chain (default: 4000).
#' @param adapt_delta Numeric in (0,1); target acceptance rate (default: 0.95).
#'   Higher values reduce divergences but slow sampling.
#' @param max_treedepth Integer; maximum tree depth for NUTS (default: 10).
#' @param seed Integer; random seed for reproducibility (default: 90210).
#' @param show_messages Logical; show Stan sampling messages (default: TRUE).
#' @param ... Additional arguments passed to CmdStanModel$sample().
#'
#' @return A dlvm_fit S3 object containing:
#'   \item{fit}{The CmdStanMCMC object from CmdStanR}
#'   \item{metadata}{Metadata from the prepared data for post-processing}
#'   \item{timing}{Named list with compilation and sampling times}
#'
#' @export
#' @examples
#' prep <- dlvm_prepare(my_data, "country", "year", "score")
#' result <- dlvm_fit(prep, chains = 4, threads_per_chain = 2)
#' dlvm_summary(result)
#'
dlvm_fit <- function(prepared_data,
                     chains = 4L,
                     parallel_chains = 4L,
                     threads_per_chain = 2L,
                     iter_warmup = 1000L,
                     iter_sampling = 4000L,
                     adapt_delta = 0.95,
                     max_treedepth = 10L,
                     seed = 90210L,
                     show_messages = TRUE,
                     ...) {

  # --- Validate input ---
  if (!inherits(prepared_data, "dlvm_prepared")) {
    stop("[dlvm] 'prepared_data' must be a dlvm_prepared object from dlvm_prepare().")
  }

  # --- Compile model ---
  t_compile_start <- proc.time()["elapsed"]
  mod <- dlvm_compile(quiet = !show_messages)
  t_compile <- proc.time()["elapsed"] - t_compile_start

  # --- Run sampler ---
  use_threads <- dlvm_has_threads()
  effective_threads <- if (use_threads) threads_per_chain else 1L

  if (show_messages) {
    if (use_threads) {
      message(sprintf(
        "[dlvm] Sampling: %d chains x %d iter (%d warmup + %d sampling), %d threads/chain",
        chains, iter_warmup + iter_sampling, iter_warmup, iter_sampling, threads_per_chain
      ))
    } else {
      message(sprintf(
        "[dlvm] Sampling: %d chains x %d iter (%d warmup + %d sampling), no threading",
        chains, iter_warmup + iter_sampling, iter_warmup, iter_sampling
      ))
    }
    message(sprintf(
      "[dlvm] Data: %d latent states, %d observations",
      prepared_data$stan_data$n_states, prepared_data$stan_data$n_obs
    ))
  }

  t_sample_start <- proc.time()["elapsed"]

  # Build sample args (conditionally include threads_per_chain)
  sample_args <- list(
    data = prepared_data$stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    show_messages = show_messages
  )
  if (use_threads) {
    sample_args$threads_per_chain <- threads_per_chain
  }

  # Merge extra args
  extra_args <- list(...)
  sample_args <- c(sample_args, extra_args)

  # Ensure CmdStan's TBB is found at runtime (not Homebrew's)
  tbb_path <- dlvm_tbb_lib_path()
  old_dyld <- Sys.getenv("DYLD_LIBRARY_PATH", unset = NA)
  if (!is.null(tbb_path)) {
    if (is.na(old_dyld) || old_dyld == "") {
      Sys.setenv(DYLD_LIBRARY_PATH = tbb_path)
    } else {
      Sys.setenv(DYLD_LIBRARY_PATH = paste(tbb_path, old_dyld, sep = ":"))
    }
  }

  fit <- tryCatch(
    do.call(mod$sample, sample_args),
    finally = {
      # Restore original DYLD_LIBRARY_PATH
      if (is.na(old_dyld)) {
        Sys.unsetenv("DYLD_LIBRARY_PATH")
      } else {
        Sys.setenv(DYLD_LIBRARY_PATH = old_dyld)
      }
    }
  )

  t_sample <- proc.time()["elapsed"] - t_sample_start

  timing <- list(
    compilation_secs = t_compile,
    sampling_secs = t_sample,
    total_secs = t_compile + t_sample
  )

  if (show_messages) {
    message(sprintf(
      "[dlvm] Done. Compilation: %s, Sampling: %s, Total: %s",
      fmt_duration(timing$compilation_secs),
      fmt_duration(timing$sampling_secs),
      fmt_duration(timing$total_secs)
    ))
  }

  result <- structure(
    list(
      fit = fit,
      metadata = prepared_data$metadata,
      timing = timing
    ),
    class = "dlvm_fit"
  )

  return(result)
}

# ============================================================================
# Print method
# ============================================================================

#' @export
print.dlvm_fit <- function(x, ...) {
  m <- x$metadata
  t <- x$timing
  cat("DLVM fit:\n")
  cat(sprintf("  Units:        %d\n", m$n_units))
  cat(sprintf("  Time range:   %.1f – %.1f\n", min(m$times), max(m$times)))
  cat(sprintf("  Latent states: %d\n", m$n_states))
  cat(sprintf("  Observations:  %d (%s mode)\n", m$n_obs, m$mode))
  cat(sprintf("  Timing:        %s (compile: %s, sample: %s)\n",
              fmt_duration(t$total_secs),
              fmt_duration(t$compilation_secs),
              fmt_duration(t$sampling_secs)))

  # Quick diagnostics
  diag <- tryCatch({
    summ <- x$fit$summary(variables = c("innov", "sigma"))
    div <- x$fit$diagnostic_summary(quiet = TRUE)
    list(
      max_rhat = max(summ$rhat, na.rm = TRUE),
      min_ess = min(summ$ess_bulk, na.rm = TRUE),
      n_divergent = sum(div$num_divergent)
    )
  }, error = function(e) NULL)

  if (!is.null(diag)) {
    cat(sprintf("  Divergences:   %d\n", diag$n_divergent))
    cat(sprintf("  Max R-hat:     %.3f\n", diag$max_rhat))
    cat(sprintf("  Min ESS:       %.0f\n", diag$min_ess))
  }
  invisible(x)
}
