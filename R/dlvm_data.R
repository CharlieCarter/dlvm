# dlvm/R/dlvm_data.R
# Data validation and preparation for the DLVM
#
# This file provides the main data preparation function that transforms
# user-supplied panel data into the Stan data list required by dlvm.stan.
#
# NOTE: All aggregation and counting uses base R to avoid data.table gforce
# segfaults observed on macOS ARM64 with certain data.table versions.
# data.table is used only for merge, sort, and shift operations.

# ============================================================================
# Main entry point
# ============================================================================

#' Prepare panel data for the Dynamic Robust Latent Variable Model
#'
#' Validates input data, constructs the latent state grid, computes time gaps,
#' builds the state adjacency structure, and returns a Stan-ready data list.
#'
#' @param data A data.frame or data.table with at least unit, time, and value columns.
#' @param unit_col Character; name of the column identifying units (e.g., "country").
#' @param time_col Character; name of the column with time values (numeric or Date).
#' @param value_col Character; name of the column with observed values.
#'   For student_t: numeric. For cumulative: integer 1..K or factor.
#'   For bernoulli: integer 0/1. For binomial: integer (successes).
#' @param se_col Character or NULL; column with known standard errors (student_t only).
#' @param trials_col Character or NULL; column with trial counts (binomial family only).
#' @param se_correction Logical; if TRUE (default), use small-sample bias correction for SE.
#' @param family Character; observation distribution family (default: "student_t").
#'   Primary values: "student_t", "cumulative", "bernoulli", "binomial".
#'   Aliases: "continuous" (→ student_t), "ordinal" (→ cumulative), "binary" (→ bernoulli).
#' @param mode Character; one of "auto", "means", or "individual".
#' @param nu_obs Numeric or NULL; degrees of freedom for observation distribution
#'   (student_t family only). If NULL (default), a mode-dependent value is used:
#'   4 for individual mode (doc-level scores can be genuine outliers) and 30 for
#'   means mode (by CLT, cell means are approximately normal, and the SEM already
#'   captures precision heterogeneity across cells).
#' @param nu_state Numeric; degrees of freedom for state innovations (default: 4).
#' @param scale_state Numeric; scale for state innovation prior (default: 4).
#' @param grainsize Integer; reduce_sum grainsize (default: 1 = automatic).
#' @param compute_gq Logical; compute generated quantities (default: FALSE).
#'
#' @return A list with components:
#'   \item{stan_data}{Named list ready to pass to CmdStanR's $sample()}
#'   \item{metadata}{List with unit labels, time values, family, and mapping info}
#' @export
dlvm_prepare <- function(data,
                         unit_col,
                         time_col,
                         value_col,
                         se_col = NULL,
                         trials_col = NULL,
                         se_correction = TRUE,
                         family = "student_t",
                         mode = c("auto", "means", "individual"),
                         nu_obs = NULL,
                         nu_state = 4,
                         scale_state = 4,
                         grainsize = 1L,
                         compute_gq = FALSE) {

  mode <- match.arg(mode)
  family <- .resolve_family(family)

  # --- Input validation ---
  .validate_inputs(data, unit_col, time_col, value_col, se_col, trials_col, family)

  # --- Build a clean data.frame with standardised columns ---
  unit_vec  <- data[[unit_col]]
  time_vec  <- data[[time_col]]
  value_vec <- data[[value_col]]

  # For cumulative family: convert factor to integer levels if needed
  if (family == "cumulative" && is.factor(value_vec)) {
    value_vec <- as.integer(value_vec)
  }

  # Convert Date/POSIXct to numeric fractional years
  if (inherits(time_vec, "Date") || inherits(time_vec, "POSIXct")) {
    message("[dlvm] Converting Date/POSIXct time column to fractional years.")
    time_vec <- as.numeric(format(time_vec, "%Y")) +
      (as.numeric(format(time_vec, "%j")) - 1) / 365.25
  }

  # Handle SE column (student_t only)
  has_se <- !is.null(se_col) && se_col %in% names(data) && family == "student_t"
  se_vec <- if (has_se) {
    sv <- data[[se_col]]
    sv[is.na(sv) | sv <= 0] <- 0
    sv
  } else {
    rep(0, length(value_vec))
  }

  # Handle trials column (binomial only)
  trials_vec <- if (family == "binomial") {
    data[[trials_col]]
  } else {
    NULL
  }

  # Build work data.frame
  df <- data.frame(
    unit  = unit_vec,
    time  = time_vec,
    value = value_vec,
    se    = se_vec,
    stringsAsFactors = FALSE
  )
  if (!is.null(trials_vec)) {
    df$trials <- trials_vec
  }

  # Remove rows where value is NA
  n_before <- nrow(df)
  df <- df[!is.na(df$value), , drop = FALSE]
  n_dropped <- n_before - nrow(df)
  if (n_dropped > 0) {
    message(sprintf("[dlvm] Dropped %d rows with NA values.", n_dropped))
  }

  # Check for empty/all-NA units
  unit_counts <- table(df$unit)
  empty_units <- names(unit_counts[unit_counts == 0])
  if (length(empty_units) > 0) {
    warning(sprintf("[dlvm] Removing %d unit(s) with no valid observations: %s",
                    length(empty_units), paste(empty_units, collapse = ", ")))
    df <- df[!df$unit %in% empty_units, , drop = FALSE]
  }

  # --- Determine mode ---
  ut_key <- paste(df$unit, df$time, sep = "|||")
  rows_per_ut <- table(ut_key)
  max_rows_per_ut <- max(rows_per_ut)

  if (mode == "auto") {
    if (max_rows_per_ut > 1) {
      mode <- "individual"
      message(sprintf("[dlvm] Auto-detected mode: 'individual' (max %d obs per unit-time).", max_rows_per_ut))
    } else {
      mode <- "means"
      message("[dlvm] Auto-detected mode: 'means' (1 obs per unit-time).")
    }
  }

  # --- Resolve mode-dependent nu_obs ---
  # For student_t family, individual doc-level scores can be genuine outliers
  # (unusual use of a term, atypical context), so nu_obs = 4 gives appropriate
  # heavy-tailed robustness. For means mode, the Central Limit Theorem implies
  # cell means are approximately normal even if individual scores are not, and
  # the SEM (se_col / sqrt(n) scaling) already captures per-cell precision
  # heterogeneity — so heavy tails are redundant and weaken informativeness.
  # nu_obs = 30 is approximately Gaussian while retaining minimal tail protection.
  if (is.null(nu_obs)) {
    if (family == "student_t") {
      nu_obs <- if (mode == "means") 30 else 4
      message(sprintf("[dlvm] Using mode-dependent nu_obs = %g (mode: %s).", nu_obs, mode))
    } else {
      nu_obs <- 4  # non-student_t families don't use nu_obs, but set a default
    }
  }

  # --- Mode restrictions by family ---
  if (mode == "means" && family == "cumulative") {
    stop("[dlvm] mode='means' is not supported for the cumulative (ordinal) family. ",
         "Ordinal observations cannot be meaningfully averaged. Use mode='individual'.")
  }

  # --- Bernoulli in means mode → auto-convert to binomial ---
  if (mode == "means" && family == "bernoulli" && max_rows_per_ut > 1) {
    message("[dlvm] Auto-converting bernoulli in means mode to binomial (summing successes and trials per cell).")
    split_idx <- split(seq_len(nrow(df)), ut_key)
    agg_list <- lapply(split_idx, function(idx) {
      data.frame(
        unit    = df$unit[idx[1]],
        time    = df$time[idx[1]],
        value   = sum(df$value[idx]),  # successes
        se      = 0,
        trials  = length(idx),         # trials
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, agg_list)
    rownames(df) <- NULL
    family <- "binomial"
    # Recalculate keys after aggregation to prevent double-aggregation
    ut_key <- paste(df$unit, df$time, sep = "|||")
    rows_per_ut <- table(ut_key)
    max_rows_per_ut <- max(rows_per_ut)
    message(sprintf("[dlvm] Converted to binomial: %d cells.", nrow(df)))
  }

  # --- Binomial aggregation: sum successes and trials per cell ---
  if (mode == "means" && family == "binomial" && max_rows_per_ut > 1) {
    ut_key <- paste(df$unit, df$time, sep = "|||")
    split_idx <- split(seq_len(nrow(df)), ut_key)
    agg_list <- lapply(split_idx, function(idx) {
      data.frame(
        unit    = df$unit[idx[1]],
        time    = df$time[idx[1]],
        value   = sum(df$value[idx]),
        se      = 0,
        trials  = sum(df$trials[idx]),
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, agg_list)
    rownames(df) <- NULL
    message(sprintf("[dlvm] Aggregated binomial data: %d cells.", nrow(df)))
  }

  # --- Student-t aggregation (existing logic) ---
  if (mode == "means" && family == "student_t" && max_rows_per_ut > 1) {
    message("[dlvm] mode='means' but multiple rows per unit-time found. Aggregating to means.")

    sigma_hat <- sd(df$value, na.rm = TRUE)
    n1_count <- sum(vapply(split(seq_len(nrow(df)), ut_key), length, integer(1)) == 1L)

    split_idx <- split(seq_len(nrow(df)), ut_key)
    agg_list <- lapply(split_idx, function(idx) {
      vals <- df$value[idx]
      ses  <- df$se[idx]
      n    <- length(idx)
      mean_val <- mean(vals, na.rm = TRUE)
      if (has_se) {
        agg_se <- sqrt(mean(ses^2, na.rm = TRUE)) / sqrt(n)
      } else {
        if (n > 1) {
          if (se_correction) {
            agg_se <- dlvm_se_cureton(vals, na.rm = TRUE)
          } else {
            agg_se <- sd(vals, na.rm = TRUE) / sqrt(n)
          }
        } else {
          agg_se <- sigma_hat
        }
      }
      data.frame(
        unit  = df$unit[idx[1]],
        time  = df$time[idx[1]],
        value = mean_val,
        se    = agg_se,
        obs_n = n,
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, agg_list)
    rownames(df) <- NULL
    if (!has_se) {
      has_se <- TRUE
      message(sprintf(
        "[dlvm] Computed SEM from within-group SD during aggregation (%d of %d cells used sigma_hat=%.4f proxy for n=1).",
        n1_count, length(split_idx), sigma_hat
      ))
    }
  }

  # --- Build latent state grid ---
  units <- sort(unique(df$unit))
  n_units <- length(units)
  all_times <- sort(unique(df$time))

  state_grid <- expand.grid(unit = units, time = all_times,
                            stringsAsFactors = FALSE)
  state_grid <- state_grid[order(state_grid$unit, state_grid$time), , drop = FALSE]
  rownames(state_grid) <- NULL
  state_grid$state_id <- seq_len(nrow(state_grid))

  # --- Compute predecessor indices and delta_t ---
  state_grid$prev_state_id <- 0L
  state_grid$delta_t <- 0.0
  for (u in units) {
    idx <- which(state_grid$unit == u)
    if (length(idx) > 1) {
      state_grid$prev_state_id[idx[-1]] <- state_grid$state_id[idx[-length(idx)]]
      state_grid$delta_t[idx[-1]] <- diff(state_grid$time[idx])
    }
  }

  # --- Map observations to states ---
  sg_key <- paste(state_grid$unit, state_grid$time, sep = "|||")
  df_key <- paste(df$unit, df$time, sep = "|||")
  df$state_id <- state_grid$state_id[match(df_key, sg_key)]

  if (any(is.na(df$state_id))) {
    stop("[dlvm] Internal error: some observations could not be mapped to latent states.")
  }

  # --- Assemble Stan data list (family-specific) ---
  n_states <- nrow(state_grid)
  n_obs <- nrow(df)

  # Shared data across all families
  stan_data <- list(
    n_states    = n_states,
    n_obs       = n_obs,
    n_units     = n_units,
    obs_to_state = as.integer(df$state_id),
    state_prev  = as.integer(state_grid$prev_state_id),
    delta_t     = state_grid$delta_t,
    nu_state    = nu_state,
    scale_state = scale_state,
    grainsize   = as.integer(grainsize),
    compute_gq  = as.integer(compute_gq)
  )

  # Family-specific data
  if (family == "student_t") {
    se_final <- if (has_se) df$se else rep(0, n_obs)
    obs_n_vec <- if ("obs_n" %in% names(df)) as.integer(df$obs_n) else rep(1L, n_obs)
    stan_data$y      <- df$value
    stan_data$has_se  <- as.integer(has_se && any(se_final > 0))
    stan_data$se      <- se_final
    stan_data$obs_n   <- obs_n_vec
    stan_data$nu_obs  <- nu_obs
  } else if (family == "cumulative") {
    K <- as.integer(max(df$value))
    stan_data$y <- as.integer(df$value)
    stan_data$K <- K
  } else if (family == "bernoulli") {
    stan_data$y <- as.integer(df$value)
  } else if (family == "binomial") {
    stan_data$y        <- as.integer(df$value)
    stan_data$n_trials <- as.integer(df$trials)
  }

  # --- Rename for user-facing output ---
  names(state_grid)[names(state_grid) == "unit"] <- unit_col
  names(state_grid)[names(state_grid) == "time"] <- time_col
  names(df)[names(df) == "unit"]  <- unit_col
  names(df)[names(df) == "time"]  <- time_col
  names(df)[names(df) == "value"] <- value_col

  # --- Build metadata for post-processing ---
  metadata <- list(
    unit_col   = unit_col,
    time_col   = time_col,
    value_col  = value_col,
    se_col     = se_col,
    trials_col = trials_col,
    family     = family,
    mode       = mode,
    n_states   = n_states,
    n_obs      = n_obs,
    n_units    = n_units,
    units      = units,
    times      = all_times,
    state_grid = state_grid,
    obs_data   = df
  )

  message(sprintf(
    "[dlvm] Prepared: %d units x %d time points = %d latent states, %d observations (family: %s, mode: %s).",
    n_units, length(all_times), n_states, n_obs, family, mode
  ))

  structure(
    list(stan_data = stan_data, metadata = metadata),
    class = "dlvm_prepared"
  )
}

# ============================================================================
# Input validation (internal)
# ============================================================================

.validate_inputs <- function(data, unit_col, time_col, value_col, se_col, trials_col, family) {
  if (!is.data.frame(data)) {
    stop("[dlvm] 'data' must be a data.frame or data.table.")
  }
  if (nrow(data) == 0) {
    stop("[dlvm] 'data' has zero rows.")
  }

  required_cols <- c(unit_col, time_col, value_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("[dlvm] Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\n  Available columns: ", paste(names(data), collapse = ", "))
  }

  if (!is.numeric(data[[time_col]]) &&
      !inherits(data[[time_col]], "Date") &&
      !inherits(data[[time_col]], "POSIXct")) {
    stop("[dlvm] Time column '", time_col, "' must be numeric, Date, or POSIXct.")
  }

  # Family-specific value column validation
  if (family == "student_t") {
    if (!is.numeric(data[[value_col]])) {
      stop("[dlvm] Value column '", value_col, "' must be numeric for family='student_t'.")
    }
    if (!is.null(se_col) && !se_col %in% names(data)) {
      stop("[dlvm] SE column '", se_col, "' not found in data.",
           "\n  Available columns: ", paste(names(data), collapse = ", "))
    }
    if (!is.null(se_col) && !is.numeric(data[[se_col]])) {
      stop("[dlvm] SE column '", se_col, "' must be numeric.")
    }

  } else if (family == "cumulative") {
    vals <- data[[value_col]]
    if (!is.integer(vals) && !is.factor(vals)) {
      # Allow numeric that is integer-valued
      if (is.numeric(vals) && all(vals == floor(vals), na.rm = TRUE)) {
        # OK — will convert
      } else {
        stop("[dlvm] Value column '", value_col,
             "' must be integer or factor for family='cumulative' (ordinal data).")
      }
    }
    numeric_vals <- if (is.factor(vals)) as.integer(vals) else as.integer(vals)
    if (min(numeric_vals, na.rm = TRUE) < 1) {
      stop("[dlvm] Ordinal values must be >= 1 (got min = ",
           min(numeric_vals, na.rm = TRUE), ").")
    }
    K <- max(numeric_vals, na.rm = TRUE)
    if (K < 2) {
      stop("[dlvm] Ordinal data must have at least 2 categories (K = ", K, ").")
    }

  } else if (family == "bernoulli") {
    vals <- data[[value_col]]
    if (!all(vals %in% c(0L, 1L, 0, 1, NA))) {
      stop("[dlvm] Value column '", value_col,
           "' must contain only 0 and 1 for family='bernoulli'.")
    }

  } else if (family == "binomial") {
    if (is.null(trials_col)) {
      stop("[dlvm] 'trials_col' is required for family='binomial'.")
    }
    if (!trials_col %in% names(data)) {
      stop("[dlvm] Trials column '", trials_col, "' not found in data.",
           "\n  Available columns: ", paste(names(data), collapse = ", "))
    }
    vals <- data[[value_col]]
    trials <- data[[trials_col]]
    if (!is.numeric(vals) || !is.numeric(trials)) {
      stop("[dlvm] Value and trials columns must be numeric for family='binomial'.")
    }
    if (any(vals > trials, na.rm = TRUE)) {
      stop("[dlvm] Some values exceed trial counts (successes > trials).")
    }
    if (any(vals < 0, na.rm = TRUE) || any(trials < 1, na.rm = TRUE)) {
      stop("[dlvm] Values must be >= 0 and trials >= 1 for family='binomial'.")
    }
  }
}

# ============================================================================
# Print method
# ============================================================================

#' @export
print.dlvm_prepared <- function(x, ...) {
  m <- x$metadata
  cat(sprintf("DLVM prepared data:\n"))
  cat(sprintf("  Family:       %s\n", m$family))
  cat(sprintf("  Units:        %d (%s)\n", m$n_units, m$unit_col))
  cat(sprintf("  Time range:   %.1f - %.1f (%s)\n",
              min(m$times), max(m$times), m$time_col))
  cat(sprintf("  Latent states: %d\n", m$n_states))
  cat(sprintf("  Observations:  %d (%s mode)\n", m$n_obs, m$mode))
  if (m$family == "student_t") {
    cat(sprintf("  Has SE:        %s\n", ifelse(x$stan_data$has_se, "yes", "no")))
  }
  if (m$family == "cumulative") {
    cat(sprintf("  Categories:    %d\n", x$stan_data$K))
  }
  invisible(x)
}
