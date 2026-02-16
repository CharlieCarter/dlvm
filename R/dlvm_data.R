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
#' @param value_col Character; name of the column with observed values (numeric).
#' @param se_col Character or NULL; name of the column with known standard errors.
#' @param se_correction Logical; if TRUE (default), use small-sample bias correction (Cureton/Holtzman) for SE.
#' @param mode Character; one of "auto", "means", or "individual".
#' @param nu_obs Numeric; degrees of freedom for observation distribution (default: 4).
#' @param nu_state Numeric; degrees of freedom for state innovations (default: 4).
#' @param scale_state Numeric; scale for state innovation prior (default: 4).
#' @param grainsize Integer; reduce_sum grainsize (default: 1 = automatic).
#'
#' @return A list with components:
#'   \item{stan_data}{Named list ready to pass to CmdStanR's $sample()}
#'   \item{metadata}{List with unit labels, time values, and mapping info}
#' @export
dlvm_prepare <- function(data,
                         unit_col,
                         time_col,
                         value_col,
                         se_col = NULL,
                         se_correction = TRUE,
                         mode = c("auto", "means", "individual"),
                         nu_obs = 4,
                         nu_state = 4,
                         scale_state = 4,
                         grainsize = 1L,
                         compute_gq = FALSE) {

  mode <- match.arg(mode)

  # --- Input validation ---
  .validate_inputs(data, unit_col, time_col, value_col, se_col)

  # --- Build a clean data.frame with standardised columns ---
  # We construct from scratch rather than using := to avoid data.table
  # internal state issues that cause gforce segfaults.
  unit_vec  <- data[[unit_col]]
  time_vec  <- data[[time_col]]
  value_vec <- data[[value_col]]

  # Convert Date/POSIXct to numeric fractional years
  if (inherits(time_vec, "Date") || inherits(time_vec, "POSIXct")) {
    message("[dlvm] Converting Date/POSIXct time column to fractional years.")
    time_vec <- as.numeric(format(time_vec, "%Y")) +
      (as.numeric(format(time_vec, "%j")) - 1) / 365.25
  }

  # Handle SE column
  has_se <- !is.null(se_col) && se_col %in% names(data)
  se_vec <- if (has_se) {
    sv <- data[[se_col]]
    sv[is.na(sv) | sv <= 0] <- 0
    sv
  } else {
    rep(0, length(value_vec))
  }

  # Build work data.frame (base R to avoid data.table internals)
  df <- data.frame(
    unit  = unit_vec,
    time  = time_vec,
    value = value_vec,
    se    = se_vec,
    stringsAsFactors = FALSE
  )

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

  if (mode == "means" && max_rows_per_ut > 1) {
    message("[dlvm] mode='means' but multiple rows per unit-time found. Aggregating to means.")
    # Aggregate using base R split/lapply
    split_idx <- split(seq_len(nrow(df)), ut_key)
    agg_list <- lapply(split_idx, function(idx) {
      vals <- df$value[idx]
      ses  <- df$se[idx]
      n    <- length(idx)
      mean_val <- mean(vals, na.rm = TRUE)
      if (has_se) {
        agg_se <- sqrt(mean(ses^2, na.rm = TRUE)) / sqrt(n)
      } else {
        # Calculate SEM from sample SD
        if (n > 1) {
          if (se_correction) {
            agg_se <- dlvm_se_cureton(vals, na.rm = TRUE)
          } else {
            agg_se <- sd(vals, na.rm = TRUE) / sqrt(n)
          }
        } else {
          agg_se <- 0
        }
      }
      data.frame(
        unit  = df$unit[idx[1]],
        time  = df$time[idx[1]],
        value = mean_val,
        se    = agg_se,
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, agg_list)
    rownames(df) <- NULL
    if (!has_se) {
      has_se <- TRUE
      message("[dlvm] Computed SEM from within-group standard deviation during aggregation.")
    }
  }

  # --- Build latent state grid ---
  units <- sort(unique(df$unit))
  n_units <- length(units)
  all_times <- sort(unique(df$time))

  # Create full state grid: all units × all times (using expand.grid, not CJ)
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
  # Create a lookup key for state_grid
  sg_key <- paste(state_grid$unit, state_grid$time, sep = "|||")
  df_key <- paste(df$unit, df$time, sep = "|||")
  df$state_id <- state_grid$state_id[match(df_key, sg_key)]

  if (any(is.na(df$state_id))) {
    stop("[dlvm] Internal error: some observations could not be mapped to latent states.")
  }

  # --- Assemble Stan data list ---
  n_states <- nrow(state_grid)
  n_obs <- nrow(df)
  se_final <- if (has_se) df$se else rep(0, n_obs)

  stan_data <- list(
    n_states = n_states,
    n_obs = n_obs,
    n_units = n_units,
    y = df$value,
    obs_to_state = as.integer(df$state_id),
    has_se = as.integer(has_se && any(se_final > 0)),
    se = se_final,
    state_prev = as.integer(state_grid$prev_state_id),
    delta_t = state_grid$delta_t,
    nu_obs = nu_obs,
    nu_state = nu_state,
    scale_state = scale_state,
    grainsize = as.integer(grainsize),
    compute_gq = as.integer(compute_gq)
  )

  # --- Rename for user-facing output ---
  names(state_grid)[names(state_grid) == "unit"] <- unit_col
  names(state_grid)[names(state_grid) == "time"] <- time_col
  names(df)[names(df) == "unit"]  <- unit_col
  names(df)[names(df) == "time"]  <- time_col
  names(df)[names(df) == "value"] <- value_col

  # --- Build metadata for post-processing ---
  metadata <- list(
    unit_col = unit_col,
    time_col = time_col,
    value_col = value_col,
    se_col = se_col,
    mode = mode,
    n_states = n_states,
    n_obs = n_obs,
    n_units = n_units,
    units = units,
    times = all_times,
    state_grid = state_grid,
    obs_data = df
  )

  message(sprintf(
    "[dlvm] Prepared: %d units x %d time points = %d latent states, %d observations (%s mode).",
    n_units, length(all_times), n_states, n_obs, mode
  ))

  structure(
    list(stan_data = stan_data, metadata = metadata),
    class = "dlvm_prepared"
  )
}

# ============================================================================
# Input validation (internal)
# ============================================================================

.validate_inputs <- function(data, unit_col, time_col, value_col, se_col) {
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

  if (!is.null(se_col) && !se_col %in% names(data)) {
    stop("[dlvm] SE column '", se_col, "' not found in data.",
         "\n  Available columns: ", paste(names(data), collapse = ", "))
  }

  if (!is.numeric(data[[value_col]])) {
    stop("[dlvm] Value column '", value_col, "' must be numeric.")
  }

  if (!is.numeric(data[[time_col]]) &&
      !inherits(data[[time_col]], "Date") &&
      !inherits(data[[time_col]], "POSIXct")) {
    stop("[dlvm] Time column '", time_col, "' must be numeric, Date, or POSIXct.")
  }

  if (!is.null(se_col) && !is.numeric(data[[se_col]])) {
    stop("[dlvm] SE column '", se_col, "' must be numeric.")
  }
}

# ============================================================================
# Print method
# ============================================================================

#' @export
print.dlvm_prepared <- function(x, ...) {
  m <- x$metadata
  cat(sprintf("DLVM prepared data:\n"))
  cat(sprintf("  Units:        %d (%s)\n", m$n_units, m$unit_col))
  cat(sprintf("  Time range:   %.1f - %.1f (%s)\n",
              min(m$times), max(m$times), m$time_col))
  cat(sprintf("  Latent states: %d\n", m$n_states))
  cat(sprintf("  Observations:  %d (%s mode)\n", m$n_obs, m$mode))
  cat(sprintf("  Has SE:        %s\n", ifelse(x$stan_data$has_se, "yes", "no")))
  invisible(x)
}
