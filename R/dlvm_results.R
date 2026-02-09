# dlvm/R/dlvm_results.R
# Results extraction, diagnostics, and visualisation
#
# This file provides functions to extract parameter estimates,
# run diagnostics, compute LOO-CV, and generate plots from a dlvm_fit object.

# ============================================================================
# Summary extraction
# ============================================================================

#' Extract latent state estimates as a data.table
#'
#' Returns the posterior median and credible intervals for each latent state,
#' merged with the unit and time labels from the original data.
#'
#' @param fit A dlvm_fit object from dlvm_fit().
#' @param prob Numeric; width of the credible interval (default: 0.95 for 95% CI).
#' @return A data.table with columns: unit, time, estimate, lower, upper, and
#'   optionally the observed values and SEs.
#' @export
dlvm_summary <- function(fit, prob = 0.95) {
  if (!inherits(fit, "dlvm_fit")) {
    stop("[dlvm] 'fit' must be a dlvm_fit object from dlvm_fit().")
  }

  alpha <- (1 - prob) / 2
  probs <- c(alpha, 0.5, 1 - alpha)

  # Extract theta draws
  theta_draws <- fit$fit$draws("theta", format = "matrix")

  # Compute quantiles
  estimates <- t(apply(theta_draws, 2, quantile, probs = probs))
  colnames(estimates) <- c("lower", "estimate", "upper")

  # Build output table from state_grid (use as plain data.frame)
  state_grid <- as.data.frame(fit$metadata$state_grid)
  result <- cbind(state_grid, as.data.frame(estimates))

  # Rename to user's original column names
  unit_col <- fit$metadata$unit_col
  time_col <- fit$metadata$time_col

  # Add observed data (means or raw, depending on mode)
  obs_data <- as.data.frame(fit$metadata$obs_data)
  if (!is.null(obs_data) && "state_id" %in% names(obs_data)) {
    value_col <- fit$metadata$value_col
    # Aggregate observations per state using base R
    obs_agg <- aggregate(
      obs_data[[value_col]],
      by = list(state_id = obs_data$state_id),
      FUN = function(v) c(mean = mean(v, na.rm = TRUE), n = length(v))
    )
    obs_agg <- data.frame(
      state_id = obs_agg$state_id,
      obs_mean = obs_agg$x[, "mean"],
      obs_n    = as.integer(obs_agg$x[, "n"])
    )

    # SE aggregation if available
    se_col_name <- fit$metadata$se_col
    if (!is.null(se_col_name) && se_col_name %in% names(obs_data)) {
      se_agg <- aggregate(
        obs_data[[se_col_name]],
        by = list(state_id = obs_data$state_id),
        FUN = mean, na.rm = TRUE
      )
      names(se_agg) <- c("state_id", "obs_se")
      obs_agg <- merge(obs_agg, se_agg, by = "state_id")
    }

    result <- merge(result, obs_agg, by = "state_id", all.x = TRUE)
  }

  # Clean up internal columns
  drop_cols <- intersect(c("prev_state_id", "prev_time"), names(result))
  if (length(drop_cols) > 0) result <- result[, !names(result) %in% drop_cols]

  # Sort by unit then time
  result <- result[order(result[[unit_col]], result[[time_col]]), ]
  rownames(result) <- NULL

  class(result) <- c("dlvm_summary", "data.frame")
  return(result)
}

# ============================================================================
# Scalar parameter summary
# ============================================================================

#' Extract scalar parameter estimates
#'
#' @param fit A dlvm_fit object.
#' @return A data.table with parameter name, estimate, CI, Rhat, and ESS.
#' @export
dlvm_parameters <- function(fit) {
  if (!inherits(fit, "dlvm_fit")) {
    stop("[dlvm] 'fit' must be a dlvm_fit object.")
  }
  summ <- fit$fit$summary(variables = c("innov", "sigma"))
  data.table::as.data.table(summ)
}

# ============================================================================
# Diagnostics
# ============================================================================

#' Run comprehensive MCMC diagnostics
#'
#' @param fit A dlvm_fit object from dlvm_fit().
#' @return A list with diagnostic information, printed to console.
#' @export
dlvm_diagnostics <- function(fit) {
  if (!inherits(fit, "dlvm_fit")) {
    stop("[dlvm] 'fit' must be a dlvm_fit object.")
  }

  # CmdStan diagnostics
  diag <- fit$fit$diagnostic_summary(quiet = TRUE)

  # Scalar parameter summaries
  scalar_summ <- fit$fit$summary(variables = c("innov", "sigma"))

  # Theta summaries (just Rhat and ESS)
  theta_summ <- fit$fit$summary(variables = "theta")
  max_rhat <- max(theta_summ$rhat, na.rm = TRUE)
  min_ess_bulk <- min(theta_summ$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(theta_summ$ess_tail, na.rm = TRUE)
  pct_rhat_ok <- mean(theta_summ$rhat < 1.01, na.rm = TRUE) * 100

  result <- list(
    n_divergent = sum(diag$num_divergent),
    n_max_treedepth = sum(diag$num_max_treedepth),
    ebfmi = diag$ebfmi,
    scalar_params = scalar_summ,
    theta_max_rhat = max_rhat,
    theta_min_ess_bulk = min_ess_bulk,
    theta_min_ess_tail = min_ess_tail,
    theta_pct_rhat_ok = pct_rhat_ok
  )

  # Print summary
  cat("=== DLVM Diagnostics ===\n")
  cat(sprintf("Divergences:       %d\n", result$n_divergent))
  cat(sprintf("Max treedepth:     %d hits\n", result$n_max_treedepth))
  cat(sprintf("E-BFMI:            %s\n", paste(round(unlist(result$ebfmi), 2), collapse = ", ")))
  cat("\n--- Scalar Parameters ---\n")
  print(scalar_summ[, c("variable", "mean", "median", "q5", "q95", "rhat", "ess_bulk")])
  cat(sprintf("\n--- Latent States (theta) ---\n"))
  cat(sprintf("Max R-hat:         %.4f\n", max_rhat))
  cat(sprintf("Min ESS (bulk):    %.0f\n", min_ess_bulk))
  cat(sprintf("Min ESS (tail):    %.0f\n", min_ess_tail))
  cat(sprintf("R-hat < 1.01:      %.1f%%\n", pct_rhat_ok))

  # Warnings
  if (result$n_divergent > 0) {
    warning(sprintf("[dlvm] %d divergent transitions. Consider increasing adapt_delta.", result$n_divergent))
  }
  if (max_rhat > 1.01) {
    warning("[dlvm] Some R-hat > 1.01. Chains may not have converged. Consider more iterations.")
  }
  if (min_ess_bulk < 400) {
    warning("[dlvm] Some ESS < 400. Consider more sampling iterations.")
  }

  invisible(result)
}

# ============================================================================
# LOO-CV
# ============================================================================

#' Compute LOO-CV using the loo package
#'
#' Uses the pointwise log-likelihood from the generated quantities block.
#'
#' @param fit A dlvm_fit object.
#' @return A loo object (from the loo package).
#' @export
dlvm_loo <- function(fit) {
  if (!inherits(fit, "dlvm_fit")) {
    stop("[dlvm] 'fit' must be a dlvm_fit object.")
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("[dlvm] Package 'loo' required. Install with: install.packages('loo')")
  }

  log_lik <- fit$fit$draws("log_lik", format = "matrix")
  loo_result <- loo::loo(log_lik)
  return(loo_result)
}

# ============================================================================
# Plotting
# ============================================================================

#' Plot latent state trajectories
#'
#' @param fit A dlvm_fit object.
#' @param units Character vector of unit labels to plot (default: all).
#' @param show_obs Logical; overlay observed data points (default: TRUE).
#' @param prob Numeric; credible interval width (default: 0.95).
#' @param ncol Integer; number of facet columns (default: auto).
#' @return A ggplot2 object.
#' @export
dlvm_plot <- function(fit, units = NULL, show_obs = TRUE, prob = 0.95, ncol = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("[dlvm] Package 'ggplot2' required.")
  }

  summ <- dlvm_summary(fit, prob = prob)
  unit_col <- fit$metadata$unit_col
  time_col <- fit$metadata$time_col

  if (!is.null(units)) {
    summ <- summ[summ[[unit_col]] %in% units, ]
  }

  # Cap facets for readability
  all_units <- unique(summ[[unit_col]])
  if (length(all_units) > 30 && is.null(units)) {
    message("[dlvm] Showing first 30 units. Pass 'units' to select specific ones.")
    summ <- summ[summ[[unit_col]] %in% all_units[1:30], ]
  }

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data[[time_col]])) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#3b82f6", alpha = 0.2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = estimate),
      color = "#1d4ed8", linewidth = 0.7
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data[[unit_col]]),
      scales = "free_y",
      ncol = ncol
    ) +
    ggplot2::labs(
      x = time_col,
      y = "Latent state (\u03b8)",
      title = "DLVM: Latent State Trajectories",
      subtitle = sprintf("%d%% credible intervals", round(prob * 100))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Overlay raw individual observations if available
  if (show_obs) {
    obs_data <- as.data.frame(fit$metadata$obs_data)
    value_col <- fit$metadata$value_col
    if (!is.null(obs_data) && value_col %in% names(obs_data)) {
      if (!is.null(units)) {
        obs_data <- obs_data[obs_data[[unit_col]] %in% units, ]
      }
      if (length(all_units) > 30 && is.null(units)) {
        obs_data <- obs_data[obs_data[[unit_col]] %in% all_units[1:30], ]
      }
      p <- p + ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(y = .data[[value_col]]),
        color = "#ef4444", size = 0.8, alpha = 0.4
      )
    } else if ("obs_mean" %in% names(summ)) {
      # Fallback to aggregated means
      obs_subset <- summ[!is.na(summ$obs_mean), ]
      p <- p + ggplot2::geom_point(
        data = obs_subset,
        ggplot2::aes(y = obs_mean),
        color = "#ef4444", size = 1, alpha = 0.6
      )
    }
  }

  return(p)
}

#' Plot posterior predictive check
#'
#' Overlays replicated data (y_rep) from the posterior predictive distribution
#' on the observed data to assess model fit.
#'
#' @param fit A dlvm_fit object.
#' @param n_rep Integer; number of posterior predictive draws to overlay (default: 50).
#' @return A ggplot2 object.
#' @export
dlvm_plot_ppc <- function(fit, n_rep = 50) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("[dlvm] Package 'ggplot2' required.")
  }

  # Extract y_rep draws
  y_rep <- fit$fit$draws("y_rep", format = "matrix")
  y_obs <- fit$metadata$obs_data[[fit$metadata$value_col]]
  n_obs <- length(y_obs)

  # Sample a subset of draws
  n_draws <- nrow(y_rep)
  draw_idx <- sample(n_draws, min(n_rep, n_draws))

  # Build long-format data for y_rep
  rep_list <- lapply(draw_idx, function(d) {
    data.frame(value = y_rep[d, ], draw = d, type = "y_rep")
  })
  rep_df <- do.call(rbind, rep_list)
  obs_df <- data.frame(value = y_obs, draw = 0, type = "y_obs")

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = rep_df,
      ggplot2::aes(x = value, group = factor(draw)),
      color = "#93c5fd", alpha = 0.3, linewidth = 0.3
    ) +
    ggplot2::geom_density(
      data = obs_df,
      ggplot2::aes(x = value),
      color = "#1d4ed8", linewidth = 1.2
    ) +
    ggplot2::labs(
      x = "Value",
      y = "Density",
      title = "Posterior Predictive Check",
      subtitle = sprintf("Dark blue = observed, light blue = %d posterior predictive draws", length(draw_idx))
    ) +
    ggplot2::theme_minimal(base_size = 12)

  return(p)
}
