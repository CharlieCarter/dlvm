# demo_showcase.R — Four DLVM Demo Scenarios
#
# Each scenario simulates realistic data for a different use case,
# fits the DLVM, and produces a ggplot visualization showing
# raw observations + latent posterior trajectory with credible intervals.
#
# Usage:
#   cd dlvm/
#   Rscript demo/demo_showcase.R           # Run all 4 scenarios
#   Rscript demo/demo_showcase.R polling   # Run only one scenario

# ============================================================================
# Setup
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
})

# Source DLVM package — use library() if installed, source() for development
if (requireNamespace("dlvm", quietly = TRUE)) {
  library(dlvm)
} else {
  for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
}

# Output directory
dir.create("demo/output", showWarnings = FALSE, recursive = TRUE)

# Parse which scenarios to run
args <- commandArgs(trailingOnly = TRUE)
run_all <- length(args) == 0
scenarios_to_run <- if (run_all) c("polling", "conflict", "financial", "social") else args

cat("=== DLVM Demo Showcase ===\n")
cat("Running scenarios:", paste(scenarios_to_run, collapse = ", "), "\n\n")

# Compile model once for all scenarios
mod <- dlvm_compile(threads = TRUE)

# ============================================================================
# Helper: enhanced plot with raw data + posterior
# ============================================================================
dlvm_demo_plot <- function(fit, raw_data, unit_col, time_col, value_col,
                            title, subtitle, source_col = NULL,
                            point_alpha = 0.35, point_size = 1.2,
                            prob = 0.95, facet_ncol = NULL) {

  summ <- dlvm_summary(fit, prob = prob)

  # Base plot: ribbon + line
  p <- ggplot(summ, aes(x = .data[[time_col]])) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "#6366f1", alpha = 0.18) +
    geom_line(aes(y = estimate),
              color = "#4338ca", linewidth = 0.9)

  # Raw data points
  if (!is.null(source_col) && source_col %in% names(raw_data)) {
    p <- p + geom_point(
      data = raw_data,
      aes(x = .data[[time_col]], y = .data[[value_col]],
          color = .data[[source_col]]),
      size = point_size, alpha = point_alpha
    ) +
    scale_color_brewer(palette = "Set2", name = source_col)
  } else {
    p <- p + geom_point(
      data = raw_data,
      aes(x = .data[[time_col]], y = .data[[value_col]]),
      color = "#ef4444", size = point_size, alpha = point_alpha
    )
  }

  # Faceting
  n_units <- length(unique(summ[[unit_col]]))
  if (n_units > 1) {
    p <- p + facet_wrap(vars(.data[[unit_col]]),
                        scales = "free_y", ncol = facet_ncol)
  }

  # Theme
  p <- p +
    labs(
      x = time_col, y = "Latent state",
      title = title, subtitle = subtitle
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 10),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      legend.position = if (!is.null(source_col)) "bottom" else "none"
    )

  return(p)
}

# ============================================================================
# Scenario 1: Election Polling
# ============================================================================
if ("polling" %in% scenarios_to_run) {
  cat("--- Scenario 1: Election Polling ---\n")

  set.seed(2024)

  # Simulate latent voter intention for 2 candidates over 120 days
  n_days <- 120
  n_pollsters <- 6
  days <- seq.Date(as.Date("2024-06-01"), by = "day", length.out = n_days)

  # Latent truth: slow drift with one "convention bounce"
  theta_true <- 48 + cumsum(rnorm(n_days, 0, 0.15))
  theta_true[40:55] <- theta_true[40:55] + seq(0, 3, length.out = 16)  # convention bounce
  theta_true[56:70] <- theta_true[56:70] + seq(3, 0.5, length.out = 15) # decay

  pollster_names <- paste0("Pollster_", LETTERS[1:n_pollsters])
  pollster_bias <- c(-1.5, 0.8, -0.3, 1.2, -0.8, 0.5)  # house effects
  pollster_noise <- c(2.0, 3.5, 1.5, 2.5, 3.0, 1.8)     # measurement noise

  # Each pollster publishes polls on random days (not every day)
  polls <- do.call(rbind, lapply(seq_len(n_pollsters), function(p) {
    # Random subset of days (more frequent closer to election)
    poll_prob <- seq(0.08, 0.25, length.out = n_days)
    poll_days <- which(runif(n_days) < poll_prob)
    if (length(poll_days) < 3) poll_days <- sample(n_days, 5)

    data.frame(
      candidate = "Harris",
      day = as.numeric(days[poll_days] - days[1]),
      date = days[poll_days],
      pct = theta_true[poll_days] + pollster_bias[p] + rnorm(length(poll_days), 0, pollster_noise[p]),
      se = pollster_noise[p] / sqrt(sample(400:1500, length(poll_days), replace = TRUE) / 100),
      pollster = pollster_names[p],
      stringsAsFactors = FALSE
    )
  }))

  cat(sprintf("  %d polls from %d pollsters over %d days\n", nrow(polls), n_pollsters, n_days))

  # Fit
  prep <- dlvm_prepare(polls, "candidate", "day", "pct", se_col = "se", mode = "individual")
  fit <- dlvm_fit(prep, chains = 2, parallel_chains = 2, threads_per_chain = 2,
                  iter_warmup = 500, iter_sampling = 500, show_messages = FALSE)

  # ggplot
  summ <- dlvm_summary(fit)
  # Map summary day back to date
  summ$date <- days[1] + summ$day

  polls_plot <- polls
  p1 <- ggplot() +
    geom_ribbon(data = summ, aes(x = date, ymin = lower, ymax = upper),
                fill = "#6366f1", alpha = 0.18) +
    geom_line(data = summ, aes(x = date, y = estimate),
              color = "#4338ca", linewidth = 1) +
    geom_point(data = polls_plot, aes(x = date, y = pct, color = pollster),
               size = 1.5, alpha = 0.5) +
    geom_line(data = data.frame(date = days, truth = theta_true),
              aes(x = date, y = truth), linetype = "dashed", color = "grey30", linewidth = 0.5) +
    scale_color_brewer(palette = "Set2", name = "Pollster") +
    labs(
      x = "Date", y = "Vote share (%)",
      title = "Election Polling: Latent Voter Intention",
      subtitle = "Points = individual polls by 6 pollsters | Line = posterior median | Ribbon = 95% CI | Dashed = true latent"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave("demo/output/01_polling.pdf", p1, width = 10, height = 6)
  ggsave("demo/output/01_polling.png", p1, width = 10, height = 6, dpi = 150)
  cat("  Saved: demo/output/01_polling.pdf\n")

  diag1 <- fit$fit$diagnostic_summary(quiet = TRUE)
  cat(sprintf("  Divergences: %d | sigma: %.2f | innov: %.3f\n",
              sum(diag1$num_divergent),
              mean(fit$fit$draws("sigma", format = "matrix")),
              mean(fit$fit$draws("innov", format = "matrix"))))
}

# ============================================================================
# Scenario 2: Conflict Event Intensity (ACLED-style)
# ============================================================================
if ("conflict" %in% scenarios_to_run) {
  cat("\n--- Scenario 2: Conflict Event Intensity ---\n")

  set.seed(2023)

  # Simulate irregularly-timed conflict events across 3 regions over 2 years
  regions <- c("Donetsk", "Kherson", "Odesa")
  n_events_per_region <- c(600, 300, 150)

  # Simulate over 1 year (365 days) — rounded to integer days
  # so events cluster into day bins with variable counts
  days_range <- 365

  events_list <- lapply(seq_along(regions), function(r) {
    n <- n_events_per_region[r]
    # Non-uniform event times (cluster around escalation periods)
    if (r == 1) {
      # Donetsk: surge early, sustained
      t_raw <- c(rbeta(n * 0.3, 2, 5), rbeta(n * 0.7, 1.5, 2)) * days_range
    } else if (r == 2) {
      # Kherson: mid-period surge (counteroffensive)
      t_raw <- c(rbeta(n * 0.4, 3, 3), rbeta(n * 0.6, 5, 2)) * days_range
    } else {
      # Odesa: sporadic throughout
      t_raw <- runif(n) * days_range
    }
    # Round to integer days — multiple events per day
    t_day <- round(t_raw)
    t_day <- pmax(1, pmin(t_day, days_range))

    data.frame(
      region = regions[r],
      day = t_day,
      # "Severity score" as noisy observation + region baseline
      severity = rnorm(n, mean = c(3.5, 2.5, 1.8)[r], sd = c(1.2, 1.0, 1.5)[r]),
      stringsAsFactors = FALSE
    )
  })

  events <- do.call(rbind, events_list)
  cat(sprintf("  %d events across %d regions over %d days\n",
              nrow(events), length(regions), days_range))

  # Fit — individual mode (many events per time-region)
  prep2 <- dlvm_prepare(events, "region", "day", "severity", mode = "individual")
  fit2 <- dlvm_fit(prep2, chains = 2, parallel_chains = 2, threads_per_chain = 2,
                   iter_warmup = 500, iter_sampling = 500, show_messages = FALSE)

  # ggplot
  summ2 <- dlvm_summary(fit2)
  summ2$date <- as.Date("2022-01-01") + summ2$day - 1
  events$date <- as.Date("2022-01-01") + events$day - 1

  p2 <- ggplot() +
    geom_point(data = events, aes(x = date, y = severity),
               color = "#ef4444", size = 0.4, alpha = 0.15) +
    geom_ribbon(data = summ2, aes(x = date, ymin = lower, ymax = upper),
                fill = "#6366f1", alpha = 0.25) +
    geom_line(data = summ2, aes(x = date, y = estimate),
              color = "#4338ca", linewidth = 0.9) +
    facet_wrap(~region, ncol = 1, scales = "free_y") +
    labs(
      x = "Date", y = "Event severity score",
      title = "Conflict Event Intensity: Latent Severity by Region",
      subtitle = "Points = individual events (irregular timestamps) | Line = latent trajectory | Ribbon = 95% CI"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 10),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )

  ggsave("demo/output/02_conflict.pdf", p2, width = 10, height = 8)
  ggsave("demo/output/02_conflict.png", p2, width = 10, height = 8, dpi = 150)
  cat("  Saved: demo/output/02_conflict.pdf\n")

  diag2 <- fit2$fit$diagnostic_summary(quiet = TRUE)
  cat(sprintf("  Divergences: %d | sigma: %.2f | innov: %.3f\n",
              sum(diag2$num_divergent),
              mean(fit2$fit$draws("sigma", format = "matrix")),
              mean(fit2$fit$draws("innov", format = "matrix"))))
}

# ============================================================================
# Scenario 3: Financial Market Sentiment
# ============================================================================
if ("financial" %in% scenarios_to_run) {
  cat("\n--- Scenario 3: Financial Market Sentiment ---\n")

  set.seed(2008)

  # Simulate multiple market indicators tracking latent "market confidence"
  # over 2007-2009 (financial crisis period)
  n_months <- 36  # 3 years
  months <- seq.Date(as.Date("2007-01-01"), by = "month", length.out = n_months)

  # Latent confidence: gradual decline → crash → partial recovery
  theta_fin <- numeric(n_months)
  theta_fin[1] <- 60  # starting confidence index
  for (t in 2:n_months) {
    # Crisis dynamics
    if (t <= 12) shock <- -0.3          # 2007: gradual worry
    else if (t <= 15) shock <- -1.5     # early 2008: Bear Stearns
    else if (t <= 21) shock <- -3.0     # mid-2008: Lehman
    else if (t <= 24) shock <- -1.0     # late 2008: TARP
    else shock <- 0.8                   # 2009: recovery
    theta_fin[t] <- theta_fin[t - 1] + shock + rnorm(1, 0, 0.5)
  }

  # 5 different market indicators, each with different precision and bias
  indicators <- data.frame(
    name = c("Consumer Confidence", "VIX (inverted)", "Credit Spreads (inv)",
             "CEO Survey", "Retail Flows"),
    bias = c(0, -5, 3, 2, -3),
    noise_sd = c(3, 8, 5, 4, 6),
    freq = c(1, 1, 1, 0.33, 0.5),         # fraction of months observed
    stringsAsFactors = FALSE
  )

  mkt <- do.call(rbind, lapply(1:nrow(indicators), function(i) {
    obs_months <- which(runif(n_months) < indicators$freq[i])
    if (length(obs_months) < 3) obs_months <- sort(sample(n_months, 5))

    data.frame(
      sector = "Market",
      month = obs_months,
      date = months[obs_months],
      value = theta_fin[obs_months] + indicators$bias[i] +
              rnorm(length(obs_months), 0, indicators$noise_sd[i]),
      se = indicators$noise_sd[i],
      indicator = indicators$name[i],
      stringsAsFactors = FALSE
    )
  }))

  cat(sprintf("  %d observations from %d indicators over %d months\n",
              nrow(mkt), nrow(indicators), n_months))

  prep3 <- dlvm_prepare(mkt, "sector", "month", "value", se_col = "se", mode = "individual")
  fit3 <- dlvm_fit(prep3, chains = 2, parallel_chains = 2, threads_per_chain = 2,
                   iter_warmup = 500, iter_sampling = 500, show_messages = FALSE)

  summ3 <- dlvm_summary(fit3)
  summ3$date <- months[summ3$month]
  mkt_plot <- mkt

  p3 <- ggplot() +
    geom_ribbon(data = summ3, aes(x = date, ymin = lower, ymax = upper),
                fill = "#6366f1", alpha = 0.18) +
    geom_line(data = summ3, aes(x = date, y = estimate),
              color = "#4338ca", linewidth = 1) +
    geom_point(data = mkt_plot, aes(x = date, y = value, color = indicator),
               size = 2, alpha = 0.6) +
    geom_line(data = data.frame(date = months, truth = theta_fin),
              aes(x = date, y = truth), linetype = "dashed", color = "grey30", linewidth = 0.5) +
    scale_color_brewer(palette = "Dark2", name = "Indicator") +
    labs(
      x = "Date", y = "Confidence index",
      title = "Financial Market Sentiment: Latent Confidence (2007–2009)",
      subtitle = "Points = 5 indicators with varying precision + bias | Line = latent | Dashed = true"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave("demo/output/03_financial.pdf", p3, width = 10, height = 6)
  ggsave("demo/output/03_financial.png", p3, width = 10, height = 6, dpi = 150)
  cat("  Saved: demo/output/03_financial.pdf\n")

  diag3 <- fit3$fit$diagnostic_summary(quiet = TRUE)
  cat(sprintf("  Divergences: %d | sigma: %.2f | innov: %.3f\n",
              sum(diag3$num_divergent),
              mean(fit3$fit$draws("sigma", format = "matrix")),
              mean(fit3$fit$draws("innov", format = "matrix"))))
}

# ============================================================================
# Scenario 4: Social Media Sentiment (multi-country, larger scale)
# ============================================================================
if ("social" %in% scenarios_to_run) {
  cat("\n--- Scenario 4: Social Media Sentiment (Multi-Country) ---\n")

  set.seed(2025)

  # 8 countries, daily posts over 6 months, variable volume
  countries <- c("USA", "UK", "Germany", "Brazil")
  n_countries <- length(countries)
  n_days <- 90

  # Generate latent sentiment per country (random walks with country-specific dynamics)
  theta_social <- matrix(0, nrow = n_days, ncol = n_countries)
  for (c in seq_len(n_countries)) {
    theta_social[1, c] <- rnorm(1, 0, 0.3)
    innov_scale <- runif(1, 0.02, 0.08)
    for (t in 2:n_days) {
      theta_social[t, c] <- theta_social[t - 1, c] + rnorm(1, 0, innov_scale)
    }
    # Add a sentiment shock for some countries
    if (c %in% c(1, 3)) {
      shock_day <- sample(30:(n_days - 20), 1)
      shock_end <- min(shock_day + 14, n_days)
      theta_social[shock_day:shock_end, c] <- theta_social[shock_day:shock_end, c] +
        seq(0, -0.8, length.out = shock_end - shock_day + 1)
    }
  }

  # Generate noisy posts — variable number per country per day
  posts <- do.call(rbind, lapply(seq_len(n_countries), function(c) {
    do.call(rbind, lapply(1:n_days, function(d) {
      # Number of posts varies (more for bigger countries)
      base_rate <- c(30, 20, 15, 20)[c]
      n_posts <- rpois(1, base_rate * runif(1, 0.3, 1.2))
      if (n_posts == 0) return(NULL)

      data.frame(
        country = countries[c],
        day = d,
        sentiment = theta_social[d, c] + rnorm(n_posts, 0, 0.5),
        stringsAsFactors = FALSE
      )
    }))
  }))

  cat(sprintf("  %d posts across %d countries over %d days\n",
              nrow(posts), n_countries, n_days))

  # Fit — individual mode
  prep4 <- dlvm_prepare(posts, "country", "day", "sentiment", mode = "individual")
  fit4 <- dlvm_fit(prep4, chains = 2, parallel_chains = 2, threads_per_chain = 2,
                   iter_warmup = 500, iter_sampling = 500, show_messages = FALSE)

  summ4 <- dlvm_summary(fit4)

  # Build truth data for overlay
  truth_df <- do.call(rbind, lapply(seq_len(n_countries), function(c) {
    data.frame(country = countries[c], day = 1:n_days,
               truth = theta_social[, c], stringsAsFactors = FALSE)
  }))

  # Subsample posts for plotting (too many points otherwise)
  posts_sub <- posts[sample(nrow(posts), min(nrow(posts), 15000)), ]

  p4 <- ggplot() +
    geom_point(data = posts_sub, aes(x = day, y = sentiment),
               color = "#94a3b8", size = 0.3, alpha = 0.08) +
    geom_ribbon(data = summ4, aes(x = day, ymin = lower, ymax = upper),
                fill = "#6366f1", alpha = 0.3) +
    geom_line(data = summ4, aes(x = day, y = estimate),
              color = "#4338ca", linewidth = 0.7) +
    geom_line(data = truth_df, aes(x = day, y = truth),
              linetype = "dashed", color = "#dc2626", linewidth = 0.4, alpha = 0.7) +
    facet_wrap(~country, ncol = 2, scales = "free_y") +
    labs(
      x = "Day", y = "Sentiment score",
      title = "Social Media Sentiment: Latent Opinion Across 8 Countries",
      subtitle = sprintf("Grey = %dk individual posts | Blue = posterior | Red dashed = true latent",
                         round(nrow(posts) / 1000))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 10),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )

  ggsave("demo/output/04_social_media.pdf", p4, width = 12, height = 7)
  ggsave("demo/output/04_social_media.png", p4, width = 12, height = 7, dpi = 150)
  cat("  Saved: demo/output/04_social_media.pdf\n")

  diag4 <- fit4$fit$diagnostic_summary(quiet = TRUE)
  cat(sprintf("  Divergences: %d | sigma: %.2f | innov: %.3f\n",
              sum(diag4$num_divergent),
              mean(fit4$fit$draws("sigma", format = "matrix")),
              mean(fit4$fit$draws("innov", format = "matrix"))))
}

cat("\n=== Demo Showcase Complete ===\n")
cat("Output saved to: demo/output/\n")
