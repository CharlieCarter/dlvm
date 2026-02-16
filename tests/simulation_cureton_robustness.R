#!/usr/bin/env Rscript
# tests/simulation_cureton_robustness.R
# Comprehensive Robustness Check & Visualization for Cureton Adjustment

suppressPackageStartupMessages({
  library(dlvm)
  library(cmdstanr)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork) # For combining plots
})

# Parameters
N_ITERS <- 5
N_UNITS <- 5
N_TIMES <- 30
INNOV_TRUE <- 0.5
SIGMA_TRUE <- 1.0

cat("=== Cureton Robustness Simulation ===\n")
cat(sprintf("Iterations: %d\nUnits: %d\nTime: %d\n", N_ITERS, N_UNITS, N_TIMES))

# Storage for results
results_df <- data.frame()
plot_data_list <- list()

# --- Simulation Loop ---
for (i in 1:N_ITERS) {
  set.seed(1000 + i)

  # 1. Generate Data with Anomalies
  thetas <- matrix(NA_real_, nrow = N_UNITS, ncol = N_TIMES)
  for (u in 1:N_UNITS) {
    thetas[u, 1] <- rnorm(1, 0, 2)
    for (t in 2:N_TIMES) {
      thetas[u, t] <- thetas[u, t - 1] + INNOV_TRUE * rnorm(1)
    }
  }

  obs_list <- list()
  anomaly_info <- list()

  for (u in 1:N_UNITS) {
    for (t in 1:N_TIMES) {
      # Inject anomaly at Unit 1, Time 15 in EVERY iteration for consistency
      is_anomaly <- (u == 1 && t == 15)

      if (is_anomaly) {
        # Mean shift = +5, N=2, tight cluster
        target <- thetas[u, t] + 5.0
        vals <- c(target, target + 0.1)
        anomaly_info <- list(u=u, t=t, val=mean(vals), true=thetas[u, t])
      } else {
        # Normal data N=3
        vals <- rnorm(3, mean = thetas[u, t], sd = SIGMA_TRUE)
      }

      obs_list[[length(obs_list) + 1]] <- data.frame(
        unit = paste0("unit_", u), time = t, value = vals, stringsAsFactors = FALSE
      )
    }
  }
  sim_data <- do.call(rbind, obs_list)

  # 2. Fit Models
  # Cureton
  prep_c <- dlvm_prepare(sim_data, "unit", "time", "value", mode = "means", se_correction = TRUE)
  fit_c <- dlvm_fit(prep_c, chains=1, threads_per_chain=1, iter_warmup=200, iter_sampling=200,
                    refresh=0, show_messages=FALSE) # Fast fit for loop

  # Standard
  prep_s <- dlvm_prepare(sim_data, "unit", "time", "value", mode = "means", se_correction = FALSE)
  fit_s <- dlvm_fit(prep_s, chains=1, threads_per_chain=1, iter_warmup=200, iter_sampling=200,
                    refresh=0, show_messages=FALSE)

  # 3. Extract Metrics
  # function to get posterior summary
  get_post <- function(fit) {
    d <- fit$fit$draws("theta", format="matrix")
    list(mean=colMeans(d), sd=apply(d, 2, sd), q05=apply(d, 2, quantile, 0.05), q95=apply(d, 2, quantile, 0.95))
  }

  post_c <- tryCatch(get_post(fit_c), error=function(e) NULL)
  post_s <- tryCatch(get_post(fit_s), error=function(e) NULL)

  if (!is.null(post_c) && !is.null(post_s)) {
    # Anomaly Index (Unit 1 is first N_TIMES states)
    idx_anom <- (anomaly_info$u - 1) * N_TIMES + anomaly_info$t

    # Global RMSE
    true_vec <- as.vector(t(thetas))
    rmse_c <- sqrt(mean((post_c$mean - true_vec)^2))
    rmse_s <- sqrt(mean((post_s$mean - true_vec)^2))

    # Store
    results_df <- rbind(results_df, data.frame(
      Iter = i,
      RMSE_Cureton = rmse_c,
      RMSE_Standard = rmse_s,
      Anom_True = anomaly_info$true,
      Anom_Obs = anomaly_info$val,
      Anom_Est_C = post_c$mean[idx_anom],
      Anom_Est_S = post_s$mean[idx_anom],
      Anom_SD_C = post_c$sd[idx_anom],
      Anom_SD_S = post_s$sd[idx_anom],
      Input_SE_C = prep_c$metadata$obs_data$se[idx_anom],
      Input_SE_S = prep_s$metadata$obs_data$se[idx_anom]
    ))

    # Save LAST iteration data for plotting
    if (i == N_ITERS) {
      plot_data <- data.frame(
        Time = 1:N_TIMES,
        True = thetas[1, ],
        Obs_Val = NA,
        Obs_SE_C = NA,
        Est_C = post_c$mean[1:N_TIMES],
        Low_C = post_c$q05[1:N_TIMES],
        High_C = post_c$q95[1:N_TIMES],
        Est_S = post_s$mean[1:N_TIMES],
        Low_S = post_s$q05[1:N_TIMES],
        High_S = post_s$q95[1:N_TIMES]
      )

      # Add obs data for unit 1
      u1_data <- prep_c$metadata$obs_data[prep_c$metadata$obs_data$unit == "unit_1", ]
      # Match times
      plot_data$Obs_Val[u1_data$time] <- u1_data$value
      plot_data$Obs_SE_C[u1_data$time] <- u1_data$se # Cureton SEs
      plot_data$Obs_SE_S_Input <- prep_s$metadata$obs_data$se[u1_data$time] # Standard SEs

      plot_data_list[[1]] <- plot_data
    }
  }

  if (i %% 5 == 0) cat(sprintf("Processed %d/%d...\n", i, N_ITERS))
}

# --- 4. summary ---
cat("\n=== Robustness Summary ===\n")
cat(sprintf("Average RMSE: Cureton %.4f vs Standard %.4f\n",
            mean(results_df$RMSE_Cureton), mean(results_df$RMSE_Standard)))

cat("Anomaly Uncertainty (Average Posterior SD at Outlier):\n")
cat(sprintf("  Cureton:  %.4f\n", mean(results_df$Anom_SD_C)))
cat(sprintf("  Standard: %.4f\n", mean(results_df$Anom_SD_S)))
cat(sprintf("  Increase: +%.1f%%\n", 100 * (mean(results_df$Anom_SD_C) - mean(results_df$Anom_SD_S))/mean(results_df$Anom_SD_S)))

# --- 5. Visualization ---
pd <- plot_data_list[[1]]

# A. Trajectory Plot
p1 <- ggplot(pd, aes(x=Time)) +
  # Truth
  geom_line(aes(y=True), linetype="dashed", color="grey50", alpha=0.8) +
  # CIs (Ribbons) - Plotted first to be behind
  geom_ribbon(aes(ymin=Low_C, ymax=High_C, fill="Cureton"), alpha=0.2) +
  geom_ribbon(aes(ymin=Low_S, ymax=High_S, fill="Standard"), alpha=0.2) +
  # Estimates
  geom_line(aes(y=Est_C, color="Cureton"), linewidth=1) +
  geom_line(aes(y=Est_S, color="Standard"), linewidth=1) +
  # Anomaly Point
  geom_point(aes(y=Obs_Val), shape=4, size=3, data=pd %>% filter(Time == 15)) +
  annotate("text", x=15, y=pd$Obs_Val[15]+1, label="Anomaly (N=2)", color="black", fontface="bold") +
  # Scales
  scale_color_manual(values=c("Cureton"="#1f77b4", "Standard"="#d62728")) +
  scale_fill_manual(values=c("Cureton"="#1f77b4", "Standard"="#d62728")) +
  labs(title = "Latent Variable Recovery & Uncertainty",
       subtitle = "Comparing Cureton correction (Blue) vs Standard SEM (Red) at an anomaly",
       y = "Latent Value", color="Method", fill="Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

# B. Uncertainty Zoom at Anomaly (Bar Plot)
anom_data <- results_df %>%
  summarise(
    Input_SE_C = mean(Input_SE_C),
    Input_SE_S = mean(Input_SE_S),
    Post_SD_C = mean(Anom_SD_C),
    Post_SD_S = mean(Anom_SD_S)
  )

zoom_df <- data.frame(
  Method = c("Standard", "Standard", "Cureton", "Cureton"),
  Type = c("Input SE", "Posterior SD", "Input SE", "Posterior SD"),
  Value = c(anom_data$Input_SE_S, anom_data$Post_SD_S, anom_data$Input_SE_C, anom_data$Post_SD_C)
)
zoom_df$Method <- factor(zoom_df$Method, levels=c("Standard", "Cureton"))

p2 <- ggplot(zoom_df, aes(x=Method, y=Value, fill=Method, alpha=Type)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  scale_fill_manual(values=c("Standard"="#d62728", "Cureton"="#1f77b4")) +
  scale_alpha_manual(values=c("Input SE"=0.4, "Posterior SD"=1.0)) +
  labs(title = "Uncertainty Propagation",
       subtitle = "Input SE vs. Posterior SD at Anomaly (Avg over 20 runs)",
       y = "Standard Deviation / Error") +
  theme_minimal() +
  theme(legend.position = "none")

# Combine
final_plot <- p1 + inset_element(p2, left=0.6, bottom=0.6, right=0.99, top=0.99)

ggsave("tests/cureton_robustness_viz.png", final_plot, width=10, height=6)
cat("Plot saved to tests/cureton_robustness_viz.png\n")
