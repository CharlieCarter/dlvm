#!/usr/bin/env Rscript
# dlvm/demo/demo_vdem.R
# Full worked example: Dynamic Robust Latent Variable Model
#
# Demonstrates the complete DLVM workflow using gapminder data
# (life expectancy by country, observed every 5 years — perfect for
# testing continuous-time/irregular spacing).
#
# Usage:
#   cd dlvm/
#   Rscript demo/demo_vdem.R

cat("============================================================\n")
cat("  DLVM Demo: Latent Trajectories from Panel Data\n")
cat("============================================================\n\n")

# --- Setup ---
dlvm_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."), mustWork = FALSE)
if (!file.exists(file.path(dlvm_root, "stan", "dlvm.stan"))) {
  dlvm_root <- normalizePath(".", mustWork = FALSE)
}
for (f in list.files(file.path(dlvm_root, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# Ensure output directory exists
output_dir <- file.path(dlvm_root, "demo", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Check dependencies ---
dlvm_check_deps(install = FALSE)

# --- Load demo data ---
if (!requireNamespace("gapminder", quietly = TRUE)) {
  stop("Install gapminder: install.packages('gapminder')")
}

data("gapminder", package = "gapminder")
gm <- as.data.frame(gapminder::gapminder)

# Select a diverse set of countries
selected <- c(
  "United States", "United Kingdom", "France", "Germany", "Japan",
  "China", "India", "Brazil", "South Africa", "Nigeria",
  "Egypt", "Australia", "Mexico", "Indonesia", "Russia",
  "Argentina", "Turkey", "Thailand", "Kenya", "Colombia",
  "Sweden", "Norway", "Chile", "Iran", "Pakistan",
  "Vietnam", "Ethiopia", "Philippines", "Poland", "Bangladesh"
)
gm_sub <- gm[gm$country %in% selected, ]

cat(sprintf("Data: %d observations, %d countries, years %d-%d\n",
            nrow(gm_sub), length(unique(gm_sub$country)),
            min(gm_sub$year), max(gm_sub$year)))
cat(sprintf("Time spacing: every 5 years (irregular — exercises delta_t)\n\n"))

# Standardise for numerical stability
gm_sub$lifeExp_z <- (gm_sub$lifeExp - mean(gm_sub$lifeExp)) / sd(gm_sub$lifeExp)

# ============================================================
# Example 1: Basic fit (means mode, no SE)
# ============================================================
cat("--- Example 1: Basic fit (means, no SE) ---\n")

prep1 <- dlvm_prepare(gm_sub, "country", "year", "lifeExp_z")
print(prep1)

t1 <- system.time({
  fit1 <- dlvm_fit(prep1,
                   chains = 4,
                   parallel_chains = 4,
                   threads_per_chain = 2,
                   iter_warmup = 1000,
                   iter_sampling = 2000,
                   show_messages = FALSE)
})

cat("\n--- Diagnostics ---\n")
diag1 <- dlvm_diagnostics(fit1)

cat("\n--- Scalar Parameters ---\n")
print(dlvm_parameters(fit1))

cat("\n--- LOO-CV ---\n")
loo1 <- dlvm_loo(fit1)
print(loo1)

# Save trajectory plot
p1 <- dlvm_plot(fit1, ncol = 5)
ggplot2::ggsave(file.path(output_dir, "trajectories_all.pdf"),
                p1, width = 16, height = 14)
cat(sprintf("\nSaved: %s\n", file.path(output_dir, "trajectories_all.pdf")))

# Select countries plot
p1_select <- dlvm_plot(fit1, units = c("United States", "China", "Nigeria", "Sweden", "Brazil"))
ggplot2::ggsave(file.path(output_dir, "trajectories_select.pdf"),
                p1_select, width = 12, height = 4)
cat(sprintf("Saved: %s\n", file.path(output_dir, "trajectories_select.pdf")))

# PPC plot
p1_ppc <- dlvm_plot_ppc(fit1)
ggplot2::ggsave(file.path(output_dir, "ppc.pdf"), p1_ppc, width = 8, height = 5)
cat(sprintf("Saved: %s\n", file.path(output_dir, "ppc.pdf")))

# ============================================================
# Example 2: With artificial measurement error (SE mode)
# ============================================================
cat("\n--- Example 2: With known standard errors ---\n")

# Add artificial SEs (larger for earlier years = less reliable data)
set.seed(42)
gm_sub$se <- pmax(0.05, 0.5 - 0.005 * (gm_sub$year - 1950) + rnorm(nrow(gm_sub), 0, 0.05))

prep2 <- dlvm_prepare(gm_sub, "country", "year", "lifeExp_z", se_col = "se")
print(prep2)

fit2 <- dlvm_fit(prep2,
                 chains = 4,
                 parallel_chains = 4,
                 threads_per_chain = 2,
                 iter_warmup = 1000,
                 iter_sampling = 2000,
                 show_messages = FALSE)

cat("\n--- Diagnostics (with SE) ---\n")
diag2 <- dlvm_diagnostics(fit2)

cat("\n--- LOO comparison: no-SE vs with-SE ---\n")
loo2 <- dlvm_loo(fit2)
loo_comp <- loo::loo_compare(list(no_se = loo1, with_se = loo2))
print(loo_comp)

# Save comparison plot
p2 <- dlvm_plot(fit2, units = c("United States", "China", "Nigeria", "Sweden", "Brazil"))
ggplot2::ggsave(file.path(output_dir, "trajectories_with_se.pdf"),
                p2, width = 12, height = 4)
cat(sprintf("Saved: %s\n", file.path(output_dir, "trajectories_with_se.pdf")))

# ============================================================
# Benchmark summary
# ============================================================
cat("\n============================================================\n")
cat("  BENCHMARK SUMMARY\n")
cat("============================================================\n")
cat(sprintf("Example 1 (no SE):   %s\n", fmt_duration(fit1$timing$total_secs)))
cat(sprintf("Example 2 (with SE): %s\n", fmt_duration(fit2$timing$total_secs)))
cat(sprintf("Data size:           %d obs, %d latent states\n",
            prep1$stan_data$n_obs, prep1$stan_data$n_states))
cat(sprintf("Chains:              4 × (1000 warmup + 2000 sampling)\n"))
cat(sprintf("Threading:           4 chains × 2 threads = 8 total\n"))

cat("\n--- Diagnostics Comparison ---\n")
cat(sprintf("%-20s  %-10s  %-10s\n", "", "No SE", "With SE"))
cat(sprintf("%-20s  %-10d  %-10d\n", "Divergences",
            diag1$n_divergent, diag2$n_divergent))
cat(sprintf("%-20s  %-10.4f  %-10.4f\n", "Max R-hat (theta)",
            diag1$theta_max_rhat, diag2$theta_max_rhat))
cat(sprintf("%-20s  %-10.0f  %-10.0f\n", "Min ESS (theta)",
            diag1$theta_min_ess_bulk, diag2$theta_min_ess_bulk))

cat("\nDone! Check demo/output/ for plots.\n")
