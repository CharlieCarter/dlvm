# dlvm — Dynamic Latent Variable Model

An R package for estimating dynamic latent variables from noisy, irregularly-spaced panel data. Uses a Bayesian random-walk state-space model with Student-t observation and innovation distributions, compiled and sampled via [CmdStan](https://mc-stan.org/cmdstanr/).

## Features

- **Irregular time spacing** — handles variable gaps between observations with automatic δt scaling
- **Variable observation counts** — multiple noisy measurements per time period per unit
- **Measurement uncertainty** — optionally incorporates known standard errors
- **Multi-threaded sampling** — within-chain parallelisation via `reduce_sum` and TBB
- **Robust distributions** — Student-t observation and innovation models for outlier resistance

## Installation

```r
# From GitHub:
devtools::install_github("charliecarter/dlvm")

# From a local clone:
devtools::install("/path/to/dlvm")
```

### Prerequisites

The package requires CmdStan (the command-line interface to Stan). If you don't have it:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan(cores = parallel::detectCores())
```

Or let the package check for you:

```r
dlvm::dlvm_check_deps(install = TRUE)
```

### First-time setup

After installing, run the one-time CmdStan configuration. This auto-detects your
platform and compiler, and writes the correct CmdStan `make/local` for your machine:

```r
library(dlvm)
dlvm_setup_cmdstan()   # auto-detects platform, configures CmdStan
dlvm_compile()         # first compile takes ~1 min, cached afterwards
```

On Linux and macOS with default Xcode tools, `dlvm_setup_cmdstan()` just enables
threading. On macOS with Homebrew Clang, it automatically detects the conflict and
configures CmdStan to use the system Apple Clang instead (see Platform Notes below).

## Quick Start

```r
library(dlvm)

# Prepare data (data.frame with unit, time, and value columns)
prep <- dlvm_prepare(my_data, unit_col = "country", time_col = "year",
                     value_col = "score", se_col = "se")

# Compile Stan model (cached after first run)
mod <- dlvm_compile()

# Fit the model
fit <- dlvm_fit(prep, chains = 4, threads_per_chain = 2)

# Extract results
summ <- dlvm_summary(fit)           # Posterior median + 95% CI
params <- dlvm_parameters(fit)      # sigma, innov estimates
diag <- dlvm_diagnostics(fit)       # Rhat, ESS, divergences

# Visualise
dlvm_plot(fit)                      # Latent trajectories
dlvm_plot_ppc(fit)                  # Posterior predictive check

# Model comparison
loo_result <- dlvm_loo(fit)
```

### Observation Families

DLVM supports four observation distribution families:

```r
# Continuous data — Student-t observation model (default)
prep <- dlvm_prepare(my_data, "country", "year", "score",
                     family = "student_t")

# Ordinal ratings (1..K, K auto-detected from data)
prep <- dlvm_prepare(my_data, "country", "year", "rating",
                     family = "cumulative")

# Individual binary outcomes (0/1)
prep <- dlvm_prepare(my_data, "country", "year", "sympathetic",
                     family = "bernoulli")

# Counts: k of n successes
prep <- dlvm_prepare(my_data, "country", "year", "n_sympathetic",
                     trials_col = "n_articles", family = "binomial")
```

Data-type aliases are also accepted: `"continuous"` → `"student_t"`, `"ordinal"` → `"cumulative"`, `"binary"` → `"bernoulli"`.

### Tuning Parameters: `nu_obs` vs `nu_state`

The model uses two separate Student-t distributions with independently controllable degrees of freedom:

| Parameter | Controls | Default | Purpose |
|---|---|---|---|
| `nu_obs` | Observation distribution tails | 4 | How much a single extreme observation can pull θ. Lower = more robust to outlier observations. |
| `nu_state` | Innovation distribution tails | 4 | How large structural breaks in the latent trajectory can be. Lower = more tolerant of abrupt regime shifts. |

Both are passed as data (not estimated). Setting either to ≥ 30 is effectively Gaussian.

> **TODO**: Systematic simulation study exploring the impact of `nu_obs` on latent state recovery
> under different outlier regimes (Student-t observation vs. Gaussian observation).

## Demos

Run the showcase from the package source directory:

```bash
cd dlvm/
Rscript demo/demo_showcase.R              # All 4 scenarios
Rscript demo/demo_showcase.R polling      # Single scenario
```

The demo generates data for four use cases:
1. **Election polling** — multi-pollster data with house effects
2. **Conflict events** — irregularly-timed events with variable density
3. **Financial sentiment** — multi-source indicator fusion
4. **Social media sentiment** — high-volume multi-country tracking

---

## Platform Notes

### Linux / HPC

Works out of the box:

```r
dlvm_setup_cmdstan()  # enables threading, no other config needed
dlvm_compile()
```

### macOS with Xcode Command Line Tools (default)

Works without issues. Apple's Clang and system `libc++` are consistent:

```r
dlvm_setup_cmdstan()  # enables threading
```

### macOS with Homebrew Clang (⚠️ common issue)

If you've installed LLVM/Clang via Homebrew (`brew install llvm`), you may
encounter **TBB linking failures** or **header search path errors** when
compiling Stan models. This happens because Homebrew's Clang uses its own
`libc++` headers (newer ABI) but links against the system's `libc++.dylib`
(older ABI).

**Automatic fix:**

```r
dlvm_setup_cmdstan(force = TRUE)
# → Detects Homebrew Clang, configures CmdStan to use /usr/bin/clang++ instead
```

**What it does:** `dlvm_setup_cmdstan()` detects Homebrew Clang in your PATH
by checking `clang++ --version` for "Homebrew". When found, it writes a
`make/local` that sets `CXX = /usr/bin/clang++` (Apple's system Clang), which
avoids the ABI mismatch entirely. It also cleans any stale pre-compiled headers.

**Manual alternative:** Edit `make/local` directly:

```bash
# Find your CmdStan path:
Rscript -e 'cat(cmdstanr::cmdstan_path())'

# Edit make/local:
nano ~/.cmdstan/cmdstan-2.38.0/make/local
```

```makefile
# Use system Apple Clang (avoids Homebrew ABI mismatch)
CXX = /usr/bin/clang++

# Enable threading
STAN_THREADS=true
```

**How to tell if you have Homebrew Clang:**

```bash
clang++ --version
# Homebrew: "Homebrew clang version 21.x.x"
# System:   "Apple clang version 16.x.x"
```

**Runtime TBB path:** If the model compiles but crashes at runtime,
`dlvm_fit()` handles this automatically by setting `DYLD_LIBRARY_PATH`
to CmdStan's bundled TBB.

### Windows

Not tested, but should work via CmdStan's MSYS2/mingw64 toolchain. TBB
issues are unlikely on Windows.

---

## Package Structure

```
dlvm/
├── DESCRIPTION          # Package metadata
├── NAMESPACE            # Exported functions
├── R/
│   ├── dlvm_data.R      # dlvm_prepare() — data validation + Stan data list
│   ├── dlvm_fit.R       # dlvm_fit() — CmdStan sampling wrapper
│   ├── dlvm_results.R   # dlvm_summary/plot/diagnostics/loo
│   ├── dlvm_utils.R     # Family resolution, compilation, path helpers
│   └── zzz.R            # Package load hooks
├── inst/stan/
│   ├── dlvm.stan            # Legacy Student-t model (backward compat)
│   ├── dlvm_student_t.stan  # Student-t observation family
│   ├── dlvm_cumulative.stan # Ordinal (ordered logistic) family
│   ├── dlvm_bernoulli.stan  # Binary (Bernoulli logit) family
│   └── dlvm_binomial.stan   # Binomial (count/proportion) family
├── demo/
│   └── demo_showcase.R  # 4 demo scenarios with ggplot visualisations
└── tests/
    ├── run_tests.R      # Test runner
    └── test_*.R         # Unit tests
```

## License

MIT
