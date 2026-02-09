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

Typically works out of the box. CmdStan's bundled TBB and the system g++ are compatible:

```r
cmdstanr::install_cmdstan()
library(dlvm)
dlvm_compile()  # Just works
```

### macOS with Xcode Command Line Tools (default)

Usually works without issues. Apple's Clang and system `libc++` are consistent:

```r
cmdstanr::install_cmdstan()
```

### macOS with Homebrew Clang (⚠️ common issue)

If you've installed LLVM/Clang via Homebrew (e.g., `brew install llvm`), you may encounter **TBB linking failures** when compiling Stan models with `reduce_sum`. This happens because Homebrew's Clang uses its own `libc++` headers (newer C++17 ABI) but links against the system's `libc++.dylib` (older ABI), causing undefined symbol errors like:

```
Undefined symbols:
  "std::__1::__hash_memory::__hash_memory()"
  "typeinfo for tbb::task"
```

#### The Fix

Edit CmdStan's `make/local` file to force Homebrew's Clang to use the **system** C++ standard library headers instead of its own:

```bash
# Find your CmdStan path:
Rscript -e 'cat(cmdstanr::cmdstan_path())'

# Edit make/local in that directory:
nano ~/.cmdstan/cmdstan-2.38.0/make/local
```

Add these lines (adjust paths for your system):

```makefile
# Force Homebrew Clang to use system libc++ headers (fixes TBB ABI mismatch)
CXXFLAGS += -nostdinc++ -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1

# Enable threading for reduce_sum
STAN_THREADS=true
```

> **Why this works:** Homebrew's Clang ships its own `libc++` headers at
> `/opt/homebrew/opt/llvm/include/c++/v1/` which define a newer ABI
> (with symbols like `__hash_memory`). But the actual `libc++.dylib`
> linked at runtime is Apple's system version, which doesn't have those
> symbols. The `-nostdinc++` flag tells Clang to skip its own headers,
> and `-isystem` points it to the system SDK headers that match the
> system library.

#### Runtime TBB Path

If the model compiles but crashes at runtime with TBB errors, you may need to ensure CmdStan's bundled TBB is loaded instead of Homebrew's. The `dlvm_fit()` function handles this automatically by setting `DYLD_LIBRARY_PATH`, but if you're calling CmdStan directly, set:

```bash
export DYLD_LIBRARY_PATH="$(Rscript -e 'cat(cmdstanr::cmdstan_path())')/stan/lib/stan_math/lib/tbb:$DYLD_LIBRARY_PATH"
```

#### How to tell if you have Homebrew Clang

```bash
which clang++
# If /opt/homebrew/opt/llvm/bin/clang++ → Homebrew
# If /usr/bin/clang++ → System (Xcode)

clang++ --version
# Homebrew: "Homebrew clang version 21.x.x"
# System:   "Apple clang version 16.x.x"
```

### Windows

Not tested, but should work via CmdStan's MSYS2/mingw64 toolchain. TBB issues are unlikely on Windows.

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
│   ├── dlvm_utils.R     # Compilation, path resolution, TBB helpers
│   └── zzz.R            # Package load hooks
├── inst/stan/
│   └── dlvm.stan        # The Stan model
├── demo/
│   └── demo_showcase.R  # 4 demo scenarios with ggplot visualisations
└── tests/
    ├── run_tests.R      # Test runner
    └── test_*.R         # Unit tests
```

## License

MIT
