# dlvm/R/dlvm_utils.R
# Internal utility functions for the DLVM package
#
# This file provides dependency management, Stan model compilation,
# and path resolution helpers.

# ============================================================================
# Family resolution
# ============================================================================

# Valid primary family names (distribution-accurate, following brms conventions)
.dlvm_families <- c("student_t", "cumulative", "bernoulli", "binomial")

# Data-type aliases accepted silently
.dlvm_family_aliases <- c(
  continuous = "student_t",
  ordinal    = "cumulative",
  binary     = "bernoulli"
)

#' Resolve and validate a family argument
#'
#' Accepts primary distribution names or data-type aliases.
#' Rejects 'gaussian' with an informative error.
#'
#' @param family Character string.
#' @return Validated canonical family name.
#' @keywords internal
.resolve_family <- function(family) {
  family <- tolower(trimws(family))

  # Reject gaussian with helpful message
  if (family == "gaussian") {
    stop("[dlvm] DLVM uses Student-t, not Gaussian, for continuous data. ",
         "Use family = 'student_t' (or the alias 'continuous').")
  }

  # Resolve aliases
  if (family %in% names(.dlvm_family_aliases)) {
    family <- .dlvm_family_aliases[[family]]
  }

  # Validate
  if (!family %in% .dlvm_families) {
    stop("[dlvm] Unknown family '", family, "'. ",
         "Valid families: ", paste(.dlvm_families, collapse = ", "),
         ". Aliases: ", paste(names(.dlvm_family_aliases), collapse = ", "), ".")
  }

  family
}

# ============================================================================
# Path resolution
# ============================================================================

#' Get path to a DLVM Stan model file
#'
#' Resolves the path to the Stan model for the given family, checking (in order):
#' 1. Installed package location (via system.file)
#' 2. Development layout: inst/stan/ relative to package source root
#' 3. Legacy layout: stan/ relative to working directory
#'
#' @param family Character; observation family name (default: "student_t").
#' @return Absolute path to the Stan model file.
#' @export
dlvm_stan_path <- function(family = "student_t") {
  family <- .resolve_family(family)
  stan_file <- paste0("dlvm_", family, ".stan")

  # 1. Installed package
  installed <- system.file("stan", stan_file, package = "dlvm")
  if (nzchar(installed)) return(installed)

  # 2. Development: inst/stan/ relative to package source root
  dev_candidates <- c(
    file.path("inst", "stan", stan_file),
    file.path(".", "inst", "stan", stan_file)
  )
  for (cand in dev_candidates) {
    if (file.exists(cand)) return(normalizePath(cand))
  }

  # 3. Legacy layout (stan/ in working directory) — student_t falls back to dlvm.stan
  if (family == "student_t") {
    legacy_candidates <- c(
      file.path("stan", "dlvm.stan"),
      file.path(".", "stan", "dlvm.stan"),
      file.path("inst", "stan", "dlvm.stan"),
      file.path(".", "inst", "stan", "dlvm.stan")
    )
    for (cand in legacy_candidates) {
      if (file.exists(cand)) return(normalizePath(cand))
    }
  }

  stop("[dlvm] Cannot locate ", stan_file, ". ",
       "If using as installed package, reinstall with devtools::install(). ",
       "If developing, run from the dlvm/ source directory.")
}

# ============================================================================
# Dependency checking
# ============================================================================

#' Check and optionally install required R packages
#'
#' @param install Logical; if TRUE, attempt to install missing packages
#' @export
#' @return Invisible TRUE if all dependencies are satisfied
dlvm_check_deps <- function(install = FALSE) {
  required <- c("cmdstanr", "data.table", "loo", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) == 0) {
    message("[dlvm] All dependencies satisfied.")
    return(invisible(TRUE))
  }

  if (!install) {
    stop(
      "[dlvm] Missing required packages: ", paste(missing, collapse = ", "), "\n",
      "  Install cmdstanr:    install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))\n",
      "  Install others:      install.packages(c('", paste(setdiff(missing, "cmdstanr"), collapse = "','"), "'))\n",
      "  Then run:            cmdstanr::install_cmdstan()\n",
      "  Or call:             dlvm_check_deps(install = TRUE)"
    )
  }

  # Install missing packages
  if ("cmdstanr" %in% missing) {
    message("[dlvm] Installing cmdstanr from stan-dev r-universe...")
    install.packages("cmdstanr",
                     repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
    missing <- setdiff(missing, "cmdstanr")
  }
  if (length(missing) > 0) {
    message("[dlvm] Installing: ", paste(missing, collapse = ", "))
    install.packages(missing)
  }

  # Check CmdStan backend
  if (!is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))) {
    message("[dlvm] CmdStan backend found at: ", cmdstanr::cmdstan_path())
  } else {
    message("[dlvm] CmdStan not found. Installing (this takes a few minutes)...")
    cmdstanr::install_cmdstan(cores = parallel::detectCores())
  }

  message("[dlvm] All dependencies installed successfully.")
  invisible(TRUE)
}

# ============================================================================
# Model compilation
# ============================================================================

# .dlvm_env is defined in zzz.R (package-level object)

#' Compile a DLVM Stan model
#'
#' Compiles the model for the specified observation family using CmdStan's
#' make system directly (which correctly handles TBB linking for reduce_sum
#' parallelisation on all platforms), then wraps the pre-compiled binary
#' in a cmdstanr CmdStanModel object.
#'
#' The compiled binary is stored alongside the .stan file and reused across
#' R sessions unless force=TRUE or the .stan file has been modified.
#' Models are cached per-family, so compiling a different family does not
#' invalidate the cache for other families.
#'
#' @param family Character; observation family (default: "student_t").
#'   Valid values: "student_t", "cumulative", "bernoulli", "binomial".
#'   Aliases: "continuous", "ordinal", "binary".
#' @param force Logical; if TRUE, recompile even if cached
#' @param quiet Logical; if TRUE, suppress compilation messages
#' @param threads Logical; if TRUE (default), compile with threading support
#' @return A CmdStanModel object
#' @export
dlvm_compile <- function(family = "student_t", force = FALSE, quiet = FALSE, threads = TRUE) {
  family <- .resolve_family(family)
  cache_key <- paste0("model_", family)

  if (!force && exists(cache_key, envir = .dlvm_env)) {
    if (!quiet) message(sprintf("[dlvm] Using cached compiled model (family: %s).", family))
    return(.dlvm_env[[cache_key]])
  }

  stan_file <- dlvm_stan_path(family)
  if (!file.exists(stan_file)) {
    stop("[dlvm] Stan model not found at: ", stan_file)
  }

  # Determine binary output path (next to the .stan file)
  exe_file <- sub("\\.stan$", "", stan_file)

  # Check if we can skip compilation (binary exists and is newer than .stan)
  if (!force && file.exists(exe_file) &&
      file.mtime(exe_file) > file.mtime(stan_file)) {
    if (!quiet) message(sprintf("[dlvm] Using existing compiled binary (family: %s).", family))
  } else {
    # Compile via CmdStan's make system, which handles TBB linking correctly
    cmdstan_dir <- cmdstanr::cmdstan_path()
    make_args <- sprintf("STANCFLAGS='--O1'")
    if (threads) {
      make_args <- paste(make_args, "STAN_THREADS=TRUE")
    }

    make_cmd <- sprintf(
      "cd '%s' && make '%s' %s",
      cmdstan_dir, exe_file, make_args
    )

    if (!quiet) {
      threading_str <- if (threads) "with threading" else "without threading"
      message(sprintf("[dlvm] Compiling Stan model (family: %s) %s via CmdStan make...", family, threading_str))
    }

    result <- system(make_cmd, intern = FALSE, ignore.stdout = quiet, ignore.stderr = quiet)

    if (result != 0) {
      stop("[dlvm] Stan model compilation failed. Check CmdStan configuration.\n",
           "  Hint: run  cmdstanr::rebuild_cmdstan()  if CmdStan was recently updated.")
    }

    if (!file.exists(exe_file)) {
      stop("[dlvm] Compilation appeared to succeed but binary not found at: ", exe_file)
    }
  }

  # Wrap the pre-compiled binary in a CmdStanModel object
  # Pass cpp_options so cmdstanr knows about threading support
  mod <- cmdstanr::cmdstan_model(
    stan_file,
    exe_file = exe_file,
    cpp_options = if (threads) list(stan_threads = TRUE) else list()
  )

  .dlvm_env[[cache_key]] <- mod
  .dlvm_env$has_threads <- threads
  threading_str <- if (threads) "with threading" else "without threading"
  if (!quiet) message(sprintf("[dlvm] Compilation successful (family: %s, %s). Binary cached for this session.", family, threading_str))
  return(mod)
}

#' Check if the compiled model has threading support
#' @return Logical
#' @export
dlvm_has_threads <- function() {

  isTRUE(.dlvm_env$has_threads)
}
# ============================================================================
# TBB / runtime helpers
# ============================================================================

#' Get path to CmdStan's TBB library directory
#' @return Absolute path to TBB lib dir, or NULL if not found
#' @export
dlvm_tbb_lib_path <- function() {
  tbb_path <- file.path(cmdstanr::cmdstan_path(), "stan", "lib", "stan_math", "lib", "tbb")
  if (dir.exists(tbb_path)) tbb_path else NULL
}

# ============================================================================
# CmdStan setup helper
# ============================================================================

#' Configure CmdStan for dlvm (auto-detects platform issues)
#'
#' Detects your platform and compiler, then writes an appropriate
#' CmdStan `make/local` configuration file. Handles known issues:
#'
#' - **Homebrew Clang on macOS**: switches to system Apple Clang
#' - **Stale C headers in /usr/local/include**: adds -isysroot to scope headers
#' - **Homebrew TBB shadowing CmdStan's bundled TBB**: pins TBB paths
#' - **Missing TBB dylibs / stale build markers**: cleans and optionally rebuilds
#'
#' On Linux and Windows, this just enables threading — no special config needed.
#'
#' @param force Logical; if TRUE, overwrite existing make/local (default: FALSE)
#' @param rebuild Logical; if TRUE, trigger a full clean rebuild of CmdStan
#'   after writing make/local. This rebuilds TBB dylibs, pre-compiled headers,
#'   and everything else from scratch. Recommended on first setup or after
#'   encountering linker errors. Takes ~2-5 minutes. (Default: FALSE)
#' @param verbose Logical; if TRUE, print diagnostic information (default: TRUE)
#' @return Invisible path to the make/local file
#' @export
dlvm_setup_cmdstan <- function(force = FALSE, rebuild = FALSE, verbose = TRUE) {
  cmdstan_dir <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) {
      stop("[dlvm] CmdStan not found. Install it first:\n",
           "  cmdstanr::install_cmdstan()")
    }
  )

  make_local <- file.path(cmdstan_dir, "make", "local")

  # If make/local exists and not forcing, show contents and return
  if (file.exists(make_local) && !force && !rebuild) {
    existing <- readLines(make_local, warn = FALSE)
    if (verbose) {
      message("[dlvm] Existing make/local found at: ", make_local)
      message("[dlvm] Contents:\n", paste("  ", existing, collapse = "\n"))
      message("[dlvm] Use dlvm_setup_cmdstan(force = TRUE) to overwrite.")
    }
    return(invisible(make_local))
  }

  # ---- Detect platform and build config ----
  os <- Sys.info()[["sysname"]]
  config_lines <- character()
  notes <- character()

  if (os == "Darwin") {
    config_lines <- c("# Auto-configured by dlvm::dlvm_setup_cmdstan()")

    # 1. Detect Homebrew Clang
    clang_version <- tryCatch(
      system2("clang++", "--version", stdout = TRUE, stderr = TRUE)[1],
      error = function(e) ""
    )
    is_homebrew_clang <- grepl("Homebrew|homebrew", clang_version)

    if (is_homebrew_clang) {
      system_clang <- "/usr/bin/clang++"
      if (!file.exists(system_clang)) {
        stop("[dlvm] Homebrew Clang detected but system Clang not found.\n",
             "  Install Xcode Command Line Tools:  xcode-select --install")
      }
      config_lines <- c(config_lines,
        paste0("# Detected Homebrew Clang: ", Sys.which("clang++")),
        "# Using system Apple Clang to avoid libc++ ABI mismatch",
        paste0("CXX = ", system_clang), ""
      )
      notes <- c(notes, paste0("Compiler: Switching to system Apple Clang ",
                                "(Homebrew detected at ", Sys.which("clang++"), ")"))
    } else {
      config_lines <- c(config_lines, "# System Apple Clang detected", "")
      notes <- c(notes, "Compiler: System Apple Clang (no issues)")
    }

    # 2. Detect stale /usr/local/include (common on machines with Homebrew history)
    #    These old headers contain C macros (isinf/isnan) that conflict with C++.
    #    We can't use -isysroot to work around this because it also breaks
    #    CmdStan's internal TBB build. The clean fix is to remove them.
    stale_math_h <- "/usr/local/include/math.h"
    if (file.exists(stale_math_h)) {
      math_content <- tryCatch(readLines(stale_math_h, n = 200, warn = FALSE),
                               error = function(e) "")
      if (any(grepl("__inline_isinff|__inline_isnanf", math_content))) {
        notes <- c(notes,
          "WARNING: Stale /usr/local/include/ detected (will cause C++ macro conflicts).",
          "  Fix:  sudo mv /usr/local/include /usr/local/include.bak")
        # Set a flag so we can stop early if this is a problem
        .dlvm_env$stale_usr_local_include <- TRUE
      }
    }

    # 3. Detect system TBB libs that shadow CmdStan's bundled TBB via -ltbb
    #    Both /usr/local/lib (old x86 leftovers) and /opt/homebrew/lib (active
    #    Homebrew TBB v12) cause linker to find the wrong TBB.
    stale_tbb_local <- list.files("/usr/local/lib", pattern = "^libtbb", full.names = TRUE)
    stale_tbb_brew  <- list.files("/opt/homebrew/lib", pattern = "^libtbb", full.names = TRUE)
    stale_tbb_libs  <- c(stale_tbb_local, stale_tbb_brew)
    if (length(stale_tbb_libs) > 0) {
      fix_cmds <- character()
      if (length(stale_tbb_local) > 0) fix_cmds <- c(fix_cmds, "sudo rm -f /usr/local/lib/libtbb*")
      if (length(stale_tbb_brew) > 0)  fix_cmds <- c(fix_cmds, "sudo rm -f /opt/homebrew/lib/libtbb*")
      notes <- c(notes,
        "WARNING: System TBB libraries will shadow CmdStan's bundled TBB and cause linker errors.",
        paste0("  Fix:  ", paste(fix_cmds, collapse = " && ")))
    }

    # 4. Detect Homebrew TBB that could shadow CmdStan's bundled TBB
    #    Only override TBB_INC (header path), NOT TBB_LIB.
    #    Setting TBB_LIB tells CmdStan's build system to treat TBB as external
    #    and skip building TBB dylibs entirely — which causes linker errors.
    has_homebrew_tbb <- dir.exists("/opt/homebrew/Cellar/tbb") ||
                        dir.exists("/opt/homebrew/include/tbb") ||
                        dir.exists("/usr/local/Cellar/tbb")

    tbb_inc_path <- file.path(cmdstan_dir, "stan", "lib", "stan_math", "lib",
                              "tbb_2020.3", "include")

    if (has_homebrew_tbb && dir.exists(tbb_inc_path)) {
      config_lines <- c(config_lines,
        "# Pin CmdStan's bundled TBB headers (Homebrew TBB detected, incompatible API)",
        "# NOTE: Only TBB_INC is set — TBB_LIB is intentionally omitted so that",
        "# CmdStan builds its own TBB dylibs from source (setting TBB_LIB would",
        "# tell the build system TBB is external and skip dylib compilation).",
        paste0("TBB_INC = ", tbb_inc_path), ""
      )
      notes <- c(notes, "TBB: Headers pinned to CmdStan bundled (Homebrew TBB bypassed)")
    }

    # 4. Threading
    config_lines <- c(config_lines,
      "# Enable threading for reduce_sum parallelisation",
      "STAN_THREADS=true"
    )
    notes <- c(notes, "Threading: enabled")

  } else if (os == "Linux") {
    config_lines <- c(
      "# Auto-configured by dlvm::dlvm_setup_cmdstan()", "",
      "# Enable threading for reduce_sum parallelisation",
      "STAN_THREADS=true"
    )
    notes <- c("Platform: Linux (no special config needed)", "Threading: enabled")
  } else {
    config_lines <- c(
      "# Auto-configured by dlvm::dlvm_setup_cmdstan()", "",
      "STAN_THREADS=true"
    )
    notes <- c(paste0("Platform: ", os), "Threading: enabled")
  }

  # ---- Write make/local ----
  writeLines(config_lines, make_local)

  # ---- Clean stale build artifacts ----
  pch_dir <- file.path(cmdstan_dir, "stan", "src", "stan", "model",
                       "model_header.hpp.gch")
  if (dir.exists(pch_dir)) {
    unlink(pch_dir, recursive = TRUE)
    notes <- c(notes, "Cleaned: stale pre-compiled headers")
  }

  tbb_dir <- file.path(cmdstan_dir, "stan", "lib", "stan_math", "lib", "tbb")
  tbb_marker <- file.path(tbb_dir, "tbb-make-check")
  tbb_dylibs <- list.files(tbb_dir, pattern = "\\.(dylib|so)$")

  if (file.exists(tbb_marker) && length(tbb_dylibs) == 0) {
    # Marker says built but no dylibs — remove marker + ancillary files
    unlink(tbb_marker)
    for (f in list.files(tbb_dir, pattern = "^(version_|tbb\\.def|tbbmalloc\\.def|tbbvars)",
                         full.names = TRUE)) {
      unlink(f)
    }
    notes <- c(notes, "Cleaned: TBB build markers (dylibs were missing)")
    if (!rebuild) {
      rebuild <- TRUE
      notes <- c(notes, "Auto-triggering rebuild (TBB dylibs need to be built)")
    }
  }

  # ---- Print config summary ----
  if (verbose) {
    message("[dlvm] CmdStan configured:")
    for (n in notes) message("  ", n)
    message("  Config: ", make_local)
  }

  # ---- Optional full rebuild ----
  if (rebuild) {
    if (verbose) {
      message("")
      message("[dlvm] Rebuilding CmdStan from scratch (2-5 minutes)...")
      message("[dlvm] This rebuilds TBB dylibs, pre-compiled headers, and all tooling.")
    }

    # Clean everything
    system2("make", c("-C", cmdstan_dir, "clean-all"),
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE)

    # Rebuild
    cores <- max(1L, parallel::detectCores())
    build_result <- system2(
      "make", c("-C", cmdstan_dir, paste0("-j", cores), "build"),
      stdout = if (verbose) "" else FALSE,
      stderr = if (verbose) "" else FALSE
    )

    if (build_result != 0) {
      stop("[dlvm] CmdStan rebuild failed (exit code ", build_result, ").\n",
           "  Check the error output above.\n",
           "  You may also try: cmdstanr::rebuild_cmdstan()")
    }

    # Verify TBB dylibs
    rebuilt_dylibs <- list.files(tbb_dir, pattern = "\\.(dylib|so)$")
    if (length(rebuilt_dylibs) > 0) {
      if (verbose) message("[dlvm] TBB dylibs built: ",
                           paste(rebuilt_dylibs, collapse = ", "))
    } else {
      warning("[dlvm] Build succeeded but TBB dylibs not found at: ", tbb_dir)
    }

    if (verbose) {
      message("[dlvm] CmdStan rebuild complete.")
      message("[dlvm] Next: run dlvm_compile() to build the Stan model.")
    }
  } else if (verbose) {
    message("")
    message("[dlvm] Next: run dlvm_compile() to build the Stan model.")
    message("[dlvm] If you get linker errors, run: dlvm_setup_cmdstan(rebuild = TRUE)")
  }

  invisible(make_local)
}

# Formatting helpers
# ============================================================================

#' Format seconds as human-readable duration
#' @param secs Numeric seconds
#' @return Character string like "2m 34s"
fmt_duration <- function(secs) {
  if (secs < 60) return(sprintf("%.1fs", secs))
  mins <- floor(secs / 60)
  remaining <- secs - mins * 60
  sprintf("%dm %.0fs", mins, remaining)
}

# ============================================================================
# Statistical helpers
# ============================================================================

#' Calculate small-sample bias corrected standard deviation (Cureton/Holtzman)
#'
#' Computes the standard deviation using the Holtzman correction factor to account
#' for bias in small samples.
#'
#' @param x Numeric vector of observations
#' @param na.rm Logical; if TRUE, remove NA values before computation
#' @return Numeric standard deviation
#' @export
dlvm_sd_cureton <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)

  if (n < 2) return(0) # Undefined or zero variance for n=0,1

  # Standard SD (with Bessel's correction n-1)
  s_bessel <- sd(x)

  # Holtzman correction factor:
  # K = (n - 1) / [ Gamma((n-1)/2) * (sqrt((n-1)/2) / Gamma(n/2)) ]^2
  # This factor K is applied to the variance sum of squares, but defined
  # relative to the uncorrected estimator?
  #
  # Let's use the user's provided implementation exactly:
  # C_N = Gamma((N-1)/2) * (sqrt((N-1)/2) / Gamma(N/2))
  # k = (N-1) / C_N^2
  # s = sqrt( sum((x-mean)^2) / k )

  c_n <- gamma((n - 1) / 2) * (sqrt((n - 1) / 2) / gamma(n / 2))
  k <- (n - 1) / c_n^2

  # Sum of squared deviations
  ssd <- sum((x - mean(x))^2)

  # Corrected SD
  sqrt(ssd / k)
}

#' Calculate small-sample bias corrected standard error of the mean
#'
#' Computes the standard error of the mean using the Cureton/Holtzman
#' corrected standard deviation.
#'
#' @param x Numeric vector of observations
#' @param na.rm Logical; if TRUE, remove NA values before computation
#' @return Numeric standard error
#' @export
dlvm_se_cureton <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return(0)

  dlvm_sd_cureton(x) / sqrt(n)
}
