# dlvm/R/dlvm_utils.R
# Internal utility functions for the DLVM package
#
# This file provides dependency management, Stan model compilation,
# and path resolution helpers.

# ============================================================================
# Path resolution
# ============================================================================

#' Get path to the Stan model file
#'
#' Resolves the path to dlvm.stan, checking (in order):
#' 1. Installed package location (via system.file)
#' 2. Development layout: inst/stan/ relative to package source root
#' 3. Legacy layout: stan/ relative to working directory
#'
#' @return Absolute path to dlvm.stan
#' @export
dlvm_stan_path <- function() {
  # 1. Installed package
  installed <- system.file("stan", "dlvm.stan", package = "dlvm")
  if (nzchar(installed)) return(installed)

  # 2. Development: inst/stan/ relative to package source root
  dev_candidates <- c(
    file.path("inst", "stan", "dlvm.stan"),
    file.path(".", "inst", "stan", "dlvm.stan")
  )
  for (cand in dev_candidates) {
    if (file.exists(cand)) return(normalizePath(cand))
  }

  # 3. Legacy layout (stan/ in working directory)
  legacy_candidates <- c(
    file.path("stan", "dlvm.stan"),
    file.path(".", "stan", "dlvm.stan")
  )
  for (cand in legacy_candidates) {
    if (file.exists(cand)) return(normalizePath(cand))
  }

  stop("[dlvm] Cannot locate dlvm.stan. ",
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

#' Compile the DLVM Stan model
#'
#' Compiles the model using CmdStan's make system directly (which correctly
#' handles TBB linking for reduce_sum parallelisation on all platforms),
#' then wraps the pre-compiled binary in a cmdstanr CmdStanModel object.
#'
#' The compiled binary is stored alongside the .stan file and reused across
#' R sessions unless force=TRUE or the .stan file has been modified.
#'
#' @param force Logical; if TRUE, recompile even if cached
#' @param quiet Logical; if TRUE, suppress compilation messages
#' @param threads Logical; if TRUE (default), compile with threading support
#' @return A CmdStanModel object
#' @export
dlvm_compile <- function(force = FALSE, quiet = FALSE, threads = TRUE) {
  if (!force && exists("model", envir = .dlvm_env)) {
    if (!quiet) message("[dlvm] Using cached compiled model.")
    return(.dlvm_env$model)
  }

  stan_file <- dlvm_stan_path()
  if (!file.exists(stan_file)) {
    stop("[dlvm] Stan model not found at: ", stan_file)
  }

  # Determine binary output path (next to the .stan file)
  exe_file <- sub("\\.stan$", "", stan_file)

  # Check if we can skip compilation (binary exists and is newer than .stan)
  if (!force && file.exists(exe_file) &&
      file.mtime(exe_file) > file.mtime(stan_file)) {
    if (!quiet) message("[dlvm] Using existing compiled binary.")
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
      message(sprintf("[dlvm] Compiling Stan model %s via CmdStan make...", threading_str))
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

  .dlvm_env$model <- mod
  .dlvm_env$has_threads <- threads
  threading_str <- if (threads) "with threading" else "without threading"
  if (!quiet) message(sprintf("[dlvm] Compilation successful (%s). Binary cached for this session.", threading_str))
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
#' CmdStan `make/local` configuration file. This is a one-time setup
#' step that ensures Stan model compilation works correctly.
#'
#' On most systems (Linux, macOS with default Xcode tools, Windows),
#' this amounts to enabling threading only. On macOS with Homebrew LLVM/Clang
#' in the PATH, it configures CmdStan to use the system Apple Clang instead,
#' avoiding a known ABI mismatch between Homebrew's libc++ headers and the
#' system's libc++ runtime library.
#'
#' @param force Logical; if TRUE, overwrite existing make/local (default: FALSE)
#' @param verbose Logical; if TRUE, print diagnostic information (default: TRUE)
#' @return Invisible path to the make/local file
#' @export
dlvm_setup_cmdstan <- function(force = FALSE, verbose = TRUE) {
  cmdstan_dir <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) {
      stop("[dlvm] CmdStan not found. Install it first:\n",
           "  cmdstanr::install_cmdstan()")
    }
  )

  make_local <- file.path(cmdstan_dir, "make", "local")

  # If make/local exists and not forcing, warn

  if (file.exists(make_local) && !force) {
    existing <- readLines(make_local, warn = FALSE)
    if (verbose) {
      message("[dlvm] Existing make/local found at: ", make_local)
      message("[dlvm] Contents:\n", paste("  ", existing, collapse = "\n"))
      message("[dlvm] Use dlvm_setup_cmdstan(force = TRUE) to overwrite.")
    }
    return(invisible(make_local))
  }

  # Detect platform
  os <- Sys.info()[["sysname"]]
  arch <- Sys.info()[["machine"]]
  config_lines <- character()
  notes <- character()

  if (os == "Darwin") {
    # Check if Homebrew Clang is the default clang++
    clang_path <- Sys.which("clang++")
    clang_version <- tryCatch(
      system2("clang++", "--version", stdout = TRUE, stderr = TRUE)[1],
      error = function(e) ""
    )
    is_homebrew <- grepl("Homebrew|homebrew", clang_version, ignore.case = TRUE)

    if (is_homebrew) {
      # Homebrew Clang detected — use system Apple Clang instead
      system_clang <- "/usr/bin/clang++"
      if (!file.exists(system_clang)) {
        stop("[dlvm] Homebrew Clang detected but system Clang not found at ",
             system_clang, ".\n",
             "  Install Xcode Command Line Tools: xcode-select --install")
      }

      config_lines <- c(
        "# Auto-configured by dlvm::dlvm_setup_cmdstan()",
        paste0("# Detected Homebrew Clang at: ", clang_path),
        "# Using system Apple Clang to avoid libc++ ABI mismatch",
        paste0("CXX = ", system_clang),
        "",
        "# Enable threading for reduce_sum parallelisation",
        "STAN_THREADS=true"
      )
      notes <- c(
        "Detected: Homebrew Clang (may cause TBB/libc++ ABI conflicts)",
        paste0("Fix: Using system Apple Clang at ", system_clang),
        "Threading: enabled"
      )
    } else {
      # System Apple Clang — just enable threading
      config_lines <- c(
        "# Auto-configured by dlvm::dlvm_setup_cmdstan()",
        "# System Apple Clang detected — no special config needed",
        "",
        "# Enable threading for reduce_sum parallelisation",
        "STAN_THREADS=true"
      )
      notes <- c(
        "Detected: System Apple Clang (no issues expected)",
        "Threading: enabled"
      )
    }
  } else if (os == "Linux") {
    # Linux — typically works out of the box
    config_lines <- c(
      "# Auto-configured by dlvm::dlvm_setup_cmdstan()",
      "",
      "# Enable threading for reduce_sum parallelisation",
      "STAN_THREADS=true"
    )
    notes <- c("Detected: Linux (no special config needed)", "Threading: enabled")
  } else {
    # Windows or other
    config_lines <- c(
      "# Auto-configured by dlvm::dlvm_setup_cmdstan()",
      "",
      "STAN_THREADS=true"
    )
    notes <- c(paste0("Detected: ", os), "Threading: enabled")
  }

  # Write make/local
  writeLines(config_lines, make_local)

  # Clean stale pre-compiled headers to ensure a fresh build
  pch_dir <- file.path(cmdstan_dir, "stan", "src", "stan", "model",
                       "model_header.hpp.gch")
  if (dir.exists(pch_dir)) {
    unlink(pch_dir, recursive = TRUE)
    notes <- c(notes, "Cleaned stale pre-compiled headers")
  }

  if (verbose) {
    message("[dlvm] CmdStan configured successfully:")
    for (n in notes) message("  ", n)
    message("  Config: ", make_local)
    message("")
    message("[dlvm] Next step: run dlvm_compile() to build the Stan model.")
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
