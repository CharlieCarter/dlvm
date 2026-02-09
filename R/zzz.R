# dlvm/R/zzz.R
# Package load/attach hooks

# Internal environment for caching compiled model across calls
.dlvm_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Nothing else needed — .dlvm_env is created at package load time
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    sprintf("dlvm %s — Dynamic Latent Variable Model",
            utils::packageVersion("dlvm"))
  )
}
