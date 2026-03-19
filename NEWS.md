# dlvm 0.2.0

## Windows Cross-Platform Fixes

* **`dlvm_compile()`**: Replaced Unix-only `system("cd ... && make ...")` with
  cross-platform `system2("make", c("-C", dir, ...))`. The old approach used
  single-quoted paths and shell `cd && make` chaining that fails under Windows
  `cmd.exe`. Also adds Windows TBB `PATH` injection during compilation.

* **`dlvm_fit()`**: TBB runtime library path injection now handles Windows
  (`PATH`) in addition to macOS (`DYLD_LIBRARY_PATH`). Previously the sampler
  could fail to find `tbb.dll` on Windows even after successful compilation.

* **`dlvm_setup_cmdstan()`**: TBB library verification pattern now includes
  `.dll` alongside `.dylib` and `.so`, preventing false "missing TBB" reports
  and unnecessary auto-rebuilds on Windows.

# dlvm 0.1.0

* Initial release.
