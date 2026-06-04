# SVMC 0.1.0

* First packaged release. The original `functions.R` is split into themed
  files under `R/` (utils, oric, identity, alignment, sv_calling,
  sv_size_modes, annotation).
* New `svmc_cohesion()` and `svmc_cohesion_screen()`: motif cohesion scoring,
  the package's headline methodology.
* Cleanup: removed superseded `annotate_SVs2/3` and three `*_prev` functions;
  `annotate_SVs4` renamed to `annotate_SVs`; duplicate `mash_commands` and
  `GC_check_R` definitions de-duplicated; inline `require()`/`library()` calls
  removed in favour of declared dependencies.
* Known gap: `oriC_scores_v2()` (called by `locate_OriC()`) is an unimplemented
  stub that errors with an explanatory message.
