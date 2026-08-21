# SVMC Consolidation Report

Generated during package assembly from 10 disparate pipeline files.

## Source files and recency ranking (highest wins on collision)

| Rank | File | Role |
|------|------|------|
| 5 (newest) | SVCM.R | Main current pipeline — owns core SV calling (stages 1-8) |
| 4 | ch4_functions.R | Chapter 4 RASR / phi exclusivity scoring |
| 4 | SVCM_result1_helpers.R | Jaccard motif grouping, states, ordination helpers |
| 3 | main.R | entry script (no functions extracted) |
| 2 | SVCM6.R, SVM7.R, i9.R | older iterations — fallback functions only |
| 1 (oldest) | functions.R | original core library — superseded except GC_check_R, compute_disparities, recursive_gmm |
| notebooks | Ch2.Rmd, SVCM.Rmd | 17 functions recovered that existed nowhere else |

## Package structure (R/ files by pipeline stage)

| File | Stage |
|------|-------|
| 00-utils.R | Setup, cache, dirs, id helpers |
| 01-download.R | Genome download / inventory / fasta staging |
| 02-oric.R | OriC locating (dnaA, GC-skew, disparities) |
| 03-identity.R | MASH identity, medoid reference selection |
| 04-alignment.R | nucmer alignment, delta read/filter/plot |
| 05-sv-calling.R | SV callers (indels, dup, transloc, RASR, sub-inversions) |
| 06-size-modes.R | KDE/GMM size-mode detection |
| 07-annotation.R | Prokka annotation, SV annotation, breakpoints |
| 08-motif-group.R | Reciprocal-overlap + Jaccard motif grouping |
| 09-rasr-phi.R | RASR classification + phi exclusivity scoring (Ch4 method) |
| 10-ordination.R | PCoA, DBSCAN, state assignment |
| 11-plots.R | Synteny panels, pathway plots |
| 12-orchestration.R | step1/step2/step3 pipeline drivers |
| 98-from-notebooks.R | Functions recovered from Rmd notebooks |
| 99-misc.R | Unassigned helpers (parsnp staging, inspection) — TRIAGE |

## Dead code EXCLUDED from package

annotate_SVs2, annotate_SVs3 (superseded by annotate_SVs = former annotate_SVs4)
delta_duplications_prev, get_pairwise_identities_prev, mash_commands_prev
rebuild_stepc_motif_events / _2 / _3 (experimental)

## KNOWN ISSUES TO RESOLVE BEFORE PUBLICATION

1. **phi vs Jaccard confirmed OK**: lineage-association score = phi (ch4_functions.R:
   phi_from_abcd -> score_event_pairs -> rasr_scoreboard). Jaccard is used SEPARATELY
   for motif grouping (08-motif-group.R). Thesis Methods are correct.

2. **Fallback functions from older files** (SVCM6.R annotation/alignment steps):
   verify step2_align_call_resumable and the svcm_step3_* CDS functions are the versions
   you actually run. They fill gaps SVCM.R doesn't cover but weren't cross-checked.

3. **R/99-misc.R needs triage**: 14 functions of uncertain relevance (mostly parsnp
   staging + inspection utilities). Decide keep/cut before release.

4. **No documentation yet**: functions have no roxygen. NAMESPACE exports everything
   via exportPattern. Add @export + @param incrementally.

5. **oriC_scores_v2**: referenced in some call sites but not found in any source file.
   Locate or recreate before the OriC stage is fully runnable.

6. **cohesion() as a named function**: the thesis describes cohesion conceptually but
   the metric is phi-based (see #1). No separate cohesion() function needs writing —
   confirmed cohesion was always describing the phi/MASH-concentration result.

## Next steps for a runnable package

1. devtools::load_all("SVMC")  — fix any parse/load errors iteratively
2. Add Imports to DESCRIPTION as load errors reveal missing packages
3. Triage 99-misc.R
4. Add roxygen docs (start with the ~15 user-facing functions)
5. Write vignette: Salmonella worked example
6. devtools::check()
