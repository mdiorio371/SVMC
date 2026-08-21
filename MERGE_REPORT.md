# SVMC Merge Report

Final consolidated package. 177 unique functions, zero duplicates.

## How the merge worked

Two package attempts were reconciled:

1. **Documented files** (utils.R, oric.R, identity.R, alignment.R, sv_calling.R,
   sv_size_modes.R, annotation.R, cohesion.R, doric.R) — an earlier, roxygen-documented
   consolidation you had already started. These WIN wherever they overlap.

2. **Full extraction** (from SVCM.R, ch4_functions.R, SVCM_result1_helpers.R, functions.R,
   SVCM6.R, SVM7.R, i9.R, and the Rmd notebooks) — complete coverage of all 12 pipeline
   stages. Fills every gap the documented files don't cover.

For 16 overlapping functions, the documented version was kept (files ending `b-...-documented.R`).
The 161 extraction-only functions were all retained. 3 new functions came from the
documented files (oriC_scores_v2, svmc_cohesion, svmc_cohesion_screen).

## File layout

Numbered files (00, 01, ...) = extraction, by pipeline stage.
Files with `b` suffix (00b, 02b, ...) = documented versions (roxygen), same stage.
13-cohesion.R, 14-doric.R = new stages from the documented set.

The `b` files and numbered files coexist without collision (dedup guaranteed).
You may later merge each `Nb` into its `N` file for tidiness, but it is not required —
R sources all files in R/ regardless of name.

## RESOLVED: the phi vs cohesion question

Earlier we believed cohesion was never a separate metric. THAT WAS WRONG.
cohesion.R contains real, documented cohesion functions:

    svmc_cohesion(motif_uid, M_all, dist_matrix)
      = 1 - (mean within-carrier MASH dist / mean overall MASH dist)

This is DISTINCT from phi (phi_from_abcd). You have TWO metrics:
  - phi        : exclusivity (2x2 co-occurrence)  -> lineage-association score
  - cohesion   : MASH-distance concentration of carriers

### ACTION REQUIRED before defense
Determine which metric produced the thesis numbers 0.978 / 0.863 / 0.437.
Run BOTH on the Salmonella data:
  svmc_cohesion_screen(M_all, dist_matrix)          # cohesion
  rasr_scoreboard(...)                               # phi
Whichever reproduces 0.978/0.863/0.437 is the one the thesis reports.
The thesis text was edited to call them "lineage-association score (max|phi|)".
If they are actually COHESION scores, that labeling is wrong and needs an erratum
or a defense-time clarification.

## Known remaining tasks (post-load)

1. Confirm `devtools::load_all()` succeeds (177 functions).
2. Resolve phi vs cohesion labeling (above) — highest priority.
3. Implement doric.R (currently a documented stub; the DoriC offset helpers).
4. Add roxygen to the ~150 undocumented functions (the `b` files show the target style).
5. Consider merging each `Nb-documented.R` into its `N` file.
6. Triage 99-misc.R (parsnp staging helpers — keep or cut).
7. NAMESPACE currently exports everything via exportPattern; refine with @export later.

## Namespace masking note

On load you will see warnings about dplyr masking S4Vectors/IRanges (first, desc,
slice, rename, union, intersect, etc.). These are harmless — dplyr wins, which is what
the code expects. If any Bioconductor-object operation errors at runtime (most likely
IRanges::slice or GenomicRanges::intersect in annotation/overlap code), fully-qualify
that specific call.
