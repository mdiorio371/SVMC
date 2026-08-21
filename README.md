# SVMC — Structural Variant Mapping and Characterization

An R package for identifying, classifying, and validating large-scale
structural variation in complete bacterial genome assemblies.

## Installation

```r
# install.packages("devtools")
devtools::install_github("YOURUSERNAME/SVMC")
```

Or load locally during development:

```r
devtools::load_all("path/to/SVMC")
```

## Pipeline overview

SVMC processes complete bacterial genomes through the following stages:

1. **Download / inventory** — sync complete assemblies from RefSeq
2. **OriC locating** — multi-marker origin scoring (dnaA, GC-skew, DnaA-box density)
3. **Identity** — MASH distances, medoid reference selection
4. **Alignment** — MUMmer/nucmer pairwise alignment, delta parsing
5. **SV calling** — indels, duplications, translocations, RASRs, sub-inversions
6. **Size modes** — KDE/GMM peak detection on SV size distributions
7. **Annotation** — Prokka CDS annotation, breakpoint context
8. **Motif grouping** — reciprocal-overlap + Jaccard clustering into recurrent motifs
9. **RASR / phi scoring** — exclusivity-based lineage-association scoring
10. **Cohesion** — MASH-distance motif cohesion (lineage-defining candidate identification)
11. **Ordination** — PCoA / DBSCAN state assignment
12. **Plots** — synteny panels, pathway diagrams

## Key functions

- `svmc_cohesion()` / `svmc_cohesion_screen()` — motif cohesion scoring
- `phi_from_abcd()` / `score_event_pairs()` / `rasr_scoreboard()` — phi exclusivity screen
- `oriC_scores_v2()` — multi-marker OriC confidence scoring
- `delta_all_SVs()` / `get_all_SVs()` — full SV calling from delta alignments
- `group_svs()` — motif grouping by reciprocal overlap
- `read_delta()` / `filter_delta()` — MUMmer delta handling

## Two association metrics

SVMC computes two complementary, independent metrics:

- **phi coefficient** (`phi_from_abcd`): exclusivity from the 2×2 motif co-occurrence
  table. This is the **lineage-association score** (max|phi| across exclusive pairs).
- **cohesion** (`svmc_cohesion`): MASH-distance concentration of a motif's carriers,
  `1 - (mean within-carrier distance / mean overall distance)`.

## Citation

D'Iorio, M. (2026). Chromosomal architecture and bacterial evolution:
Mapping and contextualizing large-scale structural variants in bacterial
phylogeny. PhD thesis, McGill University.

## License

MIT
