# SVMC

**Structural Variation Motif Cohesion for Bacterial Genomes**

SVMC is a pipeline and analysis framework for identifying, classifying, and
validating large-scale structural variation (SV) in complete bacterial genome
assemblies. It orients genomes to OriC using multi-marker scoring (dnaA,
DnaA-box density, GC-skew), aligns genomes with MUMmer, classifies five SV
classes (indels, duplications, translocations, inversions, and sub-structural
inversions), and identifies lineage-defining motifs using a cohesion-based
score validated against MASH distance clustering.

> **Status: experimental (v0.1.0).** The package installs and the functions
> run, but it is early research software. See *Known limitations* below before
> relying on it.

## Installation

SVMC depends on several Bioconductor packages, so install those first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Biostrings", "IRanges", "GenomicRanges", "S4Vectors",
  "GenomeInfoDb", "BiocGenerics", "rtracklayer"
))

# install.packages("remotes")
remotes::install_github("yourusername/SVMC")
```

The alignment and distance steps shell out to external command-line tools that
must be on your `PATH`:

- [MUMmer](https://github.com/mummer4/mummer) (`nucmer`, `delta-filter`)
- [Mash](https://github.com/marbl/Mash) (`mash`)
- [Prokka](https://github.com/tseemann/prokka) (optional, for gene annotation)

## Quick start

The headline contribution is the cohesion score, which measures how
genomically tight a motif's carriers are relative to the cohort:

```r
library(SVMC)

## toy cohort: 5 genomes, 2 motifs
ids <- paste0("g", 1:5)
M <- matrix(c(1, 1, 1, 0, 0,    # motif A carried by g1-g3
              1, 0, 0, 0, 1),   # motif B carried by g1 and g5
            nrow = 5, dimnames = list(ids, c("A", "B")))
D <- as.matrix(dist(c(g1 = 0, g2 = 0.01, g3 = 0.02, g4 = 0.5, g5 = 0.9)))

svmc_cohesion("A", M, D)              # ~0.97  (tight, lineage-defining)
svmc_cohesion("B", M, D)              # negative (carriers far apart)

svmc_cohesion_screen(M, D, min_cohesion = 0.40)   # the high-cohesion tail
```

A score near 1 marks a tightly clustered, likely lineage-defining motif; near
0 means carriers are no more similar than the cohort average; negative means
carriers are more divergent than average. A threshold of `0.40` separates the
lineage-defining tail in practice.

## Package layout

| File | Contents |
|---|---|
| `R/utils.R` | Small helpers (`%||%`, `coalesce_chr`, `file_exists_at_url`) |
| `R/oric.R` | OriC localization and orientation (dnaA, DnaA-box, GC-skew) |
| `R/identity.R` | MASH command generation and pairwise identity |
| `R/alignment.R` | nucmer commands, delta parsing/filtering/plotting, benchmarking |
| `R/sv_calling.R` | The five SV classes plus aggregation and grouping |
| `R/sv_size_modes.R` | SV size-mode decomposition (KDE / Gaussian mixtures) |
| `R/annotation.R` | Prokka GFF import and SV-to-gene annotation |
| `R/cohesion.R` | Motif cohesion scoring (`svmc_cohesion*`) |
| `R/doric.R` | DoriC integration helpers *(planned)* |

## Known limitations

- Analysis is designed for complete (closed) bacterial assemblies with >80% sequence similarity.

## License

MIT © 2026 Matthew D'Iorio. See [LICENSE.md](LICENSE.md).
