[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Status](https://img.shields.io/badge/status-work--in--progress-orange)
# SVMC <img src="figs/logo.png" align="right" width="120"/>


**SVMC** (Structural Variant Mapping and Characterization) is an R package for analyzing bacterial structural variation.  
It provides tools to locate origins of replication, align complete assemblies, parse alignments for structural variants (SVs), model SV length distributions, and annotate breakpoint contexts. 


---

## 🔧 Installation

```r
##(This package is still in development)
# install.packages("devtools")
devtools::install_github("mdiorio371/SVMC")

library(SVMC)


ncbi_table <- readRDS("ncbi_table.rds")
species_name <- "Salmonella_enterica"

# Run the SVMC workflow
SVMC(
  species    = species_name,
  ncbi_table = ncbi_table
)

```  





### The Origin of replication can be located for an individual or set of complete genome sequences
A confidence score is provided based on the three methods for locating the OriC: the GC inflection, DnaA box clusters, and the dnaA gene annotation.

assembly_dir <- "path/to/assembly"

# Load and analyze assemblies
load_assemblies(species_name, assembly_dir, n = 20)
locate_ori(assembly_dir)


![ori_location](ori.png)


### Annotation of gene-length SVs and enrichment analyses

![sv_gene](figures/sv_gene_enrichment.png)
