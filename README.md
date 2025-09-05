

# SVMC <img src="figs/logo.png" align="right" width="120"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

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


```mermaid
flowchart TD
    %% --- Preprocessing ---
    subgraph Preprocessing["Preprocessing"]
        A["Choose bacterial species<br/>with multiple assemblies"]
        B["Arrange sequences at dnaA gene"]
        C["Verify OriC location<br/>(GC skew + dnaA boxes)"]
        A --> B --> C
    end
    style Preprocessing fill:#e6f0fa,stroke:#99bbdd,stroke-width:2px

    %% --- Alignment ---
    subgraph Alignment["Alignment strategies"]
        D["Choose alignment strategy"]
        D1["All vs All"]
        D2["References vs Queries"]
        D3["One vs All"]
        E["Process alignment tables"]
        D --> D1 --> E
        D --> D2 --> E
        D --> D3 --> E
    end
    style Alignment fill:#eaf7ea,stroke:#88bb88,stroke-width:2px

    %% --- Variants ---
    subgraph Variants["Variant capture & annotation"]
        F["Capture variants"]
        F1["Indels"]
        F2["Translocations"]
        F3["Duplications"]
        F4["Inversions"]
        G["Characterize variants"]
        H["Fit Gaussian mixture model<br/>to SV length distribution"]
        I["Annotate SVs with Prokka"]
        F --> F1 --> G
        F --> F2 --> G
        F --> F3 --> G
        F --> F4 --> G
        G --> H --> I
    end
    style Variants fill:#fff2e0,stroke:#ddaa77,stroke-width:2px

    %% --- Connections ---
    C --> D
    E --> F
```


### The Origin of replication can be located for an individual or set of complete genome sequences
A confidence score is provided based on the three methods for locating the OriC: the GC inflection, DnaA box clusters, and the dnaA gene annotation.

assembly_dir <- "path/to/assembly"

# Load and analyze assemblies
load_assemblies(species_name, assembly_dir, n = 20)
locate_ori(assembly_dir)


![ori_location](ori.png)


###Annotation of gene-length SVs and enrichment snalysis

![sv_gene](sv_genes.png)
