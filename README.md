

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


### Alignment strategies summary
%% Alignment strategies schematic

```mermaid
%% Alignment strategies (vertical panels)

flowchart TB
  %% ---------- Panel A ----------
  subgraph A["A. One-vs-All"]
    direction TB
    A_R["Reference"]
    A_Q1["Query 1"]
    A_Q2["Query 2"]
    A_Q3["Query 3"]
    A_Q4["Query 4"]
    A_Qn["..."]

    A_R --> A_Q1
    A_R --> A_Q2
    A_R --> A_Q3
    A_R --> A_Q4
    A_R --> A_Qn
  end

  %% ---------- Panel B ----------
  subgraph B["B. Refs-vs-Queries"]
    direction LR

    subgraph B_L["Refs"]
      direction TB
      B_R1["Ref 1"]
      B_R2["Ref 2"]
      B_R3["Ref 3"]
      B_Rn["..."]
    end

    subgraph B_R["Queries"]
      direction TB
      B_Q1["Query 1"]
      B_Q2["Query 2"]
      B_Q3["Query 3"]
      B_Qn["..."]
    end

    %% all refs compared to all queries
    B_R1 --> B_Q1
    B_R1 --> B_Q2
    B_R1 --> B_Q3
    B_R1 --> B_Qn

    B_R2 --> B_Q1
    B_R2 --> B_Q2
    B_R2 --> B_Q3
    B_R2 --> B_Qn

    B_R3 --> B_Q1
    B_R3 --> B_Q2
    B_R3 --> B_Q3
    B_R3 --> B_Qn

    B_Rn --> B_Q1
    B_Rn --> B_Q2
    B_Rn --> B_Q3
    B_Rn --> B_Qn
  end

  %% ---------- Panel C ----------
  %% For AvA, show a clean ring + note to avoid spaghetti
  subgraph C["C. All-vs-All (every pair compared)"]
    direction TB
    C_A["Genome 1"] --- C_B["Genome 2"] --- C_C["Genome 3"]
    C_A --- C_D["Genome 4"]
    C_C --- C_E["Genome 5"]
    C_E --- C_N["..."]

    %% Minimal extra links to suggest density
    C_B --- C_D
    C_D --- C_E

    C_note["All pairs are compared; edges simplified for clarity."]
    C_A -.-> C_note
  end
```
