# SVMC 
Structural variant mapping and characterization

(This R Package is under development)

install.packages("devtools")
devtools::install_github("mdiorio371/SVMC")

library(SVMC)

ncbi_table <- readRDS("ncbi_table.rds")
species_name <- "Salmonella_enterica"

SVMC(
  species = Salmonella_enterica, 
  ncbi_table = ncbi_table
  )


  ```mermaid
flowchart TD
    subgraph Preprocessing[" "]
        A["Choose a bacterial species<br/>with multiple complete assemblies"]
        B["Load and arrange each sequence<br/>at the dnaA gene"]
        C["Verify OriC location with<br/>GC inflection and dnaA boxes"]
        A --> B --> C
    end
    style Preprocessing fill:#e6f0fa,stroke:#99bbdd,stroke-width:2px

    subgraph Alignment[" "]
        D["Choose alignment strategy"]
        D1["All-VS-all"]
        D2["References VS Queries"]
        D3["One-VS-all"]
        E["Process alignment tables"]
        D --> D1 --> E
        D --> D2 --> E
        D --> D3 --> E
    end
    style Alignment fill:#eaf7ea,stroke:#aacbaa,stroke-width:2px

    subgraph Variants[" "]
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
    style Variants fill:#fff2e0,stroke:#ddbb99,stroke-width:2px

    C --> D
    E --> F
```


### The Origin of replication can be located for an individual or set of complete genome sequences

assembly_dir <- "path/to/assembly"
load_assemblies(species_name, assembly_dir, n = 20)
locate_ori(assembly_dir)

![ori_location](ori.png)


