list.of.packages <-
    c(
        "tidyverse", "data.table", "tidyr", 
        "lubridate", "ggplot2",
        "curl", "RColorBrewer", "ape",
        "pbmcapply", "igraph", "viridis",
        "gridExtra", "gridExtra", "scales", 
        "knitr", "RCurl", "kableExtra", 
        "GenomicRanges",
        "data.tree", "Matrix", "devtools",
        "openxlsx", "FactoMineR", "ggtree", "rcartocolor",
        "Biostrings", "ggtree", "dendextend", "vegan",
        #, "msa", 
        "R.utils", "MASS", "rentrez",
        "boot",  "ggdendro",
        #"car", "modelr", "broom", "ggpubr", "cowplot",
        "ggridges", "factoextra", "cluster", "boot",
        "zoo"
    )
new.packages <- 
    list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

if (requireNamespace("BiocManager", quietly = TRUE)) { 
    install.packages("BiocManager")
}
if (!("Biostrings" %in% installed.packages())){
    BiocManager::install("Biostrings")
}
if (!("ggtree" %in% installed.packages())){
    BiocManager::install("ggtree")
}




library(tidyverse)
library(data.table)
library(Biostrings)
library(pbmcapply)
library(igraph)
library(curl)
library(RCurl)
library(ape)
library(openxlsx)
#library(msa)
library(FactoMineR)
library(factoextra)
library(cluster)
library(vegan)
library(zoo)



#colours/plots
library(rcartocolor)
library(RColorBrewer)
library(viridis)  
library(gridExtra)
library(knitr)
library(lubridate)
library(scales)
library(kableExtra)
library(Matrix)
library(kableExtra)
library(ggtree)
library(ggridges)
library(ggdendro)
library(dendextend)
library(boot)



library(R.utils)




### data loading
get_ncbi_table <- 
    function(
        db_dir = "db"
    ){
        ## type = c("refseq", "Genbank", "Both")
        if (!dir.exists(db_dir)){dir.create(db_dir)}
        ncbi_dir <- 
            "ftp://ftp.ncbi.nlm.nih.gov/genomes"
        # refseq
        refseq_file <- 
            sprintf(
                "%s/refseq/bacteria/assembly_summary.txt",
                ncbi_dir
            )
        
        # genbank genomes 
        genbank_file <- 
            sprintf(
                "%s/genbank/bacteria/assembly_summary.txt",
                ncbi_dir
            )
        refseq <- 
            fread(
                refseq_file, 
                stringsAsFactors = F, 
                quote = "", 
                select = c(1,5,6,7,8,11,12,15,16,17,18,19,20)
            ) %>% 
            filter(
                assembly_level =="Complete Genome" &
                    version_status =="latest"
            ) %>%
            `colnames<-`(c("assembly", colnames(.)[-1])) %>%
            as_tibble %>%
            mutate(db = "refseq")
        ## Access the Genbank genomes that are not in refseq
        genbank <- 
            fread(
                genbank_file, 
                stringsAsFactors = F, 
                quote = "", 
                select = c(1,5,6,7,8,11,12,15,16,17,18,19,20)
            ) %>% 
            filter(
                assembly_level =="Complete Genome" &
                    version_status =="latest"
            ) %>%
            `colnames<-`(c("assembly", colnames(.)[-1])) %>%
            as_tibble %>%
            filter(paired_asm_comp != "identical") %>%
            mutate(db = "genbank")
        
        
        
        # refseq & genbank combined  
        bac_summary <- 
            bind_rows(refseq, genbank) %>%
            mutate(
                fna_path = 
                    sprintf(
                        "%s/%s_genomic.fna.gz",
                        ftp_path, gsub(".*/", "", ftp_path)
                    )
            )
        head(bac_summary)
        # get taxonomy information
        ncbi_url <- 
            "ftp://ftp.ncbi.nlm.nih.gov"
        taxa_file <- 
            sprintf(
                "%s/pub/taxonomy/new_taxdump/new_taxdump.tar.gz",
                ncbi_url
            )
        
        taxa_con <- 
            curl_download(
                taxa_file,
                destfile = 
                    sprintf(
                        "%s/new_taxdump.tar.gz",
                        db_dir
                    )
            )
        ### extract the ranked lineage table
        untar(taxa_con, exdir = db_dir, files = "rankedlineage.dmp")
        ## read it in
        taxa_table <- 
            read_tsv(
                sprintf(
                    "%s/rankedlineage.dmp",
                    db_dir
                ), 
                col_names = 
                    c(
                        "taxid", "name", "s", "g", "f", "o", "c","p", "k", "d"
                    ),
                col_types=("i-c-c-c-c-c-c-c-c-c-")
            ) %>% 
            filter(d=="Bacteria")
        ## Combine the tables
        
        bac_taxa <- 
            left_join(bac_summary, taxa_table, by = "taxid") %>% 
            mutate(
                p2 = 
                    ifelse(
                        p=="Proteobacteria", 
                        sub("proteobacteria", "",c), 
                        p
                    )
            ) %>%
            dplyr::select(-c(version_status, assembly_level, k, d)) %>%
            mutate(
                species = 
                    unlist(
                        pbmclapply(
                            name,
                            function(x){
                                spec_out <-
                                    strsplit(x, " ")[[1]][1:2] %>% 
                                    paste(collapse = " ")
                                spec_out_2 <-
                                    strsplit(spec_out, " ")[[1]][2]
                                if (spec_out_2=="sp."){
                                    spec_out <-
                                        x
                                }
                                return(spec_out)
                            }
                        )
                    )
            ) 
        
        bt_summary <- 
            bac_taxa %>%
            count(p, g, species) %>% 
            arrange(desc(n)) %>% 
            filter(n>=10)
        
        info_tib <- 
            tibble(
                "Complete sequences" =
                    sum(bt_summary$n),
                "Species with complete_sequences >= 10"  = 
                    nrow(bt_summary),
                downlaod_date = 
                    Sys.Date()
            )
        ncbi_table_clade_count <- 
            bac_taxa %>% 
            add_count(species, name = "spec_count") %>%
            filter(spec_count>=10) %>%
            group_by(p2, species) %>%
            summarise(spec_count = n()) %>%
            ungroup %>%
            mutate(genus = sub("([A-Za-z]+).*", "\\1", species)) %>%
            group_by(p2, genus) %>%
            mutate(genus_count = sum(spec_count)) %>%
            arrange(desc(genus_count), desc(spec_count)) %>%
            mutate(species_ = sub(" ", "_", species)) %>%
            dplyr::rename(phylum = p2) %>%
            ungroup
        
        out_list <- 
            list(
                ncbi_table = 
                    bac_taxa,
                ncbi_info = 
                    info_tib,
                ncbi_cc = 
                    ncbi_table_clade_count
            )
        
        return(out_list)
    }

#### Get report summary and species summary  ####
# 2 functions
get_report_table <-
    function(species_table, core_number=7){
        
        report_details <-
            pbmclapply(
                1:nrow(species_table),
                function(i){
                    
                    ftp_path <-
                        (species_table %>% pull(ftp_path))[i]
                    if (ftp_path=="na" ){
                        report_details <- 
                            tibble(
                                asm_name = NA,
                                org_name = NA,
                                submitter = NA,
                                accession = NA, 
                                len = NA, 
                                strain = NA,
                                date_sub = NA,
                                technology = NA
                            )
                        out_tib <- 
                            report_details
                    } else {
                        rep_file <- 
                            sprintf(
                                "%s/%s_assembly_report.txt",
                                ftp_path, gsub(".*/", "", ftp_path)
                            ) 
                        if (!url.exists(rep_file)){
                            rep_file <- 
                                sprintf(
                                    "%s/%s_%s/%s_%s_assembly_report.txt",
                                    dirname(ftp_path), 
                                    pull(species_table[i,], assembly),
                                    pull(species_table[i,], asm_name),
                                    pull(species_table[i,], assembly),
                                    pull(species_table[i,], asm_name)
                                ) 
                        }
                        
                        rep_read <-
                            tryCatch(readLines(rep_file))
                        
                        strain <- 
                            rep_read[
                                grepl("# Infraspecific name:", rep_read)
                            ] %>% 
                            sub(".*strain=", "", .)
                        if (length(strain)==0){
                            strain <- NA
                        }
                        technology <- 
                            rep_read[grepl("# Sequencing", rep_read)] %>% 
                            sub(".*technology:", "", .) %>% 
                            trimws
                        if (length(technology)==0){
                            technology <- NA
                        }
                        out_tib <- 
                            tibble(
                                asm_name = 
                                    sub(".*: ", "", rep_read[1]) %>%
                                    trimws,
                                org_name = 
                                    sub(".*: ", "", rep_read[2]) %>%
                                    trimws,
                                submitter = 
                                    rep_read[grepl("Submitter", rep_read)] %>% 
                                    sub(".*Submitter:", "", .) %>% 
                                    trimws,
                                accession = 
                                    rep_read[grepl("Chromosome", rep_read)] %>% 
                                    sub(".*=\\t", "", .) %>% 
                                    sub("\t.*", "", .),
                                len = 
                                    rep_read[grepl("Chromosome", rep_read)] %>% 
                                    sub(".*Assembly\\t", "", .) %>%
                                    sub("\t.*", "", .),
                                strain = 
                                    strain,
                                date_sub = 
                                    rep_read[grepl("Date", rep_read)] %>% 
                                    sub(".*Date:", "", .) %>% 
                                    trimws,
                                technology = 
                                    technology
                            )
                    }
                    
                    return(out_tib)
                }, mc.cores = core_number
            ) %>%
            bind_rows()
        return(report_details)
    }

get_species_summary <- 
    function(
        species_table
    ){
        # species_table <- 
        #     ncbi_table %>%
        #     filter(species== sub("_", " ", species_name))
        #get_report_table(species_table[1723,])
        report_details <-
            get_report_table(species_table) %>%
            mutate(
                len = as.numeric(len),
                accession = make.unique(accession)
            ) %>% 
            distinct
        if (nrow(report_details)>nrow(species_table)){
            report_details <- 
                report_details %>%
                group_by(asm_name) %>%
                filter(!is.na(len)) %>%
                filter(len == max(len, na.rm = T)) %>% ungroup()
        }
        
        report_details <- 
            report_details %>%
            filter(
                len < mean(len, na.rm = T)+1e6 & len > mean(len, na.rm = T)-1e6
            ) 
        
        ## join the report table with the filtered table
        species_table <-   
            left_join(report_details, species_table, by = "asm_name") %>%
            mutate(len = as.numeric(len)) %>%
            filter(!duplicated(asm_name)) %>%
            dplyr::select(
                -c(
                    refseq_category, s, 
                    gbrs_paired_asm, paired_asm_comp, 
                    species_taxid
                )
            )
        return(species_table)
        
    }


## oriC locating

oriC_scores <- function(
        dbox, GC_tib, dgene, genome_len,
        window_size = 20000, step_size = 5000,
        decay_constant = 10000, gc_valley_scale = 50000,
        weights = c(box = 1, gc = 1.5, gene = 1.2),
        thin_distance = 100000
) {
    
    
    # --- Helper: Normalize with zero protection
    rescale01 <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (diff(rng) == 0) return(rep(0, length(x)))
        (x - rng[1]) / diff(rng)
    }
    
    # --- GC skew processing ---
    spline_fit <- smooth.spline(GC_tib$loc, GC_tib$GC, spar = 0.5)
    GC_tib <- GC_tib %>%
        mutate(
            smoothed = predict(spline_fit, x = loc)$y,
            curvature = predict(spline_fit, deriv = 2)$y
        )
    gc_global_min <- GC_tib$loc[which.min(GC_tib$smoothed)]
    
    # --- Identify valleys ---
    valley_locs <- 
        GC_tib %>% 
        filter(
            curvature > quantile(curvature, 0.75) & 
                smoothed < quantile(smoothed, 0.25)
            ) %>% 
        pull(loc)
    n_gc_valleys <- length(valley_locs)
    gc_confidence_boost <- ifelse(n_gc_valleys <= 3, 1.25, 1)
    
    # --- Process dnaA box signal ---
    box_density <- dbox %>%
        mutate(loc_bin = round(loc, -3)) %>%
        count(loc_bin, wt = b) %>%
        tidyr::complete(loc_bin = seq(1, genome_len, by = 1000), fill = list(n = 0)) %>%
        arrange(loc_bin) %>%
        mutate(smoothed = zoo::rollmean(n, k = 5, fill = 0))
    
    # Identify high box clusters
    box_peaks <- box_density %>% filter(smoothed > quantile(smoothed, 0.95))
    n_box_peaks <- nrow(box_peaks)
    box_confidence_boost <- ifelse(n_box_peaks <= 5, 1.2, 1)
    
    # --- Score windows ---
    starts <- seq(1, genome_len - window_size, by = step_size)
    ends <- starts + window_size - 1
    mids <- round((starts + ends) / 2)
    
    scores <- map2_dfr(starts, ends, ~{
        win <- IRanges(start = .x, end = .y)
        mid <- round((start(win) + end(win)) / 2)
        
        box_score <- approx(box_density$loc_bin, box_density$smoothed, xout = mid, rule = 2)$y
        curve_val <- approx(GC_tib$loc, GC_tib$curvature, xout = mid, rule = 2)$y
        gc_depth <- -approx(GC_tib$loc, GC_tib$smoothed, xout = mid, rule = 2)$y
        dist_to_gc_min <- min(abs(mid - gc_global_min), genome_len - abs(mid - gc_global_min))
        
        # Reduced penalty on distance to GC min
        gc_valley_score <- (0.7 * max(0, gc_depth) + 0.3 * max(0, curve_val)) * exp(-dist_to_gc_min / (gc_valley_scale * 2))
        
        gene_dist <- min(abs(mid - dgene), genome_len - abs(mid - dgene))
        gene_score <- exp(-gene_dist / decay_constant)
        
        tibble(
            start = start(win),
            end = end(win),
            mid = mid,
            box_score,
            gc_valley_score,
            gene_score
        )
    })
    
    scores <- scores %>%
        mutate(
            box_z  = rescale01(box_score),
            gc_z   = rescale01(gc_valley_score),
            gene_z = rescale01(gene_score),
            agreement_strength = box_z * gc_z * gene_z,
            ori_score = (
                box_confidence_boost * weights["box"] * box_z +
                    gc_confidence_boost * weights["gc"] * gc_z +
                    weights["gene"] * gene_z
            ) / sum(weights)
        ) %>%
        arrange(desc(ori_score))
    
    # Thin overlapping regions
    keep <- rep(TRUE, nrow(scores))
    for (i in seq_len(nrow(scores))) {
        if (!keep[i]) next
        keep <- keep & (abs(scores$mid - scores$mid[i]) > thin_distance | seq_len(nrow(scores)) <= i)
    }
    
    scores <- scores[keep, ] %>%
        mutate(
            rank = row_number(),
            confidence_tier = 
                case_when(
                    ori_score > 0.65 ~ "High",
                    ori_score > 0.3 ~ "Moderate",
                    TRUE ~ "Low"
                )
        )
    
    return(scores)
} 


plot_OriC_markers <- 
    function(
        GC_disp,
        dnaA_box_tib,
        fasta_len,
        top_regions
    ){
        
        dnaA_gene_box_GC_disp <- 
            ggplot(
                GC_disp, 
                aes(x= Loci, y = `G-C`)
            ) +
            geom_rect(
                data = top_regions,
                aes(
                    xmin = (start / 1e4)-0.5, 
                    xmax = (end / 1e4)+0.5, fill = confidence_tier
                    ),
                ymin = -Inf, ymax = Inf,
                color = "grey80", alpha = 0.6,
                inherit.aes = FALSE
            ) +
            geom_point(size = 0.1) +  
            geom_segment(
                data = dnaA_box_tib , 
                aes(
                    x = midpt/1e4, 
                    xend = midpt/1e4, 
                    y = 
                        min(GC_disp$`G-C`)-25e2, 
                    yend = 
                        min(GC_disp$`G-C`) + 
                        (b * ((max(abs(GC_disp$`G-C`))-6e3)/max(b)))-25e2
                    #1e6)  
                ), 
                color = "grey", linewidth = 1
            ) +
            geom_point(
                data = dnaA_box_tib, 
                aes(
                    x = midpt / 1e4, 
                    y = min(GC_disp$`G-C`) + 
                        (b * ((max(abs(GC_disp$`G-C`)) - 6e3) / max(b)))-25e2
                ), 
                color = "blue", size = 2,  alpha = 0.4
            ) +
            geom_line() +
            labs(
                x = "Genome sequence loci (10 Kb)", 
                y = "G-C disparity", 
                #title = ""
            ) +
            annotate(
                "segment",
                x = round(fasta_len / 2e4), xend = round(fasta_len / 2e4),  
                y = 
                    #0,
                    min(GC_disp[["G-C"]]) -25e2,
                yend = 
                    #-5e3,
                    min(GC_disp[["G-C"]]) -55e2,
                color = "black",
                linewidth = 1.5,
                arrow = 
                    arrow(
                        ends = "first", type = "closed", 
                        length = unit(0.5, "cm")
                    )
            ) +
            scale_fill_manual(values = tier_colors) +
            
            geom_text(
                data = top_regions,
                aes(x = mid / 1e4, y = max(GC_disp$`G-C`) * 0.95, label = rank),
                color = "black", fontface = "bold", size = 8,
                inherit.aes = FALSE
            ) +
            labs(
                title = 
                    paste0(
                        "OriC Confidence= ",
                        round(top_regions$ori_score[1], 2),
                        "\n (",
                        top_regions$confidence_tier[1], ")"
                    )
            ) +
            theme_classic() +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold"),
                text = element_text(size=22),
                legend.position = "none"
            )
        
        return(dnaA_gene_box_GC_disp)
        
    }

locate_OriC <- 
    function(
        fasta_file,
        disparity_file,
        fasta_len
    ){
        
        genome_seq <- 
            readDNAStringSet(
                fasta_file
            )
        
        rearranged_seq <- 
            xscat(
                subseq(genome_seq, start = round(fasta_len/2), end = fasta_len),
                subseq(genome_seq, start = 1, end = round(fasta_len/2)-1)
            )
        
        ### dnaA box check
        genome_seq <- 
            rearranged_seq
        
        dnaA_box_tib <- 
            id_dnaA_boxes(
                genome_seq
            )
        
        
        dbox <- 
            dnaA_box_tib %>%
            transmute(
                loc = midpt, 
                b
            )
        
        ### GC check
        writeXStringSet(
            rearranged_seq,
            sprintf(
                "%s/%s.txt",
                oriC_dir, 
                st2$accession[j]
            )
        )
        GC_info <- 
            GC_check(
                sprintf(
                    "%s/%s.txt",
                    oriC_dir, 
                    st2$accession[j]
                ), 
                disparity_file
            )
        
        GC_disp <- 
            fread(
                disparity_file
            ) %>%
            mutate(
                Loci = 
                    row_number()
            )
        
        GC_tib <- 
            GC_disp %>%
            transmute(
                GC = `G-C`,
                loc = Loci*1e4
            )
        
        #dbox
        #GC_tib
        dgene <- 
            round(fasta_len/2)
        
        
        ori_scores <- 
            oriC_scores(
                dbox = dbox,
                GC_tib = GC_tib,
                dgene = round(fasta_len / 2),
                genome_len = fasta_len
            ) 
        
        tier_colors <-
            c(
                High = "#B2D8B2",    # faint green
                Moderate = "#FFCC66",# amber
                Low = "#D3D3D3"      # light grey
            )
        top_regions <- 
            ori_scores %>% 
            dplyr::slice(1:3) %>%#
            mutate(
                rank_label = as.character(rank),
                fill_color = tier_colors[confidence_tier]
            )
        
        dnaA_gene_box_GC_disp <- 
            plot_OriC_markers(
                GC_disp,
                dnaA_box_tib,
                fasta_len,
                top_regions
            )
        
        out_list <- 
            list(
                top_regions = top_regions, 
                OriC_plot = dnaA_gene_box_GC_disp
            )
        return(out_list)
    }

## pairwise alignment


## Structural variant capture


## structural variant annotation





