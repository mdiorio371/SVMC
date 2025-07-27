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
# library(car)
# library(modelr)
# library(broom)
# library(ggpubr)
# library(cowplot)


library(R.utils)
# library(MASS)
# library(rentrez)
# library(GenomicRanges)
# library(data.tree)


#### NCBI download ####
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


#### Downlaod and annotate the sequences
# 2 functions
download_fna <- 
    function(species_table, fna_dir, core_number=7){
        fna_table <- 
            pbmclapply(
                1:nrow(species_table),
                function(i){
                    ## get the gff path
                    asm_name <- 
                        pull(species_table, asm_name)[i]
                    fna_path <- 
                        (species_table %>% 
                             pull(ftp_path))[i] %>%
                        sprintf(
                            "%s/%s_genomic.fna.gz",
                            ., gsub(".*/", "", .)
                        )
                    accession <- 
                        pull(species_table, accession)[i]
                    ## Set the download path
                    fna_outfile <- 
                        sprintf(
                            "%s/%s.fna.gz",
                            fna_dir, accession
                            #basename(fna_path)
                        )
                    # if it's there, don't downlaod it
                    if (file.exists(fna_outfile)){
                        fna_tib <- 
                            tibble(
                                fna = basename(fna_outfile),
                                fna_exists = T,
                                asm_name = asm_name
                            )
                    } else{ # if not, try downloading it
                        test_downlaod <- 
                            try(
                                curl_download(
                                    url = 
                                        fna_path,
                                    destfile = 
                                        fna_outfile
                                ), 
                                silent = T
                            )
                        fna_tib <- 
                            tibble(
                                fna = basename(fna_outfile),
                                fna_exists = !grepl("Error", test_downlaod),
                                asm_name = asm_name
                            )
                    }
                    return(fna_tib)
                }, mc.cores = core_number
            ) %>% 
            bind_rows()
        return(fna_table)
    }

annotate_fna <- 
    function(species_name, fna_dir, gff_dir){
        fna_files <- 
            list.files(fna_dir, full.names = T)
        prokka_coms <- 
            pbmclapply(
                1:length(fna_files), 
                function(i){
                    prokka_com <- 
                        sprintf(
                            "prokka --species '%s' --outdir %s --prefix %s %s %s",
                            gsub("_", " ", species_name),
                            gff_dir, 
                            gsub(".fna", "", basename(fna_files))[i],
                            fna_files[i],
                            "--force"
                        )
                    return(prokka_com)
                }
            ) %>% 
            unlist
        prokka_file <- 
            sprintf("%s/prokka_coms.txt", prokka_dir)
        
        prokka_out <- 
            list(
                file = prokka_file,
                comms = prokka_coms
            )
        return(prokka_out)
    }

annotate_fna2 <- 
    function(species_name, sync_dir, prokka_dir, accessions){
        sfs <- 
            list.files(sync_dir, full.names = T) 
        sfs <- 
            sfs[
            sub(".txt", "", basename(sfs)) %in% accessions
            ]
        
        accessions
        prokka_coms <- 
            pbmclapply(
                1:length(sfs), 
                function(i){
                    # prokka_check <- 
                    #     sprintf(
                    #         "if [ ! -e %s ]",
                    #         filter_out
                    #     )
                    prokka_com <- 
                        sprintf(
                            "prokka --species '%s' --outdir %s --prefix %s %s %s",
                            gsub("_", " ", species_name),
                            prokka_dir, 
                            gsub(".txt", "", basename(sfs))[i],
                            sfs[i],
                            "--force"
                        )
                    return(prokka_com)
                }
            ) %>% 
            unlist
        prokka_file <- 
            sprintf("%s/prokka_coms.txt", prokka_dir)
        
        prokka_out <- 
            list(
                file = prokka_file,
                comms = prokka_coms
            )
        return(prokka_out)
    }

annotate_ind <- 
    function(fasta_file, species_name, prokka_dir){

        prokka_com <- 
            sprintf(
                "prokka --species '%s' --outdir %s --prefix %s %s %s",
                gsub("_", " ", species_name),
                prokka_dir, 
                gsub(".txt", "", basename(fasta_file)),
                fasta_file,
                "--force"
            )
                    
        prokka_file <- 
            sprintf("%s/prokka_coms.txt", prokka_dir)
        
        prokka_out <- 
            list(
                file = prokka_file,
                comms = prokka_com
            )
    }


#### Synchronize ####

sync_dnaA <- 
    function(
        assembly_path, gff_dir, sync_dir,
        save_gene = T
    ){
        # assembly_path <- 
        #     sync_summary$asm_name
        # the fasta file
        fasta_path <- 
            sprintf(
                "%s/%s_genomic.fna.gz",
                assembly_path, 
                gsub(".*/", "", assembly_path)
            )
        
        asm_name <- 
            basename(assembly_path)
        
        #check if it exists
        fna_exists <- 
            file_exists_at_url(fasta_path)
        
        if (!fna_exists){
            out_tib <- 
                tibble(
                    asm_name = 
                        basename(assembly_path),
                    accession = 
                        genome_acc,
                    ori_pos = 
                        NA,
                    ori_strand = 
                        NA,
                    len = 
                        genome_length, 
                    out_file = NA,
                    note = "no_Fasta"
                )
        } else{
       
        
            genome_seq <- 
                readDNAStringSet(fasta_path)
            #showConnections(all = T)
            
            genome_seq <- 
                genome_seq[
                    !(grepl("plasmid", names(genome_seq),  
                            ignore.case = T))
                ]
            genome_seq <- 
                genome_seq[
                    !(grepl("phage", names(genome_seq),  
                            ignore.case = T)) |
                        grepl(
                            "phage resistant", 
                            names(genome_seq),  
                            ignore.case = T
                        )
                ]
            
            ## check chromosomes
            genome_seqs <- 
                genome_seq[order(width(genome_seq), decreasing = T),]
            
            #names(genome_seqs)
            chrom_num <- 
                length(genome_seqs)
            
            genome_length <- 
                width(genome_seq)
            
            genome_acc <- 
                strsplit(names(genome_seq), " ")[[1]][1]
            
            genome_genus <- 
                strsplit(names(genome_seq), " ")[[1]][2] 
            
            ## download the gff
            gff_path <- 
                assembly_path %>%
                sprintf(
                    "%s/%s_genomic.gff.gz",
                    ., gsub(".*/", "", .)
                )
            ## Set the download path
            gff_outfile <- 
                sprintf(
                    "%s/%s",
                    gff_dir,
                    basename(gff_path)
                )
            # if it's there, don't downlaod it
            if (file.exists(gff_outfile)){
                gff_tib <- 
                    tibble(
                        gff = basename(gff_outfile),
                        gff_exists = T,
                        asm_name = asm_name
                    )
            } else{ # if not, try downloading it
                test_downlaod <- 
                    try(
                        curl_download(
                            url = 
                                gff_path,
                            destfile = 
                                gff_outfile
                        ), 
                        silent = T
                    )
                gff_tib <- 
                    tibble(
                        gff = basename(gff_outfile),
                        gff_exists = !grepl("Error", test_downlaod),
                        asm_name = asm_name
                    )
            }
            if (gff_tib$gff_exists==F){
                out_tib <- 
                    tibble(
                        asm_name = 
                            basename(assembly_path),
                        accession = 
                            genome_acc,
                        ori_pos = 
                            NA,
                        ori_strand = 
                            NA,
                        len = 
                            genome_length, 
                        out_file = NA,
                        note = "no_GFF"
                    )
            } else{
                ## use the dnaA annotation to synchronize
                gene_tib <- 
                    dnaA_from_gff(gff_path, genome_acc)
                            #else if (j==2 & genome_genus=="Burkholderia")
                # repABC_from_gff
                
            if (is.vector(gene_tib)){
                out_tib <- 
                    tibble(
                        asm_name = 
                            basename(assembly_path),
                        accession = 
                            genome_acc,
                        ori_pos = 
                            NA,
                        ori_strand = 
                            NA,
                        len = 
                            genome_length,
                        out_file =
                            NA,
                        note = "no dnaA annotation"
                    )
                } else{
                    ## synchronize the sequence
                    gff_ori_strand <- 
                        pull(gene_tib, strand)
                    ## syncronize to the dnaA position
                    if (gff_ori_strand=="+"){
                        ori_start <- 
                            pull(gene_tib, start)
                        out_seq <- 
                            xscat(
                                subseq(
                                    genome_seq, 
                                    start = ori_start, 
                                    end = genome_length),
                                subseq(
                                    genome_seq, 
                                    start = 1, 
                                    end = (ori_start-1))
                            ) %>%
                            `names<-`(genome_acc)
                    } else {
                        gene_end <- 
                            pull(gene_tib, end) 
                        # somehow in Lacticaseibacillus paracasei, 
                        # the dnaa end was 1 nt
                        # more than the genome length...
                        if (gene_end > width(genome_seq)){
                            gene_end <- width(genome_seq)
                        }
                        ori_start <- 
                            genome_length - gene_end + 1
                        out_seq <- 
                            xscat(
                                subseq(
                                    reverseComplement(genome_seq), 
                                    start = ori_start, 
                                    end = genome_length
                                ),
                                subseq(
                                    reverseComplement(genome_seq), 
                                    start = 1, 
                                    end = (ori_start-1)
                                )
                            ) %>%
                            `names<-`(genome_acc)
                        
                    } 
                    
                    out_file <- 
                        sprintf(
                            "%s/%s.txt",
                            sync_dir, genome_acc
                        )
                    file.exists(out_file)
                    writeXStringSet(
                        out_seq, 
                        filepath = out_file
                    )
                    
                    if (save_gene){
                        gene_seq <- 
                            out_seq %>%
                            subseq(
                                1, 
                                pull(gene_tib, end) - pull(gene_tib, start)
                            )
                        gene_out_dir <- 
                            sprintf(
                                "%s/dnaA_gene",
                                sync_dir
                            )
                        if (!dir.exists(gene_out_dir)){
                            dir.create(gene_out_dir)
                        }
                        out_gene <- 
                            sprintf(
                                "%s/%s.txt",
                                gene_out_dir,
                                genome_acc
                            )
                        if(
                            length(
                                list.files(
                                    gene_out_dir
                                )
                            )<3
                        ){
                            writeXStringSet(
                                gene_seq,
                                sprintf(
                                    "%s/%s_gene.txt",
                                    gene_out_dir,
                                    genome_acc
                                )
                            )
                            
                        }
                        
                        out_tib <- 
                            tibble(
                                asm_name = 
                                    basename(assembly_path),
                                accession = 
                                    genome_acc,
                                ori_pos = 
                                    gene_tib$start[1],
                                ori_strand = 
                                    gene_tib$strand[1],
                                len = 
                                    genome_length,
                                out_file = 
                                    out_file,
                                note = NA
                            )
                        
                    }
                }
            }
        }
        
        out_tib <- 
            bind_rows(out_tib)
        return(out_tib)
        
        
    }


id_dnaA_boxes <- 
    function(
        genome_seq,
        dnaA_box = 
            DNAString("TTATCCACA")
            #DNAString("CCGTTCACA")
            #DNAString("YYMTMCRMA")
            # c(
            #     DNAStringSet("TTGTCCACA"),
            #     DNAStringSet("TTCTCCACA"),
            #     DNAStringSet(reverseComplement(DNAString("TTGTCCACA"))),
            #     DNAStringSet(reverseComplement(DNAString("TTCTCCACA")))
            # )
    ){
        
        forward_matches <-
            matchPattern(dnaA_box, unlist(genome_seq), max.mismatch = 1) %>%
            #matchPattern(dnaA_box, unlist(genome_seq), fixed = F) %>%
            as.data.frame() %>%
            rowwise() %>%
            mutate(
                midpt = round((end+start)/2),
                dbh1 = abs(midpt-1),
                dbh2 = abs(width(genome_seq)-midpt),
                min_dist = min(c(dbh1, dbh2)),
                direction = "forward"
            ) %>% 
            as_tibble()
        
        reverse_matches <- 
            matchPattern(
                dnaA_box, 
                reverseComplement(unlist(genome_seq)), max.mismatch = 1
                ) %>%
            #matchPattern(dnaA_box, reverseComplement(unlist(genome_seq)), fixed = F) %>%
            as.data.frame() %>%
            rowwise() %>%
            mutate(
                midpt = 
                    width(genome_seq) - round((end+start)/2) + 1,
                dbh1 = abs(midpt-1),
                dbh2 = abs(width(genome_seq)-midpt),
                min_dist = min(c(dbh1, dbh2)),
                direction = "rev"
            ) %>% 
            as_tibble()
        
        out_tib <- 
            bind_rows(
                forward_matches, 
                reverse_matches
            ) %>% 
            arrange(midpt) %>%
            mutate(
                prev_dist = midpt - lag(midpt, default = min(midpt)),  
                next_dist = lead(midpt, default = max(midpt)) - midpt,  
                d = prev_dist + next_dist, 
                b = 
                    (1 / d)
            )
        
        # ggplot(out_tib, aes(x = midpt/1e4, y = b)) +
        #     geom_segment(aes(x = midpt/1e4, xend = midpt/1e4, y = 0, yend = b), 
        #                  color = "gray", linewidth = 0.5) +
        #     geom_point(color = "blue") +
        #     labs(x = "Genome Position (bp)", y = "DnaA Box Density (b = 1/d)",
        #          title = "DnaA Box Clustering Across the Genome") +
        #     theme_minimal()
        # 
        # geom_segment(
        #     data = out_tib, 
        #     aes(
        #         x = midpt/1e4, 
        #         xend = midpt/1e4, 
        #         y = 0, yend = b-4e4
        #     ), 
        #     color = "gray", linewidth = 0.5
        # )
        
        # fna_dir <- 
        #     (sync_tib %>%pull(out_file))[1]
        # genome_seq <- 
        #     readDNAStringSet(fna_dir)
        # dnaA_box_tib <- 
        #     lapply(
        #         1:length(dnaA_boxes),
        #         function(k){
        #         #for (k in 1:length(dnaA_boxes))
        #             x <- 
        #                 dnaA_boxes[k][[1]]
        #             dbh <- 
        #                 matchPattern(x, genome_seq[[1]],max.mismatch = 2) %>%
        #                 as.data.frame() %>%
        #                 mutate(
        #                     midpt = round((end+start)/2),
        #                     dbh1 = abs(midpt-1),
        #                     dbh2 = abs(width(genome_seq)-midpt),
        #                     min_dist = min(c(dbh1, dbh2))
        #                     ) #%>%
                    #pull(midpt)
                # dbh1 <- 
                #     abs(dbh-1)
                # dbh2 <- 
                #     abs(width(genome_seq)-1)
                # if (dbh1<dbh2){
                #     out_dist <- 
                #         dbh1
                # } else {
                #     out_dist <- 
                #         dbh2
                # }
        #         return(dbh)
        #     }
        # ) %>%
        #     bind_rows() %>%
        #     mutate(box = as.factor(seq))
       
        return(out_tib)
        
    }

sync_gene2 <- 
    function(
        assembly_path, gff_dir, sync_dir,
        save_gene = T
    ){
        # assembly_path <- 
        #     sync_summary$asm_name
        # the fasta file
        fasta_path <- 
            sprintf(
                "%s/%s_genomic.fna.gz",
                assembly_path, 
                gsub(".*/", "", assembly_path)
            )
        
        asm_name <- 
            basename(assembly_path)
        
        genome_seq <- 
            readDNAStringSet(fasta_path)
        #showConnections(all = T)
        
        genome_seq <- 
            genome_seq[
                !(grepl("plasmid", names(genome_seq),  
                        ignore.case = T))
            ]
        genome_seq <- 
            genome_seq[
                !(grepl("phage", names(genome_seq),  
                        ignore.case = T)) |
                    grepl(
                        "phage resistant", 
                        names(genome_seq),  
                        ignore.case = T
                    )
            ]
        
        ## check chromosomes
        genome_seqs <- 
            genome_seq[order(width(genome_seq), decreasing = T),]
        
        #names(genome_seqs)
        chrom_num <- 
            length(genome_seqs)
        
        out_tibs <- 
            list()
        
        #if (chrom_num>1){
        for (j in 1:chrom_num){
            
            if (chrom_num>1){
                sync_chrom_dir <- 
                    sprintf(
                        "%s/chr_%s",
                        sync_dir, j
                    )
                if (!dir.exists(sync_chrom_dir)){
                    dir.create(sync_chrom_dir)
                }
            } else {
                sync_chrom_dir <- 
                    sync_dir
            }
            gff_chrom_dir <- 
                sprintf(
                    "%s/chr_%s",
                    gff_dir, j
                )
            if (!dir.exists(gff_chrom_dir)){
                dir.create(gff_chrom_dir)
            }
            
            genome_seq <- 
                genome_seqs[j]
            genome_length <- 
                width(genome_seq)
            
            genome_acc <- 
                strsplit(names(genome_seq), " ")[[1]][1]
            
            genome_genus <- 
                strsplit(names(genome_seq), " ")[[1]][2] 
                
            ## download the gff
            gff_path <- 
                assembly_path %>%
                sprintf(
                    "%s/%s_genomic.gff.gz",
                    ., gsub(".*/", "", .)
                )
            ## Set the download path
            gff_outfile <- 
                sprintf(
                    "%s/%s",
                    gff_chrom_dir,
                    basename(gff_path)
                )
            # if it's there, don't downlaod it
            if (file.exists(gff_outfile)){
                gff_tib <- 
                    tibble(
                        gff = basename(gff_outfile),
                        gff_exists = T,
                        asm_name = asm_name
                    )
            } else{ # if not, try downloading it
                test_downlaod <- 
                    try(
                        curl_download(
                            url = 
                                gff_path,
                            destfile = 
                                gff_outfile
                        ), 
                        silent = T
                    )
                gff_tib <- 
                    tibble(
                        gff = basename(gff_outfile),
                        gff_exists = !grepl("Error", test_downlaod),
                        asm_name = asm_name
                    )
            }
            if (gff_tib$gff_exists==F){
                out_tibs[[j]] <- 
                    tibble(
                        asm_name = 
                            basename(assembly_path),
                        accession = 
                            genome_acc,
                        chr = j,
                        ori_pos = 
                            NA,
                        ori_strand = 
                            NA,
                        len = 
                            genome_length, 
                        out_file = NA
                    )
            } else{
            
                ## use the dnaA annotation to synchronize
                if (j==1){
                    gene_tib <- 
                        dnaA_from_gff(gff_path)
                } else if ((j==2) & (genome_genus=="Vibrio")){
                    gene_tib <- 
                        RctB_from_gff(gff_path)
                } else if (j>1){
                    out_tibs[[j]] <- 
                        tibble(
                            asm_name = 
                                basename(assembly_path),
                            accession = 
                                genome_acc,
                            chr = j,
                            ori_pos = 
                                NA,
                            ori_strand = 
                                NA,
                            len = 
                                genome_length
                        )
                    next
                }
                
                #else if (j==2 & genome_genus=="Burkholderia")
                # repABC_from_gff
                
                if (is.vector(gene_tib)){
                    out_tibs[[j]] <- 
                        tibble(
                            asm_name = 
                                basename(assembly_path),
                            accession = 
                                genome_acc,
                            chr = j,
                            ori_pos = 
                                NA,
                            ori_strand = 
                                NA,
                            len = 
                                genome_length,
                            out_file =
                                NA
                        )
                } else{
                    ## synchronize the sequence
                    gff_ori_strand <- 
                        pull(gene_tib, strand)
                    ## syncronize to the dnaA position
                    if (gff_ori_strand=="+"){
                        ori_start <- 
                            pull(gene_tib, start)
                        out_seq <- 
                            xscat(
                                subseq(
                                    genome_seq, 
                                    start = ori_start, 
                                    end = genome_length),
                                subseq(
                                    genome_seq, 
                                    start = 1, 
                                    end = (ori_start-1))
                            ) %>%
                            `names<-`(genome_acc)
                    } else {
                        gene_end <- 
                            pull(gene_tib, end) 
                        # somehow in Lacticaseibacillus paracasei, 
                        # the dnaa end was 1 nt
                        # more than the genome length...
                        if (gene_end > width(genome_seq)){
                            gene_end <- width(genome_seq)
                        }
                        ori_start <- 
                            genome_length - gene_end + 1
                        out_seq <- 
                            xscat(
                                subseq(
                                    reverseComplement(genome_seq), 
                                    start = ori_start, 
                                    end = genome_length
                                ),
                                subseq(
                                    reverseComplement(genome_seq), 
                                    start = 1, 
                                    end = (ori_start-1)
                                )
                            ) %>%
                            `names<-`(genome_acc)
                        
                    } 
                    
                    out_file <- 
                        sprintf(
                            "%s/%s.txt",
                            sync_chrom_dir, genome_acc
                        )
                    file.exists(out_file)
                    writeXStringSet(
                        out_seq, 
                        filepath = out_file
                    )
                    
                    if (save_gene){
                        gene_seq <- 
                            out_seq %>%
                            subseq(
                                1, 
                                pull(gene_tib, end) - pull(gene_tib, start)
                                )
                        gene_out_dir <- 
                            sprintf(
                                "%s/chr_%s_gene",
                                sync_dir, j
                            )
                        if (!dir.exists(gene_out_dir)){
                            dir.create(gene_out_dir)
                        }
                        out_gene <- 
                            sprintf(
                                "%s/%s.txt",
                                gene_out_dir,
                                genome_acc
                            )
                        if(
                            length(
                                list.files(
                                    gene_out_dir
                                )
                                )<3
                        ){
                            writeXStringSet(
                                gene_seq,
                                sprintf(
                                    "%s/%s_gene.txt",
                                    gene_out_dir,
                                    genome_acc
                                )
                            )
                            
                        }
                    }
                    
                    
                    out_tibs[[j]] <- 
                        tibble(
                            asm_name = 
                                basename(assembly_path),
                            accession = 
                                genome_acc,
                            chr = j,
                            ori_pos = 
                                gene_tib$start[1],
                            ori_strand = 
                                gene_tib$strand[1],
                            len = 
                                genome_length,
                            out_file = 
                                out_file
                        )
                    
                }
            }
            }
            
        
        out_tib <- 
            bind_rows(out_tibs)
        return(out_tib)
        
        
    }


# 3 functions
download_gff <- 
    function(species_table, gff_dir, core_number=7){
        gff_table <- 
            pbmclapply(
                1:nrow(species_table),
                function(i){
                    ## get the gff path
                    asm_name <- 
                        pull(species_table, asm_name)[i]
                    gff_path <- 
                        (species_table %>% 
                             pull(ftp_path))[i] %>%
                        sprintf(
                            "%s/%s_genomic.gff.gz",
                            ., gsub(".*/", "", .)
                        )
                    ## Set the download path
                    gff_outfile <- 
                        sprintf(
                            "%s/%s",
                            gff_dir,
                            basename(gff_path)
                        )
                    # if it's there, don't downlaod it
                    if (file.exists(gff_outfile)){
                        gff_tib <- 
                            tibble(
                                gff = basename(gff_outfile),
                                gff_exists = T,
                                asm_name = asm_name
                            )
                    } else{ # if not, try downloading it
                        test_downlaod <- 
                            try(
                                curl_download(
                                    url = 
                                        gff_path,
                                    destfile = 
                                        gff_outfile
                                ), 
                                silent = T
                            )
                        gff_tib <- 
                            tibble(
                                gff = basename(gff_outfile),
                                gff_exists = !grepl("Error", test_downlaod),
                                asm_name = asm_name
                            )
                    }
                    return(gff_tib)
                }, mc.cores = core_number
            ) %>% 
            bind_rows()
        return(gff_table)
    }

sync_with_gff <- 
    function(
        species_table, sync_dir, 
        gff_dir, core_number=7#,
    ){
        #        species_table = species_summary
        ## check for dnaAs
        sync_report <- 
            pbmclapply(
                1:nrow(species_table),
                function(i){
                    #for (i in 1:nrow(species_table)){
                    gff_path <- 
                        sprintf(
                            "%s/%s_genomic.gff.gz",
                            gff_dir,
                            basename((species_table %>% pull(ftp_path))[i])
                        )
                    
                    ## Set the download path
                    gff_outfile <- 
                        sprintf(
                            "%s/%s",
                            gff_dir,
                            basename(gff_path)
                        )
                    #read.gff(gff_outfile)
                    genome_acc <- 
                        (species_table %>% pull(accession))[i]
                    out_file <- 
                        sprintf(
                            "%s/%s.txt",
                            sync_dir, genome_acc
                        )
                    
                    ## Check if it's been synced
                    if (file.exists(out_file)) {
                        out_string <- 
                            sprintf(
                                "%s already synced",
                                genome_acc
                            )
                        dnaA_vect <- NULL
                    } else if (!file.exists(gff_outfile)){
                        out_string <- 
                            sprintf(
                                "%s has no GFF file",
                                genome_acc
                            )
                    } else{
                        #for (i in 1:length(fna_paths)){
                        gff <- 
                            try(read.gff(gff_outfile), silent = T)
                        
                        if (any(grepl("Error", gff[[1]]))){
                            out_string <- 
                                sprintf(
                                    "%s has a missing GFF file",
                                    genome_acc
                                )
                        } else {
                            asm_path <- 
                                (species_table %>% pull(ftp_path))[i]
                            genome_file <- 
                                sprintf(
                                    "%s/%s_genomic.fna.gz",
                                    asm_path, 
                                    gsub(".*/", "", asm_path)
                                )
                            genome_seq <- 
                                readDNAStringSet(genome_file)
                            genome_seq <- 
                                genome_seq[
                                    !(grepl("plasmid", names(genome_seq),  
                                            ignore.case = T))
                                    ]
                            genome_seq <- 
                                genome_seq[
                                    !(grepl("phage", names(genome_seq),  
                                            ignore.case = T)) |
                                        grepl(
                                            "phage resistant", 
                                            names(genome_seq),  
                                            ignore.case = T
                                        )
                                    ]
                            # take the bigger sequence ie the main one
                            genome_seq <- 
                                genome_seq[
                                    width(genome_seq)==max(width(genome_seq))
                                    ]
                            genome_length <- 
                                width(genome_seq)
                            
                            ## check if dnaa is annotated
                            ## to cut space
                            product_string <- 
                                "product="
                            dnaAstring <- 
                                "DnaA"
                            replication_string <- 
                                " replication initiator "
                            dnaA_strings <- 
                                c(
                                    sprintf(
                                        "%schromosome%s%s",
                                        product_string, 
                                        replication_string, 
                                        dnaAstring
                                    ),
                                    sprintf(
                                        "%schromosome%sprotein %s",
                                        product_string,
                                        replication_string,
                                        dnaAstring
                                    )
                                )
                            
                            dnaA_vect <- 
                                gff %>% 
                                as_tibble() %>% 
                                filter(
                                    (grepl(
                                        "gene=dnaA", attributes, ignore.case = T
                                    ) |
                                        grepl(
                                            dnaA_strings[1], 
                                            attributes, ignore.case = T
                                        ) |
                                        grepl(
                                            dnaA_strings[2], 
                                            attributes, ignore.case = T
                                        ) |
                                        grepl(
                                            "product=DnaA", 
                                            attributes, ignore.case = T
                                        )
                                    ) &
                                        type == "CDS"
                                ) %>% as.data.frame()
                            
                            if (nrow(dnaA_vect)>1){
                                if (
                                    any(
                                        grepl(
                                            "gene=DnaA", 
                                            dnaA_vect$attributes, 
                                            ignore.case = T
                                        )
                                    )
                                ){
                                    dnaA_vect <- 
                                        dnaA_vect[
                                            grepl(
                                                "gene=DnaA", 
                                                dnaA_vect$attributes, 
                                                ignore.case = T
                                            ),
                                            ] %>% dplyr::slice(1)
                                } else if (all(pull(dnaA_vect, strand)=="+")){
                                    dnaA_vect <- 
                                        dnaA_vect %>%
                                        filter(start ==min(start)) %>% 
                                        dplyr::slice(1)
                                } else if (all(pull(dnaA_vect, strand)=="-")) {
                                    dnaA_vect <- 
                                        dnaA_vect %>%
                                        filter(start ==max(start))
                                } else {
                                    dnaA_vect <- 
                                        dnaA_vect %>% dplyr::slice(1)
                                }
                            }
                            
                            if (nrow(dnaA_vect)==0) {
                                out_string <- 
                                    sprintf(
                                        "%s has no dnaA annotation",
                                        genome_acc
                                    )
                            } else {
                                gff_ori_strand <- 
                                    pull(dnaA_vect, strand)
                                ## syncronize to the dnaA position
                                if (gff_ori_strand=="+"){
                                    ori_start <- 
                                        pull(dnaA_vect, start)
                                    out_seq <- 
                                        xscat(
                                            subseq(
                                                genome_seq, 
                                                start = ori_start, 
                                                end = genome_length),
                                            subseq(
                                                genome_seq, 
                                                start = 1, 
                                                end = (ori_start-1))
                                        ) %>%
                                        `names<-`(genome_acc)
                                } else {
                                    dnaA_end <- 
                                        pull(dnaA_vect, end) 
                                    # somehow in Lacticaseibacillus paracasei, 
                                    # the dnaa end was 1 nt
                                    # more than the genome length...
                                    if (dnaA_end > width(genome_seq)){
                                        dnaA_end <- width(genome_seq)
                                    }
                                    ori_start <- 
                                        genome_length - dnaA_end + 1
                                    out_seq <- 
                                        xscat(
                                            subseq(
                                                reverseComplement(genome_seq), 
                                                start = ori_start, 
                                                end = genome_length
                                            ),
                                            subseq(
                                                reverseComplement(genome_seq), 
                                                start = 1, 
                                                end = (ori_start-1)
                                            )
                                        ) %>%
                                        `names<-`(genome_acc)
                                    
                                } 
                                
                                writeXStringSet(
                                    out_seq, 
                                    filepath = out_file
                                )
                                out_string <- 
                                    sprintf(
                                        "%s synced",
                                        genome_acc
                                    )
                            }
                        }
                    }
                    
                    return(out_string)
                }, mc.cores = core_number
            ) %>% unlist
        
        ### second pass
        second_pass_index <- 
            which(
                !(
                    species_table$accession %in% 
                        (list.files(sync_dir) %>% sub(".txt", "", .))
                )
            )
        for (i in second_pass_index){
            ## get the gff path
            gff_path <- 
                sprintf(
                    "%s/%s_genomic.gff.gz",
                    gff_dir,
                    basename((species_table %>% pull(ftp_path))[i])
                )
            
            ## Set the download path
            gff_outfile <- 
                sprintf(
                    "%s/%s",
                    gff_dir,
                    basename(gff_path)
                )
            genome_acc <- 
                (species_table %>% pull(accession))[i]
            out_file <- 
                sprintf(
                    "%s/%s.txt",
                    sync_dir, genome_acc
                )
            
            ## Check if it's been synced
            if (file.exists(out_file)) {
                out_string <- 
                    sprintf(
                        "%s already synced",
                        genome_acc
                    )
            } else if (!file.exists(gff_outfile)){
                out_string <- 
                    sprintf(
                        "%s has no GFF file",
                        genome_acc
                    )
            } else{
                #for (i in 1:length(fna_paths)){
                
                gff <- 
                    try(read.gff(gff_outfile), silent = T)
                if (length(gff)==1){
                    out_string <- 
                        sprintf(
                            "%s has no GFF file",
                            genome_acc
                        )
                    next
                }
                
                asm_path <- 
                    (species_table %>% pull(ftp_path))[i]
                genome_file <- 
                    sprintf(
                        "%s/%s_genomic.fna.gz",
                        asm_path, 
                        gsub(".*/", "", asm_path)
                    )
                genome_seq <- 
                    readDNAStringSet(genome_file)
                genome_seq <- 
                    genome_seq[
                        !(grepl("plasmid", names(genome_seq),  ignore.case = T))
                        ]
                genome_seq <- 
                    genome_seq[
                        !(grepl("phage", names(genome_seq),  ignore.case = T))
                        ]
                # take the bigger sequence ie the main one
                genome_seq <- 
                    genome_seq[
                        width(genome_seq)==max(width(genome_seq))
                        ]
                genome_length <- 
                    width(genome_seq)
                ## check if dnaa is annotated
                ## check if dnaa is annotated
                ## to cut space
                product_string <- 
                    "product="
                dnaAstring <- 
                    "DnaA"
                replication_string <- 
                    " replication initiator "
                dnaA_strings <- 
                    c(
                        sprintf(
                            "%schromosome%s%s",
                            product_string, 
                            replication_string, 
                            dnaAstring
                        ),
                        sprintf(
                            "%schromosome%sprotein %s",
                            product_string,
                            replication_string,
                            dnaAstring
                        )
                    )
                
                dnaA_vect <- 
                    gff %>% 
                    as_tibble() %>% 
                    filter(
                        (grepl(
                            "gene=dnaA", attributes, ignore.case = T
                        ) |
                            grepl(
                                dnaA_strings[1], 
                                attributes, ignore.case = T
                            ) |
                            grepl(
                                dnaA_strings[2], 
                                attributes, ignore.case = T
                            ) |
                            grepl(
                                "product=DnaA", 
                                attributes, ignore.case = T
                            )
                        ) &
                            type == "CDS"
                    ) %>% as.data.frame()
                
                
                if (nrow(dnaA_vect)>1){
                    if (all(pull(dnaA_vect, strand)=="+")){
                        dnaA_vect <- 
                            dnaA_vect %>%
                            filter(start ==min(start))
                    } else if (all(pull(dnaA_vect, strand)=="-")) {
                        dnaA_vect <- 
                            dnaA_vect %>%
                            filter(start ==max(start))
                    } else {
                        dnaA_vect <- 
                            dnaA_vect %>% dplyr::slice(1)
                    }
                }
                
                if (nrow(dnaA_vect)==0) {
                    out_string <- 
                        sprintf(
                            "%s has no dnaA annotation",
                            genome_acc
                        )
                } else {
                    gff_ori_strand <- 
                        pull(dnaA_vect, strand)
                    #syncronize to the dnaA position
                    if (gff_ori_strand=="+"){
                        ori_start <- 
                            pull(dnaA_vect, start)
                        out_seq <- 
                            xscat(
                                subseq(
                                    genome_seq, 
                                    start = ori_start, end = genome_length),
                                subseq(
                                    genome_seq, 
                                    start = 1, end = (ori_start-1))
                            ) %>%
                            `names<-`(genome_acc)
                    } else {
                        ori_start <- 
                            genome_length - pull(dnaA_vect, end) + 1
                        out_seq <- 
                            xscat(
                                subseq(
                                    reverseComplement(genome_seq), 
                                    start = ori_start, end = genome_length
                                ),
                                subseq(
                                    reverseComplement(genome_seq), 
                                    start = 1, end = (ori_start-1)
                                )
                            ) %>%
                            `names<-`(genome_acc)
                        
                    }
                    writeXStringSet(
                        out_seq, out_file
                    )
                    out_string_2 <- 
                        sprintf(
                            "%s synced",
                            genome_acc
                        )
                }
            }
        }
        
        return(sync_report)
    }

get_connections <- 
    function(uri){
    
    showConnections(all = TRUE) %>%
        as.data.frame %>%
        rownames_to_column('con_id') %>%
        filter(grepl(gff_path, description)) %>%
        pull(con_id) %>%
        as.integer %>%
        map(getConnection)
    
}


close_connection <- 
    function(url_path){
    
        close(
            getConnection(
                as.integer(
                    names(
                        which(
                            showConnections(all = TRUE)[, 1]==
                                sprintf(
                                    "gzcon(%s)", gff_path
                                    )
                            )
                        )
                    )
                )
            )
        
        
}

dnaA_from_gff <- 
    function(gff_path, genome_acc){
        #
        if (!file_exists_at_url(gff_path)){ 
            dnaA_vect <- 
                sprintf(
                    "%s has no GFF annotation",
                    genome_acc
                )
        } else {
                
            gff_gzcon <- 
                gzcon(url(gff_path))
            
            gff_table <- 
                read.gff(
                    gff_gzcon
                )
            gff_gzcon %>% as.character()
            
            # connection_index <- 
            #     showConnections(all = TRUE)[, 1]
            close(gff_gzcon)
            #close_connection(gff_path)
            ##
            product_string <- 
                "product="
            dnaAstring <- 
                "DnaA"
            replication_string <- 
                " replication initiator "
            dnaA_strings <- 
                c(
                    sprintf(
                        "%schromosome%s%s",
                        product_string, 
                        replication_string, 
                        dnaAstring
                    ),
                    sprintf(
                        "%schromosome%sprotein %s",
                        product_string,
                        replication_string,
                        dnaAstring
                    )
                )
            
            dnaA_vect <- 
                gff_table %>% 
                as_tibble() %>% 
                filter(
                    (grepl(
                        "gene=dnaA", attributes, ignore.case = T
                    ) |
                        grepl(
                            dnaA_strings[1], 
                            attributes, ignore.case = T
                        ) |
                        grepl(
                            dnaA_strings[2], 
                            attributes, ignore.case = T
                        ) |
                        grepl(
                            "product=DnaA", 
                            attributes, ignore.case = T
                        )
                    ) &
                        type == "CDS"
                ) %>% as.data.frame()
            
            if (nrow(dnaA_vect)>1){
                if (
                    any(
                        grepl(
                            "gene=DnaA", 
                            dnaA_vect$attributes, 
                            ignore.case = T
                        )
                    )
                ){
                    dnaA_vect <- 
                        dnaA_vect[
                            grepl(
                                "gene=DnaA", 
                                dnaA_vect$attributes, 
                                ignore.case = T
                            ),
                            ] %>% dplyr::slice(1)
                } else if (all(pull(dnaA_vect, strand)=="+")){
                    dnaA_vect <- 
                        dnaA_vect %>%
                        filter(start ==min(start)) %>% 
                        dplyr::slice(1)
                } else if (all(pull(dnaA_vect, strand)=="-")) {
                    dnaA_vect <- 
                        dnaA_vect %>%
                        filter(start ==max(start))
                } else {
                    dnaA_vect <- 
                        dnaA_vect %>% dplyr::slice(1)
                }
            }
            
            if (nrow(dnaA_vect)==0) {
                dnaA_vect <- 
                    sprintf(
                        "%s has no dnaA annotation",
                        genome_acc
                    )
            }
            }
        
        return(dnaA_vect)
    }

## RctB 
RctB_from_gff <- 
    function(gff_path){
        #
        gff_table <- 
            read.gff(
                gzcon(url(gff_path))
            )
        ##
        product_string <- 
            "product="
        RctBstring <- 
            "RctB"
        replication_string <- 
            "replication initiator protein"
        RctB_strings <- 
            c(
                sprintf(
                    "%schromosome%s%s",
                    product_string, 
                    replication_string, 
                    RctBstring
                ),
                sprintf(
                    "%s%s %s",
                    product_string,
                    replication_string,
                    RctBstring
                )
            )
        # gff_table %>% 
        #     as_tibble() %>% 
        #     filter(grepl("replication initiator", attributes)) %>%
        #     pull(attributes)
        
        RctB_vect <- 
            gff_table %>% 
            as_tibble() %>% 
            filter(
                (grepl(
                    "gene=RctB", attributes, ignore.case = T
                ) |
                    grepl(
                        RctB_strings[1], 
                        attributes, ignore.case = T
                    ) |
                    grepl(
                        RctB_strings[2], 
                        attributes, ignore.case = T
                    ) |
                    grepl(
                        "product=RctB", 
                        attributes, ignore.case = T
                    )
                ) &
                    type == "CDS"
            ) %>% as.data.frame()
        
        if (nrow(RctB_vect)>1){
            if (
                any(
                    grepl(
                        "gene=RctB", 
                        RctB_vect$attributes, 
                        ignore.case = T
                    )
                )
            ){
                RctB_vect <- 
                    RctB_vect[
                        grepl(
                            "gene=RctB", 
                            RctB_vect$attributes, 
                            ignore.case = T
                        ),
                    ] %>% dplyr::slice(1)
            } else if (all(pull(RctB_vect, strand)=="+")){
                RctB_vect <- 
                    RctB_vect %>%
                    filter(start ==min(start)) %>% 
                    dplyr::slice(1)
            } else if (all(pull(RctB_vect, strand)=="-")) {
                RctB_vect <- 
                    RctB_vect %>%
                    filter(start ==max(start))
            } else {
                RctB_vect <- 
                    RctB_vect %>% dplyr::slice(1)
            }
        }
        
        if (nrow(RctB_vect)==0) {
            RctB_vect <- 
                sprintf(
                    "%s has no RctB annotation",
                    genome_acc
                )
        }
        return(RctB_vect)
    }

RepA_from_gff <- 
    function(gff_path){
        #
        gff_table <- 
            read.gff(
                gzcon(url(gff_path))
            )
        ##
        product_string <- 
            "product="
        repAstring <- 
            "repA"
        replication_string <- 
            "replication initiator protein"
        repA_strings <- 
            c(
                sprintf(
                    "%schromosome%s%s",
                    product_string, 
                    replication_string, 
                    repAstring
                ),
                sprintf(
                    "%s%s %s",
                    product_string,
                    replication_string,
                    repAstring
                )
            )
        gff_table %>%
            as_tibble() %>%
            filter(grepl("replication initiator", attributes)) %>%
            pull(attributes)
        
        repA_vect <- 
            gff_table %>% 
            as_tibble() %>% 
            filter(
                (grepl(
                    "gene=repA", attributes, ignore.case = T
                ) |
                    grepl(
                        repA_strings[1], 
                        attributes, ignore.case = T
                    ) |
                    grepl(
                        repA_strings[2], 
                        attributes, ignore.case = T
                    ) |
                    grepl(
                        "product=repA", 
                        attributes, ignore.case = T
                    )
                ) &
                    type == "CDS"
            ) %>% as.data.frame()
        
        if (nrow(RctB_vect)>1){
            if (
                any(
                    grepl(
                        "gene=RctB", 
                        RctB_vect$attributes, 
                        ignore.case = T
                    )
                )
            ){
                RctB_vect <- 
                    RctB_vect[
                        grepl(
                            "gene=RctB", 
                            RctB_vect$attributes, 
                            ignore.case = T
                        ),
                    ] %>% dplyr::slice(1)
            } else if (all(pull(RctB_vect, strand)=="+")){
                RctB_vect <- 
                    RctB_vect %>%
                    filter(start ==min(start)) %>% 
                    dplyr::slice(1)
            } else if (all(pull(RctB_vect, strand)=="-")) {
                RctB_vect <- 
                    RctB_vect %>%
                    filter(start ==max(start))
            } else {
                RctB_vect <- 
                    RctB_vect %>% dplyr::slice(1)
            }
        }
        
        if (nrow(RctB_vect)==0) {
            RctB_vect <- 
                sprintf(
                    "%s has no RctB annotation",
                    genome_acc
                )
        }
        return(RctB_vect)
    }



sync2 <- 
    function(
        fasta_file, gff_file, out_file
    ){
        ### reading the fasta file
        genome_seq <- 
            readDNAStringSet(fasta_file)
        genome_seq <- 
            genome_seq[
                !(grepl("plasmid", names(genome_seq),  
                        ignore.case = T))
                ]
        genome_seq <- 
            genome_seq[
                !(grepl("phage", names(genome_seq),  
                        ignore.case = T)) |
                    grepl(
                        "phage resistant", 
                        names(genome_seq),  
                        ignore.case = T
                    )
                ]
        # take the bigger sequence ie the main one
        genome_seq <- 
            genome_seq[
                width(genome_seq)==max(width(genome_seq))
                ]
        genome_length <- 
            width(genome_seq)
        ## getting the dnaA coordinates
        dnaA_tib <- 
            dnaA_from_gff(gff_path)
        
        gff_ori_strand <- 
            pull(dnaA_tib, strand)
        ## syncronize to the dnaA position
        if (gff_ori_strand=="+"){
            ori_start <- 
                pull(dnaA_tib, start)
            out_seq <- 
                xscat(
                    subseq(
                        genome_seq, 
                        start = ori_start, 
                        end = genome_length),
                    subseq(
                        genome_seq, 
                        start = 1, 
                        end = (ori_start-1))
                ) %>%
                `names<-`(genome_acc)
        } else {
            dnaA_end <- 
                pull(dnaA_tib, end) 
            # somehow in Lacticaseibacillus paracasei, 
            # the dnaa end was 1 nt
            # more than the genome length...
            if (dnaA_end > width(genome_seq)){
                dnaA_end <- width(genome_seq)
            }
            ori_start <- 
                genome_length - dnaA_end + 1
            out_seq <- 
                xscat(
                    subseq(
                        reverseComplement(genome_seq), 
                        start = ori_start, 
                        end = genome_length
                    ),
                    subseq(
                        reverseComplement(genome_seq), 
                        start = 1, 
                        end = (ori_start-1)
                    )
                ) %>%
                `names<-`(genome_acc)
            
        } 
        
        writeXStringSet(
            out_seq, 
            filepath = out_file
        )
        
        
    }


GC_check <- 
    function(fasta_file, disparity_file){
        nuc_disparities <- 
            sprintf(
                "python3 nucleotide_disparities.py %s %s",
                fasta_file,
                disparity_file
            )
        
        system(nuc_disparities)
        skew_tib <- 
            fread(
                disparity_file
            )
        acc_i <- 
            basename(
                sub(".txt", "", fasta_file)
            )
        #plot(skew_tib$`G-C`)
        ## check that the minimum skew is position 1
        if (
            which(skew_tib$`G-C`==min(skew_tib$`G-C`)) == 1
        ){
            GC_info_tib <- 
                tibble(
                    accession = acc_i,
                    GC_ori_loc = 1
                )
        } else {
            minima_coords <- 
                c(
                    1e4*(which(skew_tib$`G-C`==min(skew_tib$`G-C`))-1),
                    1e4*(which(skew_tib$`G-C`==min(skew_tib$`G-C`))+1)
                )
            disp_regions_dir <- 
                sprintf(
                    "%s/regions",
                    disparities_dir
                )
            original_ff <- 
                readDNAStringSet(fasta_file)
            if(!dir.exists(disp_regions_dir)){dir.create(disp_regions_dir)}
            ori_region_file <- 
                sprintf(
                    "%s/%s_ori.txt", 
                    disp_regions_dir, acc_i
                )
            if (minima_coords[2]>width(original_ff)){
                minima_coords[2] <- 
                    width(original_ff)
            }
            ori_region <- 
                original_ff %>%
                subseq(., minima_coords[1], minima_coords[2]) %>%
                `names<-`(c(paste0(acc_i, "_ori_region")))
            
            writeXStringSet(
                ori_region,
                ori_region_file
            )
            ori_disparity_file <- 
                sprintf(
                    "%s/%s_ori_region.csv",
                    disp_regions_dir,
                    acc_i
                )
            ## get the precise nucleotide disparities
            ori_nuc_disparities <- 
                sprintf(
                    "python3 nucleotide_disparities.py %s %s 1",
                    ori_region_file,
                    ori_disparity_file
                )
            
            system(ori_nuc_disparities)
            
            GC_min_region_idx <- 
                which(
                    (fread(ori_disparity_file) %>% pull(`G-C`))==
                        min(fread(ori_disparity_file) %>% pull(`G-C`))
                )[1]
            
            ori_location <- 
                (minima_coords[1] + GC_min_region_idx)
            
            # new_sync_seq <- 
            #     xscat(
            #         subseq(
            #             original_ff,
            #             ori_location, 
            #             width(original_ff)
            #         ),
            #         subseq(
            #             original_ff,
            #             1, ori_location-1
            #         )
            #     ) %>%
            #     `names<-`(c(acc_i))
            # # overwrite the old sequence
            # writeXStringSet(
            #     new_sync_seq,
            #     fasta_file 
            # )    
            
            GC_info_tib <- 
                tibble(
                    accession = acc_i,
                    GC_ori_loc = 
                        ori_location,
                    ori_dist = 
                        min(
                            c(
                                abs(1-ori_location),
                                abs(ori_location-width(original_ff))
                            )
                        ),
                    ori_dist_frac = 
                        round(100*(ori_dist/width(original_ff)), 4),
                    genome_len = width(original_ff)
                )
            
            ## recreate the GC skew calculation
            # nuc_disparities <- 
            #     sprintf(
            #         "python3 nucleotide_disparities.py %s %s",
            #         fasta_file,
            #         disparity_file
            #     )
            # 
            # system(nuc_disparities)
            # skew_tib <- 
            #     fread(
            #         disparity_file
            #     )
        }
        
        return(GC_info_tib)
    }


GC_skew_sync <- 
    function(fasta_file, disparity_file){
        nuc_disparities <- 
            sprintf(
                "python3 nucleotide_disparities.py %s %s",
                fasta_file,
                disparity_file
            )
        
        system(nuc_disparities)
        skew_tib <- 
            fread(
                disparity_file
            )
        acc_i <- 
            basename(
                sub(".txt", "", fasta_file)
            )
        #plot(skew_tib$`G-C`)
        ## check that the minimum skew is position 1
        if (
            which(skew_tib$`G-C`==min(skew_tib$`G-C`)) == 1
        ){
            GC_info_tib <- 
                tibble(
                    accession = acc_i,
                    GC_ori_loc = 1
                )
        } else {
            minima_coords <- 
                c(
                    1e4*(which(skew_tib$`G-C`==min(skew_tib$`G-C`))-1),
                    1e4*(which(skew_tib$`G-C`==min(skew_tib$`G-C`))+1)
                )
            disp_regions_dir <- 
                sprintf(
                    "%s/regions",
                    disparities_dir
                )
            original_ff <- 
                readDNAStringSet(fasta_file)
            if(!dir.exists(disp_regions_dir)){dir.create(disp_regions_dir)}
            ori_region_file <- 
                sprintf(
                    "%s/%s_ori.txt", 
                    disp_regions_dir, acc_i
                )
            if (minima_coords[2]>width(original_ff)){
                minima_coords[2] <- 
                    width(original_ff)
            }
            ori_region <- 
                original_ff %>%
                subseq(., minima_coords[1], minima_coords[2]) %>%
                `names<-`(c(paste0(acc_i, "_ori_region")))
            
            writeXStringSet(
                ori_region,
                ori_region_file
            )
            ori_disparity_file <- 
                sprintf(
                    "%s/%s_ori_region.csv",
                    disp_regions_dir,
                    acc_i
                )
            ## get the precise nucleotide disparities
            ori_nuc_disparities <- 
                sprintf(
                    "python3 nucleotide_disparities.py %s %s 1",
                    ori_region_file,
                    ori_disparity_file
                )
            
            system(ori_nuc_disparities)
            
            GC_min_region_idx <- 
                which(
                    (fread(ori_disparity_file) %>% pull(`G-C`))==
                        min(fread(ori_disparity_file) %>% pull(`G-C`))
                )[1]
            
            ori_location <- 
                (minima_coords[1] + GC_min_region_idx)
            
            new_sync_seq <- 
                xscat(
                    subseq(
                        original_ff,
                        ori_location, 
                        width(original_ff)
                    ),
                    subseq(
                        original_ff,
                        1, ori_location-1
                    )
                ) %>%
                `names<-`(c(acc_i))
            # overwrite the old sequence
            writeXStringSet(
                new_sync_seq,
                fasta_file 
            )    
            
            GC_info_tib <- 
                tibble(
                    accession = acc_i,
                    GC_ori_loc = 
                        ori_location
                )
            
            ## recreate the GC skew calculation
            # nuc_disparities <- 
            #     sprintf(
            #         "python3 nucleotide_disparities.py %s %s",
            #         fasta_file,
            #         disparity_file
            #     )
            # 
            # system(nuc_disparities)
            # skew_tib <- 
            #     fread(
            #         disparity_file
            #     )
        }
            
        return(GC_info_tib)
    }

synchronize <- 
    function(species_table, sync_dir, gff_dir){
        gff_table <- 
            download_gff(species_table, gff_dir)
        sync_report <-
            sync_with_gff(
                species_table,
                sync_dir, 
                gff_dir
            ) 
        
        return(sync_report)
    }

rs_sync <- 
    function(
        ref_seq,
        fasta_seq,
        delta_out,
        out_file
    ){
        
        system(
            simple_nucmer(
                ref_seq, fasta_seq, output = delta_out
                ), 
            intern = F
            )
        
        sst <- 
            read_delta(
                sprintf(
                    "%s_filtered.delta",
                    delta_out
                )
            ) 
        if (is.null(sst)){
            return(NULL)
        } else {
            
            sst <- 
                sst %>% 
                filter(rs==min(rs))
            genome_seq <- readDNAStringSet(fasta_seq)
            if (length(genome_seq)>1){
                stop(">1 sequences in fasta seq")
            }
            genome_length <- 
                width(genome_seq)
            genome_acc <- 
                names(genome_seq)
            if (sst$strand=="+"){
                ori_start <- 
                    pull(sst, qs)
                out_seq <- 
                    xscat(
                        subseq(
                            genome_seq, 
                            start = ori_start, 
                            end = genome_length),
                        subseq(
                            genome_seq, 
                            start = 1, 
                            end = (ori_start-1))
                    ) %>%
                    `names<-`(genome_acc)
            } else {
                dnaA_end <- 
                    pull(sst, qe)
                # somehow in Lacticaseibacillus paracasei, 
                # the dnaa end was 1 nt
                # more than the genome length...
                if (dnaA_end > width(genome_seq)){
                    dnaA_end <- width(genome_seq)
                }
                ori_start <- 
                    genome_length - dnaA_end + 1
                out_seq <- 
                    xscat(
                        subseq(
                            reverseComplement(genome_seq), 
                            start = ori_start, 
                            end = genome_length
                        ),
                        subseq(
                            reverseComplement(genome_seq), 
                            start = 1, 
                            end = (ori_start-1)
                        )
                    ) %>%
                    `names<-`(genome_acc)
                
            } 
            file.remove(
                sprintf("%s_filtered.delta", delta_out)
            )
            
            writeXStringSet(
                out_seq, 
                filepath = out_file
            )
            out_tib <- 
                tibble(
                    accession = 
                        names(genome_seq),
                    ori_pos = 
                        ori_start[1],
                    ori_strand = 
                        sst$strand[1],
                    len = 
                        genome_length
                )
            return(out_tib)
            
        }
            
    }


#### Pairwise identities ####
# 4 functions
# distance_matrix <- 
#     id_sld_wide_dist
get_dists <- 
    function(distance_matrix){
        diag(distance_matrix) <- NA
        sequence_averages <- rowMeans(distance_matrix, na.rm = TRUE)
        all_distances <- as.vector(distance_matrix)
        all_distances <- 
            lapply(all_distances, function(x){x[!is.na(x)]})
            #!all_distances[!is.na(all_distances)]  # Remove NA values
        overall_mean_distance <-
            mean(unlist(all_distances))
            #lapply(all_distances, mean)
            #mean(all_distances)
        difference_from_mean <- sequence_averages - overall_mean_distance
        
        distance_table <- 
            tibble(
                Sequence = names(sequence_averages), 
                Sequence_mean_distance = sequence_averages, 
                Mean_distance = overall_mean_distance,
                Difference_from_Mean = difference_from_mean
            )
        return(distance_table)
    }


mash_commands <- 
    function(genus, org, identity_dir, sync_dir){
        sketchfile <- 
            sprintf(
                "mash sketch -p 8 -o %s/mash.msh -s 10000 '%s'/*.txt ",
                identity_dir, sync_dir
            )
        mash_lines <- 
            sprintf(
                "mash dist -p 8 %s/mash.msh %s/mash.msh > %s/%s_%s_mash.txt",
                identity_dir, identity_dir, 
                identity_dir, genus, org
            )
        return(capture.output(cat(sketchfile, mash_lines, sep = "\n")))
    }

identity_sld <- 
    function(
        identity_dir, genus_name, species_name,
        species_summary, filetype = "mash",
        clust_method = "complete", sync_dir
    ){
        
        id_table <-
            fread(
                sprintf(
                    "%s/%s_%s_%s.txt", 
                    identity_dir, genus_name, species_name, filetype
                ),# stringsAsFactors = F, 
                select = c(1:3)
            ) %>%
            mutate(V1 = sub(".txt", "", basename(V1))) %>%
            mutate(V2 = sub(".txt", "", basename(V2))) %>%
            `colnames<-`(c("ref", "qry", "identity")) %>%
            mutate(identity = as.numeric(identity)) %>%
            as_tibble
        
        if (filetype=="mash"){
            id_table <- 
                id_table %>%
                mutate(
                    identity = (1-identity)*100
                )
        }
        
        id_table_wide <- 
            id_table %>%
            dplyr::select(c(ref, qry, identity)) %>%
            spread(., qry, identity) %>%
            as.data.frame() %>%
            column_to_rownames("ref")
        
        ## for clustering purposes
        id_dist <- 
            100-as.matrix(id_table_wide) 
        id_dist[upper.tri(id_dist, diag = T)]<-NA
        
        id_dist <- 
            as.dist(id_dist)
        id_clust <- 
            hclust(id_dist, method = clust_method)
        clust_order <-
            id_clust$labels[id_clust$order]
        
        id_table_long <- 
            id_table_wide %>%
            rownames_to_column("ref") %>%
            gather(key = qry,value = "identity", na.rm = FALSE, -c(ref)) %>%
            group_by(ref) %>% 
            mutate(identity_mean = mean(identity)) %>%
            arrange(desc(identity_mean)) %>%
            ungroup
        
        seq_lens <- 
            species_summary %>% 
            dplyr::select(c(accession, len))
        
        # if (any(is.na(pull(seq_lens, len)))){
        #     
        # }
        if (
            any(
                !(
                    (pull(id_table_long, ref) %>% unique) %in% 
                    pull(seq_lens, accession)
                )
            )
        ){
            accs <- 
                (pull(id_table_long, ref) %>% unique)[
                    !(
                        (pull(id_table_long, ref) %>% unique) %in% 
                            pull(seq_lens, accession)
                    )
                    ]
            accs_tib <- 
                lapply(accs, function(x){
                    #for (x in accs){
                    tibble(
                        accession = x,
                        len = 
                            width(
                                readDNAStringSet(
                                    sprintf("%s/%s.txt", sync_dir, x)
                                )
                            )
                    )
                }) %>% bind_rows
            seq_lens <- 
                seq_lens %>%
                bind_rows(accs_tib)
        }
        
        
        
        id_sld <- 
            id_table_long[,c(1:3)] %>%
            left_join(
                ., rename(seq_lens, accession = "ref"), 
                by = "ref") %>%
            rename(len = "reflen") %>%
            left_join(
                ., rename(seq_lens, accession = "qry"), 
                by = "qry") %>%
            rename(len = "qrylen") %>%
            mutate(
                ld = abs(reflen-qrylen),
                sld = 100*(abs(reflen-qrylen)/mean(c(reflen, qrylen))),
                identity_sld = identity-sld
            ) %>%
            mutate(
                ref = factor(ref, levels = clust_order),
                qry = factor(qry, levels = clust_order)
            )
        return(id_sld)
    }

get_pairwise_identities <- 
    function(
        sync_dir, identity_dir, species_table,
        genus_name, species_name
    ){
        # run mash if it hasn't been run
        if (!
            file.exists(
                sprintf(
                    "%s/%s_%s_mash.txt", 
                    identity_dir, genus_name, species_name
                )
            )
        ){
            mash_bash <- 
                mash_commands(genus_name, species_name, identity_dir, sync_dir)
            system(mash_bash[1], ignore.stderr = T)
            system(mash_bash[2])
        }
        mash_tib <-
            identity_sld(
                identity_dir, genus_name, species_name, species_table, 
                filetype = "mash", 
                sync_dir = sync_dir
            )
        
        #Recreate the wide matrices
        id_wide <- 
            mash_tib %>%
            dplyr::select(c(ref, qry, identity)) %>%
            spread(., qry, identity) %>%
            as.data.frame() %>%
            column_to_rownames("ref") %>%
            as.matrix()
        
        id_wide_dist <- 
            100 - id_wide
        
        id_sld_wide <- 
            mash_tib %>%
            dplyr::select(c(ref, qry, identity_sld)) %>%
            spread(., qry, identity_sld) %>%
            as.data.frame() %>%
            column_to_rownames("ref")
        
        id_sld_wide_dist <- 
            100 - id_sld_wide
        
        
        distance_matrix <- 
            id_wide_dist
        
        id_dists <- 
            get_dists(distance_matrix) %>%
            arrange(Sequence_mean_distance)
        
        id_sld_dists <-
            get_dists(id_sld_wide_dist) %>%
            arrange(Sequence_mean_distance)
        
        id_dist_lt <- 
            as.dist(distance_matrix)
        
        seq_order <- 
            colnames(distance_matrix)[
                hclust( id_dist_lt, method = "ward.D" )$order
            ]
            
        #
        mash_plot <- 
            ggplot(
                data = mash_tib, 
                aes(
                    x= factor(ref, seq_order),
                    y= factor(qry, seq_order),
                    fill = identity#_sld#identity
                )
            ) + 
            geom_tile() + 
            theme_classic() +
            theme(
                #axis.text.x = element_text(angle=45, hjust = 1),
                panel.background=element_rect(fill="white"),
                plot.title = element_text(size=20,face="bold", hjust = 0.5),
                axis.title = element_blank(),
                axis.text = element_blank()
            ) + 
            ggtitle(
                sprintf(
                    "%s \n Sequence pairwise MASH similarity", 
                    sub("_", " ", species_name)
                )
            ) + 
            theme(
                legend.title=element_text( size=18),
                axis.text=element_blank(),
                    #element_text(size=14), 
                legend.text = element_text(size=18)
            ) + 
            scale_fill_viridis(
                option = "B", limits = c(95,100)
            ) 
        
        mash_list <- 
            list(
                mash_tib = mash_tib,
                mash_plot = mash_plot,
                id_dists = id_dists,
                id_sld_dists = id_sld_dists
            )
        return(mash_list)
    }
#### Clustering
## 8 functions
# Alignment, 2 functions
generate_nucmer_commands <- 
    function(genome_matrix, delta_dir, sync_dir, nuc_remove = T){
        nuc_commands <- 
            apply(
                genome_matrix,
                1,
                function(x){
                    
                    nuc_params <-
                        "--mum --maxgap=500 --mincluster=100"
                    nucmer <-
                        sprintf(
                            "nucmer %s --prefix=%s/%s_v_%s",
                            nuc_params,
                            delta_dir,
                            x[1], x[2]
                        )
                    nucmer_command <- 
                        sprintf(
                            "then %s %s/%s.txt %s/%s.txt",
                            nucmer,
                            sync_dir, x[1], 
                            sync_dir, x[2]
                        )
                    filter_out <- 
                        sprintf(
                            "%s/%s_v_%s_filtered.delta",
                            delta_dir,
                            x[1], 
                            x[2]
                        )
                    
                    nuc_check <- 
                        sprintf(
                            "if [ ! -e %s ]",
                            filter_out
                        )
                    
                    nucmer_filter <-
                        sprintf(
                            "delta-filter -r -q %s/%s_v_%s.delta > %s",
                            delta_dir,
                            x[1], x[2],
                            filter_out
                        )
                    ## to save space, remove the original delta file
                    if (nuc_remove){
                        nuc_remove <- 
                            sprintf(
                                "rm %s/%s_v_%s.delta",
                                delta_dir,
                                x[1], x[2]
                            )
                    } else {
                        unfiltered_dir <- 
                            sprintf(
                                "%s/unfiltered",
                                delta_dir
                            )
                        if (!dir.exists(unfiltered_dir)){
                            dir.create(unfiltered_dir)
                            }
                        
                        nuc_remove <- 
                            sprintf(
                                "mv %s/%s_v_%s.delta %s/%s_v_%s.delta",
                                delta_dir,
                                x[1], x[2],
                                unfiltered_dir,
                                x[1], x[2]
                            )
                    }
                    
                    
                    return(
                        capture.output(
                            cat(
                                nuc_check,
                                nucmer_command, 
                                nucmer_filter, 
                                nuc_remove,
                                "fi",
                                sep = "; "
                            )
                        )
                    )
                    
                }
            )
        return(nuc_commands)
    }

## for 1 vs all
nucmer_ref_v_all <- 
    function(
        most_related, accession_list, sync_dir, delta_dir
    ){
        alignment_matrix <- 
            tibble(
                ref = most_related,
                qry = accession_list[accession_list != most_related]
            ) %>% 
            as.matrix()
        # nuc_params <-
        #     "--mum --maxgap=500 --mincluster=100"
        delta_dir_rva <- 
            sprintf(
                "%s/ref_v_all",
                delta_dir
            )
        if (
            !dir.exists(delta_dir_rva)
        ){
            dir.create(delta_dir_rva, recursive = T)
        }
        nuc_commands <- 
            generate_nucmer_commands(
                genome_matrix = alignment_matrix, 
                delta_dir =  delta_dir_rva, 
                sync_dir = sync_dir 
            )
        pbmclapply(nuc_commands, system, mc.cores = 7)
    }


## for all vs all 
nucmer_all_v_all <- 
    function(
        accession_list, sync_dir, delta_dir
    ){
        alignment_matrix <- 
            t(combn(accession_list, 2))
            
        nuc_params <-
            "--mum --maxgap=500 --mincluster=100"
        delta_dir_ava <- 
            sprintf(
                "%s/all_v_all",
                delta_dir
            )
        if (
            !dir.exists(delta_dir_ava)
        ){
            dir.create(delta_dir_ava, recursive = T)
        }
        nuc_commands <- 
            generate_nucmer_commands(
                genome_matrix = alignment_matrix, 
                delta_dir =  delta_dir_ava, 
                sync_dir = sync_dir 
            )
        pbmclapply(nuc_commands, system, mc.cores = 7)
    }


## show-snps

show_snps <- 
    function(
        delta_file, snp_file
    ){
        ss_command <- 
            sprintf(
                "show-snps -CTHlr %s > %s",
                delta_file, snp_file
            )
        system(ss_command)
        snp_table <- 
            fread(snp_file) %>% 
            `colnames<-`(
                c(
                    "ref_pos",
                    "ref_snp", "qry_snp",
                    "qry_pos", "BUFF", "DIST",
                    "reflen", "qrylen", "direction", "TAG", "ref", "qry"
                )
            ) %>% 
            dplyr::select(
                c(
                    "ref_pos", "qry_pos",
                    "ref_snp", "qry_snp",
                    "reflen", "qrylen",  
                    "ref", "qry"
                )
            )
        return(snp_table)
    }

show_coords <- 
    function(
        delta_file, coord_file
    ){
        sc_command <- 
            sprintf(
                "show-coords -dbcrTlH %s > %s",
                #"show-coords -B %s > %s",
                delta_file, coord_file
            )
        system(sc_command)
        coord_table <- 
            fread(coord_file) %>% 
            `colnames<-`(
                c(
                    "R1", "R2", "Q1", "Q2",
                    "LEN_1","LEN_2","LEN_R", 
                    "LEN_Q", "COV_R", "COV_Q", 
                    "REF", "QRY"
                )
            ) %>%
            rowwise %>%
            mutate(
                strand = 
                    ifelse(
                        Q1<Q2, 
                        "+", "-"
                    )
            )
        return(coord_table)
    }

dna_diff <- 
    function(
        delta_file, dnadiff_file
    ){
        dnadiff_file
        
        dd_command <- 
            sprintf(
                "dnadiff -d %s -p %s",
                delta_file,
                dnadiff_file
            )
        system(dd_command)
        readLines(
            sprintf(
                "%s.report",
                dnadiff_file 
            )
        )
        
    }

# dnadiff_file <- 
#     "../data/processing/test/dnadiff/tes1"

        
        
        
# 
# coord_tib <-
#     show_coords(
#         delta_file,
#         coord_file
#     ) %>%
#     `colnames<-`(
#         c(
#             "rs", "re", "qs", "qe",
#             "rcov", "qcov",  "rlen", "qlen",
#             "rcov_pc", "qcov_pc", "rid", "qid", "strand"
#         )
#     ) %>%
#     as_tibble
# 
# 
# 
# delta_table %>%
#     filter(!(rs %in% coord_tib$rs))

# 
# plot_delta(coord_tib)
# delta_table %>% colnames(
# )
# Reading delta files, 3 functions
read_delta <- 
    function(delta_path){
        delta_lines <- 
            (readLines(delta_path) %>% strsplit(., " "))[-c(1,2)]
        if (length(delta_lines )==0) {
            return(NULL)
        } else {
            
            id_lines <- 
                delta_lines[lengths(delta_lines)==4][[1]]
            delta_alignments <- 
                delta_lines[lengths(delta_lines)==7] %>%
                unlist %>%
                matrix(., ncol = 7, byrow = T) 
            
            delta_alignments <- 
                delta_alignments %>%
                apply(., 2, as.numeric) %>%
                matrix(., ncol = 7) %>%
                `colnames<-`(c("rs", "re", "qs", "qe", "error", "e2", "zero")) %>% 
                as_tibble %>%
                dplyr::select(c(1:5)) %>%
                mutate(
                    strand = ifelse(qe-qs > 0, '+', '-'),
                    rid = strsplit(id_lines, " ")[[1]] %>% sub(">", "", .),
                    qid = strsplit(id_lines, " ")[[2]],
                    rlen = strsplit(id_lines, " ")[[3]] %>% as.numeric,
                    qlen = strsplit(id_lines, " ")[[4]] %>% as.numeric,
                    rcov = abs(re-rs+1),
                    qcov = abs(qe-qs+1),
                    perc_error = (round(100*error/ pmax(rcov, qcov), 2))
                ) %>% 
                rowwise() %>%
                mutate(
                    meanlen = 
                        ceiling(mean(c(rcov, qcov))),
                    refmid = 
                        (rs+re)/2,
                    qrymid = 
                        (qs+qe)/2,
                    expected_offset_ref = 
                        (refmid * ((qlen / rlen) - 1)) / sqrt(2),
                    expected_offset_qry = 
                        (qrymid * ((qlen / rlen) - 1)) / sqrt(2),
                    mean_offset = 
                        mean(
                            c(
                                expected_offset_ref, 
                                expected_offset_qry
                            )
                        ),
                    fwd_dist = 
                        abs(
                            abs(refmid - qrymid) / sqrt(2) - mean_offset
                            ),
                    rev_dist = 
                        abs(
                            ((refmid + qrymid) - mean(c(rlen, qlen))) / 
                                sqrt(2)
                        ) - mean_offset,
                    X_dist = 
                        round(
                            max(c(min(c(fwd_dist, rev_dist)), 0))
                        ),
                    # X_dist = 
                    #     ifelse(
                    #         strand == "+",
                    #         abs(refmid - qrymid) / sqrt(2),
                    #         abs(
                    #             (refmid + qrymid) - mean(c(rlen, qlen))) / 
                    #             sqrt(2)
                    #         ) - mean_offset,
                    # X_dist = 
                    #     ifelse(
                    #         X_dist<0,
                    #         0, X_dist
                    #     ),
                    X_dist_weight = 
                        ceiling(
                            round(X_dist * (meanlen/mean(c(rlen, qlen))), 2)
                            )
                ) %>%
                ungroup %>%
                arrange(rs) %>%
                dplyr::select(-c(
                    refmid, qrymid, expected_offset_ref,
                    expected_offset_qry, mean_offset
                    ))
        
            return(delta_alignments)
        }
    }
###
filter_delta <- 
    function(
        delta_table, contig_summary = FALSE, 
        maxgap = 1e4, minlen = 1e4, X_dist_diff = 5e4
    ){
        if (is.character(delta_table)){
            delta_table <- read_delta(delta_table)
        }
        
        out_tibble <- 
            delta_table %>%
            group_by(strand) %>% 
            mutate(
                qry_gapsize = 
                    qs - lag(qe, default = qs[1]),
                ref_gapsize = 
                    rs - lag(re, default = rs[1]),
                X_diff = 
                    X_dist-lag(X_dist, default = 0)
            ) %>%
            ungroup() %>%
            mutate(
                qry_gapsize = 
                    ifelse(
                        strand=="+",
                        qry_gapsize,
                        qry_gapsize*-1
                    ),
                qry_gaps_up = 
                    ifelse(
                        (qry_gapsize) < maxgap,
                        0,1
                    ),
                qry_gaps_up = 
                    cumsum(qry_gaps_up),
                qry_gaps_down = 
                    ifelse(
                        (qry_gapsize) > -maxgap,
                        0,1
                    ),
                qry_gaps_down = 
                    cumsum(qry_gaps_down),
                ref_gaps = 
                    ifelse(
                        abs(ref_gapsize) < maxgap,
                        0,1
                    ),
                ref_gaps = 
                    cumsum(ref_gaps),
                XDD = ## new
                    X_dist-lag(X_dist, default = 0),
                XDF = 
                    ifelse(
                        abs(XDD) < X_dist_diff,
                        0,1
                    ),
                XDF_diff = 
                    cumsum(XDF),
                new_contigs = ## new
                    rleid(
                        strand, qry_gaps_up, qry_gaps_down, ref_gaps, XDF_diff
                    )
            ) %>%
            group_by(new_contigs, strand) %>% 
            dplyr::summarise(
                "X_dist" = mean(X_dist),
                "meanlen" = sum(meanlen),
                "rs" = min(rs),
                "re" = max(re),
                "qs" = 
                    unique(ifelse(
                        strand=="+",
                        min(qs),
                        max(qs)
                    )),
                "qe" = 
                    unique(ifelse(
                        strand=="+",
                        max(qe),
                        min(qe)
                    )),
                "rid" = unique(rid),
                "qid" = unique(qid),
                "slope" = (qe-qs)/(re-rs),
                "rlen" = unique(rlen),
                "qlen" = unique(qlen),
                .groups = "keep"
            ) %>%
            ungroup %>%
            filter(
                #Each contig is greater than the minlen
                (meanlen > minlen) |
                    # or is the first
                    (new_contigs==1) |
                    # or is the last
                    (new_contigs == max(new_contigs))
            ) %>%
            mutate(
                qry_gapsize = 
                    qs - lag(qe, default = qs[1]),
                ref_gapsize = 
                    rs - lag(re, default = rs[1]),
                XDD = 
                    X_dist-lag(X_dist, default = 0),
                XDF = 
                    ifelse(
                        abs(XDD) < X_dist_diff,
                        0,1
                    ),
                XDF_diff = 
                    cumsum(XDF),
                new_contigs = 
                    rleid(strand, XDF_diff)
            ) %>%
            ungroup() %>%
            group_by(new_contigs, strand) %>% 
            dplyr::summarise(
                #"X_dist" =  mean(X_dist),
                "rs" = min(rs),
                "re" = max(re),
                "qs" = 
                    unique(ifelse(
                        strand=="+",
                        min(qs),
                        max(qs)
                    )),
                "qe" = 
                    unique(ifelse(
                        strand=="+",
                        max(qe),
                        min(qe)
                    )),
                "rid" = unique(rid),
                "qid" = unique(qid),
                "slope" = (qe-qs)/(re-rs),
                "rlen" = unique(rlen),
                "qlen" = unique(qlen), 
                "meanlen" = mean(c(abs(qe-qs), abs(re-rs))),
                .groups = "keep"
            ) %>% 
            ungroup %>%
            rowwise() %>%
            mutate(
                X_dist = 
                    ifelse(
                        strand=="+",
                        mean(
                            c(
                                abs((rs-qs)/sqrt(2)), # start position
                                abs((re-qe)/sqrt(2))
                            )
                        ),
                        mean(
                            c(
                                abs((qs + rs - qlen)/sqrt(2)), # start position
                                abs((qe + re - qlen)/sqrt(2))
                            )
                        )
                    )
            ) %>% 
            ungroup
        
        if (contig_summary){
            out_tibble <- 
                out_tibble %>%
                mutate(
                    c_contig = rleid(strand)
                ) %>%
                group_by(c_contig, strand) %>%
                dplyr::summarise(
                    "rs" = min(rs),
                    "re" = max(re),
                    "qs" = 
                        unique(ifelse(
                            strand=="+",
                            min(qs),
                            max(qs)
                        )),
                    "qe" = 
                        unique(ifelse(
                            strand=="+",
                            max(qe),
                            min(qe)
                        )
                        ),
                    "rlen" = unique(rlen),
                    "qlen" = unique(qlen),
                    "slope" = mean(slope),
                    "X_dist" = mean(X_dist),
                    "rid" = unique(rid),
                    "qid" = unique(qid),
                    .groups = "keep"
                ) %>%
                ungroup()
        }
        
        return(out_tibble)
    }

filter_combined_delta <- 
    function(
        dd, 
        minimum_alignment_percentage, minimum_inversion_length
    ){

        #dd =sprintf("%s/ref_v_all",  delta_dir)
        combined_delta <- 
            pbmclapply(
                1:length(list.files(dd, pattern = ".delta")),
                function(i){
                    dfo <- 
                        read_delta(
                            list.files(dd, full.names = T)[i]
                        ) 
                    if (!is.null(dfo)){
                        dfo <- 
                            filter_delta(dfo)
                    }
                    return(dfo)
                }, mc.cores = 6
            ) %>% bind_rows
        
        # for (ii in 1:length(list.files(dd, pattern = ".delta"))){
        #     read_delta(
        #         list.files(dd, full.names = T)[ii]
        #     ) %>% filter_delta() 
        # }
        
        
        # read_delta(list.files(dd, pattern = ".delta", full.names = T)[1])
        # read_delta(list.files(dd, pattern = ".delta", full.names = T)[2])
        
        low_coverage <- 
            combined_delta %>%
            group_by(rid, qid) %>% 
            summarise(
                cov_length = sum(meanlen),
                cov_prop = 100*(cov_length/mean(c(unique(rlen), unique(qlen)))),
                rlen = unique(rlen),
                qlen = unique(qlen), 
                .groups = "keep"
            ) %>% 
            arrange(cov_prop)
        
        low_cov_accs <- 
            low_coverage %>% 
            filter(cov_prop<minimum_alignment_percentage) %>% 
            pull(qid)
        
        ## remove from the table
        filtered_combined_delta <- 
            combined_delta %>%
            group_by(rid, qid) %>% 
            filter(
                ## reomve the small inversions
                (meanlen > minimum_inversion_length) | 
                    (rs ==1) |
                    (re==max(re)),
                ## remove the low coverage queries
                !(qid %in% low_cov_accs)
            ) %>%
            mutate(
                segment = rleid(strand),
                segments = length(rle(strand)$lengths),
                invs = sum(rle(strand)$values=="-")
            ) %>% 
            mutate(
                ref_al_len =re-rs+1
            ) %>%
            dplyr::select(
                rid, qid,
                strand, X_dist, rs, re, qs, qe, qid, rlen, qlen, 
                meanlen, segment, segments, invs
            )
        
        
        return(filtered_combined_delta)
    }

plot_delta <- 
    function(
        delta_table, gtitle = "", 
        xlb = NULL, ylb = NULL
    ){
        
        if (is.character(delta_table)){
            delta_table <- read_delta(delta_table)
        }
        myColors <- brewer.pal(5,"Set1")[c(1,2)]
        names(myColors) <- c("-", "+")
        colScale <- scale_colour_manual(name = "grp",values = myColors)
        
        p <- 
            ggplot(
                delta_table %>% 
                    mutate(
                        rid = sub("NZ_", "", rid),
                        qid = sub("NZ_", "", qid)
                    ), 
                aes(
                    x=rs, xend=re, y=qs, yend=qe, 
                    colour= factor(strand, levels = c("-", "+"))
                ) 
            )+
            geom_segment(alpha=1, linewidth = 1.5) + 
            theme_classic() + 
            theme(
                legend.position= "none", 
                legend.justification=c(1,0),
                #axis.text=element_text(size=14),
                plot.title = element_text(hjust = 0.35),
                axis.title = element_text(size=12)
            ) +
            geom_abline(
                intercept = 
                    mean(c(delta_table$rlen[1], delta_table$qlen[1])), 
                slope = -1,
                color="darkgrey",
                linetype="dashed", linewidth=.5) +
            geom_abline(
                intercept = 0, slope = 1,
                color="darkgrey",
                linetype="dashed", linewidth=.5) +
            colScale +
            scale_x_continuous(
                limits=c(1, delta_table$rlen[1]), expand = c(0,0)
            ) +
            scale_y_continuous(
                limits=c(1, delta_table$qlen[1]), expand = c(0,0)
            )
        
        if (is.null(xlb) & is.null(ylb)){
            p <- 
                p +             
                xlab(sub("NZ_", "", unique(pull(delta_table, rid)))) + 
                ylab(sub("NZ_", "", unique(pull(delta_table, qid)))) 
        } else{
            p <-
                p +
                xlab(xlb) + 
                ylab(ylb) 
        }
        return(p)
    }
# clustering, 3 functions
merge_overlaps <- 
    function(
        irange_in, min_overlap = 0.95, edge_proximity = 1e4,
        silent = T
    ){
        ## returns a data frame with grouped rows
        
        groups_i <- list()
        overlap_rows <- list()
        
        if (!silent){
            pb <- txtProgressBar(min = 1, max = length(irange_in), initial = 0)
        }
        for (i in 1:length(irange_in)) {
            #     # Skip interivals that are already assigned to a group
            if (i %in% unlist(groups_i)){
                groups_i[[i]] <- NA
                next
            } 
            
            # Calculate the overlap length threshold for the current interval
            
            overlap_length_threshold <- min_overlap * width(irange_in[i])
            
            # Search for overlapping intervals 
            subj_overlaps_1 <- 
                subjectHits(
                    IRanges::findOverlaps(
                        irange_in[i], irange_in, 
                        minoverlap = overlap_length_threshold
                    )
                )
            
            subj_overlaps_1 <- subj_overlaps_1[subj_overlaps_1!=i]
            
            if (length(subj_overlaps_1)==0){
                groups_i[[i]] <- 
                    NA
                #print(c(i))
                #next
            } else if (length(subj_overlaps_1)>0){
                hit_vector <- 
                    c()
                for (j in subj_overlaps_1){
                    overlap_length_threshold_j <- 0.95 * width(irange_in[j])
                    
                    subj_overlaps_2 <- 
                        subjectHits(
                            IRanges::findOverlaps(
                                irange_in[i], irange_in[j], 
                                minoverlap = overlap_length_threshold_j
                            )
                        )
                    starts_distance <- 
                        abs(
                            start(irange_in[i]) - 
                                start(irange_in[j])
                        )
                    ends_distance <- 
                        abs(
                            end(irange_in[i]) - 
                                end(irange_in[j])
                        )
                    
                    breakpoint_proximity <- 
                        (starts_distance < edge_proximity) &
                        (ends_distance < edge_proximity) 
                    
                    
                    if (
                        length(subj_overlaps_2)==0 | 
                        i>=j | 
                        (!breakpoint_proximity)
                    ){
                        hit_vector <- 
                            append(hit_vector, NA)
                        
                    } else {
                        hit_vector <- 
                            append(hit_vector, j)
                    }
                    groups_i[[i]] <- 
                        hit_vector
                    groups_i[[i]] <- 
                        unique(groups_i[[i]])
                }
            } 
            
            groups_i[[i]] <- 
                groups_i[[i]][!is.na(groups_i[[i]])]
            if (length(groups_i[[i]]) >0){
                overlap_rows[[i]] <- 
                    c(i, groups_i[[i]])
            } else{
                overlap_rows[[i]] <- 
                    NA
            }
            
            if(!silent){setTxtProgressBar(pb,i)}
        }
        
        or4 <-
            tibble(
                idx = unlist(overlap_rows), 
                cluster_i = rep(
                    seq(length(overlap_rows)), lengths(overlap_rows)
                    )
            ) 
        
        out_table <- 
            irange_in %>% 
            as.data.frame() %>%
            as_tibble %>% 
            mutate(original_idx = row_number()) %>%
            mutate(
                idx = 1:n()
                #overlap_row_group = overlap_row_group
            ) %>% 
            left_join(
                ., or4, by = "idx"
            ) %>% 
            mutate(
                unordered_overlap_cluster = 
                    ifelse(
                        is.na(cluster_i), idx, cluster_i
                    )
            ) %>% 
            #dplyr::add_count(c(cluster_i, idx), name = "oc_count") 
            dplyr::add_count(
                unordered_overlap_cluster, name = "cluster_count"
                ) %>% 
            arrange(desc(cluster_count), desc(unordered_overlap_cluster)) %>% 
            mutate(
                overlap_cluster = 
                    rleid(
                        unordered_overlap_cluster
                    )
            ) %>%
            dplyr::select(
                c(
                    original_idx, start, end, width, 
                    overlap_cluster, cluster_count
                    )
            )
        
        
        return(out_table)
        
    }




cluster_segments <- 
    function(inner_segment_table){
        
        inner_segment_table_summary <- 
            inner_segment_table %>%
            group_by(qid) %>%
            arrange(rs) %>% 
            summarise(
                strand = unique(strand),
                X_dist_mean_prop = 
                    round(100*(mean(X_dist)/(mean(c(rlen, qlen)))), 2),
                X_dist_check = 
                    ifelse(
                        is.nan(round(mean(diff(X_dist)))),
                        0, round(mean(diff(X_dist)))
                    ),
                rs = min(rs),
                re = max(re),
                inv_segments = n(),
                invs = unique(invs),
            )
        segment_ranges <- 
            IRanges(
                start = inner_segment_table_summary$rs,
                end = inner_segment_table_summary$re
            )
        
        clustered_segments <- 
            merge_overlaps(segment_ranges)
        
        clustered_segment_summary <-
            inner_segment_table_summary %>% 
            mutate(original_idx = row_number()) %>%
            left_join(., clustered_segments, by = "original_idx") %>%
            group_by(overlap_cluster) %>%
            summarise(
                accessions = paste(qid, collapse = ", "),
                cluster_count = unique(cluster_count),
                start = mean(start),
                end = mean(end),
                width = (end-start)+1,
                strand = unique(strand),
                mean_X_dist_perc = mean(X_dist_mean_prop)
            )
        return(clustered_segment_summary)
        
    }


delta_diversity <- 
    function(delta_file){
        delta_table <- 
            read_delta(delta_file) 
        fd <- 
            filter_delta(delta_table)
        plot_delta(fd)
        
        xdist_summary <- 
            summary(pull(delta_table, X_dist))
        
        d2 <- 
            delta_table %>%
            rowwise() %>%
            mutate(
                qs2 = min(c(qs, qe)),
                qe2 = max(c(qs, qe))
            ) %>% 
            ungroup
        
        ref_ranges <- 
            IRanges(
                start = 
                    delta_table$rs,
                end = 
                    delta_table$re
            ) %>% 
            reduce() 
        
        qry_ranges <- 
            IRanges(
                start = 
                    d2$qs2,
                end = 
                    d2$qe2
            ) %>% 
            reduce()
        
        ref_gaps <- 
            ref_ranges %>%
            gaps(start = 1, end = delta_table$rlen[1]) 
        qry_gaps <- 
            qry_ranges %>%
            gaps(start = 1, end = delta_table$qlen[1]) 
        
        
        
    }

nucmer_diversity_metrics <- 
    function(
        delta_file,
        diff_threshold = 5e4
    ){
        delta_table <- 
            read_delta(delta_file) 
        fd <- 
            filter_delta(delta_table)
        
        delta_identity <- 
            (sum(delta_table$meanlen)-sum(delta_table$error))/
            mean(c(delta_table$rlen, delta_table$qlen))
        ## Gaps > 5e3
        ref_gaps <- 
            IRanges(
                start = 
                    delta_table$rs,
                end = 
                    delta_table$re
            ) %>% 
            reduce() %>%
            gaps(start = 1, end = delta_table$rlen[1]) 
        d2 <- 
            delta_table %>%
            rowwise() %>%
            mutate(
                qs2 = min(c(qs, qe)),
                qe2 = max(c(qs, qe))
            )
        
        qry_gaps <- 
            IRanges(
                start = 
                    d2$qs2,
                end = 
                    d2$qe2
            ) %>% 
            reduce() %>%
            gaps(start = 1, end = delta_table$qlen[1])
        ref_gap_sum <- 
            sum(width(ref_gaps))
        qry_gap_sum <- 
            sum(width(qry_gaps))
        large_ref_gaps <- 
            sum(ref_gaps>diff_threshold)
        large_qry_gaps <- 
            sum(qry_gaps>diff_threshold)
        ref_qry_seg <- 
            c(
                IRanges(
                    start = 
                        delta_table$rs,
                    end = 
                        delta_table$re
                ),
                IRanges(
                    start = 
                        d2$qs2,
                    end = 
                        d2$qe2
                )
            ) %>%
            reduce()
        
        # segs <- 
        #     ref_qry_seg %>% reduce() %>%as.data.frame() %>% as_tibble() %>%
        #         mutate(
        #             mp = (end+start)/2,
        #             mp2 = 
        #                 mean(c(delta_table$rlen[1], delta_table$qlen[1]))-mp
        #             )
        # pairwise_sums <- 
        #     outer(segs$mp, segs$mp, "+")
        # pairwise_sums >= (mean(c(delta_table$rlen[1], delta_table$qlen[1]))-5e4) & 
        #     pairwise_sums <= (mean(c(delta_table$rlen[1], delta_table$qlen[1]))+5e4) 
        
        
        ref_qry_gaps <- 
            ref_qry_seg %>% 
            gaps(
                start = 1, end = max(c(delta_table$rlen[1],delta_table$qlen[1]))
            )%>% reduce() %>%as.data.frame() %>% as_tibble() %>%
            mutate(
                mp = (end+start)/2,
                mp2 = 
                    mean(c(delta_table$rlen[1], delta_table$qlen[1]))-mp
            ) 
        
        # ref_qry_gaps %>%
        #     as.data.frame()
        # 
        # ref_gaps %>% as.data.frame() %>% arrange(desc(width)) %>% head
        # qry_gaps %>% as.data.frame() %>% arrange(desc(width))
        
        large_gaps <- 
            ref_qry_seg %>% 
            reduce() %>%
            gaps(
                start = 1, end = max(c(delta_table$rlen[1],delta_table$qlen[1]))
            ) %>%
            as.data.frame() %>%
            as_tibble %>%
            filter(width>diff_threshold) 
        
        large_gaps_count <- 
            large_gaps %>%
            nrow()
        large_gaps_nt <- 
            large_gaps %>%
            pull(width) %>%
            sum()
        
        
        ## Translocations
        large_translocations <-
            delta_table %>%
            filter(
                X_dist >diff_threshold, 
                meanlen>diff_threshold,
                strand=="+"
            ) 
        
        large_translocations_count <- 
            large_translocations %>% 
            nrow()
        large_translocations_nt <- 
            large_translocations %>%
            pull(meanlen) %>%
            sum()
        
        ## Length difference
        length_difference <-
            abs(delta_table$rlen[1] - delta_table$qlen[1])
        
        ##inversions
        large_ns_inversions <- 
            delta_table %>% 
            filter(
                strand=="-",
                meanlen>diff_threshold,
                X_dist>diff_threshold 
            ) 
        
        large_ns_inversions_count <- 
            large_ns_inversions %>% 
            nrow()
        large_ns_inversions_nt <- 
            large_ns_inversions %>%
            pull(meanlen) %>%
            sum()
        
        large_symmetric_inversions <- 
            delta_table %>% 
            filter(
                strand=="-",
                meanlen>diff_threshold,
                X_dist<diff_threshold
            ) 
        
        large_symmetric_inversions_count <- 
            large_symmetric_inversions %>% 
            nrow()
        large_symmetric_inversions_nt <- 
            large_symmetric_inversions %>%
            pull(meanlen) %>%
            sum()
        
        delta_summary <- 
            tibble(
                ref = 
                    delta_table$qid[1],
                qry = 
                    delta_table$rid[1],
                delta_identity = 
                    delta_identity,
                large_gap_nts = 
                    large_gaps_nt,
                large_translocations_nt = 
                    large_translocations_nt,
                length_difference = 
                    length_difference,
                large_ns_inversions_nt = 
                    large_ns_inversions_nt,
                large_symmetric_inversions_nt = 
                    large_symmetric_inversions_nt
            )
        
        return(delta_summary)
    }


cluster_genomes <- 
    function(
        filtered_combined_delta,
        minimum_inversion_length = 5e4,
        minimum_alignment_coverage = 90
    ){
        
        ## the first cluster is colinear (group 0)
        cluster_0 <- 
            filtered_combined_delta %>%
            group_by(qid) %>% 
            mutate(segments = length(rle(strand)$lengths)) %>% 
            filter(max(segments)==1) %>%
            summarise(al_cov = sum(meanlen), rlen = unique(rlen)) %>%
            rowwise %>%
            mutate(cov = al_cov/rlen) %>% 
            ungroup %>% 
            arrange(cov)
        
        cluster_others <-
            filtered_combined_delta %>%
            group_by(qid) %>% 
            mutate(
                segment = rleid(strand),
                segments = length(rle(strand)$lengths),
                invs = sum(rle(strand)$values=="-")
            ) %>% 
            mutate(
                ref_al_len =re-rs+1
            ) %>%
            filter(invs >0) %>% 
            dplyr::select(
                strand, X_dist, rs, re, qs, qe, 
                qid, rlen, qlen, meanlen, segment, 
                segments, invs
            )
        
        ## if there are no inversions
        if (nrow(cluster_others)==0){
            colinear_table <- 
                tibble(
                    cluster_id = 0,
                    accessions = 
                        c(
                            filtered_combined_delta$rid[1],
                            cluster_0$qid
                        ),
                    grouped_sequences = length(accessions),
                    inversions = 0
                )
            
            clusters_and_coords <-
                list(
                    cluster_summary = 
                        colinear_table,
                    cluster_coords = 
                        NULL
                )
            
        } else {
            
            #cluster_out_table <- list()
            cluster_multiple_inversions <- 
                filtered_combined_delta %>%
                filter(invs>0) %>% 
                ungroup 
            
            multiple_inversion_types <- 
                cluster_multiple_inversions %>%
                count(invs) %>% 
                pull(invs)
            # filtered_combined_delta %>% filter(qid=="CP049377.1") %>% plot_delta
            # cluster_multiple_inversions %>% filter(invs==5)
            
            inner_segments <- 
                lapply(
                    multiple_inversion_types,
                    function(j){
                        segment_list <- 
                            seq(from = 2, to = j*2, by = 1)
                        return(segment_list)
                    }
                )
            
            clustered_inversion_table <- 
                list()
            clustered_ids <- 
                list()
            inversion_set <- 
                list()
            for (i in multiple_inversion_types){
                for (j in inner_segments[[
                    which(i==multiple_inversion_types)
                    ]]){
                    inv_seg <- 
                        cluster_multiple_inversions %>%
                        filter(
                            (invs==i) &
                                (segment==j)
                        )
                    inversion_set[[
                        #which(inner_segments[[i]]==j)
                        #which(i==multiple_inversion_types)
                        j
                        ]] <- 
                        cluster_segments(inv_seg) %>%
                        mutate(
                            inversion_group = i,
                            inner_segment = j
                        )
                }
                ## this table has the coordinates of each inversion
                clustered_inversion_table[[i]] <- 
                    bind_rows(inversion_set)
                
                ## this shows which group each accession is in based on 
                ## sharing all inversions
                clustered_ids[[i]] <- 
                    clustered_inversion_table[[i]] %>% 
                    dplyr::select(
                        c(overlap_cluster, accessions, inner_segment)
                        ) %>% 
                    group_by(overlap_cluster, inner_segment) %>%
                    separate_rows(accessions, sep = ", ") %>% 
                    ungroup %>%
                    group_by(accessions) %>% 
                    summarise(
                        inversion_clusters = 
                            paste(overlap_cluster, collapse = ", "),
                        segments = 
                            paste(inner_segment, collapse = ", "),
                    ) %>% 
                    arrange(inversion_clusters) %>% 
                    group_by(inversion_clusters) %>% 
                    add_count(
                        inversion_clusters, 
                        name = "grouped_sequences", sort = T
                    ) %>%
                    mutate(
                        inversions = i
                    )
                
                ## reset the inversion set
                inversion_set <- 
                    list()
            }
            
            combined_clustered_inversion_table <- 
                clustered_ids %>% 
                bind_rows %>% 
                arrange(desc(grouped_sequences)) %>%
                ungroup %>%
                mutate(cluster_id = rleid(inversion_clusters)) %>%
                dplyr::select(
                    c(
                        cluster_id, accessions, grouped_sequences, inversions
                    )
                )
            
            ## add the colinear genomes
            colinear_table <- 
                tibble(
                    cluster_id = 0,
                    accessions = 
                        c(
                            filtered_combined_delta$rid[1],
                            cluster_0$qid
                        ),
                    grouped_sequences = length(accessions),
                    inversions = 0
                )
            final_clusters_summary <- 
                bind_rows(
                    colinear_table,
                    combined_clustered_inversion_table
                )
            
            
            ## get the inversions
            inversion_coords <- 
                bind_rows(clustered_inversion_table)[,-1] %>%
                mutate(cluster = rleid(accessions))
            
            clusters_and_coords <-
                list(
                    cluster_summary = 
                        final_clusters_summary,
                    cluster_coords = 
                        inversion_coords
                )
        }
        return(clusters_and_coords)
        
    }

cluster_with_identity <- 
    function(
        clustered_genomes,
        mash_tib
    ){
        clust_identity <- 
            suppressWarnings(
                inner_join(
                    mash_tib,
                    (clustered_genomes[,c(1:2,4)] %>% 
                         `colnames<-`(c("ref_clust", "ref", "ref_invs_to_0"))),
                    by = "ref"
                ) %>%
                    inner_join(
                        ., 
                        (clustered_genomes[,c(1:2,4)] %>% 
                             `colnames<-`(
                                 c("qry_clust", "qry", "qry_invs_to_0"))
                        ), 
                        by = "qry"
                    ) %>%
                    mutate(
                        ref_clust =
                            factor(
                                ref_clust,
                                levels =
                                    clustered_genomes$cluster_id %>% unique
                            ),
                        qry_clust =
                            factor(
                                qry_clust,
                                levels =  
                                    clustered_genomes$cluster_id %>% unique
                            )
                    )
            )
        
        clust_identity_diff <- 
            clust_identity %>%
            filter(
                ref_clust != qry_clust
            ) %>%
            group_by(ref_clust, qry_clust) %>%
            filter(
                identity_sld==max(identity_sld) | identity ==max(identity)
            ) %>% 
            mutate(
                clusts = 
                    paste(
                        sort(
                            c(
                                as.character(ref_clust)[1], 
                                as.character(qry_clust)[1]
                            )
                        ), 
                        collapse = ""
                    )
            ) %>%
            ungroup %>%
            group_by(clusts) %>%
            filter(identity_sld==max(identity_sld)) %>% 
            arrange(clusts) %>% 
            filter(!duplicated(clusts))
        
        clust_identity_list <- 
            list(
                clust_identity = clust_identity,
                clust_identity_diff = clust_identity_diff
            )
        return(clust_identity_list)
    }



#### SV identification

Chromosome_level_SVs <- 
    function(
            delta_file
            ){
        fd <- 
            filter_delta(delta_file) %>%
            mutate(
                segment = rleid(strand),
                segments = length(rle(strand)$lengths),
                invs = sum(rle(strand)$values=="-")
            )
        
        fd2 <- 
            fd %>%
            filter( ## not the first or last alignment, 
                # rs !=min(rs), re !=max(re), 
                # qs !=1
                segment != min(segment),
                segment != max(segment)
            )%>%
            filter(
                X_dist<1e6,
                meanlen>=5e4
                )
        events <- 
            list()
        
        if (nrow(fd2)==0){
            
            events1 <- 
                fd %>%
                filter( ## not the first or last alignment, 
                    # rs !=min(rs), re !=max(re), 
                    # qs !=1
                    segment != min(segment),
                    segment != max(segment)
                ) 
            
            if (nrow(events1)==0){
                events <- 
                    tibble()
            } else {
                events <- 
                    events1 %>%
                    rowwise() %>%
                    mutate(
                        rsp = 
                            ((100*rs)/rlen) %>% round(.),
                        rep = 
                            ((100*re)/rlen) %>% round(.),
                        qsp = 
                            ((100*qs)/qlen) %>% round(.),
                        qep = 
                            ((100*qe)/qlen) %>% round(.),
                        rmp = round((rs+re)/2),
                        qmp = round((qs+qe)/2),
                        # r_midprop = 
                        #     round(100*rmp/rlen),
                        # q_midprop = 
                        #     round(100*qmp/qlen),
                        # meansize = 
                        #     round((100*meanlen)/mean(c(rlen, qlen))),
                        # rmprop = 
                        #     round(100*rmp/rlen),
                        # qmprop = 
                        #     round(100*qmp/qlen),
                        midpoint = 
                            mean(c(rmp,qmp)),
                            #mean(c(rmprop, qmprop)),
                        segment = as.character(segment),
                        ms =
                            ifelse(
                                strand=="+",
                                mean(c(rsp, qsp)),
                                mean(c(rsp, qep))
                            ),
                        me = 
                            ifelse(
                                strand=="+",
                                mean(c(rep, qep)),
                                mean(c(rep, qsp))
                            )
                    ) %>%
                    dplyr::select(
                        rid, qid,# rsp, rep, qsp, qep, 
                        strand, segment,
                        X_dist, meanlen, 
                        midpoint, ms, me, rlen, qlen,
                        rs, re, qs, qe
                        #r_midprop, q_midprop
                    ) %>%
                    ungroup()
                
            }
                
        } else{
            fd3 <- 
                fd %>%
                filter( ## not the first or last alignment, 
                    segment != min(segment),
                    segment != max(segment)
                )%>%
                filter(
                    X_dist<1e6,
                    meanlen >=5e4
                    ) %>%
                group_by(segment) %>%
                summarise(
                    strand = unique(strand),
                    rs = min(rs), 
                    re = max(re),
                    qs = ifelse(strand=="+", min(qs), max(qs)),
                    qe = ifelse(strand=="+", max(qe), min(qe)),
                    X_dist = round(mean(X_dist)),
                    rlen = unique(rlen),
                    qlen = unique(qlen),
                    meanlen = round(sum(meanlen)),
                    qid = unique(qid),
                    rid = unique(rid)
                )
            
            segdiff <- 
                fd3 %>%
                mutate(segdiff = segment-lag(segment, default = 1)) %>%
                pull(segdiff)
            
            if (any(segdiff>1)) {
                events <- 
                    fd3%>%
                    rowwise() %>%
                    mutate(
                        rsp = 
                            ((100*rs)/rlen) %>% round(.),
                        rep = 
                            ((100*re)/rlen) %>% round(.),
                        qsp = 
                            ((100*qs)/qlen) %>% round(.),
                        qep = 
                            ((100*qe)/qlen) %>% round(.),
                        rmp = round((rs+re)/2),
                        qmp = round((qs+qe)/2),
                        # r_midprop = 
                        #     round(100*rmp/rlen),
                        # q_midprop = 
                        #     round(100*qmp/qlen),
                        # meansize = 
                        #     round((100*meanlen)/mean(c(rlen, qlen))),
                        # rmprop = 
                        #     round(100*rmp/rlen),
                        # qmprop = 
                        #     round(100*qmp/qlen),
                        midpoint = 
                            mean(c(rmp,qmp)),
                        #mean(c(rmprop, qmprop)),
                        segment = as.character(segment),
                        ms =
                            ifelse(
                                strand=="+",
                                mean(c(rsp, qsp)),
                                mean(c(rsp, qep))
                            ),
                        me = 
                            ifelse(
                                strand=="+",
                                mean(c(rep, qep)),
                                mean(c(rep, qsp))
                            )
                    ) %>%
                    dplyr::select(
                        rid, qid,# rsp, rep, qsp, qep, 
                        strand, segment,
                        X_dist, meanlen, 
                        midpoint, ms, me, rlen, qlen,
                        rs, re, qs, qe
                        #r_midprop, q_midprop
                    ) %>%
                    ungroup()
                
            } else{
            
                fd2_props <- 
                    fd3 %>%
                    rowwise() %>%
                    mutate(
                        rsp = 
                            ((100*rs)/rlen) %>% round(.),
                        rep = 
                            ((100*re)/rlen) %>% round(.),
                        qsp = 
                            ((100*qs)/qlen) %>% round(.),
                        qep = 
                            ((100*qe)/qlen) %>% round(.),
                        rmp = round((rs+re)/2),
                        qmp = round((qs+qe)/2),
                        rmprop = round(100*rmp/rlen),
                        qmprop = round(100*qmp/qlen),
                        #meansize = round((100*meanlen)/mean(c(rlen, qlen))),
                        midpoint = 
                            mean(c(rmp,qmp)),
                            #mean(rmprop, qmprop)
                    ) %>%
                    dplyr::select(
                        rid, qid, rsp, rep, qsp, qep, strand, segment,
                        X_dist, meanlen, midpoint, rlen, qlen,
                        rs, re, qs, qe
                    ) %>%
                    ungroup() %>%
                    filter(
                        strand=="-" |
                            (
                                (lag(strand, default = "+")=="-") & 
                                    (lead(strand, default = "+")=="-")
                            )
                    )
                
                inner_seg <- 
                    median(unique(fd2_props$segment))
                
                if (nrow(fd2_props)==1){
                    middle_seg_row <- 
                        fd2_props %>%
                        filter(
                            segment==median(segment)
                        ) %>%
                        mutate(
                            ms =
                                ifelse(
                                    strand=="+",
                                    mean(c(rsp, qsp)),
                                    mean(c(rsp, qep))
                                ),
                            me = 
                                ifelse(
                                    strand=="+",
                                    mean(c(rep, qep)),
                                    mean(c(rep, qsp))
                                ),
                            segment = as.character(segment)
                        ) %>%
                        dplyr::select(
                            rid, qid, strand, segment, 
                            X_dist, meanlen, midpoint, ms, me, rlen, qlen,
                            rs, re, qs, qe
                        )
                    events <- 
                        middle_seg_row
                } else {
                    
                    middle_seg_row <- 
                        fd2_props %>%
                        filter(
                            segment==median(segment)
                        ) %>%
                        mutate(
                            ms =
                                ifelse(
                                    strand=="+",
                                    mean(c(rsp, qsp)),
                                    mean(c(rsp, qep))
                                ),
                            me = 
                                ifelse(
                                    strand=="+",
                                    mean(c(rep, qep)),
                                    mean(c(rep, qsp))
                                ),
                            segment = as.character(segment)
                        ) %>%
                        dplyr::select(
                            rid, qid, strand, segment, X_dist, 
                            meanlen, midpoint, ms, me, rlen, qlen,
                            rs, re, qs, qe
                        )
                    
                    outer_segs <- 
                        fd2_props %>%
                        filter(segment != median(segment)) 
                    
                    possible_outers <- 
                        nrow(outer_segs) %/% 2
                    outer_events <- 
                        list()
                    for (j in 1:possible_outers){
                        outer_events[[j]] <- 
                            fd2_props %>% 
                            filter(segment %in% c(inner_seg-j, inner_seg+j)) %>%
                            summarise(
                                rid = unique(rid),
                                qid = unique(qid),
                                strand= unique(strand),
                                segment = paste(segment, collapse = ","),
                                X_dist = mean(X_dist),
                                ms =
                                    ifelse(
                                        strand=="+",
                                        mean(c(min(rsp), min(qsp))),
                                        mean(c(min(rsp), min(qep)))
                                    ),
                                me =
                                    ifelse(
                                        strand=="+",
                                        mean(c(max(rep), max(qep))),
                                        mean(c(max(rep), max(qsp)))
                                    ),
                                rs = 
                                    ifelse(
                                        strand=="+",
                                        min(re),
                                        min(rs)
                                    ),
                                re = 
                                    ifelse(
                                        strand=="+",
                                        max(rs),
                                        max(re)
                                    ),
                                qs = 
                                    ifelse(
                                        strand=="+",
                                        min(qe),
                                        min(qs)
                                    ),
                                qe = 
                                    ifelse(
                                        strand=="+",
                                        max(qs),
                                        max(qe)
                                    ),
                                
                                meanlen = 
                                    me-ms,
                                midpoint = 
                                    (ms+me)/2,
                                rlen = unique(rlen),
                                qlen = unique(qlen)
                            )
                    }
                    inv_summary <- 
                        bind_rows(
                            middle_seg_row,
                            bind_rows(outer_events)
                        )
                    
                    events <- 
                        inv_summary 
                }
                
        }
    }
        
        return(events)
    }




### final SV detection
assign_domains <- 
    function(midpoints, genome_length) {
        
        pct <- (midpoints / genome_length) * 100
        
        case_when(
            pct > 87.5 | pct <= 12.5   ~ "O",
            pct > 12.5  & pct <= 37.5  ~ "R",
            pct > 37.5  & pct <= 62.5  ~ "T",
            pct > 62.5  & pct <= 87.5  ~ "L",
            TRUE ~ NA_character_
        )
    }

delta_indels <- 
    function(
        delta_table, minlen = 50
    ){
        if (is.character(delta_table)){
            delta_table <- 
                read_delta(delta_table)
        }
        
        d2 <- 
            delta_table %>%
            rowwise() %>%
            mutate(
                qs2 = min(c(qs, qe)),
                qe2 = max(c(qs, qe))
            ) %>% 
            ungroup
        
        ## gap insertion positions
        # insertions
        ref_ranges <- 
            IRanges(
                start = 
                    d2$rs,
                end = 
                    d2$re
            )
        qry_ranges <- 
            IRanges(
                start = 
                    d2$qs2,
                end = 
                    d2$qe2
            ) 
        qry_gr <- 
            GRanges(
                seqnames = d2$qid,
                ranges   = IRanges(start = d2$qs2, end = d2$qe2),
                ref_break = d2$re   
        )
        qry_gaps <- 
            reduce(qry_ranges) %>%
            gaps(start = 1, end = d2$qlen[1]) 
        if (length(qry_gaps)==0){
            insertions_df <- NULL
        } else {
            
            ins_gr <- 
                qry_gaps %>%
                subset(start > 1 & end < d2$qlen[1]) %>%
                GRanges(seqnames = d2$qid[1], strand = "*")
            
            pre_hits <- 
                precede(ins_gr, qry_gr)
            
            insertions_df <- 
                as.data.frame(ins_gr) %>%
                as_tibble() %>%
                mutate(
                    ref_insertion_point = qry_gr$ref_break[pre_hits],
                    variant = "Insertion",
                    rid = d2$rid[1],
                    qid = d2$qid[1]
                ) %>%
                dplyr::select(
                    c(rid, qid, start, end, width, ref_insertion_point, variant)
                )
        }
        #insertions_df %>% colnames()
        ### deletions in the query are gaps in the reference
        
        deletions <- 
            ref_ranges %>%
            reduce() %>%
            gaps(start = 1, end = d2$rlen[1]) %>%
            subset(start > 1 & end < d2$rlen[1])
        deletions_df <- 
            deletions %>%
            as.data.frame() %>% 
            mutate(
                rid = d2$rid[1],
                qid = d2$qid[1],
                variant = "Deletion"
            ) %>%
            as_tibble()
        
        
        all_indels <- 
            bind_rows(
                insertions_df, 
                deletions_df
            ) %>%
            filter(width>=minlen)
        
        return(all_indels)
    }



delta_substructural_inversions <- 
    function(
        delta_table,
        delta_SR
    ){
        if (is.character(delta_table)){
            delta_table <- 
                read_delta(delta_table)
        }
        
        inv_list <- 
            list()
        if (
            (nrow(delta_SR)>0) 
        ){
            for (j in 1:nrow(delta_SR)){
                
                if (j>1){
                    dtj <- 
                        delta_table %>%
                        filter(
                            (re <= (delta_SR$rs[j-1])) |
                                (rs >= (delta_SR$re[j-1]))
                        )
                } else{
                    dtj <- 
                        delta_table
                }
                CSV <- 
                    delta_SR[j,]
                
                inv_tib <- 
                    dtj %>% 
                    filter(
                        (rs >= (CSV$rs[1])) & (rs <= (CSV$re[1])),
                        strand!=CSV$strand[1]
                    ) %>%
                    mutate(
                        variant = 
                            ifelse(
                                X_dist < 1e6,
                                "Ori inversion",
                                "local inversion"
                            )
                    ) %>%
                    dplyr::select(
                        c(
                            rs, re, qs, qe, meanlen, variant, rid, qid
                        )
                    ) %>%
                    mutate(
                        variant = as.character(variant)
                    )
                
                inv_list[[j]] <- 
                    inv_tib
            }
            
            invs_outer <- 
                delta_table %>%
                filter(
                    (re <= (delta_SR$rs[j])) |
                        (rs >= (delta_SR$re[j])),
                    strand=="-"
                ) %>% 
                mutate(
                    variant = 
                        ifelse(
                            X_dist < 1e6,
                            "Ori inversion",
                            "local inversion"
                        )
                ) %>%
                dplyr::select(
                    c(
                        rs, re, qs, qe, meanlen, variant, rid, qid
                    )
                ) %>%
                mutate(
                    variant = as.character(variant)
                )
            
            inv_list[[j+1]] <- 
                invs_outer
        } else {
            inv_list <- 
                delta_table %>%
                filter(
                    strand=="-"
                ) %>%
                mutate(
                    variant = 
                        ifelse(
                            X_dist < 1e6,
                            "Ori inversion",
                            "local inversion"
                        )
                ) %>%
                dplyr::select(
                    c(
                        rs, re, qs, qe, meanlen, variant, rid, qid
                    )
                ) %>%
                mutate(
                    variant = as.character(variant)
                )
            
        }
        
        invs <- 
            bind_rows(
                inv_list
            ) 
        return(invs)
    }





delta_structural_rearrangements <- 
    function(
        delta_table
    ){
        
        if (is.character(delta_table)){
            delta_table <- 
                read_delta(delta_table)
        }
        
        d2 <- 
            delta_table %>%
            rowwise() %>%
            mutate(
                qs2 = min(c(qs, qe)),
                qe2 = max(c(qs, qe))
            ) %>% 
            ungroup()
        
        ## segment alignments
        fd <- 
            #filter_delta(delta_file) %>%
            filter_delta(d2) %>%
            mutate(
                segment = rleid(strand),
                # segments = length(rle(strand)$lengths),
                # invs = sum(rle(strand)$values=="-")
            ) %>%
            filter(
                segment != min(segment), 
                segment != max(segment), 
                meanlen >= 1e4
            ) 
        
        events <- 
            list()
        
        if (nrow(fd)==0){
            
            events <- 
                tibble()
            
        } else{
            fd2 <- 
                fd %>%
                mutate(
                    qs2 = min(c(qs, qe)),
                    qe2 = max(c(qs, qe))
                ) %>%
                ungroup()
            fd3 <- 
                fd2 %>%
                group_by(segment) %>%
                summarise(
                    rid = unique(rid),
                    qid = unique(qid),
                    strand   = unique(strand),
                    rs       = min(rs), re = max(re),
                    qs       = if_else(strand=="+", min(qs2), max(qs2)),
                    qe       = if_else(strand=="+", max(qe2), min(qe2)),
                    X_dist   = round(mean(X_dist)),
                    meanlen  = sum(meanlen),
                    rlen     = unique(rlen), qlen = unique(qlen),
                    rmp = round((rs + re) / 2),
                    qmp = round((qs + qe) / 2),
                    mid = round((rmp + qmp) / 2),
                    .groups  = "drop"
                ) 
            
            
            segdiff <- 
                fd3 %>%
                mutate(segdiff = segment-lag(segment, default = 1)) %>%
                pull(segdiff)
            
            if (any(segdiff>1)) {
                events <- 
                    fd3# %>%
                # dplyr::select(
                #     rid, qid,# rsp, rep, qsp, qep, 
                #     strand, segment,
                #     X_dist, meanlen, 
                #     midpoint, ms, me, rlen, qlen,
                #     rs, re, qs, qe
                #     #r_midprop, q_midprop
                # ) %>%
                # ungroup()
                
            } else{
                
                fd2_props <- 
                    fd2 %>%
                    group_by(segment) %>%
                    summarise(
                        rid = unique(rid),
                        qid = unique(qid),
                        strand   = unique(strand),
                        rs       = min(rs), re = max(re),
                        qs       = if_else(strand=="+", min(qs2), max(qs2)),
                        qe       = if_else(strand=="+", max(qe2), min(qe2)),
                        X_dist   = round(mean(X_dist)),
                        meanlen  = sum(meanlen),
                        rlen     = unique(rlen), qlen = unique(qlen),
                        rmp = round((rs + re) / 2),
                        qmp = round((qs + qe) / 2),
                        mid = round((rmp + qmp) / 2),
                        .groups  = "drop"
                    )  %>%
                    ungroup() %>%
                    arrange(segment) %>%
                    filter(
                        strand=="-" |
                            (
                                (lag(strand, default = "+")=="-") & 
                                    (lead(strand, default = "+")=="-")
                            )
                    )
                seg_ids <- 
                    sort(unique(fd2_props$segment))
                inner_seg <- 
                    seg_ids[ ceiling(length(seg_ids)/2) ]
                # inner_seg <- 
                #     median(unique(fd2_props$segment))
                
                if (nrow(fd2_props)==1){
                    middle_seg_row <- 
                        fd2_props %>%
                        filter(
                            segment==inner_seg #median(segment)
                        ) %>%
                        mutate(
                            ms =
                                round(mean(c(rs, qs))),
                            me = 
                                round(mean(c(re, qe))),
                            segment = as.character(segment)
                        )
                    events <- 
                        middle_seg_row
                } else {
                    
                    middle_seg_row <- 
                        fd2_props %>%
                        filter(
                            segment==inner_seg #median(segment)
                        ) %>%
                        mutate(
                            ms =
                                round(mean(c(rs, qs))),
                            me = 
                                round(mean(c(re, qe))),
                            segment = as.character(segment)
                        )
                    
                    outer_segs <- 
                        fd2_props %>%
                        filter(segment !=inner_seg) #median(segment)) 
                    
                    possible_outers <- 
                        nrow(outer_segs) %/% 2
                    
                    outer_events <- 
                        list()
                    for (j in 1:possible_outers){
                        outer_events[[j]] <- 
                            fd2_props %>% 
                            filter(
                                segment %in% c(inner_seg-j, inner_seg+j)
                            ) %>%
                            mutate(
                                expected_nested_mid = 
                                    mean(c(mid, lag(mid)), na.rm = TRUE),
                                mid_diff = 
                                    abs(expected_nested_mid -
                                            middle_seg_row$mid),
                                nested_group = 
                                    ifelse(
                                        mid_diff <
                                            5e5,
                                        c(1,1),
                                        c(1,2)
                                    )
                            ) %>%
                            group_by(nested_group) %>%
                            summarise(
                                rid = unique(rid),
                                qid = unique(qid),
                                strand= unique(strand),
                                segment = paste(segment, collapse = ","),
                                X_dist = mean(X_dist),
                                ms =
                                    mean(c(min(rs), min(qs))),
                                me =
                                    mean(c(max(re), max(qe))),
                                rs = 
                                    min(rs),
                                re = 
                                    max(re),
                                qs = 
                                    min(qs),
                                qe = 
                                    max(qe),
                                meanlen = 
                                    me-ms,
                                midpoint = 
                                    (ms+me)/2,
                                rlen = unique(rlen),
                                qlen = unique(qlen)
                            )
                    }
                    inv_summary <- 
                        bind_rows(
                            middle_seg_row,
                            bind_rows(outer_events)
                        )
                    
                    events <- 
                        inv_summary 
                }
                
            }
        }
        
        return(events)
    }




delta_duplications <- 
    function(
        delta_table,
        unfiltered_delta_table,
        minoverlap    = 50,
        max_tandem_gap = 10000,
        min_prop = 0.9
        ) {
    
        if (is.character(delta_table)) {
            delta_table <- read_delta(delta_table)
        }
        if (is.character(unfiltered_delta_table)) {
            unfiltered_delta_table <- read_delta(unfiltered_delta_table)
        }
        
        
        ufd <- 
            unfiltered_delta_table %>%
            mutate(
                qs2 = pmin(qs, qe),
                qe2 = pmax(qs, qe)
            ) %>%
            ungroup()
        
        stopifnot(all(c("qs","qe","rs","re") %in% names(ufd)))
    
        ref_len <- unique(ufd$rlen)
        stopifnot(length(ref_len) == 1)
        
        primary_gr <- 
            GRanges(
                seqnames = delta_table$rid,
                ranges = 
                    IRanges(
                        start = delta_table$rs, 
                        end = delta_table$re)
                )
        
        
        ufd_qry_gr <- 
            GRanges(
                seqnames   = ufd$qid,
                ranges     = IRanges(start = ufd$qs2, end = ufd$qe2),
                strand     = ufd$strand,
                rs         = ufd$rs,      # reference start
                re         = ufd$re,      # reference end
                qry_start  = ufd$qs2,
                qry_end    = ufd$qe2,
                qry_name   = ufd$qid,
                ref_name   = ufd$rid
            )
        
        find_dups <- 
            function(gr, source_label) {
                dup_summary <- 
                    as.data.frame(findOverlaps(gr, gr, minoverlap = minoverlap))
                dup_summary <- 
                    dup_summary[dup_summary$queryHits != dup_summary$subjectHits, ]
                
                starts   <- start(gr)
                ends     <- end(gr)
                meta     <- mcols(gr)
                strands  <- as.character(strand(gr))
                
                dup_summary %>%
                    transmute(
                        row1 = queryHits,
                        row2 = subjectHits,
                        pair_id = 
                            paste0(pmin(row1, row2), "_", pmax(row1, row2)),
                        
                        # full alignment coordinates
                        orig_ref_start = meta$rs[row1],
                        orig_ref_end   = meta$re[row1],
                        dup_ref_start  = meta$rs[row2],
                        dup_ref_end    = meta$re[row2],
                        
                        # query-space overlap length
                        dup_len = pmax(
                            0,
                            pmin(ends[row1], ends[row2]) -
                                pmax(starts[row1], starts[row2]) + 1
                        ),
                        
                        # projected duplicated region in reference space
                        orig_ref_dup_start = round(
                            meta$rs[row1] +
                                ((pmax(starts[row1], starts[row2]) - 
                                      starts[row1]) /
                                     (ends[row1] - starts[row1])) * 
                                (meta$re[row1] - meta$rs[row1])
                        ),
                        
                        orig_ref_dup_end = orig_ref_dup_start + dup_len - 1,
                        
                        # reference separation between copies
                        dup_gap = case_when(
                            dup_ref_start > 
                                orig_ref_end ~ dup_ref_start - orig_ref_end,
                            orig_ref_start > 
                                dup_ref_end ~ orig_ref_start - dup_ref_end,
                            TRUE ~ 0
                        ),
                        
                        duplication_type = if_else(dup_gap <= max_tandem_gap,
                                                   "tandem", "dispersed"),
                        
                        orientation = if_else(strands[row1] == strands[row2],
                                              "same", "inverted"),
                        
                        ref_name = meta$ref_name[row1],
                        qry_name = meta$qry_name[row1],
                        source   = source_label
                    ) %>%
                    filter(
                        !(dup_ref_start >= 
                              orig_ref_start & dup_ref_end <= orig_ref_end) &
                            !(orig_ref_start >= 
                                  dup_ref_start & orig_ref_end <= dup_ref_end)
                    ) %>%
                    distinct(pair_id, .keep_all = TRUE)
        }
        
        
        
        ref_dups <- 
            find_dups(ufd_qry_gr, "ref")
        
        prim_dup_hits <- 
            findOverlaps(
                    IRanges(
                        ref_dups$dup_ref_start,
                        ref_dups$dup_ref_end
                        ), 
                    ranges(primary_gr), 
                    type = "within"
                    )
        
        
        if (length(prim_dup_hits) > 0) {
            ref_dups <- ref_dups[-queryHits(prim_dup_hits), ]
        }
        
        if (nrow(ref_dups) == 0){
            return(
                tibble(
                    ref_name = delta_table$rid[1],
                    qry_name = delta_table$qid[1]
                )
                )
        } else {
            
            all_dups <- 
                ref_dups %>%
                mutate(
                    orig_mid      = round((orig_ref_start + orig_ref_end) / 2),
                    dup_mid       = round((dup_ref_start  + dup_ref_end)  / 2),
                    mid_domain    = assign_domains(orig_mid, ref_len),
                    refmid_domain = assign_domains(dup_mid,  ref_len)
                ) %>%
                group_by(pair_id) %>%
                slice_max(order_by = dup_len, n = 1, with_ties = FALSE) %>%
                ungroup() %>%
                transmute(
                    rid = ref_name, qid = qry_name,
                    orig_ref_dup_start,
                    orig_ref_dup_end,
                    orig_ref_start, orig_ref_end, orig_mid,    mid_domain,
                    dup_ref_start,  dup_ref_end,  dup_mid,     refmid_domain,
                    duplication_type, orientation,
                    dup_len, dup_gap,
                    domains       = paste(mid_domain, refmid_domain, sep = "_"),
                    orig_ref_len  = orig_ref_end  - orig_ref_start + 1,
                    dup_ref_len   = dup_ref_end   - dup_ref_start  + 1,
                    prop_orig     = dup_len / orig_ref_len,
                    prop_dup      = dup_len / dup_ref_len,
                    prop_min      = dup_len / pmin(orig_ref_len, dup_ref_len)
                )
            
            
            ad2 <- 
                all_dups %>%
                    ungroup %>%
                    filter(
                        prop_min > min_prop
                        #duplication_type=="dispersed"
                        ) %>%
                    dplyr::select(
                        c(
                            rid, qid,
                            orig_ref_dup_start, orig_ref_dup_end,
                            orig_ref_start, orig_ref_end, 
                            dup_ref_start, dup_ref_end, dup_len, dup_gap,
                            domains,
                            prop_orig, duplication_type
                            )
                        )
            # collapse by similarity of the projected duplicated reference region
            
            # GRanges based on **projected** duplicated region
            dup_proj_gr <- GRanges(
                ad2$rid,
                IRanges(ad2$orig_ref_dup_start, ad2$orig_ref_dup_end)
            )
            
            # GRanges for the destination duplication locations
            dup_target_gr <- GRanges(
                ad2$rid,
                IRanges(ad2$dup_ref_start, ad2$dup_ref_end)
            )
            
            # find reciprocal overlaps ≥ 90% on both sides
            ov1 <- as.data.frame(findOverlaps(dup_proj_gr, dup_proj_gr, minoverlap = 1))
            ov2 <- as.data.frame(findOverlaps(dup_target_gr, dup_target_gr, minoverlap = 1))
            
            good1 <- ov1 %>%
                filter(
                    width(pintersect(dup_proj_gr[queryHits], dup_proj_gr[subjectHits])) >=
                        min_prop * pmin(width(dup_proj_gr[queryHits]), width(dup_proj_gr[subjectHits]))
                )
            
            good2 <- ov2 %>%
                filter(
                    width(pintersect(dup_target_gr[queryHits], dup_target_gr[subjectHits])) >=
                        min_prop * pmin(width(dup_target_gr[queryHits]), width(dup_target_gr[subjectHits]))
                )
            
            common <- inner_join(good1, good2, by = c("queryHits", "subjectHits"))
            
            if (nrow(common)==0){
                return(
                    tibble(
                        ref_name = delta_table$rid[1],
                        qry_name = delta_table$qid[1]
                    )
                )
            } else {
                edges <- common %>%
                    filter(queryHits < subjectHits) %>%
                    select(queryHits, subjectHits) %>%
                    as.matrix()
                
                nr <- nrow(ad2)
                edges_df <- common %>%
                    filter(queryHits < subjectHits) %>%
                    transmute(from = queryHits, to = subjectHits)
                
                g <- graph_from_data_frame(d = edges_df, directed = FALSE,
                                           vertices = data.frame(name = 1:nr))
                
                cl <- components(g)$membership
                
                
                dups2 <- ad2 %>%
                    mutate(cluster = cl) %>%
                    group_by(cluster) %>%
                    summarize(
                        ref_name       = unique(rid),
                        qry_name       = unique(qid),
                        
                        dup_s = min(orig_ref_dup_start),
                        dup_e = max(orig_ref_dup_end),
                        
                        # other context
                        orig_ref_start = min(orig_ref_start),
                        orig_ref_end   = max(orig_ref_end),
                        dup_ref_start  = min(dup_ref_start),
                        dup_ref_end    = max(dup_ref_end),
                        
                        # recompute length
                        dup_len        = dup_e - dup_s + 1,
                        dup_gap        = min(dup_gap),
                        domains        = unique(domains),
                        duplication_type = unique(duplication_type),
                        collapsed      = n(),
                        .groups        = "drop"
                    )
                
                
                
                return(dups2)
            }
            
            
        }
}




delta_duplications_prev2 <- 
    function(delta_table, minoverlap = 50, max_tandem_gap = 10000){
        if (is.character(delta_table)){
            delta_table <- 
                read_delta(delta_table)
        }
        
        ufd <- 
            delta_table %>%
            mutate(
                qs2 = pmin(qs, qe),
                qe2 = pmax(qs, qe)
            ) %>% 
            ungroup()
        ref_len <- 
            (ufd %>% pull(rlen))[1]
        qry_len <- 
            (ufd %>% pull(qlen))[1]
        stopifnot(length(unique(ufd$rlen)) == 1,
                  length(unique(ufd$qlen)) == 1)
        
        ## ufd is the unfiltered delta 
        # ufd <- 
        #     d2
        
        ufd_ref_gr <- 
            GRanges(
                seqnames = ufd$rid,
                ranges = IRanges(start = ufd$rs, end = ufd$re),
                strand = ufd$strand,
                qry_start = ufd$qs2,
                qry_end = ufd$qe2,
                qry_name = ufd$qid,
                ref_name  = ufd$rid,   
                rs = ufd$rs,
                re = ufd$re
            )
        
        ufd_qry_gr <- 
            GRanges(
                seqnames = ufd$qid,
                ranges = IRanges(start = ufd$qs2, end = ufd$qe2),
                strand = ufd$strand,
                qry_start = ufd$qs2,
                qry_end = ufd$qe2,
                qry_name = ufd$qid,
                ref_name  = ufd$rid,   
                rs = ufd$rs,
                re = ufd$re
            )
        
        find_dups <- 
            function(gr, source_label) {
                gr_hits <- 
                    IRanges::findOverlaps(gr, gr, minoverlap = minoverlap)
                gr_df   <- 
                    as.data.frame(gr_hits) %>% filter(queryHits != subjectHits)
            
            # pull out the vectorized starts and ends once
                starts <- start(gr)
                ends   <- end(gr)
                meta   <- mcols(gr)
                strands<- as.character(strand(gr))
                
                
                
            dups_tib <- 
                gr_df %>%
                #rowwise() %>%
                mutate(
                    row1        = queryHits,
                    row2        = subjectHits,
                    pair_id = 
                        paste0(
                            pmin(queryHits, subjectHits),
                            "_",
                            pmax(queryHits, subjectHits)
                        ),
                        #paste(sort(c(queryHits, subjectHits)), collapse = "_"),
                    ## ORIGINAL copy on ref:
                    orig_ref_start  = meta$rs[row1],       
                    orig_ref_end    = meta$re[row1],       
                    ## DUPLICATED copy on ref:
                    dup_ref_start   = meta$rs[row2],       
                    dup_ref_end     = meta$re[row2],      
                    
                    a1          = starts[queryHits],
                    a2          = ends[queryHits],
                    b1          = starts[subjectHits],
                    b2          = ends[subjectHits],
                    
                    dup_len     = pmax(0, pmin(a2, b2) - pmax(a1, b1) + 1),
                    
                    ref_name = meta$ref_name[row1],
                    qry_name = as.character(meta$qry_name[row1]),
                    qry_start   = meta$qry_start[queryHits],
                    qry_end     = meta$qry_end[queryHits],
                    
                    dup_start   = meta$qry_start[subjectHits],
                    dup_end     = meta$qry_end[subjectHits],
                    
                    strand1     = as.character(strand(gr)[queryHits]),
                    strand2     = as.character(strand(gr)[subjectHits]),
                    dup_gap = 
                        dplyr::case_when(
                            dup_ref_start > orig_ref_end  ~
                                dup_ref_start - orig_ref_end,  # B is to the right
                            orig_ref_start > dup_ref_end  ~ 
                                orig_ref_start - dup_ref_end,  # B is to the left
                            TRUE                          ~ 0                          
                        ),
                    duplication_type = 
                        if_else(
                            dup_gap <= 
                                max_tandem_gap, "tandem", "dispersed"
                            ),
                    orientation    = if_else(strands[row1] == strands[row2],
                                             "same", "inverted"),
                    source      = source_label
                ) %>%
                ungroup() %>%
                distinct(pair_id, .keep_all = TRUE) %>%
                dplyr::select(
                    pair_id, row1, row2,
                    orig_ref_start, orig_ref_end,       
                    dup_ref_start,  dup_ref_end,        
                    ref_name,# ref_start, ref_end,
                    qry_name, #qry_start, qry_end,
                    dup_start, dup_end,
                    dup_len, dup_gap, duplication_type,
                    source, orientation
                )
            return(dups_tib)
            }
        
        ref_dups <- 
            find_dups(ufd_ref_gr, "ref") %>% 
            filter(
                !(dup_ref_start >= orig_ref_start & dup_ref_end <= orig_ref_end)
            ) 
        # qry_dups <- 
        #     find_dups(ufd_qry_gr, "qry") %>%
        #     filter(
        #         !(dup_ref_start >= orig_ref_start & dup_ref_end <= orig_ref_end)
        #     ) 
        
        
        ref_gr <- 
            GenomicRanges::GRanges(
                seqnames = ref_dups$ref_name,
                ranges   = IRanges::IRanges(
                start = ref_dups$dup_ref_start,
                end   = ref_dups$dup_ref_end
            )
        )
        qry_gr <- 
            GenomicRanges::GRanges(
                seqnames = qry_dups$qry_name,
                ranges   = IRanges::IRanges(
                start = qry_dups$dup_ref_start,
                end   = qry_dups$dup_ref_end
            )
        )
        
        ov <- 
            GenomicRanges::findOverlaps(
                ranges(ref_gr), ranges(qry_gr), 
                minoverlap = 1
                )
        shared_ref_idx <- unique(S4Vectors::queryHits(ov))
        shared_qry_idx <- unique(S4Vectors::subjectHits(ov))
        
        
        all_dups <- 
            bind_rows(
              ref_dups,
              qry_dups
              )%>%
            mutate(
                orig_mid      = round((orig_ref_start + orig_ref_end) / 2),
                dup_mid       = round((dup_ref_start  + dup_ref_end)  / 2),
                
                mid_domain    = assign_domains(orig_mid, ref_len), 
                refmid_domain = assign_domains(dup_mid,  ref_len)           
                ) %>%
            transmute(
                rq = paste(ref_name, qry_name, sep = "_"),
                pair_id, source,
                ref_name, qry_name,
                orig_ref_start, orig_ref_end, orig_mid, mid_domain,
                dup_ref_start, dup_ref_end, dup_mid, refmid_domain,
                duplication_type, orientation,
                dup_len, dup_gap,
                domains = paste(mid_domain, refmid_domain, sep ="_")
            )  %>% 
            ungroup() %>%
            group_by(pair_id) %>%
            slice_max(order_by = dup_len, n = 1, with_ties = FALSE) %>% 
            ungroup() %>%
            mutate(idx = row_number()) %>%
            mutate(shared = idx %in% c(shared_ref_idx, shared_qry_idx)) %>%
            dplyr::select(-idx)  %>%
            mutate(
                orig_ref_len = orig_ref_end  - orig_ref_start + 1,
                dup_ref_len  = dup_ref_end   - dup_ref_start  + 1,
                prop_orig = dup_len / orig_ref_len,
                prop_dup  = dup_len / dup_ref_len,
                prop_min  = dup_len / pmin(orig_ref_len, dup_ref_len)
            )
        
        
        return(all_dups)
        
    }






delta_duplications_prev <- 
    function(delta_table){
        
        if (is.character(delta_table)){
            delta_table <- 
                read_delta(delta_table)
        }
        
        d2 <- 
            delta_table %>%
            rowwise() %>%
            mutate(
                qs2 = min(c(qs, qe)),
                qe2 = max(c(qs, qe))
            ) %>% 
            ungroup()
        ref_len <- 
            (d2 %>% pull(rlen))[1]
        qry_len <- 
            (d2 %>% pull(qlen))[1]
        
        ## ufd is the unfiltered delta 
        ufd <- 
            d2
        
        ufd_ref_gr <- 
            GRanges(
                seqnames = ufd$rid,
                ranges = IRanges(start = ufd$rs, end = ufd$re),
                strand = ufd$strand,
                qry_start = ufd$qs2,
                qry_end = ufd$qe2,
                qry_name = ufd$qid
            )
        
        ufd_qry_gr <- 
            GRanges(
                seqnames = ufd$qid,
                ranges = IRanges(start = ufd$qs2, end = ufd$qe2),
                strand = ufd$strand,
                qry_start = ufd$qs2,
                qry_end = ufd$qe2,
                qry_name = ufd$rid
            )
        
        qry_dups <- 
            IRanges::findOverlaps(ufd_ref_gr, ufd_ref_gr, minoverlap = 50) %>%
            as.data.frame() %>%
            filter(queryHits != subjectHits) %>%
            rowwise() %>%
            mutate(
                pair_id     = paste(sort(c(queryHits, subjectHits)), collapse = "_"),
                a1          = start(ufd_ref_gr)[queryHits],
                a2          = end(ufd_ref_gr)[queryHits],
                b1          = start(ufd_ref_gr)[subjectHits],
                b2          = end(ufd_ref_gr)[subjectHits],
                dup_len     = pmax(0, pmin(a2, b2) - pmax(a1, b1) + 1),
                ref_name    = as.character(seqnames(ufd_ref_gr)[queryHits]),
                ref_start   = a1,
                ref_end     = a2,
                qry_name    = mcols(ufd_ref_gr)$qry_name[queryHits],
                qry_start   = mcols(ufd_ref_gr)$qry_start[queryHits],
                qry_end     = mcols(ufd_ref_gr)$qry_end[queryHits],
                dup_start   = mcols(ufd_ref_gr)$qry_start[subjectHits],
                dup_end     = mcols(ufd_ref_gr)$qry_end[subjectHits],
                dup_dist    = abs(dup_start - qry_start),
                strand1     = as.character(strand(ufd_ref_gr)[queryHits]),
                strand2     = as.character(strand(ufd_ref_gr)[subjectHits]),
                orientation = if_else(strand1 == strand2, "same", "inverted"),
                duplication_type = 
                    if_else(dup_dist <= 10000, "tandem", "dispersed"),
                source      = "qry"
            ) %>%
            ungroup() %>%
            distinct(pair_id, .keep_all = TRUE) %>%
            dplyr::select(
                pair_id, ref_name, ref_start, ref_end,
                qry_name, qry_start, qry_end,
                dup_start, dup_end,
                dup_len, dup_dist, duplication_type, source, orientation
            )
        
        ref_dups <- 
            IRanges::findOverlaps(ufd_qry_gr, ufd_qry_gr, minoverlap = 50) %>%
            as.data.frame() %>%
            filter(queryHits != subjectHits) %>%
            rowwise() %>%
            mutate(
                pair_id     = paste(sort(c(queryHits, subjectHits)), collapse = "_"),
                a1          = start(ufd_qry_gr)[queryHits],
                a2          = end(ufd_qry_gr)[queryHits],
                b1          = start(ufd_qry_gr)[subjectHits],
                b2          = end(ufd_qry_gr)[subjectHits],
                dup_len     = pmax(0, pmin(a2, b2) - pmax(a1, b1) + 1),
                ref_name    = as.character(seqnames(ufd_qry_gr)[queryHits]),
                ref_start   = a1,
                ref_end     = a2,
                qry_name    = mcols(ufd_qry_gr)$qry_name[queryHits],
                qry_start   = mcols(ufd_qry_gr)$qry_start[queryHits],
                qry_end     = mcols(ufd_qry_gr)$qry_end[queryHits],
                dup_start   = mcols(ufd_qry_gr)$qry_start[subjectHits],
                dup_end     = mcols(ufd_qry_gr)$qry_end[subjectHits],
                dup_dist    = abs(dup_start - qry_start),
                strand1     = as.character(strand(ufd_qry_gr)[queryHits]),
                strand2     = as.character(strand(ufd_qry_gr)[subjectHits]),
                orientation = if_else(strand1 == strand2, "same", "inverted"),
                duplication_type = if_else(dup_dist <= 10000, "tandem", "dispersed"),
                source      = "ref"
            ) %>%
            ungroup() %>%
            distinct(pair_id, .keep_all = TRUE) %>%
            dplyr::select(
                pair_id, ref_name, ref_start, ref_end,
                qry_name, qry_start, qry_end,
                dup_start, dup_end,
                dup_len, dup_dist, duplication_type, source, orientation
            )
        
        
        out_tib <- 
            bind_rows(ref_dups, qry_dups) %>% 
            mutate(
                ref_midpoint = round((ref_start + ref_end) / 2),
                qry_midpoint_1 = round((qry_start + qry_end) / 2),
                qry_midpoint_2 = round((dup_start + dup_end) / 2),
                ref_domain   = assign_domains(ref_midpoint, ref_len),
                qry_domain_1 = assign_domains(qry_midpoint_1, qry_len),
                qry_domain_2   = assign_domains(qry_midpoint_2, qry_len),
                rq_replichores = 
                    paste(
                        ref_domain, 
                        qry_domain_2,
                        sep = "-"
                    ),
                qq_replichores = 
                    paste(
                        #ref_domain, 
                        #qry_domain, 
                        qry_domain_1, 
                        qry_domain_2,
                        sep = "-"
                        )
            )
        
        return(out_tib)
        
    }



delta_translocations <- 
    function(
        delta_table,
        minlen = 50
    ){
        
        if (is.character(delta_table)){
            delta_table <- 
                read_delta(delta_table)
        }
        
        d2 <- 
            delta_table %>%
            rowwise() %>%
            mutate(
                qs2 = min(c(qs, qe)),
                qe2 = max(c(qs, qe))
            ) %>% 
            ungroup()
        
        
        
        tlocs <- 
            delta_table %>%
            filter(
                !is.na(X_dist),
                X_dist >quantile(X_dist, 0.95)
            ) %>%
            mutate(
                #clust = clust_i,
                refmid = round((rs + re) / 2),
                qrymid = round((qs + qe) / 2),
                ref_domain = assign_domains(refmid, rlen),
                qry_domain = assign_domains(qrymid, qlen),
                replichores = paste(ref_domain, qry_domain, sep = "-")
            ) %>%
            rowwise() %>%
            mutate(
                start = 
                    rs,
                #round(mean(c(rs, qe))),
                end = 
                    re,
                #round(mean(c(re, qs))),
                width =  abs(start-end+1),
                variant = "Translocation"
            )  %>%
            ungroup() %>%
            dplyr::select(
                c(
                    start, end, width, refmid, qrymid, 
                    rid, qid, strand, variant,replichores,
                    X_dist
                )
            )
        
        return(tlocs)
    }


#####

delta_SVs <- 
    function(
        delta_table,
        chrom_SVs
    ){
        if (is.character(delta_table)){
            delta_table <- 
                read_delta(delta_table)
        }
        
        d2 <- 
            delta_table %>%
            rowwise() %>%
            mutate(
                qs2 = min(c(qs, qe)),
                qe2 = max(c(qs, qe))
            ) %>% 
            ungroup
        
        ## gap insertion positions
        # insertions
        ref_ranges <- 
            IRanges(
                start = 
                    d2$rs,
                end = 
                    d2$re
            )#
        qry_ranges <- 
            IRanges(
                start = 
                    d2$qs2,
                end = 
                    d2$qe2
            ) 
        
        ### deletions in the query are gaps in the reference
        
        deletions <- 
            ref_ranges %>%
            reduce() %>%
            gaps(start = 1, end = d2$rlen[1]) 
        
        deletions_df <- 
            deletions %>%
            as.data.frame() %>% 
            mutate(
                acc = d2$qid[1]#,
                #clust = clust_i
            ) %>%
            mutate(
                variant = "Deletion"
            )
        
        ## Insertions in the reference are gaps in the query 
        
        insertions <- 
            qry_ranges %>%
            reduce() %>%
            gaps(start = 1, end = d2$qlen[1])
        
        insertions_df <- 
            insertions %>%
            as.data.frame() %>% 
            mutate(
                acc = d2$qid[1]#,
                #clust = clust_i
            ) %>%
            mutate(
                variant = "Insertion"
            )
        
        qlength <- 
            d2$qlen[1]
        
        rlength <- 
            d2$rlen[1]
        
        ref_covered <- 
            reduce(ref_ranges) %>% 
            width() %>% sum
        
        ref_uncovered <- 
            rlength-ref_covered
        
        qry_covered <- 
            reduce(qry_ranges) %>% 
            width() %>% sum
        
        qry_uncovered <- 
            qlength-qry_covered
        
        error_perc <- 
            sum(
                d2$error
            )/
            mean(c(rlength, qlength))
        
        length_dif_perc <- 
            abs(rlength-qlength)/
            mean(c(rlength, qlength))
        
        delta_identity <- 
            (mean(c(ref_covered, qry_covered))/
                 mean(c(rlength, qlength))) -
            error_perc - length_dif_perc
        
        ## overlaps 
        qols <-
            findOverlaps(qry_ranges, type = "any")
        qunique <- 
            queryHits(qols) != subjectHits(qols)
        
        qor <- 
            pintersect(
                qry_ranges[queryHits(qols)][qunique], 
                qry_ranges[subjectHits(qols)][qunique]
            ) %>%
            unique
        # 
        # qor_tib <- 
        #     qor %>%
        #     as.data.frame()
       
        ## invs
        
        inv_list <- 
            list()
        if (
            (nrow(chrom_SVs)>0) 
            ){
            for (j in 1:nrow(chrom_SVs)){
                
                if (j>1){
                    dtj <- 
                        delta_table %>%
                        filter(
                            (re <= (chrom_SVs$rs[j-1])) |
                                (rs >= (chrom_SVs$re[j-1]))
                        )
                } else{
                    dtj <- 
                        delta_table
                }
                CSV <- 
                    chrom_SVs[j,]
                
                inv_tib <- 
                    dtj %>% 
                        filter(
                            (rs >= (CSV$rs[1])) & (rs <= (CSV$re[1])),
                            strand!=CSV$strand[1]
                            ) %>%
                    mutate(
                        acc = qid,
                        start = rs, 
                        end = re,
                        width = abs(start-end+1),
                        variant = 
                            ifelse(
                                X_dist < 1e6,
                                "Ori inversion",
                                "local inversion"
                            )
                    ) %>%
                    dplyr::select(
                        c(
                            start, end, width, acc, variant
                        )
                    ) %>%
                    mutate(
                        variant = as.character(variant)
                    )
                    
                inv_list[[j]] <- 
                    inv_tib
            }
            
            invs_outer <- 
                delta_table %>%
                filter(
                    (re <= (chrom_SVs$rs[j])) |
                        (rs >= (chrom_SVs$re[j])),
                    strand=="-"
                ) %>% 
                mutate(
                        acc = qid,
                        start = rs, 
                        end = re,
                        width = abs(start-end+1),
                        variant = 
                            ifelse(
                                X_dist < 1e6,
                                "Ori inversion",
                                "local inversion"
                            )
                    ) %>%
                dplyr::select(
                    c(
                        start, end, width, acc,  variant
                    )
                ) %>%
                mutate(
                    variant = as.character(variant)
                )
            
            inv_list[[j+1]] <- 
                invs_outer
        } else {
            inv_list <- 
                delta_table %>%
                filter(
                    strand=="-"
                ) %>%
                mutate(
                    acc = qid,
                    start = rs, 
                    end = re,
                    width = abs(start-end+1),
                    variant = 
                        ifelse(
                            X_dist < 1e6,
                            "Ori inversion",
                            "local inversion"
                        )
                ) %>%
                dplyr::select(
                    c(
                        start, end, width, acc,  variant
                    )
                ) %>%
                mutate(
                    variant = as.character(variant)
                )
                
        }
        
        invs <- 
            bind_rows(
                inv_list
            ) 
        
        mean_x_dist <- 
            mean(d2$X_dist)
    
        
        tlocs <- 
            d2 %>% 
            filter(
                !is.na(X_dist),
                X_dist > (mean(X_dist) + 1 * sd(X_dist))
            ) %>%
            mutate(
                #clust = clust_i,
                acc = qid
            )%>%
            rowwise() %>%
            mutate(
                start = 
                    rs,
                    #round(mean(c(rs, qe))),
                end = 
                    re,
                    #round(mean(c(re, qs))),
                width =  abs(start-end+1),
                variant = "Translocation"
            )  %>%
            ungroup() %>%
            dplyr::select(
                c(
                    start, end, width, acc,  variant
                    #clust
                    )
                )

        if (
            nrow(invs)==0 &
            nrow(gaps_df)==0 &
            nrow(ins_df)==0 &
            nrow(tlocs)==0 
        ){
            all_events <- 
                tibble(
                    start = NA,
                    end = NA,
                    width = NA,
                    acc = qrys[i],
                    #clust = clust_i,
                    variant = NA
                ) %>%
                mutate(
                    delta_identity = delta_identity
                )
        } else {
            
            all_events <- 
                bind_rows(
                    invs, 
                    deletions_df,
                    insertions_df,
                    tlocs
                ) %>%
                arrange(start) %>%
                mutate(
                    delta_identity = delta_identity
                )
        }
        
        return(all_events)
        
    }



# Phylogenetic networking ####
get_potential_adjacents <- 
    function(
        clust_identity_diff
    ){
        
        clust_identity_diff_adjacent <- 
            clust_identity_diff %>% 
            mutate(adjacency = abs(ref_invs_to_0-qry_invs_to_0)) %>%
            filter(
                adjacency == 1
            )
        
        cida_out <- 
            clust_identity_diff_adjacent %>%
            filter(
                ref_clust !=0,
                qry_clust !=0
            )
        return(cida_out)
    }

mummer_all_v_all <- 
    function(
        alignment_matrix, delta_dir, sync_dir
    ){
        
        nuc_params <-
            "--mum --maxgap=500 --mincluster=100"
        delta_dir_ava <- 
            sprintf(
                "%s/all_v_all",
                delta_dir
            )
        if (
            !dir.exists(delta_dir_ava)
        ){
            dir.create(delta_dir_ava, recursive = T)
        }
        nuc_commands <- 
            generate_nucmer_commands(alignment_matrix, delta_dir_ava, sync_dir)
        pbmclapply(nuc_commands, system, mc.cores = 7)
    }

get_adjaceny <- 
    function(
        adj_delta_dir, cluster_identity_diff
    ){
        
        filtered_combined_delta_adj <- 
            filter_combined_delta(
                adj_delta_dir, 
                minimum_alignment_percentage, 
                minimum_inversion_length
            )
        inversion_counts_adj <- 
            filtered_combined_delta_adj %>%
            summarise(
                inversions = unique(invs)
            ) %>%
            ungroup %>% 
            distinct 
        
        #cluster_identity_diff <- cluster_identity$clust_identity_diff
        zero_adjacent <- 
            cluster_identity_diff %>% 
            ungroup %>%
            filter(
                ref_clust==0 | qry_clust==0,
                ref_invs_to_0==1 | qry_invs_to_0==1
            ) %>% 
            mutate(
                inversions=1,
                ref_clust = as.numeric(as.character(ref_clust)),
                qry_clust = as.numeric(as.character(qry_clust))
            ) %>%
            dplyr::select(
                ref, qry, inversions, ref_clust, qry_clust
            ) %>%
            `colnames<-`(
                c(
                    "rid", "qid", "inversions", "ref_clust", "qry_clust"
                )
            )
        
        adjacent_clusters <-
            bind_rows(
                zero_adjacent,
                inner_join(
                    (inversion_counts_adj %>% filter(inversions==1)),
                    clustered_genomes %>%
                        dplyr::select(cluster_id, accessions) %>%
                        `colnames<-`(c("ref_clust", "rid")),
                    by = "rid"
                ) %>%
                    inner_join(
                        ., 
                        clustered_genomes %>%
                            dplyr::select(cluster_id, accessions) %>%
                            `colnames<-`(c("qry_clust", "qid")),
                        by = "qid"
                    )
            )
        return(adjacent_clusters)
    }

# between cluster variants
between_cluster_variants <- 
    function(
        representative_table, 
        filtered_combined_delta
    ){
        
        events <- 
            list()
        
        for (i in 1:nrow(representative_table)){
            ####
            qi <- 
                pull(representative_table, accessions)[i]
            # the delta file filtered
            fdi <- 
                filtered_combined_delta %>% 
                ungroup %>%
                filter(qid==qi) 
            # plot_delta(fdi)
            fd2 <- 
                fdi %>%
                filter( ## not the first or last alignment, 
                    # rs !=min(rs), re !=max(re), 
                    # qs !=1
                    segment != min(segment),
                    segment != max(segment)
                )%>%
                filter(X_dist<1e6)
            
            if (nrow(fd2)==0){
                events[[i]] <- 
                    fdi %>%
                    filter( ## not the first or last alignment, 
                        # rs !=min(rs), re !=max(re), 
                        # qs !=1
                        segment != min(segment),
                        segment != max(segment)
                    ) %>%
                    rowwise() %>%
                    mutate(
                        rsp = 
                            ((100*rs)/rlen) %>% round(.),
                        rep = 
                            ((100*re)/rlen) %>% round(.),
                        qsp = 
                            ((100*qs)/qlen) %>% round(.),
                        qep = 
                            ((100*qe)/qlen) %>% round(.),
                        rmp = round((rs+re)/2),
                        qmp = round((qs+qe)/2),
                        rmprop = round(100*rmp/rlen),
                        qmprop = round(100*qmp/qlen),
                        meansize = round((100*meanlen)/mean(c(rlen, qlen))),
                        midpoint = mean(rmprop, qmprop),
                        segment = as.character(segment)
                    ) %>%
                    dplyr::select(
                        qid, rsp, rep, qsp, qep, 
                        strand, segment,
                        X_dist, meansize, midpoint
                    ) %>%
                    ungroup()
                
            } else{
                fd3 <- 
                    fdi %>%
                    filter( ## not the first or last alignment, 
                        # rs !=min(rs), re !=max(re), 
                        # qs !=1
                        segment != min(segment),
                        segment != max(segment)
                    )%>%
                    filter(X_dist<1e6) %>%
                    group_by(segment) %>%
                    summarise(
                        strand = unique(strand),
                        rs = min(rs), 
                        re = max(re),
                        qs = ifelse(strand=="+", min(qs), max(qs)),
                        qe = ifelse(strand=="+", max(qe), min(qe)),
                        X_dist = round(mean(X_dist)),
                        rlen = unique(rlen),
                        qlen = unique(qlen),
                        meanlen = round(sum(meanlen)),
                        qid = unique(qid)
                    )
                
                segdiff <- 
                    fd3 %>%
                    mutate(segdiff = segment-lag(segment, default = 1)) %>%
                    pull(segdiff)
                
                if (any(segdiff>1)) next
                
                fd2_props <- 
                    fd3 %>%
                    rowwise() %>%
                    mutate(
                        rsp = 
                            ((100*rs)/rlen) %>% round(.),
                        rep = 
                            ((100*re)/rlen) %>% round(.),
                        qsp = 
                            ((100*qs)/qlen) %>% round(.),
                        qep = 
                            ((100*qe)/qlen) %>% round(.),
                        rmp = round((rs+re)/2),
                        qmp = round((qs+qe)/2),
                        rmprop = round(100*rmp/rlen),
                        qmprop = round(100*qmp/qlen),
                        meansize = round((100*meanlen)/mean(c(rlen, qlen))),
                        midpoint = mean(rmprop, qmprop)
                    ) %>%
                    dplyr::select(
                        qid, rsp, rep, qsp, qep, strand, segment,
                        X_dist, meansize, midpoint
                    ) %>%
                    ungroup() %>%
                    filter(
                        strand=="-" |
                            (
                                (lag(strand, default = "+")=="-") & 
                                    (lead(strand, default = "+")=="-")
                            )
                    )
                
                inner_seg <- 
                    median(unique(fd2_props$segment))
                
                if (nrow(fd2_props)==1){
                    middle_seg_row <- 
                        fd2_props %>%
                        filter(
                            segment==median(segment)
                        ) %>%
                        mutate(
                            ms =
                                ifelse(
                                    strand=="+",
                                    mean(c(rsp, qsp)),
                                    mean(c(rsp, qep))
                                ),
                            me = 
                                ifelse(
                                    strand=="+",
                                    mean(c(rep, qep)),
                                    mean(c(rep, qsp))
                                ),
                            segment = as.character(segment)
                        ) %>%
                        dplyr::select(
                            qid, strand, segment, 
                            X_dist, meansize, midpoint,
                            ms, me
                        )
                    events[[i]] <- 
                        middle_seg_row
                } else {
                    
                    middle_seg_row <- 
                        fd2_props %>%
                        filter(
                            segment==median(segment)
                        ) %>%
                        mutate(
                            ms =
                                ifelse(
                                    strand=="+",
                                    mean(c(rsp, qsp)),
                                    mean(c(rsp, qep))
                                ),
                            me = 
                                ifelse(
                                    strand=="+",
                                    mean(c(rep, qep)),
                                    mean(c(rep, qsp))
                                ),
                            segment = as.character(segment)
                        ) %>%
                        dplyr::select(
                            qid, strand, segment, 
                            X_dist, meansize, midpoint, 
                            ms, me
                        )
                    outer_segs <- 
                        fd2_props %>%
                        filter(segment != median(segment)) 
                    
                    possible_outers <- 
                        nrow(outer_segs) %/% 2
                    outer_events <- 
                        list()
                    for (j in 1:possible_outers){
                        outer_events[[j]] <- 
                            fd2_props %>% 
                            filter(segment %in% c(inner_seg-j, inner_seg+j)) %>%
                            summarise(
                                qid = unique(qid),
                                strand= unique(strand),
                                segment = paste(segment, collapse = ","),
                                X_dist = mean(X_dist),
                                ms =
                                    ifelse(
                                        strand=="+",
                                        mean(c(min(rsp), min(qsp))),
                                        mean(c(min(rsp), min(qep)))
                                    ),
                                me =
                                    ifelse(
                                        strand=="+",
                                        mean(c(max(rep), max(qep))),
                                        mean(c(max(rep), max(qsp)))
                                    ),
                                meansize = 
                                    me-ms,
                                midpoint = 
                                    (ms+me)/2
                            )
                    }
                    inv_summary <- 
                        bind_rows(
                            middle_seg_row,
                            bind_rows(outer_events)
                        )
                    events[[i]] <- 
                        inv_summary
                }
                
            }
            #plot_delta(fdi)
            #print(i)
        }
        
        evd <- 
            bind_rows(events) %>% 
            arrange(meansize)
        
        ### all inversions included
        evd2 <- 
            left_join(
                evd,
                (representative_table %>% dplyr::rename(qid = accessions)),
                by = "qid"
            ) 
        return(evd2)
        
    }


# additional functions
simple_nucmer <- 
    function(
        ref, qry, output, nuc_remove = F
    ){
        nuc_params <-
            "--mum --maxgap=500 --mincluster=100"
        nucmer <-
            sprintf(
                "nucmer %s --prefix=%s",
                nuc_params,
                output
            )
        nucmer_command <- 
            sprintf(
                "%s %s %s",
                nucmer,
                ref, 
                qry
            )
        filter_out <- 
            sprintf(
                "%s_filtered.delta",
                output
            )
        # if (nuc_check){
        #     nuc_check <- 
        #         sprintf(
        #             "if [ ! -e %s ] then ",
        #             filter_out
        #         )
        # } else {
        #     nuc_check <- NULL
        # }
        
        
        nucmer_filter <-
            sprintf(
                "delta-filter -r -q %s.delta > %s",
                output,
                filter_out
            )
        ## to save space, remove the original delta file
        if (nuc_remove){
            nuc_remove <- 
                sprintf(
                    "rm %s/%s_v_%s.delta",
                    delta_dir,
                    x[1], x[2]
                )
        } else {
            unfiltered_dir <- 
                sprintf(
                    "%s/unfiltered",
                    delta_dir
                )
            if (!dir.exists(unfiltered_dir)){
                dir.create(unfiltered_dir)
            }
            
            nuc_remove <- 
                sprintf(
                    "mv %s/%s_v_%s.delta %s/%s_v_%s.delta",
                    delta_dir,
                    x[1], x[2],
                    unfiltered_dir,
                    x[1], x[2]
                )
        }
        
        return(
            capture.output(
                cat(
                    #nuc_check,
                    nucmer_command, 
                    nucmer_filter, 
                    nuc_remove,
                    sep = "; "
                )
            )
        )
    }


file_exists_at_url <- 
    function(url) {
        con <- tryCatch({
            suppressWarnings({
                con <- url(url, "r")  # Try to open the URL
                close(con)            # Close the connection if it was successful
                TRUE                  # Return TRUE if the file exists
            })
        }, error = function(e) {
            FALSE  # Return FALSE if there's an error (file doesn't exist)
        })
        return(con)
    }


# end ####

