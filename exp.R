## expediency



saved_results <- 
    list()

ssn <- 
    list()

length_summary_tables <- 
    list()


species_name <- 
    top_species$Species_[5]
#(species_list %>% pull(species2))[i]

#set up directories
genus_name <- 
    strsplit(species_name, "_")[[1]][1]

clade_name <- 
    filter(top_species, Species_==species_name) %>%
    pull(Phylum) %>% gsub(" ", "_", .)

species_variable_names <-
    c(
        "fna_dir",
        "sync_dir", "gff_dir", "identity_dir", "mummer_dir",
        "oriC_dir",
        "delta_dir", "r_vars_dir", "results_dir", 
        "repeats_dir", "invs_dir",
        "disparities_dir",# "c01_dir",
        "subset_dir",
        "breakpoint_dir", "prokka_dir",
        "breakpoint_comparisons_dir", "snp_dir",
        "kmer_dir", "bd_dir", "plot_dir", "roary_dir",
        "sync2_dir", "test_dir", "IS_dir"
    )

species_dirs <- 
    gsub("_dir", "", species_variable_names)

for (var_name in species_variable_names){
    dir_name <- 
        gsub("_dir", "", var_name)
    full_dir_name <- 
        sprintf(
            "../data/processing/%s/%s/%s/%s", 
            dir_name, clade_name, 
            genus_name, species_name
        )
    
    assign(
        var_name, 
        full_dir_name
    )
    if (!dir.exists(full_dir_name)){dir.create(full_dir_name, recursive = T)}
}

species_table <- 
    ncbi_table %>%
    filter(species== sub("_", " ", species_name))

spec_count <- 
    nrow(species_table)

if(spec_count<250){
    subsample_idx <- 
        1:spec_count
} else {
    set.seed(123)
    subsample_idx <- 
        sample(1:spec_count, 200)
}

st_ss <- 
    species_table %>%
    dplyr::slice(subsample_idx)

the_rest_idx <- 
    species_table %>%
    filter(!row_number() %in% subsample_idx)

species_summary <- #species_table
    get_species_summary(st_ss) 

ssn[[i]] <- 
    tibble(
        spec_ = species_name,
        stn = nrow(species_table),
        ssn = nrow(species_summary)
    )
saveRDS(
    species_summary, 
    sprintf(
        "%s/species_summary.rds",
        r_vars_dir
    )
)



ASV %>%
    count(spec, variant, variant_specific) %>% as.data.frame()



### Finding OriC and rearranging

ftp_paths <- 
    species_summary$ftp_path

gff_paths <- 
    ftp_paths %>%
    sprintf(
        "%s/%s_genomic.gff.gz",
        ., gsub(".*/", "", .)
    )

sync_tib <- 
    pbmclapply(
        1:nrow(species_summary),
        function(ii){
            #for (ii in nrow(species_summary):1){
            assembly_path <- 
                species_summary$ftp_path[ii]
            
            #species_summary$accession[ii]
            sync_summary <- 
                sync_dnaA(
                    assembly_path = 
                        assembly_path,
                    sync_dir = 
                        sync_dir, 
                    gff_dir = 
                        gff_dir
                ) %>%
                mutate(
                    idx = ii
                )
            return(sync_summary)
        }, mc.cores = 4
    ) %>% bind_rows()

st2 <- 
    sync_tib %>%
    filter(is.na(note))


if (nrow(sync_tib)==nrow(st2)){
    "All genomes rearranged"
} else{
    sprintf(
        "%i/200 sequences had issues",
        sum(!is.na(sync_tib$note))
    )
}


### identities

saveRDS(
    sync_tib,
    sprintf(
        "%s/sync_tib.rds",
        results_dir
    )
)

sync_tib <-
    readRDS(
        sprintf(
            "%s/sync_tib.rds",
            results_dir
            )
        )


#species_summary %>% filter(accession =="NZ_CP062406.1")

sdc <- 
    sprintf(
        "%s/ss_final",
        sync_dir
    )
if(!dir.exists(sdc)){dir.create(sdc)}

lapply(
    sprintf(
        "%s/%s.txt",
        sync_dir,
        st2$accession
        ),
    function(x){
        file.copy(
            from = x, 
            to= 
                sprintf(
                    "%s/%s",
                    sdc, basename(x)
                )
        )
    }
    
    
)



# species_identities <- 
#     readRDS(
#         sprintf(
#             "%s/species_identities.rds",
#             results_dir
#         )
#     )

# files <- list.files(path = sdc, pattern = "\\.txt\\.txt$", full.names = TRUE)
# file.rename(files, sub("\\.txt$", "", files))


#}
list.files(sdc)
## Identities
species_identities <- 
    get_pairwise_identities(
        sync_dir = sdc, 
        identity_dir = identity_dir, 
        species_table = st2, 
        genus_name, species_name
    )


species_identities <- 
    sprintf(
        "%s/species_identities.rds",
        results_dir
    ) %>% readRDS()

#species_identities$mash_plot
saveRDS(
    species_identities, 
    sprintf(
        "%s/species_identities.rds",
        results_dir
    )
)

sampled_id_table <- 
    species_identities$id_dists %>%
    arrange(desc(Sequence_mean_distance))


accessions <- 
    sync_tib %>%
    filter(is.na(note)) %>%
    pull(accession)

most_related <- 
    sampled_id_table %>%
    arrange(Sequence_mean_distance) %>%
    filter(Sequence %in% accessions) %>%
    dplyr::slice(1) %>% 
    pull(Sequence)

top_3 <- 
    sampled_id_table %>%
    arrange(Sequence_mean_distance) %>%
    filter(Sequence %in% accessions) %>%
    dplyr::slice(1:3) %>% 
    pull(Sequence)


###nucmer

nucmer_commands <- 
    nucmer_ref_v_all(
        most_related, 
        accessions,
        delta_dir =  delta_dir, 
        sync_dir = sdc
    )

minimum_alignment_percentage <- 
    90
minimum_inversion_length <- 
    5e4
filtered_combined_delta <- 
    filter_combined_delta(
        sprintf("%s/ref_v_all", delta_dir), 
        minimum_alignment_percentage, 
        minimum_inversion_length
    ) 


fcd <- 
    sprintf(
        "%s/fcd.rds",
        results_dir
    ) %>% readRDS()

fcd <- 
    filtered_combined_delta %>%
    filter(!(qid %in% likely_missassembled)) %>%
    filter(qid %in% accessions)

saveRDS(
    fcd, 
    sprintf(
        "%s/fcd.rds",
        results_dir
    )
)



clustered_genomes_list <- 
    cluster_genomes(fcd)

saveRDS(
    clustered_genomes_list, 
    sprintf(
        "%s/clustered_genomes_list.rds",
        results_dir
    )
)


delta_identities <- 
    pbmclapply(
        1:length((list.files(sprintf("%s/ref_v_all", delta_dir)))),
        function(i){
            #for (i in 1:length((list.files(sprintf("%s/ref_v_all", delta_dir))))){
            ff <-     
                list.files(
                    sprintf("%s/ref_v_all", delta_dir), 
                    full.names = T
                )[i]
            if (file.size(ff)!=0 & grepl("filtered", ff)) {
                delta_tib <- 
                    read_delta(
                        list.files(
                            sprintf("%s/ref_v_all", delta_dir), 
                            full.names = T
                        )[i]
                    ) %>%
                    rowwise() %>%
                    mutate(
                        qs2 = ifelse(qe<qs, qe, qs),
                        qe2 = ifelse(qe<qs, qs, qe),
                    ) %>% ungroup
                ref_ranges <- 
                    IRanges(
                        start = delta_tib$rs,
                        end = delta_tib$re
                    ) %>% reduce
                
                length_dif_prop <- 
                    abs(
                        delta_tib$rlen[1] -
                            delta_tib$qlen[1]
                    )/
                    (max(
                        c(
                            delta_tib$rlen[1], 
                            delta_tib$qlen[1])
                    ))
                
                nucmer_id <- 
                    sum(
                        delta_tib$meanlen
                    )/
                    (max(
                        c(
                            delta_tib$rlen[1], 
                            delta_tib$qlen[1])
                    ))
                
                out_tib <- 
                    tibble(
                        accession = 
                            delta_tib$qid[1],
                        identity = 
                            nucmer_id-length_dif_prop
                    )
            } else {
                out_tib <- 
                    tibble(
                        accession = 
                            basename(ff),
                        identity = NA
                    )
            }
            
            
            
            return(out_tib)
            
        }, mc.cores = 4
    ) %>% bind_rows() %>%
    filter(accession %in% accessions) 



# all_delta_tibs <-
#     pbmclapply(
#         1:length((list.files(sprintf("%s/ref_v_all", delta_dir)))),
#         function(i){
#             delta_tib <- 
#                 read_delta(
#                     list.files(
#                         sprintf("%s/ref_v_all", delta_dir), 
#                         full.names = T
#                     )[i]
#                 ) %>%
#                 rowwise() %>%
#                 mutate(
#                     qs2 = ifelse(qe<qs, qe, qs),
#                     qe2 = ifelse(qe<qs, qs, qe),
#                 ) %>% ungroup
#         }
#     )
# 
# saveRDS(
#     all_delta_tibs, 
#     sprintf(
#         "%s/all_delta_tibs.rds",
#         r_vars_dir
#     )
# )


dt2 <- 
    right_join(
        delta_identities,
        species_identities$id_dists[c(1,2)] %>%
            `colnames<-`(c("accession", "MASH Identity"))
    )

delta_identities$accession

clustered_genomes <- 
    clustered_genomes_list$cluster_summary %>%
    left_join(
        ., 
        st2[,c("accession", "len")] %>% 
            `colnames<-`(c("accessions", "len")),
        by = "accessions"
    ) 




## All structural variants versus identity

qrys <-
    fcd$qid %>% unique
accessions[!(accessions==most_related)]

best_identity <- 
    inner_join(
        #delta_identities %>% `colnames<-`(c("accessions", "identity")),
        delta_identities %>% `colnames<-`(c("accessions", "identity")),
        clustered_genomes[,-5],
        by = "accessions"
    ) %>%
    # mutate(
    #     identity = 
    #         ifelse(identity>1,  1, identity),
    #     identity = 
    #         ifelse(accessions==most_related,  1.01, identity),
    # ) %>%
    arrange(cluster_id, desc(identity))


bis <- 
    best_identity %>% 
    group_by(cluster_id) %>%
    arrange(desc(identity)) %>% 
    filter(!is.na(cluster_id)) %>%
    dplyr::slice(1)



### ref-qry matrix and plots




#reference accessions
set.seed(123)
sampled_refs <- 
    c(
        best_identity %>% 
            filter(
                grouped_sequences==max(grouped_sequences)
            ) %>% 
            sample_n(3) %>%
            pull(accessions), ### 3 from the main cluster
        best_identity %>%
            filter(
                grouped_sequences!=max(grouped_sequences)
            ) %>% 
            sample_n(2) %>%
            pull(accessions)
    )


best_identity %>% group_by(cluster_id) %>% summarise(accessions = n())


clustered_genomes %>% count(cluster_id)    
#clustered_genomes %>% filter(cluster_id %in% 0:5)
### qry sampled
sampled_qrys <- 
    clustered_genomes[,c(1,2)] %>% 
    filter(
        !(accessions %in% sampled_refs),
        cluster_id %in% 0:5
    ) %>%
    group_by(cluster_id) %>% 
    dplyr::slice(1:2) %>%
    ungroup %>%
    sample_n(5) %>%
    pull(accessions)

##
ref_qry_pairs <- 
    expand_grid(ref = sampled_refs, qry = sampled_qrys)


ref_qry_pairs
#expand_grid(ref = c(1:5), qry = c(1:5))


delta_dir2 <- 
    sprintf(
        "%s/dd2",
        delta_dir
    )

if (!dir.exists(delta_dir2)) {dir.create(delta_dir2)}
nuc_comms <- 
    generate_nucmer_commands(
        ref_qry_pairs, delta_dir = delta_dir2, sync_dir = sync_dir
    )


delta_dir3 <- 
    sprintf(
        "%s/dd3",
        delta_dir
    )

if (!dir.exists(delta_dir3)) {dir.create(delta_dir3)}
nuc_comms <- 
    generate_nucmer_commands(
        ref_qry_pairs, delta_dir = delta_dir3, sync_dir = sync_dir, 
        nuc_remove = F
    )

pbmclapply(
    nuc_comms,
    system
)


bis_1 <- 
    bis[-1,]

if (nrow(bis_1)==0) next
rq_plots <- 
    lapply(
        1:12,#length(list.files(delta_dir2)),
        function(j){
            dfj <- 
                list.files(delta_dir3, full.names = T)[j]
            # if (j==5){
            #     pj <- 
            #         plot_delta(dfj)  +
            #         xlab("Cluster 0")
            # } else {
            pj <- 
                plot_delta(dfj)  +
                xlab("")
            #}
            pjo <- 
                pj +
                theme(
                    axis.text = element_text(size = 18),
                    axis.title = element_text(size = 22)
                )
            return(pjo)
        }
    )



do.call(
    "grid.arrange",
    c(
        rq_plots,
        ncol=4
    )
)



#reference accessions
set.seed(123)
sampled_refs <- 
    c(
        best_identity %>% 
            filter(
                grouped_sequences==max(grouped_sequences)
            ) %>% 
            sample_n(3) %>%
            pull(accessions), ### 3 from the main cluster
        best_identity %>%
            filter(
                grouped_sequences!=max(grouped_sequences)
            ) %>% 
            sample_n(2) %>%
            pull(accessions)
    )


best_identity %>% group_by(cluster_id) %>% summarise(accessions = n())


clustered_genomes %>% count(cluster_id)    
#clustered_genomes %>% filter(cluster_id %in% 0:5)
### qry sampled
sampled_qrys <- 
    clustered_genomes[,c(1,2)] %>% 
    filter(
        !(accessions %in% sampled_refs),
        cluster_id %in% 0:5
    ) %>%
    group_by(cluster_id) %>% 
    dplyr::slice(1:2) %>%
    ungroup %>%
    sample_n(5) %>%
    pull(accessions)

##
ref_qry_pairs <- 
    expand_grid(ref = sampled_refs, qry = sampled_qrys)


ref_qry_pairs
#expand_grid(ref = c(1:5), qry = c(1:5))


delta_dir2 <- 
    sprintf(
        "%s/dd2",
        delta_dir
    )

if (!dir.exists(delta_dir2)) {dir.create(delta_dir2)}
nuc_comms <- 
    generate_nucmer_commands(
        ref_qry_pairs, delta_dir = delta_dir2, sync_dir = sync_dir
    )


delta_dir3 <- 
    sprintf(
        "%s/dd2",
        delta_dir
    )

if (!dir.exists(delta_dir3)) {dir.create(delta_dir3)}
nuc_comms <- 
    generate_nucmer_commands(
        ref_qry_pairs, delta_dir = delta_dir2, sync_dir = sync_dir, 
        nuc_remove = F
    )

pbmclapply(
    nuc_comms,
    system
)


bis_1 <- 
    bis[-1,]

if (nrow(bis_1)==0) next
rq_plots <- 
    lapply(
        1:21,#length(list.files(delta_dir2)),
        function(j){
            dfj <- 
                list.files(delta_dir2, full.names = T)[j]
            # if (j==5){
            #     pj <- 
            #         plot_delta(dfj)  +
            #         xlab("Cluster 0")
            # } else {
            pj <- 
                plot_delta(dfj)  +
                xlab("")
            #}
            pjo <- 
                pj +
                theme(
                    axis.text = element_text(size = 18),
                    axis.title = element_text(size = 22)
                )
            return(pjo)
        }
    )



do.call(
    "grid.arrange",
    c(
        rq_plots,
        ncol=7
    )
)


### the prokka

prokka_accs <- 
    c(sampled_refs, sampled_qrys)

prokka_dir2 <- 
    "../data/prokka2"
pd3 <- 
    sprintf(
        "../data/prokka3/%s",
        species_name
    )

pd4 <- 
    sprintf(
        "../data/prokka4/%s",
        species_name
    )

if (!dir.exists(pd3)){dir.create(pd3, recursive = T)}
if (!dir.exists(pd4)){dir.create(pd4, recursive = T)}
prok_com <- 
    annotate_fna2(
        species_name = 
            species_name,
        sync_dir = 
            sync_dir, 
        prokka_dir = 
            pd3, 
        accessions = 
            prokka_accs
    )



pbmclapply(
    1:length(prok_com$comms),
    function(i){
        system(
            prok_com$comms[i]
            )
    }, mc.cores = 1
)


list.files(pd3)
list.files(pd4)

lapply(
    1:length(prokka_accs),
    function(i){
        file.copy(
            from = 
                list.files(
                    pd3,
                    pattern = sprintf("%s.gff", prokka_accs[i]), 
                    full.names = T
                ),
            to=
                sprintf(
                    "%s/%s.gff",
                    pd4, prokka_accs[i]
                )
            
        )
    }
    
)


prokka_list <- 
    pbmclapply(
        1:length(prokka_accs),
        function(i){
        #for (i in 1:length(prokka_accs)){
            gff_i <- 
                read.delim(
                    sprintf(
                        '%s/%s.gff', pd3, prokka_accs[i]
                    ), header=F, comment.char = "#")  %>%
                as_tibble() %>%
                `colnames<-`(
                    c(
                        "accession", "source", "type", "start", 
                        "end", "score", "strand", "phase", "attributes")
                ) %>%
                filter(!is.na(start)) %>%
                mutate(
                    width = end-start+1,
                    gene = str_extract(attributes, "(?<=gene=)[^;]+"),
                    product = str_extract(attributes, "(?<=product=)[^;]+")
                ) %>%
                dplyr::select(-c(attributes, source))
            return(gff_i)
        }
    )


prokka_list %>%
    lapply(., head)
all_prokka <- 
    bind_rows(prokka_list)

all_prokka <- 
    sprintf(
        "results/%s_all_prokka.rds",
        species_name
    ) %>% readRDS()

saveRDS(
    all_prokka,
    sprintf(
        "results/%s_all_prokka.rds",
        species_name
    )
)



clustered_genomes <- 
    clustered_genomes %>%
    rowwise() %>%
    mutate(
        len = 
            ifelse(
                is.na(len),
                width(
                    readDNAStringSet(
                        sprintf("%s/%s.txt", sync_dir, accessions)
                    )
                ),
                len
            )
    )

gene_regions <- 
    lapply(
        1:length(prokka_list),
        function(i){
            #for (i in 1:length(prokka_list)){
            gff_i <- 
                prokka_list[[i]]
            
            gff_ranges <- 
                IRanges(
                    start = 
                        gff_i$start,
                    end = 
                        gff_i$end
                )
            genes_merged <- 
                reduce(gff_ranges)
            
            gene_coverage <- 
                sum(width(genes_merged))
            
            ref_len <- 
                clustered_genomes %>%
                distinct(accessions, .keep_all = T) %>%
                filter(accessions ==gff_i$accession[1]) %>%
                pull(len)
            
            intergenic_coverage <- 
                ref_len - gene_coverage
            
            genome_coverage_df <- 
                tibble(
                    region = c("Genic", "Intergenic"),
                    bp = c(gene_coverage, intergenic_coverage)
                ) %>%
                mutate(
                    proportion = bp / sum(bp),
                    acc= 
                        gff_i %>% pull(accession) %>% unique
                )
            return(genome_coverage_df)
            
        }
    )

gene_regions_tib <- 
    gene_regions %>% bind_rows


saveRDS(
    gene_regions_tib,
    sprintf(
        "%s/gene_regions_tib.rds",
        r_vars_dir
    )
)

df3 <- 
    list.files(delta_dir3, full.names = T)


df_tib <- 
    ref_qry_pairs %>%
    mutate(
        expected_basename = str_c(ref, "_v_", qry, "_filtered.delta")
    ) %>% 
    inner_join(
        ., 
        tibble(
            expected_basename = basename(df3),
            fullname = df3
        ), 
        by = "expected_basename"
    ) %>%
    mutate(
        unfiltered_fullname = 
            sprintf(
                "%s/unfiltered/%s.delta",
                dirname(df3[1]),
                str_c(ref, "_v_", qry)
            )
    )

df_inner <- 
    df_tib$fullname

row_i <- 
    df_tib[i,]
delta_i <- 
    row_i %>%
    pull(fullname)
delta_i %>% read_delta()
## SV calling

all_indels <- 
    pbmclapply(
        1:nrow(df_tib),
        function(i){
            #for (i in 1:nrow(df_tib)){
            row_i <- 
                df_tib[i,]
            delta_i <- 
                row_i %>%
                pull(fullname)
            delta_indels_i <- 
                delta_indels(delta_i) %>%
                mutate(idx = i)
            
            if (nrow(delta_indels_i)==0){
                delta_indels_i <- 
                    tibble(
                        rid = df_tib$ref[i],
                        qid = df_tib$qry[i],
                        idx = i
                    )
                
            }
            
            #print(i)}
            return(delta_indels_i)
        }
    ) %>% 
    bind_rows %>%
    mutate(
        bin = 
            cut(
                width, 
                breaks = 
                    c(
                        50, 500, 1000, 
                        2000, 10000, Inf
                    ), 
                right = FALSE),
        species = species_name
    ) %>%
    rowwise() %>%
    mutate(rq = paste(rid, qid, sep ="_")) %>%
    group_by(bin, rq) %>%
    mutate(
        SV_count = n()
    )

all_indels %>% filter(is.na(bin))

#all_indels %>% ungroup %>% count(idx) %>% pull(n) %>% boxplot
#all_indels %>% pull((idx)) %>% unique

saveRDS(
    all_indels,
    sprintf(
        "results/%s_all_indels.rds",
        species_name
    )
)

all_dups <- 
    pbmclapply(
        1:nrow(df_tib),
        function(i){
            #for (i in 1:nrow(df_tib)){
            row_i <- 
                df_tib[i,]
            delta_i <- 
                row_i %>%
                pull(fullname)
            
            ufd_i <- 
                row_i %>%
                pull(unfiltered_fullname)
            
            delta_dups_i <- 
                delta_duplications(
                    delta_table = delta_i, 
                    unfiltered_delta_table = ufd_i
                ) 
            
            dd_check <- 
                ncol(delta_dups_i)
            
            delta_dups_i <- 
                delta_dups_i %>%
                mutate(
                    nr = 
                        ifelse(
                            dd_check==2,
                            0,
                            n()
                        ),
                    idx = i
                )
            
            if (nrow(delta_dups_i)==0){
                delta_dups_i <- 
                    tibble(
                        rid = df_tib$ref[i],
                        qid = df_tib$qry[i],
                        idx = i
                    )
                
               # print(i)}
                
            }
            return(delta_dups_i)
        }, mc.cores = 4
    ) %>% 
    bind_rows %>%
    mutate(
        bin = 
            cut(
                dup_len, 
                breaks = 
                    c(
                        50, 500, 1000, 
                        2000, 10000, Inf
                    ), 
                right = FALSE),
        species = species_name
    ) 
# group_by(bin, rq) %>%
# mutate(
#     SV_count = n()
# ) 


all_dups %>% 
    dplyr::select(c(domains, duplication_type, dup_len, dup_gap))
all_dups %>% ungroup %>% count(idx) %>% pull(n) %>% boxplot

all_dups$dup_len %>%summary

boxplot(
    all_indels %>% ungroup %>% count(idx, sort = T) %>% pull(n),
    all_dups %>% ungroup %>% count(idx) %>% pull(n),
    all_tlocs %>% ungroup %>% count(idx) %>% pull(n)    
)





all_13s <- 
    bind_rows(
        all_dups %>%
            transmute(
                s= dup_s, e = dup_e, 
                rid = ref_name,qid= qry_name, len = dup_len,
                variant = "dup"
            ),
        all_tlocs %>% ungroup %>%
            transmute(
                s= start, e = end, rid = rid,qid= qid, len = width,
                variant = "tloc"
            ),
        substruct_invs %>% ungroup %>%
            transmute(
                s= rs, e = re, rid = rid,qid= qid, len = meanlen,
                variant = "inv"
            )
    ) %>%
    filter(!is.na(s))%>%
    mutate(
        idx = row_number()
    )

saveRDS(
    all_dups,
    sprintf(
        "results/%s_all_dups.rds",
        species_name
    )
)


all_dups %>% count(duplication_type)


all_dups %>% ungroup %>% count(bin)
all_dups %>% ungroup %>%
    count(replichores)


all_tlocs <- 
    pbmclapply(
        1:nrow(df_tib),
        function(i){
            row_i <- 
                df_tib[i,]
            delta_i <- 
                row_i %>%
                pull(fullname)
            delta_tlocs_i <- 
                delta_translocations(delta_i) %>%
                mutate(idx = i)
            
            
            
            return(delta_tlocs_i)
        }
    ) %>% 
    bind_rows%>% 
    mutate(
        bin = 
            cut(
                width, 
                breaks = 
                    c(
                        50, 500, 1000, 
                        2000, 10000, Inf
                    ), 
                right = FALSE),
        species = species_name
    ) %>%
    rowwise() %>%
    mutate(rq = paste(rid, qid, sep ="_")) %>%
    group_by(bin, rq) %>%
    mutate(
        SV_count = n()
    ) %>%
    mutate(
        replichores = as.character(replichores),
        replichores = sapply(
            strsplit(replichores, "-"),
            function(x) paste(sort(x), collapse = "-")
        ),
        replichores = factor(replichores)           
    ) %>% ungroup()


saveRDS(
    all_tlocs,
    sprintf(
        "results/%s_all_tlocs.rds",
        species_name
    )
)


boxplot(
    #all_indels %>% ungroup %>% pull(width),
    all_dups %>% ungroup %>% filter(dup_len<5e3) %>% pull(dup_len),
    all_tlocs %>% ungroup %>% filter(width<5e3) %>% pull(width),
    substruct_invs%>% ungroup %>% filter(width<5e3) %>% pull(width)
)


all_tlocs %>% filter(width <1e4) %>% pull(width) %>% density %>% plot

top_species

### inversions
all_SRs <- 
    pbmclapply(
        1:nrow(df_tib),
        function(i){
            #for (i in 1:nrow(df_tib)){
            row_i <- 
                df_tib[i,]
            delta_i <- 
                row_i %>%
                pull(fullname)
            delta_SR_i <- 
                delta_structural_inversions(delta_i) %>%
                mutate(idx = i)
            # delta_SR <- 
            #     delta_SR_i
            delta_subSR_i <- 
                delta_substructural_inversions(delta_i, delta_SR_i)
            
            
            ssr_check <- 
                nrow(delta_subSR_i)
            
            delta_subSR_i <- 
                delta_subSR_i %>%
                mutate(
                    nr = 
                        ifelse(
                            ssr_check==2,
                            0,
                            n()
                        ),
                    idx = i
                )
            
            #plot_delta(delta_i)
            
            both_types <- 
                list(
                    "structural_invs" = 
                        delta_SR_i,
                    "substructural_invs" =
                        delta_subSR_i
                )
            #print(i)}
            
            return(both_types)
        }
    )


struct_invs <- 
    lapply(
        all_SRs,
        function(x){
            x[[1]]
        }
    ) %>% 
    bind_rows()


substruct_invs %>% filter(width <1e4) %>% pull(width) %>% density %>% plot

substruct_invs$meanlen[which.max(density(substruct_invs$meanlen)$y)]

substruct_invs <- 
    lapply(
        all_SRs,
        function(x){
            x[[2]]
        }
    ) %>% 
    bind_rows() %>%
    mutate(
        width= meanlen,
        bin = 
            cut(
                width, 
                breaks = 
                    c(
                        50, 500, 1000, 
                        2000, 10000, Inf
                    ), 
                right = FALSE),
        species = species_name
    ) %>%
    rowwise() %>%
    mutate(rq = paste(rid, qid, sep ="_")) %>%
    group_by(bin, rq) %>%
    mutate(
        SV_count = n()
    )
substruct_invs$meanlen[which.max(density(substruct_invs$meanlen)$y)]

saveRDS(
    all_SRs,
    sprintf(
        "results/%s_all_SRs.rds",
        species_name
    )
)

list.files("results")

small_invs <- 
    substruct_invs %>% 
    ungroup %>%
    group_by(rq) %>%
    count(variant)

substruct_invs %>% 
    ungroup %>%
    group_by(bin) %>%
    count(variant)




all_indels

# install.packages("pracma")
# library(pracma)
# x <- density(all_indels$width)
# x_peaks <- findpeaks(x$y)
# x$x[x_peaks[,2]]

all_SVs <- 
    bind_rows(
        all_indels %>%
            ungroup %>%
            transmute(
                rid,qid, width, 
                start, end,
                refmid = round((start+end)/2),
                #indel_type = variant, 
                bin, 
                variant_specific = variant,
                variant = "Indel",
                species, idx
            ),
        all_dups %>%
            transmute(
                rid = ref_name, 
                qid = qry_name,
                width = dup_len,
                start = dup_s, end = dup_e,
                refmid = round((start+end)/2),
                bin, 
                variant_specific = 
                    ifelse(
                        duplication_type =="tandem",
                        duplication_type,
                        paste(
                            duplication_type, domains, sep = " "
                        )
                    ),
                variant = "Duplication",
                species , idx
            ),
        all_tlocs %>%
            ungroup %>%
            transmute(
                rid, qid, 
                width,
                start, end,
                refmid,
                bin, variant,
                variant_specific =
                    replichores,
                species, idx
            ),
        substruct_invs %>%
            ungroup() %>%
            transmute(
                rid, qid, 
                width = meanlen,
                start = rs, end = re, 
                refmid = round((rs+re)/2),
                bin, 
                variant_specific = 
                    variant,
                variant = "Inversion",
                species, idx
            )
    ) 


saveRDS(
    all_SVs,
    sprintf(
        "results/%s_all_SVs.rds",
        species_name
    )
)







occ <- 
    all_SVs %>%
    mutate(
        #v2 = if_else(grepl("inversion", variant), "Inversion", variant),
        # reverse the stacking order so Indel is on top
        v2 = factor(variant, levels = c("Translocation",
                                   "Inversion",
                                   "Duplication",
                                   "Indel")),
        bin = factor(bin,
                     levels = c("[50,500)",
                                "[500,1e+03)",
                                "[1e+03,2e+03)",
                                "[2e+03,1e+04)",
                                "[1e+04,Inf)"))
    ) %>%
    count(bin, v2, name = "Occurrence")

totals <- occ %>%
    group_by(bin) %>%
    summarise(total = sum(Occurrence), .groups = "drop")

len_plot <-
    ggplot() +
        # stacked segments with grey85 outlines BETWEEN them
        geom_col(
            data    = occ,
            aes(x = bin, y = Occurrence, fill = v2),
            width   = 0.7,
            color   = "grey65",
            size    = 0.3
        ) +
        # outer black outline
        geom_col(
            data  = totals,
            aes(x = bin, y = total),
            fill  = NA,
            color = "black",
            size  = 0.6,
            width = 0.7
        ) +
        # total count labels
        geom_text(
            data        = totals,
            aes(x = bin, y = total, label = total),
            vjust       = -0.5,
            size        = 4
        ) +
        # new palette & no legend title
        scale_fill_viridis_d(
            option = "C",
            begin  = 0.2,
            end    = 0.8,
            direction = -1,
            name = NULL
        ) +
        scale_x_discrete(
            labels = c(
                "[50,500)"     = "50–500",
                "[500,1e+03)"  = "500–1,000",
                "[1e+03,2e+03)"= "1,000–2,000",
                "[2e+03,1e+04)"= "2,000–10,000",
                "[1e+04,Inf)"  = ">10,000"
            )
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(
            x = "\n SV length range (bp)",
            y = "SV Occurrance"
        ) +
        theme_classic(base_size =18) +
        theme(
            legend.position = "top"
        )
ggsave(
    sprintf(
        "results/%s_all_SVs.png",
        species_name
    ), plot = len_plot, width = 8, height = 6
    )

len_plot

spec_proks <- 
    list.files(
        "results", pattern = "all_prokka.rds", full.names = T
        )


all_proks <- 
    lapply(
        spec_proks,
        function(x){
            readRDS(x) %>%
                mutate(
                    species = 
                        basename(x) %>%
                        sub("_all_prokka.rds", "", .)
                )
        }
    ) %>% 
    bind_rows() %>%
    filter(!is.na(start), type =="CDS")



all_proks %>% count(type)

saveRDS(
    all_prokka,
    sprintf(
        "results/%s_all_prokka.rds",
        species_name
    )
)


list.files("results")


readRDS("results/Bordetella_pertussis_all_SRs.rds")

asv_specs <- 
    list.files("results", pattern = "all_SVs.rds", full.names = T)

dups_specs <- 
    list.files("results", pattern = "all_dups.rds", full.names = T)

asv_specs[2] %>% readRDS()
asv_specs[1] %>% readRDS()

aa_dups <- 
    lapply(
        dups_specs,
        readRDS
    ) %>% bind_rows()


aa_dups %>%count(species, collapsed) %>% filter(!is.na(collapsed))




all_dups %>% count(collapsed) 


ASV <- 
    lapply(
        asv_specs,
        readRDS
    ) %>% bind_rows() %>%
    mutate(
        s2 = sub( "_", " ", species),
        s2 = factor(s2, levels = top_species$Species[5:1]),
        spec = 
            recode(
                species,
                Salmonella_enterica        = "S. enterica",
                Staphylococcus_aureus      = "S. aureus",
                Mycobacterium_tuberculosis = "M. tuberculosis",
                Bordetella_pertussis       = "B. pertussis",
                Helicobacter_pylori        = "H. pylori"
        )
    )








ASV2 <- 
    ASV %>%
    mutate(
        spec = 
            factor(
                spec,
                levels = 
                    c(
                        "S. enterica",
                        "S. aureus",
                        "M. tuberculosis",
                        "B. pertussis",
                        "H. pylori"
                        )
                )
        )

sv_stats <- 
    ASV2 %>%
    # 1) flag each call as one of the four types
    mutate(
        is_indel        = variant =="Indel",
        is_translocation = variant == "Translocation",
        is_duplication   = variant == "Duplication",
        is_inversion     = variant == "Inversion"
        ) %>%
    group_by(spec, rid, qid) %>%
    summarise(
        SVs_per_comp   = n(),
        indels         = sum(is_indel),
        translocations = sum(is_translocation),
        duplications   = sum(is_duplication),
        inversions     = sum(is_inversion),
        .groups = "drop"
    ) %>%
    
    # 3) roll up to species‐level
    group_by(spec) %>%
    summarise(
        total_SVs            = sum(SVs_per_comp),
        avg_SVs_per_comp     = round(mean(SVs_per_comp), 1),
        pct_indels           = round(sum(indels)         / total_SVs * 100, 1),
        pct_translocations   = round(sum(translocations) / total_SVs * 100, 1),
        pct_duplications     = round(sum(duplications)   / total_SVs * 100, 1),
        pct_inversions       = round(sum(inversions)     / total_SVs * 100, 1),
        .groups = "drop"
    ) %>%
    arrange(desc(total_SVs))

sv_stats_pretty <- 
    sv_stats %>%
    dplyr::rename(
        Species                     = spec,
        `Total SVs`                 = total_SVs,
        `Mean SVs per<br>comparison`= avg_SVs_per_comp,
        Indels                      = pct_indels,
        Translocations              = pct_translocations,
        Duplications                = pct_duplications,
        Inversions                  = pct_inversions
    ) %>%
    mutate(
        # wrap species with cell_spec for italics
        Species = cell_spec(Species, italic = TRUE, escape = FALSE)
    )

# 2. Render with kableExtra
sv_stats_pretty %>%
    kbl(
        booktabs  = TRUE,
        escape    = FALSE,               # allow <br> and HTML markup
        align     = "lrrrrrr",
        linesep   = ""
    ) %>%
    # grouped header
    add_header_above(c(" " = 3, "Structural variant (%)" = 4), bold = TRUE) %>%
    
    # styling options
    kable_styling(
        bootstrap_options = c("condensed", "bordered", "hover"),
        full_width        = FALSE,
        position          = "center",
        font_size         = 12
    ) %>%
    # shade header row and bold it
    row_spec(0, bold = TRUE, background = "#D3D3D3") %>%
    # thick line under the group header
    row_spec(0, extra_css = "border-bottom: 2px solid #333;") %>%
    # vertical rules: left of Indels and right of Inversions
    column_spec(3, border_right  = TRUE) %>%
    column_spec(7, border_right = TRUE)




# 2.Join in the total per species and compute %
sv_type_pct <- sv_type_counts %>%
    left_join(sv_summary %>% select(species, total_SVs), by = "species") %>%
    mutate(
        pct = n_type / total_SVs,            # fraction
        pct_label = percent(pct, accuracy=0.1)  # nicely formatted
    ) %>%
    select(species, variant, n_type, total_SVs, pct, pct_label)



sv_summary <- 
    ASV %>%
    group_by(species) %>%
    mutate(
        SVs = n()
    ) %>%
    rowwise %>%
    mutate(
        rq = paste(rid, qid, sep = "_v_")
    ) %>%
    group_by(species, rq) %>%
    summarise(
        SV_per_comp = n(),
        mean_SV_per_comp = 
            SV_per_comp
    )
    
    group_by(species) %>%
    summarise(
        comparisons      = n_distinct(paste(rid,qid)), 
        total_SVs        = n(),
        mean_SVs         = mean(n()),               
        median_SVs       = median(n()),
        sd_SVs           = sd(n()),
        pct_duplication  = sum(variant=="duplication")/total_SVs,
        pct_indel        = sum(variant %in% c("insertion","deletion"))/total_SVs,
        pct_inversion    = sum(variant=="inversion")/total_SVs,
        pct_translocation= sum(variant=="translocation")/total_SVs
    ) %>%
    ungroup()

sv_summary %>% 
    mutate(across(starts_with("pct_"), ~ scales::percent(.x, accuracy=1))) %>%
    kbl(
        booktabs     = TRUE,
        caption      = "Summary of structural-variant load per species",
        align        = c("l", "r", "r", "r", "r", "r", rep("r",4)),
        col.names    = c(
            "Species", 
            "# Comparisons", 
            "Total SVs", 
            "Mean SVs/Comp.", 
            "Median SVs/Comp.", 
            "SD SVs/Comp.",
            "Dup", "Indel", "Inv", "Trans"
        )
    ) %>%
    add_header_above(c(
        " " = 6,    # first six cols are standalone
        "Proportion by SV type" = 4
    )) %>%
    kable_styling(full_width = FALSE)
# 4. Display ----------------------------------------------------------------

abundance_fig
ASV %>% count(species_abbrev, variant, variant_specific)


ASV %>%
    summary


all_proks

ASV %>% filter(s2=="Bordetella pertussis", variant=="Inversion")


ASV_simple <- 
    ASV %>%
    transmute(
        rid ,
        start, end, variant, variant_specific, width,
        species = s2, refmid,
        spec
        )

)

ASV_indels %>% count(species)
indel_plots <- 
    list()


indel_plots <- list()
gmm_summary <- list()

species_vec <- unique(ASV_simple$spec)

recursive_gmm <- 
    function(
        x, max_depth = 5, depth = 0, min_n = 100, 
        G_range = 1:6, proportion_cutoff = 0.1, spread_ratio_cutoff = 1.0
    ) {
        if (length(x) < min_n || depth >= max_depth) {
            return(tibble())
        }
        
        fit <- Mclust(x, G = G_range, modelNames = "V")
        means <- fit$parameters$mean
        sds   <- sqrt(fit$parameters$variance$sigmasq)
        props <- fit$parameters$pro
        
        summary_table <- tibble(
            depth = depth,
            submode = paste0("D", depth, ".", seq_along(means)),
            mean_log10 = round(means, 3),
            mean_bp = round(10^means),
            sd_log10 = round(sds, 3),
            sd_bp = round((10^(means + sds)) - (10^means)),
            proportion = round(props, 3),
            n = length(x)
        )
        
        nested_tables <- list(summary_table)
        
        for (i in seq_along(means)) {
            mean_log10 <- means[i]
            sd_log10 <- sds[i]
            prop <- props[i]
            mean_bp <- 10^mean_log10
            sd_bp <- (10^(mean_log10 + sd_log10)) - mean_bp
            
            if (prop > proportion_cutoff && (sd_bp / mean_bp) > spread_ratio_cutoff) {
                idx <- abs(x - mean_log10) < 2 * sd_log10
                if (sum(idx) >= min_n) {
                    nested_result <- recursive_gmm(
                        x[idx], max_depth = max_depth, depth = depth + 1,
                        min_n = min_n, G_range = 1:3
                    )
                    nested_tables <- append(nested_tables, list(nested_result))
                }
            }
        }
        
        bind_rows(nested_tables)
    }



gmm_tree <- 
    recursive_gmm(spec_indels$logw)
summary_table <- 
    flatten_recursive_gmm(gmm_tree)


plot_recursive_gmm_density <- function(x, gmm_table) {
    ggplot(data.frame(x = x), aes(x = x)) +
        geom_histogram(aes(y = ..density..), bins = 100, fill = "grey80", color = "white") +
        geom_density(color = "black", size = 0.6, adjust = 1.2) +
        geom_vline(data = gmm_table, aes(xintercept = mean_log10, color = as.factor(depth)), 
                   linetype = "dashed", size = 0.8) +
        scale_color_brewer(palette = "Dark2", name = "Depth") +
        labs(
            x = "log10(Indel length in bp)", 
            y = "Density",
            title = "Nested GMM Modes"
        ) +
        theme_minimal(base_size = 14)



for (i in seq_along(species_vec)) {
    spec_i <- species_vec[i]
    
    if (spec_i == "B. pertussis") next
    
    spec_indels <- 
        ASV_indels %>%
        filter(spec == spec_i, variant == "Indel") %>%
        mutate(logw = log10(width))
    
    gmm_tree <- recursive_gmm(spec_indels$width)
    summary_table <- flatten_recursive_gmm(gmm_tree)
    
    
    recursive_gmm(spec_indels$width, max_depth = 2)
    
}


nested_tables <- gmm_summary[grepl("_nested$", names(gmm_summary))]
final_nested_summary <- bind_rows(nested_tables)



ASV_simple %>% filter(spec=="S. enterica", width>=600, width<=800) %>% count(rid, variant_specific)


all_proks %>%
    filter(species=="Salmonella_enterica", grepl("acid stress", product)) 


all_proks %>%
    pull()


readDNAStringSet(
    sprintf(
        "../data/processing/sync/Pseudomonadota/Salmonella/Salmonella_enterica/%s.txt", 
        (ASV_simple %>% filter(spec=="S. enterica", width>=600, width<=800) %>% pull(rid))[2]
        )
) %>%
    subseq(
        start = 
            (ASV_simple %>% filter(spec=="S. enterica", width>=600, width<=800) %>% pull(start))[2],
        end = 
            (ASV_simple %>% filter(spec=="S. enterica", width>=600, width<=800) %>% pull(end))[2]
    ) %>% as.character()


x_log <- log10(spec_indels$width)
gmm_tree <- recursive_gmm(x_log)

# Plot
plot_recursive_gmm_density(x_log, gmm_tree)


bind_rows(gmm_summary)

    
    # indel_plots[[i]] <- 
    #     spec_indels %>%
    #     mutate(log_width = log10(width)) %>%
    #     ggplot(aes(x = log_width)) +
    #     geom_density(fill = "grey80", color = "black") +
    #     facet_wrap(~ spec, scales = "free_y") +
    #     labs(
    #         x = "log10(Indel Length, bp)",
    #         y = "Density",
    #         title = "Indel Length Distributions Across Species"
    #     ) +
    #     theme_classic(base_size = 14)
    
    
}

do.call(grid.arrange, c(indel_plots, ncol = 2))


ASV_indels <- 
    ASV_simple %>%
    filter(variant=="Indel")


##
ASV_indels %>% count(spec)

##3

hp_indels <- 
    ASV_indels %>%
    filter(spec =="H. pylori")


rm(HP_indels)
log_width <- log10(hp_indels$width)

# 3. Fit Gaussian mixture model
gmm_fit <- Mclust(log_width)

# 4. View model summary
summary(gmm_fit)


hp_indels$log_width <- log_width
hp_indels$cluster <- factor(gmm_fit$classification)








##########

log_width <- log10(ASV_indels %>% filter(spec == "H. pylori") %>% pull(width))

gmm_bic <- mclustBIC(log_width, G = 1:9) 
plot(gmm_bic)
best_model <- Mclust(log_width, x = gmm_bic)
best_model$G

gmm_constrained <- Mclust(log_width, G = 1:4) 


hp_modes <- 10^gmm_constrained$parameters$mean
names(hp_modes) <- paste0("Mode_", seq_along(hp_modes))

hp_modes
gmm_fixed4 <- Mclust(log_width, G = 4, modelNames = "V")
means <- gmm_fixed4$parameters$mean
sds   <- sqrt(gmm_fixed4$parameters$variance$sigmasq)
props <- gmm_fixed4$parameters$pro

# Create a grid of x values
x_vals <- seq(min(log_width), max(log_width), length.out = 1000)

# Construct density of each Gaussian component
dens_components <- sapply(1:4, function(i) {
    props[i] * dnorm(x_vals, mean = means[i], sd = sds[i])
})

# Total mixture density
mixture_density <- rowSums(dens_components)

# Make a tidy data frame for plotting
plot_df <- data.frame(
    x = x_vals,
    total = mixture_density,
    component1 = dens_components[,1],
    component2 = dens_components[,2],
    component3 = dens_components[,3],
    component4 = dens_components[,4]
)

# Plot
ggplot(plot_df, aes(x = x)) +
    geom_line(aes(y = total), size = 1.2, color = "black") +
    geom_line(aes(y = component1), color = "#D55E00", linetype = "dashed") +
    geom_line(aes(y = component2), color = "#0072B2", linetype = "dashed") +
    geom_line(aes(y = component3), color = "#009E73", linetype = "dashed") +
    geom_line(aes(y = component4), color = "#CC79A7", linetype = "dashed") +
    labs(x = "log10(Indel Length, bp)", y = "Density", title = "GMM Components (4 Modes)") +
    theme_classic()


10**gmm_fixed4$parameters$mean
### modelling each distribution

ASV_simple_gs <- 
    ASV_simple %>%
    filter(!is.na(width)) %>%
    count(spec, variant) %>%
    arrange(n)



## local peak detection
assign_kde_peaks <- function(df, core_thresh = 0.1, shoulder_thresh = 0.3) {
    if (nrow(df) < 5 || !"width" %in% names(df)) return(NULL)
    
    df <- df %>%
        filter(!is.na(width)) %>%
        mutate(log_length = log10(width + 1))
    
    kde <- density(df$log_length, bw = "nrd0")
    peaks <- findpeaks(kde$y, minpeakdistance = 5)
    if (is.null(peaks)) return(NULL)
    
    peak_pos <- kde$x[peaks[, 2]]
    peak_ids <- paste0("Peak_", seq_along(peak_pos))
    
    df %>%
        rowwise() %>%
        mutate(
            nearest = which.min(abs(log_length - peak_pos)),
            mode_group = peak_ids[nearest],
            mode_mean_log = peak_pos[nearest],
            mode_mean_bp  = 10^mode_mean_log,
            log_distance  = abs(log_length - mode_mean_log),
            zone = case_when(
                log_distance <= core_thresh ~ "core",
                log_distance <= shoulder_thresh ~ "shoulder",
                TRUE ~ "outlier"
            ),
            n_modes = length(peak_pos)
        ) %>%
        ungroup()
}


annotate_sv_modes <- function(df, width_col = "width", n_components = 1:6,
                              core_thresh = 1, shoulder_thresh = 2) {
    df <- df %>% filter(!is.na(.data[[width_col]]))
    if (nrow(df) < 50) return(NULL)
    
    log_lengths <- log10(df[[width_col]] + 1)
    gmm <- Mclust(log_lengths, G = n_components)
    pred <- predict(gmm)
    
    df %>%
        mutate(
            log_length       = log_lengths,
            mode_group       = paste0("Comp_", pred$classification),
            mode_mean_log    = gmm$parameters$mean[pred$classification],
            mode_mean_bp     = 10^mode_mean_log,
            log_distance     = abs(log_length - mode_mean_log),
            zone             = case_when(
                log_distance <= core_thresh ~ "core",
                log_distance <= shoulder_thresh ~ "shoulder",
                TRUE ~ "outlier"
            ),
            n_modes = length(gmm$parameters$mean)
        )
}



ap_cds <- 
    all_prokka %>%
    filter(type=="CDS") %>% 
    transmute(
        rid = accession, 
        width, start
    ) %>%
    inner_join(
        ., 
        ASV_simple %>% 
            dplyr::select(c(rid, spec)), 
        by = "rid", relationship = 'many-to-many'
    ) %>%
    distinct()


ap_cds %>% count(spec)




ASV_simple %>% count(spec)
all_prokka %>% filter(type !="CDS") %>%  
    group_by(type) %>% 
    summarise(
        ml = mean(width)
    )
ASV_v3 <- 
    ASV %>%
    transmute(
        accession = rid,
        start, end, width,
        variant, 
        vs = variant_specific,
        species, spec, refmid
    )


duplication specific


ASV %>% count(species, variant, variant_specific) %>% arrange(desc(n))

duplications


SV_list <-
    list()
for (i in seq_len(nrow(ASV_simple_gs))) {
    spec_i <- ASV_simple_gs$spec[i]
    var_i  <- ASV_simple_gs$variant[i]
    
    all_proks <- 
        spec_proks[
            grepl(
                strsplit(spec_i, " " )[[1]][2],
                spec_proks
                )
                ] %>%
        readRDS() 

    
    gr_genes <- 
        GRanges(
            seqnames = all_proks$accession,
            ranges   = IRanges(start = all_proks$start, end = all_proks$end),
            gene     = all_proks$gene,
            product  = all_proks$product
        )
    asv_i3 <- 
        ASV_v3 %>% 
        filter(spec ==spec_i, variant==var_i, !is.na(start))
    
    
    gr_asv <- 
        GRanges(
            seqnames = asv_i3$accession,
            ranges   = IRanges(start = asv_i3$start, end = asv_i3$end),
            sv_width = asv_i3$width,                 
            variant  = asv_i3$variant,
            species  = asv_i3$species,
            refmid   = asv_i3$refmid
        )
    
    all_proks$accession %>% unique
    all_proks$accession %>% unique
    
    SVs_ann_i <- 
        annotate_SVs(
            SV_granges = gr_asv,
            gene_granges = gr_genes 
        )
    
    SVs_ann_i %>% count(gene)
    
    
    spec_var %>% count(variant_specific)
    
    spec_var <- 
        ASV_simple %>%
        filter(spec == spec_i, variant == var_i, !is.na(width))
    
    
    hp_genomes <- 
        unique(spec_var$rid)
    
    hp_genes <- all_proks %>%
        filter(accession %in% hp_genomes)
    
    
    spec_var$width %>% density %>% plot
    
    n <- nrow(spec_var)
    kde_result <- 
        assign_kde_peaks(spec_var)
    
    
    kde_result %>% count(mode_mean_bp)
    
    ggplot(kde_result, aes(x = log10(width + 1), fill = zone)) +
        geom_density(alpha = 0.5) +
        labs(title = paste0(spec_i, " ", var_i, " (n = ", n, ")"),
             x = "log10(SV length)", y = "Density") +
        scale_fill_manual(values = c("core" = "steelblue", "shoulder" = "orange", "outlier" = "grey50")) +
        theme_classic()
    
    if (n >= 300) {
        tib_i <- annotate_sv_modes(spec_var, n_components = 1:6)
    } else if (n >= 100) {
        tib_i <- annotate_sv_modes(spec_var, n_components = 1:3)
    } else {
        tib_i <- 
            assign_kde_peaks(spec_var)
    }
    
    if (!is.null(tib_i)) {
        SV_list[[i]] <- tib_i
    }
}

# objs <- ls(envir = .GlobalEnv)
# sizes <- sapply(objs, function(x) object.size(get(x)))
# sizes_df <- data.frame(
#     object = objs,
#     size_MB = as.numeric(sizes) / (1024^2)
# ) %>%
#     arrange(desc(size_MB))
# 
# print(head(sizes_df, 20))

SV_char <- bind_rows(SV_list)


ASV %>%
    filter(species %in% c("Bordetella_pertussis", ))



SV_char %>%
    count(spec, variant) %>%
    ggplot(aes(x = variant, y = fct_rev(spec), fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), color = "black", size = 3) +
    scale_fill_gradient(low = "grey90", high = "steelblue") +
    labs(title = "SV Counts by Species and Variant Type", fill = "SV Count") +
    theme_minimal() +
    theme(axis.title = element_blank())


SV_char %>%
    distinct(spec, variant, n_modes) %>%
    ggplot(aes(x = variant, y = fct_rev(spec), fill = n_modes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_modes), color = "black", size = 3) +
    scale_fill_gradient(low = "grey90", high = "darkred") +
    labs(title = "# Modes Detected per Group", fill = "n Modes") +
    theme_minimal() +
    theme(axis.title = element_blank())
SV_char %>%
    count(spec, variant, zone) %>%
    group_by(spec, variant) %>%
    mutate(prop = n / sum(n)) %>%
    ggplot(aes(x = variant, y = fct_rev(spec), fill = zone, alpha = prop)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("core" = "steelblue", "shoulder" = "orange", "outlier" = "grey50")) +
    labs(title = "Zone Composition of SV Length Distributions") +
    theme_minimal() +
    theme(axis.title = element_blank())


mode_summary <- 
    SV_char %>%
    group_by(spec, variant, mode_group, n_modes) %>%
    summarise(
        Mode_size      = n(),
        Core_count     = sum(zone == "core"),
        Shoulder_count = sum(zone == "shoulder"),
        Outlier_count  = sum(zone == "outlier"),
        Core_pct       = Core_count / Mode_size,
        Median_length  = median(width, na.rm = TRUE),
        Mean_length    = mean(width, na.rm = TRUE),
        Mode_mean_bp   = dplyr::first(mode_mean_bp),
        .groups = "drop"
    )
ggplot(SV_char, aes(x = log10(width + 1), fill = zone)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ paste(spec, variant), scales = "free_y") +
    labs(x = "log10(SV length)", y = "Density", title = "SV length zones by group") +
    scale_fill_manual(values = c("core" = "steelblue", "shoulder" = "orange", "outlier" = "grey50")) +
    theme_classic()


SV_char %>%
    count(spec, variant, zone) %>%
    group_by(spec, variant) %>%
    mutate(prop = n / sum(n)) %>%
    ggplot(aes(x = interaction(spec, variant), y = prop, fill = zone)) +
    geom_col() +
    coord_flip() +
    labs(x = "Species / Variant", y = "Proportion", fill = "Zone",
         title = "Proportion of SVs in Core / Shoulder / Outlier") +
    scale_fill_manual(values = c("core" = "steelblue", "shoulder" = "orange", "outlier" = "grey50")) +
    theme_classic()

SV_char %>%
    distinct(spec, variant, n_modes) %>%
    ggplot(aes(x = interaction(spec, variant), y = n_modes)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(x = "Species / Variant", y = "# of Modes Detected",
         title = "Modal Complexity of SV Length Distributions") +
    theme_classic()

SV_char %>%
    group_by(spec, variant, zone) %>%
    summarise(median_bp = median(width), .groups = "drop") %>%
    ggplot(aes(x = zone, y = median_bp, fill = zone)) +
    geom_boxplot() +
    facet_wrap(~ interaction(spec, variant), scales = "free_y") +
    labs(y = "Median SV Length (bp)", x = "Zone", title = "SV Length by Zone") +
    scale_fill_manual(values = c("core" = "steelblue", "shoulder" = "orange", "outlier" = "grey50")) +
  theme_classic()



bind_rows(SV_list) %>% count(spec, variant, mode_mean_bp)



### annotation version 2










annotate_SVs


all_proks

asv_s2 <- 
    ASV %>%
    transmute(
        accession = rid,
        start, end, width,
        variant, 
        vs = variant_specific,
        species, refmid
        ) %>%
    filter(!is.na(start))
    filter(width>800, width<2000)


    all_proks <- 
        all_proks %>% filter(!is.na(start), type=="CDS")
    
    
    
    all_proks %>% 
        filter(
            #species=="Bordetella_pertussis") | 
        species=="Helicobacter_pylori")
gr_genes <- 
    GRanges(
        seqnames = all_proks$accession,
        ranges   = IRanges(start = all_proks$start, end = all_proks$end),
        gene     = all_proks$gene,
        product  = all_proks$product
        )


gr_asv <- 
    GRanges(
        seqnames = asv_s2$accession,
        ranges   = IRanges(start = asv_s2$start, end = asv_s2$end),
        sv_width = asv_s2$width,                 
        variant  = asv_s2$variant,
        species  = asv_s2$species,
        refmid   = asv_s2$refmid
    )
SVs_ann <- 
    annotate_SVs(
        SV_granges = gr_asv,
        gene_granges = gr_genes 
    )

SVs_ann %>%
    count(species, variant, func_group)


SVs_ann %>% 

SVs_ann %>% filter(func_group=="Other") %>% 
    count(product, sort = T) %>% dplyr::slice(1:20) %>% pull(product)

SVs_ann %>%
    count(species, variant, func_group) %>%
    # order species by total SVs (so the largest bar is leftmost)
    mutate(
        species = fct_reorder(species, n, .fun = sum),
        species_abbrev = 
            recode(
                species,
                Bordetella_pertussis       = "B. pertussis",
                Helicobacter_pylori        = "H. pylori",
                Mycobacterium_tuberculosis = "M. tuberculosis",
                Salmonella_enterica        = "S. enterica",
                Staphylococcus_aureus      = "S. aureus"
            )
        ) %>%
    ggplot(aes(x = species_abbrev, y = n, fill = func_group)) +
    geom_col(color = "black", width = 0.7) +
    facet_wrap(~ variant, scales = "free_y", nrow = 1) +
    # scale_fill_manual(
    #     "Functional\nGroup", 
    #     values = c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
    #     ) +
    scale_fill_manual(
        "Functional\nGroup", 
        values = 
            c(
                "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                "#A65628", "orangered", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
                "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77",
                "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#A6761D", "#666666"
    )) +
    #scale_fill_brewer("Functional\nGroup", palette = "Set3") +
    labs(
        x     = "Species",
        y     = "Number of SVs",
        title = "Gene annoations within each SV group"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title       = element_text(hjust = 0.5, face = "bold"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.text       = element_text(face = "bold", size = 16),
        legend.position  = "right",
    )

df_lumped <- SVs_ann %>%
    count(species, variant, func_group) %>%        # get raw counts
    group_by(species, variant) %>%                 # by facet
    mutate(
        total    = sum(n),
        prop     = n / total,
        func_grp = if_else(prop < 0.05,              # threshold: 5%
                           "Other",
                           as.character(func_group))
    ) %>%
    ungroup() %>%
    # 2) re‐summarize so all the “Other” bits collapse together
    group_by(species, variant, func_grp) %>%
    summarize(n = sum(n), .groups = "drop") %>%
    mutate(func_grp = fct_relevel(func_grp, "Other", after = Inf))


SVs_ann %>%
    filter(func_group=="Other") %>%
    count(product, sort = T)


# 3) Plot exactly as before, but with the new func_grp
ggplot(
    df_lumped,
       aes(x = fct_reorder(species, n, .fun = sum),
           y = n,
           fill = func_grp)) +
    geom_col(color = "black", width = 0.7) +
    facet_wrap(~ variant, scales = "free_y", nrow = 1) +
    #scale_fill_brewer("Functional\nGroup", palette = "Set3")
    scale_fill_manual(values = brewer.pal(12, "Paired")) +
    #scale_fill_viridis_d("Functional\nGroup",option = "C")
    labs(
        title = "Distribution of SV Functional Groups by Species and Variant Type",
        x     = "",
        y     = "Number of SVs"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        axis.text.x     = element_text(angle = 45, hjust = 1),
        strip.text      = element_text(face = "bold", size = 12),
        legend.position = "right"
    )



# generate a 12‑colour palette
set.seed(123)
my_pal12 <- distinctColorPalette(12)


annotate_SVs <- 
    function(SV_granges, gene_granges, overlap_perc = 90) {
    # 0) make sure both GRanges share the same seqlevels
        common_seq <- intersect(seqlevels(SV_granges), seqlevels(gene_granges))
        SVg   <- keepSeqlevels(SV_granges, common_seq, pruning.mode="coarse")
        Geneg <- keepSeqlevels(gene_granges, common_seq, pruning.mode="coarse")
        
        # 1) find all overlaps
        hits <- findOverlaps(SVg, Geneg)
        qh   <- queryHits(hits)
        sh   <- subjectHits(hits)
        
        if (length(hits)) {
            # 2) compute bp overlap and percent of SV
            ovl_gr     <- pintersect(SVg[qh], Geneg[sh])
            overlap_bp <- width(ovl_gr)
            sv_len     <- width(SVg)[qh]
            pct_overlap<- overlap_bp / sv_len * 100
            
            # 3) filter by threshold
            keep   <- pct_overlap >= overlap_perc
            qh2    <- qh[keep]
            sh2    <- sh[keep]
            ov2    <- overlap_bp[keep]
            pct2   <- pct_overlap[keep]
            
            matched <- tibble(
                # gene side
                gene_accession = as.character(seqnames(Geneg))[sh2],
                gene           = Geneg$gene[sh2],
                product        = Geneg$product[sh2],
                # SV side
                sv_accession   = as.character(seqnames(SVg))[qh2],
                species        = SVg$species[qh2],
                variant        = SVg$variant[qh2],
                sv_start       = start(SVg)[qh2],
                sv_end         = end(SVg)[qh2],
                sv_width       = sv_len[keep],
                refmid         = SVg$refmid[qh2],
                overlap_bp     = ov2,
                pct_overlap    = pct2
            )
        } else {
            matched <- tibble()  # no overlaps 
        }
        
        # 4) build the “non‐overlapping” SVs
        all_idx   <- seq_along(SVg)
        nohit_idx <- setdiff(all_idx, unique(qh))
        if (length(nohit_idx)==0) {
            nonmatched <- tibble(
                gene_accession = NA_character_,
                gene           = NA_character_,
                product        = NA_character_,
                sv_accession   = as.character(seqnames(SVg))[nohit_idx],
                species        = SVg$species[nohit_idx],
                variant        = SVg$variant[nohit_idx],
                sv_start       = start(SVg)[nohit_idx],
                sv_end         = end(SVg)[nohit_idx],
                sv_width       = width(SVg)[nohit_idx],
                refmid         = SVg$refmid[nohit_idx],
                overlap_bp     = 0L,
                pct_overlap    = 0
            )
        } else {
            nonmatched <- tibble()
        }
        
        # 5) bind and return
        bind_rows(matched, nonmatched) %>%
            arrange(sv_accession, sv_start) %>%
            mutate(
                func_group = case_when(
                    grepl(
                        "transposase",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Transposase",
                    grepl(
                        "hypothetical",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Hypothetical protein",
                    grepl(
                        "putative",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Putative protein",
                    grepl(
                        "^23S ribosomal RNA",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "rRNA",
                    grepl(
                        "ribosomal RNA",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "rRNA",
                    grepl(
                        "adhesin",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Adhesin",
                    grepl(
                        "flagellin|flagellar",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Flagellar proteins",
                    grepl(
                        "restriction enzyme|EcoKI specificity",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Restriction‑modification",
                    grepl(
                        "phosphomanno|phosphogluco|malic enzyme|amidohydrolase|hydrolase|mutase",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Metabolic enzyme",
                    grepl(
                        "chaperone",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Chaperone",
                    grepl(
                        "competence protein",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Competence",
                    grepl(
                        "toxin|antitoxin",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Toxin‑antitoxin system",
                    grepl(
                        "phage|capsid|tail",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Phage‑associated",
                    grepl(
                        "ESAT-6|Esx|type VII|ESX-",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "ESX secretion family",
                    grepl(
                        "ABC transporter",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "ABC transporter",
                    grepl(
                        "transporter|efflux pump|uptake",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Transporter",
                    grepl(
                        "kinase",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Signaling",
                    grepl(
                        "dehydrogenase|ligase|reductase|synthase|oxygenase",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Enzyme",
                    grepl(
                        "ribosomal protein",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Ribosomal protein",
                    grepl(
                        "catalase-peroxidase",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Oxidative stress",
                    grepl(
                        "partition protein|Smc",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Chromosome partition",
                    grepl(
                        "coagulase",
                        product,
                        ignore.case = TRUE
                    ) ~
                        "Coagulase",
                    gene %in% c("rho", "nusA") ~
                        "Regulatory",
                    gene %in% c("rpfA", "rpfE") ~
                        "Dormancy/Resuscitation",
                    grepl(
                        "repeat(-|_)?containing|repeat", product, 
                        ignore.case=TRUE) 
                    ~ "Adhesin",
                    grepl("ase$", product, ignore.case=TRUE)  
                    ~ "Enzyme",
                    is.na(product) & is.na(gene) ~
                        NA_character_,
                    TRUE ~
                        "Other"
                )
            )
        
        
}




annotate_SVs2 <- function(SV_granges, gene_granges, overlap_perc = 90) {
    # 0) Align seqlevels
    common_seq <- intersect(seqlevels(SV_granges), seqlevels(gene_granges))
    SVg   <- keepSeqlevels(SV_granges, common_seq, pruning.mode = "coarse")
    Geneg <- keepSeqlevels(gene_granges, common_seq, pruning.mode = "coarse")
    
    # 1) Find overlaps
    hits <- findOverlaps(SVg, Geneg)
    qh   <- queryHits(hits)
    sh   <- subjectHits(hits)
    
    sv_len <- width(SVg)
    gene_len <- width(Geneg)
    
    all_idx <- seq_along(SVg)
    annotated <- vector("list", length(SVg))
    
    for (i in all_idx) {
        overlaps <- which(qh == i)
        if (length(overlaps) == 0) {
            annotated[[i]] <- tibble(
                gene_accession = NA_character_,
                gene           = NA_character_,
                product        = NA_character_,
                sv_accession   = as.character(seqnames(SVg)[i]),
                species        = SVg$species[i],
                variant        = SVg$variant[i],
                sv_start       = start(SVg)[i],
                sv_end         = end(SVg)[i],
                sv_width       = sv_len[i],
                refmid         = SVg$refmid[i],
                overlap_bp     = 0L,
                pct_overlap    = 0,
                func_group     = "intergenic"
            )
        } else {
            rows <- lapply(overlaps, function(j) {
                ovl <- pintersect(SVg[qh[j]], Geneg[sh[j]])
                ov_bp <- width(ovl)
                pct_sv <- 100 * ov_bp / sv_len[qh[j]]
                pct_gene <- 100 * ov_bp / gene_len[sh[j]]
                
                g <- Geneg$gene[sh[j]]
                p <- Geneg$product[sh[j]]
                
                func_group <- if (pct_sv >= overlap_perc && pct_gene >= overlap_perc) {
                    case_when(
                        grepl("transposase", p, ignore.case = TRUE) ~ "Transposase",
                        grepl("hypothetical", p, ignore.case = TRUE) ~ "Hypothetical protein",
                        grepl("putative", p, ignore.case = TRUE) ~ "Putative protein",
                        grepl("^23S ribosomal RNA", p, ignore.case = TRUE) ~ "rRNA",
                        grepl("ribosomal RNA", p, ignore.case = TRUE) ~ "rRNA",
                        grepl("adhesin", p, ignore.case = TRUE) ~ "Adhesin",
                        grepl("flagellin|flagellar", p, ignore.case = TRUE) ~ "Flagellar proteins",
                        grepl("restriction enzyme|EcoKI specificity", p, ignore.case = TRUE) ~ "Restriction‑modification",
                        grepl("phosphomanno|phosphogluco|malic enzyme|amidohydrolase|hydrolase|mutase", p, ignore.case = TRUE) ~ "Metabolic enzyme",
                        grepl("chaperone", p, ignore.case = TRUE) ~ "Chaperone",
                        grepl("competence protein", p, ignore.case = TRUE) ~ "Competence",
                        grepl("toxin|antitoxin", p, ignore.case = TRUE) ~ "Toxin‑antitoxin system",
                        grepl("phage|capsid|tail", p, ignore.case = TRUE) ~ "Phage‑associated",
                        grepl("ESAT-6|Esx|type VII|ESX-", p, ignore.case = TRUE) ~ "ESX secretion family",
                        grepl("ABC transporter", p, ignore.case = TRUE) ~ "ABC transporter",
                        grepl("transporter|efflux pump|uptake", p, ignore.case = TRUE) ~ "Transporter",
                        grepl("kinase", p, ignore.case = TRUE) ~ "Signaling",
                        grepl("dehydrogenase|ligase|reductase|synthase|oxygenase", p, ignore.case = TRUE) ~ "Enzyme",
                        grepl("ribosomal protein", p, ignore.case = TRUE) ~ "Ribosomal protein",
                        grepl("catalase-peroxidase", p, ignore.case = TRUE) ~ "Oxidative stress",
                        grepl("partition protein|Smc", p, ignore.case = TRUE) ~ "Chromosome partition",
                        grepl("coagulase", p, ignore.case = TRUE) ~ "Coagulase",
                        g %in% c("rho", "nusA") ~ "Regulatory",
                        g %in% c("rpfA", "rpfE") ~ "Dormancy/Resuscitation",
                        grepl("repeat(-|_)?containing|repeat", p, ignore.case = TRUE) ~ "Adhesin",
                        grepl("ase$", p, ignore.case = TRUE) ~ "Enzyme",
                        is.na(p) & is.na(g) ~ NA_character_,
                        TRUE ~ "Other"
                    )
                } else {
                    "partially genic"
                }
                
                tibble(
                    gene_accession = as.character(seqnames(Geneg)[sh[j]]),
                    gene           = if (func_group == "partially genic") NA_character_ else g,
                    product        = if (func_group == "partially genic") NA_character_ else p,
                    sv_accession   = as.character(seqnames(SVg)[qh[j]]),
                    species        = SVg$species[qh[j]],
                    variant        = SVg$variant[qh[j]],
                    sv_start       = start(SVg)[qh[j]],
                    sv_end         = end(SVg)[qh[j]],
                    sv_width       = sv_len[qh[j]],
                    refmid         = SVg$refmid[qh[j]],
                    overlap_bp     = ov_bp,
                    pct_overlap    = pct_sv,
                    func_group     = func_group
                )
            })
            annotated[[i]] <- bind_rows(rows) %>% slice_max(overlap_bp, with_ties = FALSE)
        }
    }
    
    bind_rows(annotated)
}








ASV_simple2 <- 
    ASV_simple %>%
    mutate(
        variant = factor(variant, levels = c(
            "Indel","Duplication","Inversion","Translocation"
        ))
    ) %>% filter(!is.na(width), width<5e4)


ggplot(data = ASV_simple2 %>% filter(species == "Helicobacter pylori", variant == "Indel"),
       aes(x = log10(width + 1))) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    labs(x = "Log10 SV Length", y = "Density", title = "Indel length distribution in H. pylori") +
    theme_classic()

hp_indels <-  
    ASV_simple2 %>% filter(species == "Helicobacter pylori", variant == "Indel")%>%
    mutate(log_length = log10(width + 1))

gmm_hp <- Mclust(hp_indels$log_length, G = 1:6) 
gmm_table <- tibble(
    Component         = seq_along(gmm_hp$parameters$mean),
    Mean_log10_length = gmm_hp$parameters$mean,
    SD_log10_length   = sqrt(gmm_hp$parameters$variance$sigmasq),
    Mixing_proportion = gmm_hp$parameters$pro
) %>%
    mutate(
        Mean_length_bp = 10^Mean_log10_length,
        SD_length_bp   = 10^(Mean_log10_length + SD_log10_length) - 10^Mean_log10_length
    )


hp_pred <- predict(gmm_hp)

ggplot(hp_indels, aes(x = log_length)) +
    geom_histogram(aes(y = ..density..), bins = 80, fill = "grey85", color = "white") +
    stat_function(
        fun = function(x) {
            rowSums(sapply(1:length(gmm_hp$parameters$mean), function(k) {
                dnorm(x,
                      mean = gmm_hp$parameters$mean[k],
                      sd   = sqrt(gmm_hp$parameters$variance$sigmasq)[k]) * gmm_hp$parameters$pro[k]
            }))
        },
        color = "blue", size = 1.2
    ) +
    labs(title = "GMM Fit for H. pylori Indel Lengths",
         x = "log10(SV Length in bp)", y = "Density") +
    theme_classic()


kde <- density(hp_indels$log_length, bw = "nrd0")
peaks <- findpeaks(kde$y, minpeakdistance = 10)  # Tune minpeakdistance if needed
peak_positions <- kde$x[peaks[, 2]]


plot(kde, main = "KDE of H. pylori Indels (log10 length)")
points(peak_positions, peaks[,1], col = "red", pch = 19)

hp_indels <- 
    hp_indels %>%
    mutate(
        assigned_comp = hp_pred$classification,
        comp_mean     = gmm_hp$parameters$mean[assigned_comp],
        comp_sd       = sqrt(gmm_hp$parameters$variance$sigmasq)[assigned_comp],
        z_score       = (log_length - comp_mean) / comp_sd,
        zone          = case_when(
            abs(z_score) < 1 ~ "core",
            abs(z_score) < 2 ~ "shoulder",
            TRUE             ~ "outlier"
        )
    )

hp_indels


core_summary <- hp_indels %>%
    filter(zone == "core") %>%
    count(assigned_comp, name = "Core_SVs") %>%
    right_join(
        hp_indels %>% count(assigned_comp, name = "Total_SVs"),
        by = "assigned_comp"
    ) %>%
    mutate(
        Proportion_core = Core_SVs / Total_SVs,
        Mean_log10_length = gmm_hp$parameters$mean[assigned_comp],
        SD_log10_length   = sqrt(gmm_hp$parameters$variance$sigmasq)[assigned_comp],
        Mean_length_bp = 10^Mean_log10_length,
        SD_length_bp   = 10^Mean_log10_length * sqrt(10^(SD_log10_length^2) - 1)
    ) %>%
    arrange(Mean_length_bp)

core_summary %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    kable()

log_lengths <- log10(hp_indels$width + 1)
gmm <- Mclust(log_lengths, G = 1:10)  # try up to 6 components



kde <- density(hp_indels$log_length, bw = "nrd0")
peaks <- findpeaks(kde$y, minpeakdistance = 10)  # Tune minpeakdistance if needed
peak_positions <- kde$x[peaks[, 2]]

pred <- predict(gmm)

assigned_mode <- pred$classification
mode_means <- gmm$parameters$mean[assigned_mode]
mode_sds   <- sqrt(gmm$parameters$variance$sigmasq)[assigned_mode]

z_scores <- (log_lengths - mode_means) / mode_sds

within_1sd_counts <- table(assigned_mode[abs(z_scores) <= 1])
total_by_mode     <- table(assigned_mode)


proximity_summary <- tibble(
    Component = as.integer(names(total_by_mode)),
    Total_assigned     = as.integer(total_by_mode),
    Within_1SD         = as.integer(within_1sd_counts[names(total_by_mode)]),
    Proportion_within  = Within_1SD / Total_assigned
)


means <- gmm$parameters$mean
sds   <- sqrt(gmm$parameters$variance$sigmasq)
props <- gmm$parameters$pro

# Component assignment and z-scores
log_lengths <- log10(width + 1)
pred <- predict(gmm)
assigned_mode <- pred$classification
assigned_mean <- means[assigned_mode]
assigned_sd   <- sds[assigned_mode]
z_scores <- (log_lengths - assigned_mean) / assigned_sd

# Counts
within_1sd_counts <- table(assigned_mode[abs(z_scores) <= 1])
total_by_mode     <- table(assigned_mode)

# Combine everything
library(tibble)

combined_summary <- tibble(
    Component         = seq_along(means),
    Mixing_proportion = props,
    Mean_log10_length = means,
    SD_log10_length   = sds,
    Mean_length_bp    = 10^means,
    SD_length_bp      = 10^(means + sds) - 10^means,
    Total_assigned    = as.integer(total_by_mode[as.character(seq_along(means))]),
    Within_1SD        = as.integer(within_1sd_counts[as.character(seq_along(means))]),
    Proportion_within = Within_1SD / Total_assigned
)


tibble(
    Component = seq_along(gmm$parameters$mean),
    Mean_log10_length = gmm$parameters$mean,
    SD_log10_length   = sqrt(gmm$parameters$variance$sigmasq),
    Mixing_proportion = gmm$parameters$pro
) %>%
    mutate(
        Mean_length_bp = 10^Mean_log10_length,
        SD_length_bp   = 10^(Mean_log10_length + SD_log10_length) - 10^Mean_log10_length
    )


log_length_df <- data.frame(
    log_length = log10(hp_indels$width + 1),
    cluster = factor(pred$classification)
)

ggplot(log_length_df) +
    geom_density(aes(x = log_length), fill = "grey90") +
    geom_density(
        data = subset(log_length_df, cluster == 5),
        aes(x = log_length), fill = "red", alpha = 0.5
    ) +
    labs(
        title = "Component 5 SVs",
        x = "log10(SV length)",
        y = "Density"
    ) +
    theme_classic()


ggplot(data.frame(
    log_length = log10(width + 1),
    cluster = factor(pred$classification)
)) +
    geom_density(aes(x = log_length), fill = "grey90") +
    geom_density(data = subset(data.frame(
        log_length = log10(width + 1),
        cluster = factor(pred$classification)
    ), cluster == 5), aes(x = log_length), fill = "red", alpha = 0.5) +
    labs(title = "Component 5 SVs", x = "log10(SV length)", y = "Density")
plot(gmm, what = "BIC")
summary(gmm)
plot(gmm, what = "density")  #


asv_groups <- 
    ASV_simple2 %>% count(species, variant)
    
ASV_simple2 %>% add_count(species, variant)

ASV_simple2 %>% filter(species=="Helicobacter pylori", variant=="Indel") %>% pull(width) %>% boxplot



ASV_simple2 %>% count(species, variant)

make_panel <- function(var, show_y = FALSE) {
    ggplot(filter(ASV_simple2, variant == var),
           aes(x = species, y = width, fill = species)) +
        geom_violin(
            trim  = FALSE,
            scale = "width",
            alpha = 0.6,
            color = "darkgrey",   # dark grey border
            size  = 0.5
        ) +
        geom_boxplot(
            width        = 0.2,
            outlier.size = 1,
            color        = "black",
            alpha        = 0.8
        ) +
        scale_y_log10(
            breaks = c(100, 1e3, 1e4, 1e5, 1e6),
            labels = c("100","1k","10k","100k","1M")
        ) +
        scale_fill_brewer(palette = "Dark2") +  # darker, no yellow
        coord_flip() +
        labs(title = var) +
        theme_minimal(base_size = 18) +
        theme(
            plot.title       = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title     = element_blank(),       
            axis.text.x      = element_text(size = 18),
            axis.text.y      = if (show_y) element_text(size = 18) else element_blank(),
            axis.ticks.y     = if (show_y) element_line() else element_blank(),
            panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_blank(),
            legend.position  = "none"
        )
}

# 3) Build each panel (only the first shows species labels)
p1 <- make_panel("Indel",    show_y = TRUE)
p2 <- make_panel("Duplication")
p3 <- make_panel("Inversion")
p4 <- make_panel("Translocation")

# 4) Stitch them in one row, collect any legends (there are none), 
#    then add a single x‐axis label and a title
combined <- p1 + p2 + p3 + p4 +
    plot_layout(ncol = 4)

# Add global title + a single “x‑axis” label via caption
final_plot <- combined +
    plot_annotation(
        title   = "",
        caption = "Length (bp)",
        theme = theme(
            plot.title   = element_text(size = 22, hjust = 0.5),
            plot.caption = element_text(size = 18, hjust = 0.5, margin = ggplot2::margin(t = 10))
        )
    ) &
    theme(
        plot.title       = element_text(size = 20, hjust = 0.5),
        plot.caption     = element_text(size = 24, hjust = 0.5, margin = ggplot2::margin(t = 10)),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text       = element_text(face = "bold", size = 20),
        axis.text        = element_text(size = 20),
        axis.title.x     = element_blank(),    # we’re using caption instead
        legend.position  = "none"
    )





final_plot



ASV_simple2 %>%
    group_by(species, variant) %>%
    summarise(n = n()) %>% arrange(n)





#library(diptest) 

##

get_modes <- function(x, dip_p_cut = 0.05, G = 1:5) {
    # dip.test: H0 = unimodal.  p > cut ⇒ unimodal ⇒ use density().
    if (dip.test(x)$p.value > dip_p_cut) {
        d <- density(x)
        d$x[which.max(d$y)]
    } else {
        # multimodal ⇒ fit a Gaussian mixture, letting Mclust choose #components
        mc <- Mclust(x, G = G, verbose = FALSE)
        sort(mc$parameters$mean)
    }
}



modes_by_grp <- 
    ASV_simple2 %>%
    filter(!is.na(width)) %>%
    group_by(species, variant) %>%
    summarize(
        SVs = n(),
        modes = list(get_modes(width)),
        n_modes = length(modes[[1]]),
        is_normal = if (n() > 30 & n_modes == 1) {
            test_data <- if (n() > 5000) sample(width, 5000) else width
            shapiro.test(test_data)$p.value > 0.05
        } else {
            F
        },
        .groups = "drop"
    )

# If you want one row per mode (long format):
modes_long <- 
    modes_by_grp %>%
    unnest_longer(modes) %>%
    dplyr::rename(mode = modes)



modes_by_grp %>%
    mutate(modes = sapply(modes, function(xs) paste(round(xs,1), collapse = ", ")))%>%
    kable(
        format = "html",
        caption = "Modes of SV‑width distributions by species and variant"
    ) %>%
    kable_styling(full_width = FALSE)


ASV_simple2 %>%filter(!is.na(width)) %>% filter(species=="Helicobacter pylori") %>% pull(width) %>% density() %>% plot

ASV_assigned <- ASV_simple2 %>%
    filter(!is.na(width), width<1e4) %>%
    group_by(species, variant) %>%
    group_modify(~{
        x <- .x$width
        if (dip.test(x)$p.value > 0.01) {
            # unimodal ⇒ everyone in one cluster
            .x %>% mutate(
                mode_cluster = 1L,
                mode_value   = density(x)$x[which.max(density(x)$y)]
            )
        } else {
            # multimodal ⇒ fit GMM, assign cluster
            mc <- Mclust(x, G = 1:5, verbose = FALSE)
            cl <- predict(mc, x)$classification
            .x %>%
                mutate(
                    mode_cluster = cl,
                    mode_value   = mc$parameters$mean[cl]
                )
        }
    }) %>%
    ungroup()


get_modes_peaks <- function(x,
                            adj = 0.5,            # smaller = less smoothing
                            min_prom = 0.05) {    # relative to max density
    x <- na.omit(x)
    d <- density(x, adjust = adj)
    # find all peaks whose height ≥ min_prom * max(d$y)
    peaks_info <- findpeaks(
        d$y,
        minpeakheight = max(d$y) * min_prom
    )
    if (is.null(peaks_info)) {
        # no little peaks? fallback to the global max
        return(d$x[which.max(d$y)])
    }
    # column 2 of peaks_info are the indices in d$y
    peaks_idx <- peaks_info[,2]
    sort(d$x[peaks_idx])
}



ASV_clustered2 <- ASV_simple2 %>%
    filter(!is.na(width)) %>%
    group_by(species, variant) %>%
    group_modify(~{
        x <- .x$width
        
        # Decide whether to fit a mixture
        #   - Always for Indels
        #   - Or if dip.test p ≤ 0.01 (more stringent)
        force_gmm <- as.character(.y$variant) == "Indel"
        multi     <- dip.test(x)$p.value <= 0.01
        
        if (force_gmm || multi) {
            mc <- Mclust(x, G = 1:10, verbose = FALSE)
            pr <- predict(mc, x)
            
            .x %>% mutate(
                mode_cluster  = pr$classification,
                mode_value    = mc$parameters$mean[mode_cluster],
                posterior_max = apply(pr$z, 1, max)
            )
        } else {
            dens <- density(x, adjust = 0.5)
            peak <- dens$x[which.max(dens$y)]
            .x %>% mutate(
                mode_cluster  = 1L,
                mode_value    = peak,
                posterior_max = 1
            )
        }
    }) %>%
    ungroup()

# Inspect:
ASV_clustered2 %>% 
    count(species, variant, mode_cluster) %>% 
    print(n = Inf)


ASV_clustered2 %>% filter(mode_value>800, mode_value<2000)




###

all_prokka

























library(circlize)



df <- 
    ASV_clustered2 %>%
    filter(species=="Helicobacter pylori", variant=="Indel", mode_cluster==6)

clusters      <- sort(unique(df$mode_cluster))
genome_length <- max(df$refmid)

circos.clear()
# little gaps between sectors
circos.par(start.degree = 90, gap.degree = rep(2, length(clusters)))

# initialize one sector per cluster, all spanning 0 → genome_length
circos.initialize(
    factors = as.character(clusters),
    xlim    = matrix(
        c(rep(0, length(clusters)), rep(genome_length, length(clusters))),
        ncol=2
    )
)

# draw a single track; panel.fun is called once per sector
circos.trackPlotRegion(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
        cl <- CELL_META$sector.index
        pos <- df %>% filter(mode_cluster == as.integer(cl)) %>% pull(refmid)
        # plot all points for this cluster at y = 0.5
        circos.points(pos, rep(0.5, length(pos)), pch = 16, cex = 0.4)
        # optionally add the mode line
        mode_val <- unique(df$mode_value[df$mode_cluster == as.integer(cl)])
        circos.lines(c(mode_val, mode_val), c(0,1), lwd=2, col="red")
    }
)

circos.clear()


# example: detect modes for indels
indel_modes2 <- ASV_simple2 %>%
    filter(variant=="Indel", !is.na(width)) %>%
    pull(width) %>%
    get_modes_peaks(adj = 0.3, min_prom = 0.03)



ggplot(ASV_simple2, aes(x = species, y = width, fill = species)) +
    geom_violin(
        trim  = FALSE,
        scale = "width",
        alpha = 0.5,
        color = NA
    ) +
    geom_boxplot(
        width        = 0.2,
        outlier.size = 1,
        color        = "black",
        alpha        = 0.8
    ) +
    scale_y_log10(
        breaks = c(100, 1e3, 1e4, 1e5, 1e6),
        labels = c("100","1k","10k","100k","1M")
    ) +
    scale_fill_brewer(palette = "Set2") +
    coord_flip() +
    
    # **facet_wrap across variant** on top
    facet_wrap(~ variant, nrow = 1, scales = "free_y", strip.position = "top") +
    
    labs(
        title = "SV Width Distributions by Species & Variant",
        x     = NULL,
        y     = "Width (bp; log scale)"
    ) +
    
    # **remove legend** and style
    theme_minimal(base_size = 14) +
    theme(
        plot.title             = element_text(hjust = 0.5),
        strip.background       = element_rect(fill = "grey90", color = NA),
        strip.text             = element_text(face = "bold"),
        legend.position        = "none",
        panel.grid.major.x     = element_blank(),
        panel.grid.minor       = element_blank(),
        axis.text.y            = element_text(size = 11)
    )


ASV_simple %>%
    mutate(
        species = factor(species, levels = c(
            "Bordetella_pertussis",
            "Helicobacter_pylori",
            "Salmonella_enterica",
            "Staphylococcus_aureus",
            "Mycobacterium_tuberculosis"
        )),
        variant = factor(variant, levels = c(
            "Indel","Duplication","Inversion","Translocation"
        )),
        log_width = log10(width)
    ) %>%
    ggplot(aes(x = species, y = log_width, fill = variant, color = variant)) +
    
    # half‑violins
    stat_halfeye(
        aes(fill = variant),
        position = position_nudge(x = -0.2),
        adjust   = 0.5,
        width    = 0.6,
        justification = -0.2,
        .width   = 0      # draw only the density, no interval shading
    ) +
    
    # point + interval (median and 50% CI by default)
    stat_pointinterval(
        aes(color = variant),
        position = position_nudge(x = 0.2),
        point_size    = 1.2,
        interval_size = 1.2,
        .width        = 0.50  # 50% interval; try c(.50, .95) if you want two
    ) +
    
    # flip & scales
    coord_flip() +
    scale_y_continuous(
        name   = "Width (bp; log10 scale)",
        breaks = log10(c(100, 1e3, 1e4, 1e5, 1e6)),
        labels = c("100","1k","10k","100k","1M")
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    
    # facets by variant
    facet_wrap(~ variant, ncol = 2, scales = "free_x") +
    
    labs(
        title = "Raincloud SV Width Distributions by Species & Variant",
        x     = NULL, 
        y     = NULL,
        fill  = "SV Type",
        color = "SV Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title       = element_text(hjust = 0.5),
        legend.position  = "bottom",
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.major.x = element_blank()
    )



ASV_simple %>%
    mutate(
        species = factor(
            species,
            levels = c(
                "Bordetella_pertussis",
                "Helicobacter_pylori",
                "Salmonella_enterica",
                "Staphylococcus_aureus",
                "Mycobacterium_tuberculosis"
            )
        ),
        variant = factor(
            variant,
            levels = c("Indel","Duplication","Inversion","Translocation")
        ),
        log_width = log10(width)
    ) %>%
    ggplot(aes(x = log_width, y = species, fill = species)) +
    geom_density_ridges(
        scale = 1.0,
        rel_min_height = 0.01,
        alpha = 0.8,
        color = NA
    ) +
    facet_wrap(~ variant, ncol = 2, scales = "free_x") +
    scale_x_continuous(
        name = "Width (bp; log10 scale)",
        breaks = log10(c(100,1e3,1e4,1e5,1e6)),
        labels = c("100","1k","10k","100k","1M")
    ) +
    theme_minimal(base_size = 14) +
    theme(
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        legend.position = "none"
    ) +
    ggtitle("SV Width Distributions by Species — ridgelines per SV type")

ASV_simple2 <- ASV_simple %>%
    mutate(
        species = factor(species, levels = c(
            top_species$Species_[1:5]
        )),
        variant = factor(variant, levels = c(
            "Indel", "Duplication", "Inversion", "Translocation"
        ))
    )

ggplot(ASV_simple2, aes(x = species, y = width, fill = variant)) +
    # violins
    geom_violin(
        position = position_dodge(width = 0.9),
        width    = 0.8,
        trim     = FALSE,
        scale    = "width",
        alpha    = 0.3,
        color    = NA
    ) +
    # boxplots
    geom_boxplot(
        position    = position_dodge(width = 0.9),
        width       = 0.2,
        outlier.size= 1,
        color       = "black",
        alpha       = 0.8
    ) +
    scale_y_log10(
        breaks = c(100, 1e3, 1e4, 1e5, 1e6),
        labels = c("100","1k","10k","100k","1M")
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    coord_flip() +
    labs(
        title = "SV Width Distributions by Species and Variant Type",
        x = NULL,
        y = "Width (bp; log scale)",
        fill = "SV Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title         = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(size = 11),
        legend.position    = "right"
    )


ggplot(ASV_simple, aes(x = species, y = width, fill = variant)) +
    # violins in back
    geom_violin(
        position = position_dodge(width = 0.85),
        width = 0.8,
        alpha = 0.4,
        trim = TRUE
    ) +
    # boxplots on top
    geom_boxplot(
        position = position_dodge(width = 0.85),
        width = 0.2,
        outlier.size = 1,
        alpha = 0.8,
        color = "black"
    ) +
    scale_y_log10(
        breaks = c(100, 1e3, 1e4, 1e5, 1e6),
        labels = c("100","1k","10k","100k","1M")
    ) +
    scale_fill_brewer(palette = "Set2") +
    labs(
        title = "SV Width Distributions by Species and Variant Type",
        x = "Species",
        y = "Width (bp; log scale)",
        fill = "SV Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title             = element_text(hjust = 0.5),
        axis.text.x            = element_text(angle = 30, hjust = 1),
        panel.grid.major.y     = element_line(color = "lightgrey"),
        panel.grid.major.x     = element_blank(),
        panel.grid.minor       = element_blank()
    )
###




###

ggplot(ASV_simple, aes(x = variant, y = width, fill = variant)) +
    geom_boxplot(outlier.size = 1, alpha = 0.8) +
    facet_wrap(~ species, scales = "free_y") +
    scale_y_log10() +    # if widths span orders of magnitude
    labs(
        x = "SV Type",
        y = "Width (bp)",
        title = "SV Width Distributions by Variant and Species"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title   = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 30, hjust = 1),
        legend.position = "none"
    )

ggplot(ASV_simple, aes(x = width, y = variant, fill = variant)) +
    geom_density_ridges(
        scale = 1.1,
        bandwidth = 0.05,    # tweak for peak sharpness
        alpha = 0.8, 
        color = "white"
    ) +
    facet_wrap(~ species, scales = "free_x") +
    scale_x_log10() +
    labs(
        x = "SV Width (bp, log10)",
        y = "SV Type",
        title = "SV Width Densities by Variant and Species"
    ) +
    theme_ridges(font_size = 13) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
    )

ASV_simple %>%
    ggplot(aes(x = variant, y = width, fill = variant)) +
    geom_boxplot(alpha = 0.8) +
    facet_wrap(~ species, scales = "free_y") +
    scale_y_log10() +
    labs(
        title = "SV Width Distributions by Variant and Species",
        x = "SV Type", y = "Width (bp, log scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none"
    )

ggplot(ASV_simple, aes(x = variant, y = width, fill = variant)) +
    geom_violin(trim = TRUE, alpha = 0.8) +
    facet_wrap(~ species, scales = "free_y") +
    scale_y_log10() +
    labs(
        x = "SV Type",
        y = "Width (bp)",
        title = "SV Width Distributions by Variant and Species"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title   = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 30, hjust = 1),
        legend.position = "none"
    )

indel_df <- 
    ASV %>% filter(variant == "Indel") %>%
    filter(!is.na(width), width <5e4) %>%
    dplyr::select(width, variant, species) #%>%
    #filter(!species=="Bordetella_pertussis")
    #filter(species=="Bordetella_pertussis") %>% pull(width)

indel_df <- 
    indel_df %>% mutate(log_len = log10(width))
ggplot(indel_df, aes(x = log_len, y = species, fill = species)) +
    geom_density_ridges(
        scale = 1,
        rel_min_height = 0.01, # tail
        bandwidth = 0.03,
        alpha = 0.8,
        color = "white"
    ) +
    scale_fill_cyclical(
        values = RColorBrewer::brewer.pal(5, "Set2"),
        name = "Species"
    ) +
    labs(
        title = "Indel Length Distributions by Species (log10 scale)",
        x     = "log10(Indel length, bp)",
        y     = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
    )

other_df <- 
    ASV %>% filter(variant != "Indel") %>%
    filter(!is.na(width), width <5e4) %>%
    dplyr::select(width, variant, species) #%>%
#filter(!species=="Bordetella_pertussis")
#filter(species=="Bordetella_pertussis") %>% pull(width)

other_df <- 
    other_df %>% mutate(log_len = log10(width))
ggplot(other_df, aes(x = log_len, y = species, fill = species)) +
    geom_density_ridges(
        scale = 1,
        rel_min_height = 0.01, # tail
        bandwidth = 0.03,
        alpha = 0.8,
        color = "white"
    ) +
    scale_fill_cyclical(
        values = RColorBrewer::brewer.pal(5, "Set2"),
        name = "Species"
    ) +
    labs(
        title = "Indel Length Distributions by Species (log10 scale)",
        x     = "log10(Indel length, bp)",
        y     = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
    )


indel_df %>% filter(species=="Bordetella_pertussis") %>%
    pull(width) %>% density %>% plot
indel_df %>% filter(species=="Mycobacterium_tuberculosis", width<1e4) %>%
    pull(width) %>% density %>% plot

specs <- unique(indel_df$species)
pairs <- combn(specs, 2, simplify = FALSE)






kruskal.test(log_len ~ species, data = indel_df)

ks_results <- 
    map_dfr(pairs, ~{
    s1 <- filter(indel_df, species == .x[1])$log_len
    s2 <- filter(indel_df, species == .x[2])$log_len
    tst <- ks.test(s1, s2)
    tibble(
        sp1 = .x[1], sp2 = .x[2],
        statistic = tst$statistic,
        p.value   = tst$p.value
    )
})




ASV %>%
    count(species, idx) %>%
    group_by(species) %>%
    summarise(mean_hits = mean(n))


ASV %>%
    pull(species) %>%
    unique()

ASV_1

ASV %>% count(variant)
sv_df <- 
    asv_i

plot_mix_by_species <- 
    function(
        sv_df
        ){
    # clean & transform
    lengths     <- 
        sv_df$width[sv_df$width > 0 & !is.na(sv_df$width)]
    log_lengths <- log10(lengths)
    
    # fit mixture model
    mod   <- Mclust(log_lengths)
    means <- mod$parameters$mean
    sds   <- sqrt(mod$parameters$variance$sigmasq)
    wts   <- mod$parameters$pro
    K     <- length(wts)
    
    # grid for density
    x_vals <- seq(min(log_lengths), max(log_lengths), length.out = 500)
    
    # total mixture density
    total_den <- sapply(x_vals, function(x)
        sum(wts * dnorm(x, mean = means, sd = sds))
    )
    total_df <- tibble(x = x_vals, density = total_den)
    
    # per‐component densities
    comps_df <- map_dfr(seq_len(K), function(i){
        tibble(
            x       = x_vals,
            density = wts[i] * dnorm(x_vals, mean = means[i], sd = sds[i]),
            component = factor(i)
        )
    })
    
    cut_pts_log <- 
        (means[-1] + means[-length(means)])/2
    
    cut_pts_bp  <- 
        10^cut_pts_log
    
    
    ggplot(comps_df, aes(x = x, y = density)) +
        geom_line(color = "darkred", size = 0.8) +
        facet_wrap(~ component, ncol = 3, scales = "free_y") +
        labs(
            x = "log10(SV length)",
            y = "Density"
        ) +
        theme_classic(base_size = 12) +
        theme(
            strip.background = element_rect(fill = "grey90", color = NA),
            strip.text = element_text(face = "bold"),
            panel.grid.major = element_line(color = "grey95")
        )
    
    # build plot
    ggplot() +
        # geom_histogram(aes(x = log_lengths, y = ..density..),
        #                bins = 60, fill = "grey85", color = "black", alpha = 0.5) +
        # geom_line(data = total_df, aes(x = x, y = density),
        #           color = "steelblue", size = 1) +
        geom_line(data = comps_df, aes(x = x, y = density, linetype = component),
                  color = "darkred", size = 0.8) +
        labs(
            #title = species,
            x     = "log10(SV length)",
            y     = "Density",
            linetype = "Comp."
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5))
}


ASV %>% filter(species==ASV$species[1]) %>%
    plot_mix_by_species()
aarounts<- 
    ASV %>%group_by(s2) %>%
    count(variant)

asv_counts %>%
    mutate(
        species = fct_reorder(s2, -n, .fun = sum),  # order by total SVs
        variant = factor(variant, levels = c("Indel", "Duplication", "Inversion", "Translocation"))
    ) %>%
    ggplot(aes(x = species, y = n, fill = variant)) +
    geom_col(position = "dodge", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    #scale_y_log10() +
    #facet_wrap(~ variant, scales = "free_y", ncol = 4) +
    labs(
        title = "Structural Variant Counts by Species and Type",
        x = "",
        y = "SV Count",
        fill = "SV Type"
    ) +
    theme_classic(base_size = 18) +
    theme(
        plot.title = element_text(hjust = 0.5),
        #axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )



asv_counts %>%
    mutate(
        species = fct_reorder(s2, -n, .fun = sum),
        variant = factor(variant,
                         levels = c("Indel","Duplication","Inversion","Translocation")),
        panel = if_else(variant == "Indel",
                        "(a) Indels",
                        "(b) Other variants")
    ) %>%
    ggplot(aes(x = species, y = n, fill = variant)) +
    geom_col(position = "dodge", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    facet_wrap(~ panel,
               scales = "free_y",
               ncol   = 2,
               strip.position = "top") +
    labs(
        title = "Structural Variant Counts by Species and Type",
        x     = NULL,
        y     = "SV Count",
        fill  = "SV Type"
    ) +
    theme_classic(base_size = 18) +
    theme(
        plot.title       = element_text(hjust = 0.5),
        strip.text       = element_text(face = "bold", size = 16),
        axis.text.x      = element_text(angle = 30, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )


library(ggforce)
p <- 
    ggplot(asv_counts, aes(s2, n, fill = variant)) +
    geom_col(position = "dodge") +
    labs(
        x     = "Species",
        y     = "SV count",
        fill  = "SV Type",
        title = "SV Counts with Zoomed‐Inset"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p + facet_zoom(ylim = c(0, 500), split = TRUE)
asv_i <- 
    readRDS(asv_specs[4])

library(mclust)
mixmod <- 
    Mclust(
        asv_i %>% 
            filter(!is.na(width)) %>%
            mutate(lg10 = log10(width)) %>%
            pull(lg10)
        )

# Fit mixture model


# Extract model parameters
means <- mixmod$parameters$mean
sds <- sqrt(mixmod$parameters$variance$sigmasq)
weights <- mixmod$parameters$pro
n_components <- length(weights)

components_df <- map_dfr(1:n_components, function(i) {
    tibble(
        x = x_vals,
        density = dnorm(x_vals, mean = means[i], sd = sds[i]) * weights[i],
        component = paste0("Component ", i)
    )
})

# Create a grid of x-values across the data range
x_vals <- seq(min(sv_log), max(sv_log), length.out = 500)

# Manually compute total mixture density
total_density <- rep(0, length(x_vals))
for (i in 1:n_components) {
    total_density <- total_density + weights[i] * dnorm(x_vals, mean = means[i], sd = sds[i])
}



ggplot() +
    # geom_histogram(aes(x = sv_log, y = ..density..),
    #                bins = 60, fill = "grey85", color = "black", alpha = 0.5) +
    # geom_line(data = plot_df, aes(x = x, y = density), 
    #           color = "steelblue", size = 1.2) +
    geom_line(data = components_df, aes(x = x, y = density, color = component), 
              linetype = "dashed", size = 0.9, show.legend = FALSE) +
    geom_vline(data = annot_df, aes(xintercept = x), 
               linetype = "dotted", color = "grey40") +
    geom_text_repel(data = annot_df, aes(x = x, y = 0, label = label), 
                    nudge_y = 0.08, size = 4, direction = "y") +
    labs(
        title = "Mixture Model of log10(SV lengths)",
        x = "log10(SV length)",
        y = "Density"
    ) +
    scale_color_manual(values = brewer.pal(n = length(unique(components_df$component)), name = "Set1")) +
    theme_minimal(base_size = 14)

# Build data frame for ggplot
plot_df <- data.frame(x = x_vals, density = total_density)


components_df$facet <- components_df$component

# For total line (add it across all facets)
total_df <- plot_df %>%
    mutate(component = "Total", facet = "Total")

combined_df <- bind_rows(components_df, total_df)



annot_df <- tibble(
    x = means,
    label = paste0(round(10^means), " bp")
)

ggplot() +
    geom_histogram(aes(x = sv_log, y = ..density..),
                   bins = 60, fill = "grey85", color = "black", alpha = 0.5) +
    geom_line(data = plot_df, aes(x = x, y = density), color = "steelblue", size = 1.2) +
    geom_line(data = components_df, aes(x = x, y = density, color = component), 
              linetype = "dashed", size = 0.9) +
    geom_vline(data = annot_df, aes(xintercept = x), linetype = "dotted", color = "grey40") +
    geom_text_repel(data = annot_df, aes(x = x, y = 0, label = label), nudge_y = 0.1, size = 4) +
    labs(
        title = "Mixture Model with Component Peaks",
        x = "log10(SV length)",
        y = "Density"
    ) +
    theme_minimal(base_size = 14)

ggplot() +
    geom_histogram(aes(x = sv_log, y = ..density..),
                   bins = 60, fill = "grey85", color = "black", alpha = 0.5) +
    geom_line(data = combined_df, aes(x = x, y = density, color = component),
              size = 1) +
    facet_wrap(~ facet, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Set1"), "black")) +
    labs(
        title = "Mixture Components of log10(SV lengths)",
        x = "log10(SV length)",
        y = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")


ggplot() +
    geom_histogram(aes(x = sv_log, y = ..density..), 
                   bins = 60, fill = "grey80", color = "black", alpha = 0.6) +
    geom_line(data = plot_df, aes(x = x, y = density), color = "steelblue", size = 1.2) +
    labs(
        title = "Mixture Model Fit on log10(SV lengths)",
        x = "log10(SV length)",
        y = "Density"
    ) +
    theme_minimal(base_size = 14)

ggplot() +
    geom_histogram(aes(x = sv_log, y = ..density..),
                   bins = 60, fill = "grey90", color = "black", alpha = 0.5) +
    geom_line(data = plot_df, aes(x = x, y = density), 
              color = "steelblue", size = 1.2) +
    geom_line(data = components_df, aes(x = x, y = density, color = component), 
              linetype = "dashed", size = 0.9) +
    labs(
        title = "Mixture Model Fit with Individual Components",
        x = "log10(SV length)",
        y = "Density",
        color = "Mixture Component"
    ) +
    theme_minimal(base_size = 14)






