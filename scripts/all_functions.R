#####  Loading and arranging #####

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


dnaA_from_gff <- 
  function(gff_path, genome_acc){
    #
    gff_gzcon <- 
      gzcon(url(gff_path))
    
    gff_table <- 
      read.gff(
        gff_gzcon
      )
    #gff_gzcon %>% as.character()
    
    close(gff_gzcon)
    
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
          ]# %>% dplyr::slice(1)
        # if (nrow(dnaA_vect)>1){
        #   stop(">1 dnaA")
        }
      if (all(pull(dnaA_vect, strand)=="+")){
        dnaA_vect <-
          dnaA_vect %>%
          filter(start ==min(start)) %>%
          dplyr::slice(1)
      } else if (all(pull(dnaA_vect, strand)=="-")) {
        dnaA_vect <-
          dnaA_vect %>%
          filter(start ==max(start))%>% dplyr::slice(1)
      } else {
        dnaA_vect <-
          dnaA_vect %>% dplyr::slice(1)
      }
    }
    
    if (nrow(dnaA_vect)==0) {
      dnaA_vect <- 
        sprintf(
          "no dnaA annotation",
        )
    }
    
    return(dnaA_vect)
  }

rearrange_to_position <- 
  function(
    seq_file, position, pos_strand
  ){
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
    
  }


dnaA_sync <- 
  function(fasta_path, dnaA_table, out_dir){
    
    genome_seq <- 
      readDNAStringSet(fasta_path)
    
    genome_seqs <- 
      genome_seq[order(width(genome_seq), decreasing = T),]
    
    #showConnections(all = T)
    
    genome_seq <- 
      genome_seq[order(width(genome_seq), decreasing = T),][1]

    
    #names(genome_seqs)
    chrom_num <- 
      length(genome_seqs)
    
    genome_length <- 
      width(genome_seq)
    
    genome_acc <- 
      (names(genome_seq) %>%
         strsplit(., " "))[[1]][1]
    
    ## synchronize the sequen
    gff_ori_strand <- 
      pull(dnaA_table, strand)
    ## syncronize to the dnaA position
    if (gff_ori_strand=="+"){
      ori_start <- 
        pull(dnaA_table, start)
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
        pull(dnaA_table, end) 
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
        out_dir, genome_acc
      )
    
    writeXStringSet(
      out_seq, 
      filepath = out_file
    )
    
    out_tib <- 
      tibble(
        asm_name = 
          basename(dirname(fasta_path)),
        accession = 
          genome_acc,
        ori_pos = 
          dnaA_table$start[1],
        ori_strand = 
          dnaA_table$strand[1],
        len = 
          genome_length,
        out_file = 
          out_file,
        note = NA
      )
    
    
    return(out_tib)
    
  }







##### identities #####
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
  function(species_name, identity_dir, sync_dir){
    sketchfile <- 
      sprintf(
        "mash sketch -p 8 -o %s/mash.msh -s 10000 '%s'/*.txt ",
        identity_dir, sync_dir
      )
    mash_lines <- 
      sprintf(
        "mash dist -p 8 %s/mash.msh %s/mash.msh > %s/%s_mash.txt",
        identity_dir, identity_dir, 
        identity_dir, species_name
      )
    return(capture.output(cat(sketchfile, mash_lines, sep = "\n")))
  }

get_pairwise_identities <- 
  function(
    sync_dir, id_dir,
    species_name, clust_method = "complete"
  ){
    
    mash_file <- 
      sprintf(
        "%s/%s_mash.txt", 
        id_dir, species_name
      )
    # run mash if it hasn't been run
    if (!
        file.exists(
          mash_file
        )
    ){
      mash_bash <- 
        mash_commands(species_name, id_dir, sync_dir)
      system(mash_bash[1], ignore.stderr = T)
      system(mash_bash[2])
    }
    
    
    id_table <-
      fread(
        mash_file,
        select = c(1:3)
      ) %>%
      mutate(V1 = sub(".txt", "", basename(V1))) %>%
      mutate(V2 = sub(".txt", "", basename(V2))) %>%
      `colnames<-`(c("ref", "qry", "identity")) %>%
      mutate(identity = as.numeric(identity)) %>%
      as_tibble %>%
      mutate(
        identity = (1-identity)*100
      )

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
        
        #Recreate the wide matrices
    id_wide <- 
      id_table %>%
      dplyr::select(c(ref, qry, identity)) %>%
      spread(., qry, identity) %>%
      as.data.frame() %>%
      column_to_rownames("ref") %>%
      as.matrix()
    
    id_wide_dist <- 
      100 - id_wide
    
    id_sld_wide <- 
      id_table %>%
      dplyr::select(c(ref, qry, identity)) %>%
      spread(., qry, identity) %>%
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
        data = id_table_long, 
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
        mash_tib = id_table,
        mash_plot = mash_plot,
        id_dists = id_dists
      )
    return(mash_list)
  }

#### creating delta files with MUMmer 

simple_nucmer <- function(
    ref, qry,
    alignment_dir,              # directory for outputs
    output,                     # prefix name (no extension)
    minmatch   = 32,
    nuc_remove = FALSE,         # TRUE = rm unfiltered .delta, FALSE = mv to /unfiltered
    check_nuc  = TRUE           # skip if filtered delta already exists
) {
  # local dirs/paths
  unfiltered_dir <- file.path(alignment_dir, "unfiltered")
  dir.create(unfiltered_dir, showWarnings = FALSE, recursive = TRUE)
  
  prefix <- file.path(alignment_dir, output)
  delta  <- sprintf('%s.delta', prefix)
  filt   <- file.path(alignment_dir, sprintf('%s_filtered.delta', output))
  unf    <- file.path(unfiltered_dir, sprintf('%s.delta', output))
  
  # same params as generate_nucmer_commands (note: --maxmatch)
  nuc_params  <- sprintf('--maxmatch --maxgap=500 --mincluster=100 --minmatch=%s', minmatch)
  mkdir_cmd   <- sprintf('mkdir -p "%s" "%s"', alignment_dir, unfiltered_dir)
  nucmer_cmd  <- sprintf('nucmer %s --prefix="%s" "%s" "%s"', nuc_params, prefix, ref, qry)
  filter_cmd  <- sprintf('delta-filter -r -q "%s" > "%s"', delta, filt)
  post_cmd    <- if (isTRUE(nuc_remove)) sprintf('rm -f "%s"', delta) else sprintf('mv "%s" "%s"', delta, unf)
  
  core <- paste(mkdir_cmd, nucmer_cmd, filter_cmd, post_cmd, sep = " && ")
  
  if (isTRUE(check_nuc)) {
    sprintf('if [ ! -e "%s" ]; then %s ; fi', filt, core)
  } else {
    core
  }
}


####
generate_nucmer_commands <- function(
    genome_matrix, delta_dir, sync_dir,
    minmatch = 32,
    nuc_remove = FALSE,
    check_nuc = TRUE
) {
  # Also create unfiltered dir on the R side (harmless if it already exists)
  unfiltered_dir <- file.path(delta_dir, "unfiltered")
  dir.create(unfiltered_dir, showWarnings = FALSE, recursive = TRUE)
  
  cmds <- apply(genome_matrix, 1, function(x) {
    ref <- x[1]; qry <- x[2]
    
    prefix <- file.path(delta_dir, sprintf("%s_v_%s", ref, qry))
    delta  <- sprintf('%s.delta', prefix)
    filt   <- file.path(delta_dir, sprintf("%s_v_%s_filtered.delta", ref, qry))
    unf    <- file.path(unfiltered_dir, sprintf("%s_v_%s.delta", ref, qry))
    
    nuc_params <- sprintf('--maxmatch --maxgap=500 --mincluster=100 --minmatch=%s', minmatch)
    mkdir_cmd  <- sprintf('mkdir -p "%s" "%s"', delta_dir, unfiltered_dir)
    nucmer_cmd <- sprintf('nucmer %s --prefix="%s" "%s/%s.txt" "%s/%s.txt"',
                          nuc_params, prefix, sync_dir, ref, sync_dir, qry)
    filter_cmd <- sprintf('delta-filter -r -q "%s" > "%s"', delta, filt)
    post_cmd   <- if (isTRUE(nuc_remove)) sprintf('rm -f "%s"', delta) else sprintf('mv "%s" "%s"', delta, unf)
    
    # core pipeline: no leading/trailing semicolons
    core <- paste(mkdir_cmd, nucmer_cmd, filter_cmd, post_cmd, sep = " && ")
    
    if (isTRUE(check_nuc)) {
      sprintf('if [ ! -e "%s" ]; then %s ; fi', filt, core)
    } else {
      core
    }
  })
  
  unname(cmds)
}


### parsing delta files ####
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

#filter_delta(delta_table) %>% plot_delta

filter_delta <- 
  function(
    delta_table, 
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
          rleid(
            strand, qry_gaps_up, qry_gaps_down, ref_gaps
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
        (meanlen > minlen) | 
          (new_contigs==1) | 
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
        XDF_diff2   = cumsum(abs(XDD) >= X_dist_diff), 
        new_contigs = rleid(strand, XDF_diff2) 
      ) %>% 
      ungroup() %>% 
      group_by(new_contigs, strand) %>% 
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
                abs((rs-qs)/sqrt(2)), 
                abs((re-qe)/sqrt(2))
              )
            ), 
            mean(
              c(
                abs((qs + rs - qlen)/sqrt(2)), 
                abs((qe + re - qlen)/sqrt(2))
              )
            )
          )
      ) %>% 
      ungroup
    return(out_tibble)
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

##### Capturing variants #####

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


delta_indels <- function(
    delta_table,
    minlen = NULL,              # NULL = no length filter
    drop_terminal_gaps = TRUE,
    size_bands = c(10, 50),     # thresholds for labeling
    add_size_class = TRUE
){
  if (is.character(delta_table)) delta_table <- read_delta(delta_table)
  
  dt <- delta_table %>% dplyr::mutate(qs2 = pmin(qs, qe), qe2 = pmax(qs, qe))
  qlen <- dt$qlen[1]; rlen <- dt$rlen[1]
  rid1 <- as.character(dt$rid[1]); qid1 <- as.character(dt$qid[1])
  
  # -------- INSERTIONS (query gaps) --------
  qry_cov  <- IRanges::reduce(IRanges::IRanges(dt$qs2, dt$qe2))
  qry_gaps <- IRanges::gaps(qry_cov, start = 1, end = qlen)
  if (drop_terminal_gaps && length(qry_gaps) > 0L) {
    qry_gaps <- qry_gaps[
      IRanges::start(qry_gaps) > 1 &
        IRanges::end(qry_gaps)   < qlen
    ]
  }
  
  insertions_df <- if (length(qry_gaps) > 0L) {
    qry_gr <- GenomicRanges::GRanges(
      seqnames = dt$qid,
      ranges   = IRanges::IRanges(dt$qs2, dt$qe2),
      strand   = dt$strand,
      ref_break = ifelse(dt$strand == "+", dt$re, dt$rs)  # orientation-aware
    )
    ins_gr  <- GenomicRanges::GRanges(qid1, qry_gaps)
    prev_ix <- GenomicRanges::follow(ins_gr, qry_gr)      # block ending before the gap
    ref_pt  <- S4Vectors::mcols(qry_gr)$ref_break[prev_ix]
    
    tibble::tibble(
      rid = rid1, qid = qid1,
      start = as.integer(IRanges::start(ins_gr)),
      end   = as.integer(IRanges::end(ins_gr)),
      width = as.integer(IRanges::width(ins_gr)),
      ref_insertion_point = as.numeric(ref_pt),
      variant = "Insertion"
    )
  } else tibble::tibble()
  
  # -------- DELETIONS (reference gaps) -----
  ref_cov  <- IRanges::reduce(IRanges::IRanges(dt$rs, dt$re))
  ref_gaps <- IRanges::gaps(ref_cov, start = 1, end = rlen)
  if (drop_terminal_gaps && length(ref_gaps) > 0L) {
    ref_gaps <- ref_gaps[
      IRanges::start(ref_gaps) > 1 &
        IRanges::end(ref_gaps)   < rlen
    ]
  }
  
  deletions_df <- if (length(ref_gaps) > 0L) {
    tibble::tibble(
      rid = rid1, qid = qid1,
      start = as.integer(IRanges::start(ref_gaps)),
      end   = as.integer(IRanges::end(ref_gaps)),
      width = as.integer(IRanges::width(ref_gaps)),
      ref_insertion_point = NA_real_,
      variant = "Deletion"
    )
  } else tibble::tibble()
  
  # -------- Combine, guard empty schema ----
  out <- dplyr::bind_rows(insertions_df, deletions_df)
  
  if (nrow(out) == 0L || ncol(out) == 0L) {
    out <- tibble::tibble(
      rid = character(), qid = character(),
      start = integer(),  end = integer(),  width = integer(),
      ref_insertion_point = numeric(),
      variant = character()
    )
    if (isTRUE(add_size_class)) {
      out <- tibble::add_column(out, size_class = character())
    }
    return(out)  # nothing to arrange
  }
  
  if (!is.null(minlen)) out <- dplyr::filter(out, width >= minlen)
  
  if (add_size_class) {
    labs <- c(paste0("<", size_bands[1]),
              paste0(size_bands[1], "–", size_bands[2]-1),
              paste0("≥", size_bands[2]))
    out <- out %>%
      dplyr::mutate(size_class = cut(
        width, breaks = c(-Inf, size_bands, Inf), labels = labs, right = FALSE))
  }
  
  dplyr::arrange(out, start)
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
          qry_name = delta_table$qid[1],
          dup_s = NA_integer_,
          dup_e = NA_integer_,
          orig_ref_start = NA_integer_,
          orig_ref_end = NA_integer_,
          dup_ref_start = NA_integer_,
          dup_ref_end = NA_integer_,
          dup_len = NA_integer_,
          dup_gap = NA_integer_,
          domains = NA_character_,
          duplication_type = NA_character_,
          collapsed = 0,
          duplication_found = FALSE
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
            qry_name = delta_table$qid[1],
            dup_s = NA_integer_,
            dup_e = NA_integer_,
            orig_ref_start = NA_integer_,
            orig_ref_end = NA_integer_,
            dup_ref_start = NA_integer_,
            dup_ref_end = NA_integer_,
            dup_len = NA_integer_,
            dup_gap = NA_integer_,
            domains = NA_character_,
            duplication_type = NA_character_,
            collapsed = 0,
            duplication_found = FALSE
          )
        )
        #)
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


delta_structural_rearrangements <- 
  function(delta_table, localized_threshold = 1e6) {
  if (is.character(delta_table)) {
    delta_table <- read_delta(delta_table)
  }
  fd <- delta_table %>%
    filter_delta() %>%
    rowwise() %>%
    mutate(
      qs2 = min(qs, qe),
      qe2 = max(qs, qe)
    ) %>%
    ungroup()
  
  fd2 <- fd %>%
    arrange(rs) %>%
    mutate(idx = row_number())
  
  # Find all inversion indices (strand == "-")
  inv_idx <- fd2$idx[fd2$strand == "-"]
  
  # Classify as localized or symmetric based on X_dist
  localized_idx <- fd2$idx[fd2$strand == "-" & fd2$X_dist > localized_threshold]
  symmetric_idx <- fd2$idx[fd2$strand == "-" & fd2$X_dist <= localized_threshold]
  
  events <- list()
  
  # Localized inversion(s)
  if (length(localized_idx) > 0) {
    events$localized <- fd2 %>%
      filter(idx %in% localized_idx) %>%
      summarise(
        variant_specific = "localized",
        rs      = min(rs),
        re      = max(re),
        qs      = min(qs2),
        qe      = max(qe2),
        contigs = paste(new_contigs, collapse = ","),
        strand  = unique(strand),
        X_dist  = mean(X_dist),
        width   = (re - rs + 1),
        rid     = unique(rid),
        qid     = unique(qid)
      )
  }
  
  # Symmetric inversion(s)
  if (length(symmetric_idx) > 0) {
    events$symmetric <- fd2 %>%
      filter(idx %in% symmetric_idx) %>%
      summarise(
        variant_specific = "symmetric",
        rs      = min(rs),
        re      = max(re),
        qs      = min(qs2),
        qe      = max(qe2),
        contigs = paste(new_contigs, collapse = ","),
        strand  = unique(strand),
        X_dist  = mean(X_dist),
        width   = (re - rs + 1),
        rid     = unique(rid),
        qid     = unique(qid)
      )
  }
  
  # Optionally add other event types here (nested, outer, etc.) if you wish
  
  events_df <- bind_rows(events)
  
  # If nothing found, return "none"
  if (nrow(events_df) == 0) {
    return(
      tibble(
        variant_specific = "none",
        rs      = fd2$rs[1],
        re      = fd2$re[1],
        qs      = fd2$qs2[1],
        qe      = fd2$qe2[1],
        contigs = as.character(fd2$new_contigs[1]),
        strand  = fd2$strand[1],
        X_dist  = fd2$X_dist[1],
        width   = fd2$re[1] - fd2$rs[1] + 1,
        rid     = fd2$rid[1],
        qid     = fd2$qid[1]
      )
    )
  }
  return(events_df)
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



delta_all_SVs <- 
  function(delta_file, unfiltered_delta_file, species_name){
    
    ### Indels
    delta_indels_i <- 
      delta_indels(delta_file) %>%
      transmute(
        rid, qid, width, 
        start, end,
        refmid = round((start + end) / 2),
        variant_specific = variant,
        variant = "Indel",
        species = species_name
      )
    
    if (nrow(delta_indels_i) == 0){
      delta_indels_i <- tibble(
        rid = character(),
        qid = character(),
        width = numeric(),
        start = numeric(),
        end = numeric(),
        refmid = numeric(),
        variant_specific = character(),
        variant = character(),
        species = character()
      )
    }
    
    ### Duplications
    delta_dups_i <- 
      delta_duplications(
        delta_table = delta_file, 
        unfiltered_delta_table = unfiltered_delta_file
      ) %>% 
      transmute(
        rid = ref_name, qid = qry_name,
        start = dup_ref_start,
        end = dup_ref_end,
        width = abs(end - start),
        refmid = round((start + end) / 2),
        variant_specific = duplication_type,
        variant = "Duplication",
        species = species_name
      )
    
    if (nrow(delta_dups_i) == 0){
      delta_dups_i <- delta_indels_i[0,]
    }
    
    ### Translocations
    delta_tlocs_i <- delta_translocations(delta_file) %>%
      transmute(
        rid, qid, width,
        start, end,
        refmid = round((start + end) / 2),
        variant_specific = replichores,
        variant = "Translocation",
        species = species_name
      )
    
    if (nrow(delta_tlocs_i) == 0){
      delta_tlocs_i <- delta_indels_i[0,]
    }
    
    ### Structural rearrangements
    delta_SR_temp <- delta_structural_rearrangements(delta_file)
    
    delta_subSR_i <- delta_substructural_inversions(delta_file, delta_SR_temp) %>%
      ungroup() %>%
      transmute(
        rid, qid,
        width = meanlen,
        start = rs, end = re,
        refmid = round((rs + re) / 2),
        variant_specific = variant,
        variant = "Inversion",
        species = species_name
      )
    
    if (nrow(delta_subSR_i) == 0){
      delta_subSR_i <- delta_indels_i[0,]
    }
    
    delta_SR_i <- delta_SR_temp %>%
      transmute(
        rid, qid,
        width,
        start = rs,
        end = re,
        refmid = round((rs + re) / 2),
        variant = "Structural rearrangement",
        variant_specific,
        species = species_name
      )
    
    if (nrow(delta_SR_i) == 0){
      delta_SR_i <- delta_indels_i[0,]
    }
    
    ### Combine all SVs
    all_SVs <- bind_rows(
      delta_indels_i,
      delta_dups_i,
      delta_tlocs_i,
      delta_subSR_i,
      delta_SR_i
    ) %>% arrange(desc(width))
    
    return(all_SVs)
  }


benchmark_alignments <- 
  function(genome_matrix, delta_dir, sync_dir, strategy_name, mc.cores = 6) {
  

  message("Starting benchmarking: ", strategy_name)

  # Precompute all commands
  cmds <- 
    generate_nucmer_commands(
      genome_matrix, delta_dir, sync_dir, nuc_remove = F, check_nuc = F
      )

  results <- 
    pbmclapply(
      1:nrow(genome_matrix),
      function(i) {
        ref <- genome_matrix[i,1]
        qry <- genome_matrix[i,2]
        pair_id <- paste(ref, qry, sep = "_vs_")
        start_time <- Sys.time()
  
        res <- tryCatch({
          cmd <- cmds[i]
          exit_code <- system(cmd)
          end_time <- Sys.time()
          tibble(
            strategy = strategy_name,
            pair = pair_id,
            ref = ref,
            qry = qry,
            runtime_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
            exit_code = exit_code
          )
        }, error = function(e) {
          message("Error on pair: ", pair_id, " — ", e$message)
          tibble(
            strategy = strategy_name,
            pair = pair_id,
            ref = ref,
            qry = qry,
            runtime_sec = NA_real_,
            exit_code = NA_integer_
          )
        })
  
        res
    },
    mc.cores = mc.cores
  )

  # Filter out any unexpected NULLs
  results <- results[!sapply(results, is.null)]

  # Force all elements to be tibbles
  results <- lapply(results, function(x) {
    if (!inherits(x, "data.frame")) {
      tibble(
        strategy = strategy_name,
        pair = NA_character_,
        ref = NA_character_,
        qry = NA_character_,
        runtime_sec = NA_real_,
        exit_code = NA_integer_
      )
    } else {
      x
    }
  })

  bind_rows(results)
}

#library(pracma)







get_all_SVs <- 
  function(alignment_directory, species_name) {
    
    all_delta_filt_files <- 
      list.files(
        alignment_directory,
        full.names = TRUE, pattern = "filtered.delta"
      )
    
    all_delta_unfilt_files <- 
      list.files(
        sprintf("%s/unfiltered", alignment_directory),
        full.names = TRUE
      )
    
    filt_basenames  <- basename(all_delta_filt_files)
    filt_core_names <- sub("_filtered\\.delta$", "", filt_basenames)
    
    unfilt_basenames  <- basename(all_delta_unfilt_files)
    unfilt_core_names <- sub("\\.delta$", "", unfilt_basenames)
    
    unfilt_lookup <- setNames(all_delta_unfilt_files, unfilt_core_names)
    
    n <- length(all_delta_filt_files)
    all_res <- vector("list", n)
    
    pb <- txtProgressBar(min = 1, max = n, style = 3)
    
    for (i in 1:n) {
      core_name <- filt_core_names[i]
      filtered_delta <- all_delta_filt_files[i]
      
      if (!core_name %in% names(unfilt_lookup)) {
        warning(
          sprintf(
            "No unfiltered file found for %s, skipping...", 
            filtered_delta
          )
        )
        next
      }
      unfiltered_delta <- unfilt_lookup[[core_name]]
      
      all_res[[i]] <- 
        delta_all_SVs(
          delta_file = filtered_delta, 
          unfiltered_delta_file = unfiltered_delta,
          species_name = species_name
        )
      
      setTxtProgressBar(pb, i)
    }
    
    close(pb)
    
    return(bind_rows(Filter(Negate(is.null), all_res)))
  }



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



recursive_gmm <- function(
    x,
    max_depth = 5,
    depth = 0,
    min_n = 100,
    G_range = 1:6,
    proportion_cutoff = 0.1,
    spread_ratio_cutoff = 1.0,
    modelNames = "V"
) {
  # ---- input hygiene ---------------------------------------------------------
  x <- x[is.finite(x)]
  if (!is.numeric(x)) stop("`x` must be a numeric vector of log10-lengths.", call. = FALSE)
  if (length(x) < min_n || depth >= max_depth) {
    return(tibble::tibble())
  }
  
  # ---- fit mixture -----------------------------------------------------------
  fit <- try(mclust::Mclust(x, G = G_range, modelNames = modelNames), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$parameters$mean)) {
    return(tibble::tibble())
  }
  
  means <- as.numeric(fit$parameters$mean)
  # for spherical models, 'sigmasq' is a scalar per component
  sigsq <- fit$parameters$variance$sigmasq
  sds   <- sqrt(as.numeric(sigsq))
  props <- as.numeric(fit$parameters$pro)
  
  # guard lengths
  k <- length(means)
  if (any(lengths(list(sds, props)) != k)) {
    return(tibble::tibble())
  }
  
  mean_bp <- 10^means
  sd_bp   <- 10^(means + sds) - mean_bp  # +1σ offset in bp (for banding)
  
  summary_table <- tibble::tibble(
    depth      = depth,
    submode    = paste0("D", depth, ".", seq_len(k)),
    mean_log10 = means,
    mean_bp    = mean_bp,
    sd_log10   = sds,
    sd_bp      = sd_bp,
    proportion = props,
    n          = length(x)
  )
  
  # ---- decide which components to recurse into -------------------------------
  nested_tables <- list(summary_table)
  
  for (i in seq_len(k)) {
    prop_i <- props[i]
    spread_ratio_i <- ifelse(mean_bp[i] > 0, sd_bp[i] / mean_bp[i], 0)
    
    if (is.finite(prop_i) && is.finite(spread_ratio_i) &&
        prop_i > proportion_cutoff && spread_ratio_i > spread_ratio_cutoff) {
      
      idx <- abs(x - means[i]) < 2 * sds[i]  # ±2 SD window in log10 space
      if (sum(idx) >= min_n) {
        nested_result <- recursive_gmm(
          x = x[idx],
          max_depth = max_depth,
          depth = depth + 1,
          min_n = min_n,
          G_range = if (length(G_range) > 0) pmax(1, pmin(3, G_range)) else 1:3,
          proportion_cutoff = proportion_cutoff,
          spread_ratio_cutoff = spread_ratio_cutoff,
          modelNames = modelNames
        )
        if (nrow(nested_result) > 0) nested_tables <- append(nested_tables, list(nested_result))
      }
    }
  }
  
  dplyr::bind_rows(nested_tables)
}



### SV annotation:

### prokka gff handling
library("rtracklayer")

`%||%` <- function(a,b) ifelse(is.null(a), b, a)



coalesce_chr <- function(...) {
  x <- list(...)
  out <- x[[1]]
  for (i in seq_along(x)) out <- dplyr::coalesce(out, x[[i]])
  out
}

read_prokka_gff2 <- function(path,
                             keep_types = c("gene","CDS","rRNA","tRNA","tmRNA")) {
  gr <- rtracklayer::import(path)                 # auto-detects GFF3
  if (!is.null(keep_types)) {
    has_type <- "type" %in% names(mcols(gr))
    if (has_type) gr <- gr[mcols(gr)$type %in% keep_types]
  }
  
  df <- as_tibble(as.data.frame(gr))
  # Ensure expected columns exist
  for (nm in c("ID","Parent","gene","Name","locus_tag","product","phase","type","source")) {
    if (!nm %in% names(df)) df[[nm]] <- NA
  }
  
  features <- df %>%
    transmute(
      seqid     = as.character(seqnames),
      source    = as.character(source),
      type      = as.character(type),
      start     = as.integer(start),
      end       = as.integer(end),
      strand    = as.character(strand),
      phase     = suppressWarnings(as.integer(phase)),
      ID        = as.character(ID),
      Parent    = as.character(Parent),
      # prefer gene, else Name, else locus_tag/product as a last resort
      gene      = coalesce_chr(as.character(gene), as.character(Name)),
      locus_tag = as.character(locus_tag),
      product   = as.character(product)
    )
  
  # --- collapse CDS fragments to one row per locus_tag (gene-level) ---
  cds_collapsed <- features %>%
    filter(type == "CDS") %>%
    group_by(seqid, locus_tag) %>%
    summarise(
      start     = min(start, na.rm = TRUE),
      end       = max(end, na.rm = TRUE),
      strand    = dplyr::first(na.omit(strand)),
      n_parts   = dplyr::n(),
      gene      = dplyr::first(na.omit(gene)),
      product   = dplyr::first(na.omit(product)),
      .groups   = "drop"
    ) %>%
    mutate(length = end - start + 1L)
  
  list(features = features, cds = cds_collapsed)
}









annotate_SVs2 <- function(SV_granges, gene_granges, overlap_perc = 80) {
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



annotate_SVs3 <- function(
    SV_granges, gene_granges,
    sv_thresh   = 80,
    gene_thresh = 80,
    mode        = c("both","sv","gene","either","jaccard"),
    jaccard_thresh = 0.5,
    keep        = c("best","all")   # best = highest bp overlap per SV
) {
  mode <- match.arg(mode)
  keep <- match.arg(keep)
  
  # --- 0) Align seqlevels
  common_seq <- intersect(seqlevels(SV_granges), seqlevels(gene_granges))
  SVg   <- keepSeqlevels(SV_granges, common_seq, pruning.mode = "coarse")
  Geneg <- keepSeqlevels(gene_granges, common_seq, pruning.mode = "coarse")
  
  # Safely access mcols (fallback to NA if missing)
  get_col <- function(gr, nm) if (nm %in% names(mcols(gr))) mcols(gr)[[nm]] else rep(NA_character_, length(gr))
  sv_species <- get_col(SVg, "species")
  sv_variant <- get_col(SVg, "variant")
  sv_refmid  <- get_col(SVg, "refmid")
  gene_gene  <- get_col(Geneg, "gene")
  gene_prod  <- get_col(Geneg, "product")
  
  # --- 1) Overlaps
  hits <- findOverlaps(SVg, Geneg, ignore.strand = TRUE)
  if (length(hits) == 0L) {
    # All SVs are intergenic
    out <- tibble::tibble(
      gene_accession = NA_character_,
      gene           = NA_character_,
      product        = NA_character_,
      sv_accession   = as.character(GenomeInfoDb::seqnames(SVg)),
      species        = sv_species,
      variant        = sv_variant,
      sv_start       = BiocGenerics::start(SVg),
      sv_end         = BiocGenerics::end(SVg),
      sv_width       = IRanges::width(SVg),
      refmid         = sv_refmid,
      overlap_bp     = 0L,
      pct_overlap_sv   = 0,
      pct_overlap_gene = 0,
      jaccard          = 0,
      func_group       = "intergenic"
    )
    return(out)
  }
  
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  
  # precompute lengths
  sv_w   <- IRanges::width(SVg)
  gene_w <- IRanges::width(Geneg)
  
  # overlap widths
  ovl     <- IRanges::pintersect(SVg[qh], Geneg[sh])
  ov_bp   <- IRanges::width(ovl)
  pct_sv  <- 100 * ov_bp / sv_w[qh]
  pct_gene<- 100 * ov_bp / gene_w[sh]
  jacc    <- ov_bp / (sv_w[qh] + gene_w[sh] - ov_bp)
  
  # --- 2) Mode-specific filter (decides "fully genic" vs "partially genic")
  fully_genic <- switch(
    mode,
    "both"   = (pct_sv >= sv_thresh) & (pct_gene >= gene_thresh),
    "sv"     = (pct_sv >= sv_thresh),
    "gene"   = (pct_gene >= gene_thresh),
    "either" = (pct_sv >= sv_thresh) | (pct_gene >= gene_thresh),
    "jaccard"= (jacc >= jaccard_thresh)
  )
  
  # --- 3) Functional classifier
  classify_func <- function(g, p, fully) {
    if (!fully) return("partially genic")
    if (is.na(p) && is.na(g)) return(NA_character_)
    # NB: order matters (specific -> general)
    if (grepl("transposase", p, TRUE))                                      return("Transposase")
    if (grepl("hypothetical", p, TRUE))                                     return("Hypothetical protein")
    if (grepl("\\bputative\\b", p, TRUE))                                   return("Putative protein")
    if (grepl("^23S\\s+ribosomal\\s+RNA", p, TRUE))                         return("rRNA")
    if (grepl("ribosomal\\s+RNA", p, TRUE))                                 return("rRNA")
    if (grepl("\\badhesin\\b", p, TRUE))                                    return("Adhesin")
    if (grepl("flagellin|flagellar", p, TRUE))                              return("Flagellar proteins")
    if (grepl("restriction enzyme|EcoKI specificity", p, TRUE))             return("Restriction‑modification")
    if (grepl("phosphomanno|phosphogluco|malic enzyme|amidohydrolase|\\bhydrolase\\b|\\bmutase\\b", p, TRUE))
      return("Metabolic enzyme")
    if (grepl("\\bchaperone\\b", p, TRUE))                                  return("Chaperone")
    if (grepl("\\bcompetence protein\\b", p, TRUE))                         return("Competence")
    if (grepl("\\b(toxin|antitoxin)\\b", p, TRUE))                          return("Toxin‑antitoxin system")
    if (grepl("\\bphage\\b|\\bcapsid\\b|tail( fiber)?\\b", p, TRUE))        return("Phage‑associated")
    if (grepl("ESAT-6|\\bEsx\\b|type\\s*VII|\\bESX-", p, TRUE))             return("ESX secretion family")
    if (grepl("\\bABC transporter\\b", p, TRUE))                             return("ABC transporter")
    if (grepl("\\btransporter\\b|efflux pump|\\buptake\\b", p, TRUE))       return("Transporter")
    if (grepl("\\bkinase\\b", p, TRUE))                                     return("Signaling")
    if (grepl("dehydrogenase|\\bligase\\b|reductase|synthase|oxygenase", p, TRUE))
      return("Enzyme")
    if (grepl("ribosomal protein", p, TRUE))                                return("Ribosomal protein")
    if (grepl("catalase-peroxidase", p, TRUE))                              return("Oxidative stress")
    if (grepl("partition protein|\\bSmc\\b", p, TRUE))                      return("Chromosome partition")
    if (grepl("\\bcoagulase\\b", p, TRUE))                                  return("Coagulase")
    if (!is.na(g) && g %in% c("rho","nusA"))                                return("Regulatory")
    if (!is.na(g) && g %in% c("rpfA","rpfE"))                               return("Dormancy/Resuscitation")
    if (grepl("repeat(-|_)?containing|\\brepeat\\b", p, TRUE))              return("Adhesin")
    if (grepl("\\b\\w+ase\\b", p, TRUE))                                    return("Enzyme")
    "Other"
  }
  
  func_group <- mapply(classify_func, gene_gene[sh], gene_prod[sh], fully_genic, USE.NAMES = FALSE)
  
  # For partially genic rows, you might prefer to KEEP gene/product;
  # if you want them blanked, set them to NA where !fully_genic:
  gene_out <- ifelse(fully_genic, gene_gene[sh], NA_character_)
  prod_out <- ifelse(fully_genic, gene_prod[sh], NA_character_)
  
  overl_df <- tibble::tibble(
    gene_accession   = as.character(GenomeInfoDb::seqnames(Geneg))[sh],
    gene             = gene_out,
    product          = prod_out,
    sv_accession     = as.character(GenomeInfoDb::seqnames(SVg))[qh],
    species          = sv_species[qh],
    variant          = sv_variant[qh],
    sv_start         = BiocGenerics::start(SVg)[qh],
    sv_end           = BiocGenerics::end(SVg)[qh],
    sv_width         = sv_w[qh],
    refmid           = sv_refmid[qh],
    overlap_bp       = as.integer(ov_bp),
    pct_overlap_sv   = as.numeric(pct_sv),
    pct_overlap_gene = as.numeric(pct_gene),
    jaccard          = as.numeric(jacc),
    func_group       = func_group
  )
  
  # --- 4) Add intergenic rows for SVs with no hits
  hit_sv_idx <- unique(qh)
  all_sv_idx <- seq_along(SVg)
  nohit_idx  <- setdiff(all_sv_idx, hit_sv_idx)
  
  intergenic_df <- if (length(nohit_idx)) {
    tibble::tibble(
      gene_accession   = NA_character_,
      gene             = NA_character_,
      product          = NA_character_,
      sv_accession     = as.character(GenomeInfoDb::seqnames(SVg))[nohit_idx],
      species          = sv_species[nohit_idx],
      variant          = sv_variant[nohit_idx],
      sv_start         = BiocGenerics::start(SVg)[nohit_idx],
      sv_end           = BiocGenerics::end(SVg)[nohit_idx],
      sv_width         = sv_w[nohit_idx],
      refmid           = sv_refmid[nohit_idx],
      overlap_bp       = 0L,
      pct_overlap_sv   = 0,
      pct_overlap_gene = 0,
      jaccard          = 0,
      func_group       = "intergenic"
    )
  } else {
    tibble::tibble()
  }
  
  all_df <- dplyr::bind_rows(overl_df, intergenic_df)
  
  # --- 5) Best vs all overlaps
  if (keep == "best") {
    # keep the overlap with max bp per SV (ties broken arbitrarily)
    all_df <- all_df |>
      dplyr::group_by(sv_accession, sv_start, sv_end, variant) |>
      dplyr::slice_max(order_by = overlap_bp, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::arrange(sv_accession, sv_start)
  } else {
    all_df <- dplyr::arrange(all_df, sv_accession, dplyr::desc(overlap_bp))
  }
  
  all_df
}


# Improved annotate_SVs4 function with fixes and enhancements







annotate_SVs4 <- 
  function(
    SV_granges, 
    gene_granges,
    sv_thresh = 80, 
    gene_thresh = 80,
    mode = c("both","sv","gene","either","jaccard"),
    jaccard_thresh = 0.5,
    keep = c("best","all"),
    bp_pad = 0L,
    inv_like = c("inversion","local inversion","rearrangement","translocation","RAI"),
    verbose = FALSE) {
  
  mode <- match.arg(mode)
  keep <- match.arg(keep)
  
  # Input validation
  if (!is(SV_granges, "GRanges")) stop("SV_granges must be a GRanges object")
  if (!is(gene_granges, "GRanges")) stop("gene_granges must be a GRanges object")
  if (length(SV_granges) == 0) return(tibble())
  if (length(gene_granges) == 0) {
    warning("No genes provided - all SVs will be marked as intergenic")
    return(tibble(
      hit_type = "none",
      sv_accession = as.character(seqnames(SV_granges)),
      sv_start = start(SV_granges),
      sv_end = end(SV_granges),
      sv_width = width(SV_granges),
      func_group = "intergenic"
    ))
  }
  
  # --- 0) Align seqlevels safely
  common_seq <- intersect(seqlevels(SV_granges), seqlevels(gene_granges))
  
  if (length(common_seq) == 0) {
    warning("No common sequence names between SVs and genes")
    return(tibble())
  }
  
  SVg <- keepSeqlevels(SV_granges, common_seq, pruning.mode = "coarse")
  Geneg <- keepSeqlevels(gene_granges, common_seq, pruning.mode = "coarse")
  
  if (verbose) {
    cat("Processing", length(SVg), "SVs against", length(Geneg), "genes\n")
    cat("Common sequences:", length(common_seq), "\n")
  }
  
  # Helper function with better error handling
  get_col <- function(gr, nm) {
    if (nm %in% names(mcols(gr))) {
      mcols(gr)[[nm]]
    } else {
      rep(NA_character_, length(gr))
    }
  }
  
  # Extract metadata
  sv_species <- get_col(SVg, "species")
  sv_variant <- {
    v <- get_col(SVg, "variant")
    if (all(is.na(v))) get_col(SVg, "variant_specific") else v
  }
  sv_refmid <- get_col(SVg, "refmid")
  gene_gene <- get_col(Geneg, "gene")
  gene_prod <- get_col(Geneg, "product")
  
  # Precompute widths
  sv_w <- width(SVg)
  gene_w <- width(Geneg)
  
  # --- 1A) SPAN overlaps
  span_hits <- findOverlaps(SVg, Geneg, ignore.strand = TRUE)
  
  if (length(span_hits) > 0) {
    qh <- queryHits(span_hits)
    sh <- subjectHits(span_hits)
    span_ov <- pintersect(SVg[qh], Geneg[sh])
    span_bp <- width(span_ov)
    pct_sv <- 100 * span_bp / sv_w[qh]
    pct_gene <- 100 * span_bp / gene_w[sh]
    jacc <- span_bp / (sv_w[qh] + gene_w[sh] - span_bp)
    
    fully_genic <- switch(
      mode,
      "both"    = (pct_sv >= sv_thresh) & (pct_gene >= gene_thresh),
      "sv"      = (pct_sv >= sv_thresh),
      "gene"    = (pct_gene >= gene_thresh),
      "either"  = (pct_sv >= sv_thresh) | (pct_gene >= gene_thresh),
      "jaccard" = (jacc >= jaccard_thresh)
    )
    
    span_df <- tibble(
      hit_type         = "span",
      gene_accession   = as.character(seqnames(Geneg))[sh],
      gene             = ifelse(fully_genic, gene_gene[sh], NA_character_),
      product          = ifelse(fully_genic, gene_prod[sh], NA_character_),
      sv_accession     = as.character(seqnames(SVg))[qh],
      species          = sv_species[qh],
      variant          = sv_variant[qh],
      sv_start         = start(SVg)[qh],
      sv_end           = end(SVg)[qh],
      sv_width         = sv_w[qh],
      refmid           = sv_refmid[qh],
      overlap_bp       = as.integer(span_bp),
      pct_overlap_sv   = as.numeric(pct_sv),
      pct_overlap_gene = as.numeric(pct_gene),
      jaccard          = as.numeric(jacc),
      fully_genic      = fully_genic
    )
  } else {
    span_df <- tibble()
  }
  
  # --- 1B) BREAKPOINT overlaps
  if (bp_pad < 0L) bp_pad <- 0L
  
  # Improved breakpoint range creation
  if (bp_pad > 0 || length(grep(paste(tolower(inv_like), collapse="|"), 
                                tolower(sv_variant), ignore.case=TRUE)) > 0) {
    
    # Create breakpoint ranges
    bp_starts <- GRanges(
      seqnames = seqnames(SVg),
      ranges = IRanges(start = pmax(1, start(SVg) - bp_pad),
                       end = start(SVg) + bp_pad)
    )
    bp_ends <- GRanges(
      seqnames = seqnames(SVg),
      ranges = IRanges(start = pmax(1, end(SVg) - bp_pad),
                       end = end(SVg) + bp_pad)
    )
    
    # Combine and add SV index
    BPg <- c(bp_starts, bp_ends)
    mcols(BPg)$sv_idx <- c(seq_along(SVg), seq_along(SVg))
    mcols(BPg)$bp_side <- c(rep("start", length(SVg)), rep("end", length(SVg)))
    
    bp_hits <- findOverlaps(BPg, Geneg, ignore.strand = TRUE)
    
    if (length(bp_hits) > 0) {
      bq_idx <- mcols(BPg)$sv_idx[queryHits(bp_hits)]
      bs <- subjectHits(bp_hits)
      bp_side <- mcols(BPg)$bp_side[queryHits(bp_hits)]
      
      bp_ov <- pintersect(BPg[queryHits(bp_hits)], Geneg[bs])
      bp_bp <- width(bp_ov)
      
      bp_df <- tibble(
        hit_type         = "breakpoint",
        gene_accession   = as.character(seqnames(Geneg))[bs],
        gene             = gene_gene[bs],
        product          = gene_prod[bs],
        sv_accession     = as.character(seqnames(SVg))[bq_idx],
        species          = sv_species[bq_idx],
        variant          = sv_variant[bq_idx],
        sv_start         = start(SVg)[bq_idx],
        sv_end           = end(SVg)[bq_idx],
        sv_width         = sv_w[bq_idx],
        refmid           = sv_refmid[bq_idx],
        overlap_bp       = as.integer(bp_bp),
        pct_overlap_sv   = NA_real_,
        pct_overlap_gene = 100 * as.numeric(bp_bp) / gene_w[bs],
        jaccard          = NA_real_,
        fully_genic      = TRUE,
        bp_side          = bp_side  # Added: which breakpoint
      )
    } else {
      bp_df <- tibble()
    }
  } else {
    bp_df <- tibble()
  }
  
  # --- 2) Enhanced functional classifier
  classify_func <- function(g, p, fully) {
    if (!isTRUE(fully)) return("partially genic")
    if (is.na(p) && is.na(g)) return(NA_character_)
    
    # Convert to lowercase for matching
    p_lower <- tolower(as.character(p))
    g_lower <- tolower(as.character(g))
    
    # Transposable elements (highest priority)
    if (grepl("transpos|\\bis\\d+|\\btn\\d+|insertion sequence", p_lower)) return("Transposase")
    
    # Phage-related
    if (grepl("phage|prophage|capsid|tail|terminase|integrase|portal", p_lower)) return("Phage-associated")
    
    # Hypothetical/Unknown
    if (grepl("hypothetical|uncharacterized|unknown function", p_lower)) return("Hypothetical protein")
    if (grepl("\\bputative\\b|probable|predicted", p_lower)) return("Putative protein")
    
    # RNA genes
    if (grepl("ribosomal\\s+rna|\\brrna\\b|23s|16s|5s", p_lower)) return("rRNA")
    if (grepl("\\btrna\\b|transfer rna", p_lower)) return("tRNA")
    
    # Specific systems
    if (grepl("esx|esat-6|type vii secretion", p_lower)) return("ESX secretion family")
    if (grepl("toxin|antitoxin|\\bta\\b", p_lower)) return("Toxin-antitoxin system")
    if (grepl("restriction|modification|methylase", p_lower)) return("Restriction-modification")
    
    # Transporters
    if (grepl("abc transporter|atp-binding cassette", p_lower)) return("ABC transporter")
    if (grepl("transporter|permease|efflux|uptake|channel", p_lower)) return("Transporter")
    
    # Enzymes (more specific)
    if (grepl("polymerase|helicase|recombinase|topoisomerase", p_lower)) return("DNA metabolism")
    if (grepl("dehydrogenase|reductase|oxidase|catalase|peroxidase", p_lower)) return("Oxidoreductase")
    if (grepl("synthase|synthetase|ligase", p_lower)) return("Synthase/Ligase")
    if (grepl("hydrolase|peptidase|protease|nuclease", p_lower)) return("Hydrolase")
    if (grepl("kinase|phosphatase", p_lower)) return("Signaling")
    if (grepl("transferase|aminotransferase", p_lower)) return("Transferase")
    
    # Structural/Functional
    if (grepl("ribosomal protein|\\brps|\\brpl", p_lower)) return("Ribosomal protein")
    if (grepl("chaperone|heat shock|\\bhsp|\\bdnak", p_lower)) return("Chaperone")
    if (grepl("flagell|pili|fimbr|adhesin", p_lower)) return("Surface structure")
    
    # Regulatory
    if (grepl("regulator|repressor|activator|\\bfur\\b|\\blexr", p_lower)) return("Regulatory")
    
    # Gene-specific classifications
    if (!is.na(g)) {
      if (g_lower %in% c("rho", "nusa", "nusg")) return("Transcription")
      if (g_lower %in% c("rpfa", "rpfb", "rpfc", "rpfd", "rpfe")) return("Resuscitation")
      if (grepl("^pe\\d+|^ppe\\d+", g_lower)) return("PE/PPE family")
    }
    
    # Catch-all for other enzymes
    if (grepl("\\base$|\\base\\b", p_lower)) return("Enzyme")
    
    "Other"
  }
  
  # Combine span and breakpoint results
  both_df <- bind_rows(span_df, bp_df)
  
  if (nrow(both_df) > 0) {
    # Variant-aware prioritization
    inv_mask <- !is.na(both_df$variant) & 
      tolower(both_df$variant) %in% tolower(inv_like)
    
    both_df$hit_priority <- case_when(
      inv_mask & both_df$hit_type == "breakpoint" ~ 3L,
      both_df$hit_type == "span" & both_df$fully_genic ~ 2L,
      TRUE ~ 1L
    )
    
    # Add functional classification
    both_df$func_group <- mapply(classify_func, 
                                 both_df$gene, 
                                 both_df$product, 
                                 both_df$fully_genic, 
                                 USE.NAMES = FALSE)
  }
  
  # --- 3) Add intergenic rows for SVs with no hits
  if (nrow(both_df) > 0) {
    sv_key <- paste(as.character(seqnames(SVg)), start(SVg), end(SVg), sep = "_")
    hit_key <- paste(both_df$sv_accession, both_df$sv_start, both_df$sv_end, sep = "_")
    nohit_idx <- which(!sv_key %in% hit_key)
  } else {
    nohit_idx <- seq_along(SVg)
  }
  
  if (length(nohit_idx) > 0) {
    intergenic_df <- tibble(
      hit_type         = "none",
      gene_accession   = NA_character_,
      gene             = NA_character_,
      product          = NA_character_,
      sv_accession     = as.character(seqnames(SVg))[nohit_idx],
      species          = sv_species[nohit_idx],
      variant          = sv_variant[nohit_idx],
      sv_start         = start(SVg)[nohit_idx],
      sv_end           = end(SVg)[nohit_idx],
      sv_width         = sv_w[nohit_idx],
      refmid           = sv_refmid[nohit_idx],
      overlap_bp       = 0L,
      pct_overlap_sv   = 0,
      pct_overlap_gene = 0,
      jaccard          = 0,
      fully_genic      = FALSE,
      hit_priority     = 0L,
      func_group       = "intergenic"
    )
    
    out <- bind_rows(both_df, intergenic_df)
  } else {
    out <- both_df
  }
  
  # --- 4) Keep best or all hits
  if (keep == "best" && nrow(out) > 0) {
    out <- out %>%
      arrange(desc(hit_priority), desc(overlap_bp)) %>%
      group_by(sv_accession, sv_start, sv_end) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      arrange(sv_accession, sv_start)
  } else if (nrow(out) > 0) {
    out <- arrange(out, sv_accession, desc(hit_priority), desc(overlap_bp))
  }
  
  if (verbose && nrow(out) > 0) {
    cat("\nAnnotation summary:\n")
    print(table(out$hit_type))
    cat("\nFunctional groups:\n")
    print(table(out$func_group))
  }
  
  return(out)
}






