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
        if (nrow(dnaA_vect)>1){
          stop(">1 dnaA")
        }
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
          "no dnaA annotation",
        )
    }
    
    return(dnaA_vect)
  }


generate_nucmer_commands <- 
  function(genome_matrix, delta_dir, sync_dir, nuc_remove = F){
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
    #   genome_seq[
    #     !(grepl("plasmid", names(genome_seq),  
    #             ignore.case = T))
    #   ]
    # genome_seq <- 
    #   genome_seq[
    #     !(grepl("phage", names(genome_seq),  
    #             ignore.case = T)) |
    #       grepl(
    #         "phage resistant", 
    #         names(genome_seq),  
    #         ignore.case = T
    #       )
    #   ]
    # 
    # ## check chromosomes
    # genome_seqs <- 
    #   genome_seq[order(width(genome_seq), decreasing = T),]
    
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
    sync_dir, identity_dir, species_table,
    species_name
  ){
    # run mash if it hasn't been run
    if (!
        file.exists(
          sprintf(
            "%s/%s_mash.txt", 
            identity_dir, species_name
          )
        )
    ){
      mash_bash <- 
        mash_commands(species_name, identity_dir, sync_dir)
      system(mash_bash[1], ignore.stderr = T)
      system(mash_bash[2])
    }
    
    
    id_table <-
      fread(
        sprintf(
          "%s/%s_%s.txt", 
          identity_dir, species_name, filetype
        ),# stringsAsFactors = F, 
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
simple_nucmer <- 
  function(
    ref, qry, alignment_dir, output, nuc_remove = F
  ){
    nuc_params <-
      "--mum --maxgap=500 --mincluster=100"
    nucmer <-
      sprintf(
        "nucmer %s --prefix=%s/%s",
        nuc_params,
        alignment_dir,
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
        "%s/%s_filtered.delta",
        alignment_dir,
        output
      )

    nucmer_filter <-
      sprintf(
        "delta-filter -r -q %s/%s.delta > %s",
        alignment_dir,
        output,
        filter_out
      )
    ## to save space, remove the original delta file
    if (nuc_remove){
      nuc_remove <- 
        sprintf(
          "rm %s/%s.delta",
          alignment_dir,
          output
        )
    } else {
      unfiltered_dir <- 
        sprintf(
          "%s/unfiltered",
          alignment_dir
        )
      if (!dir.exists(unfiltered_dir)){
        dir.create(unfiltered_dir)
      }
      
      nuc_remove <- 
        sprintf(
          "mv %s/%s.delta %s/%s.delta",
          alignment_dir,
          output,
          unfiltered_dir,
          output
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

filter_delta <- function(delta_table,
                         maxgap      = 1e4,
                         minlen      = 1e4,
                         X_dist_diff = 5e4) {
  # you still need data.table for rleid()
  require(dplyr)
  require(data.table)
  
  delta_table %>%
    # 1) order & group per ref–query strand
    arrange(rid, qid, strand, rs) %>%
    group_by(rid, qid, strand) %>%
    
    # 2) compute the three gap‐metrics and flag real breaks
    mutate(
      # for "+" strands, gap = qs - previous qe;
      # for "-" we mirror it so it still measures “distance along the query”
      qry_gap = if_else(strand == "+",
                        qs - lag(qe,  default = qs[1]),
                        lag(qe,  default = qe[1]) - qs),
      
      ref_gap = rs - lag(re, default = rs[1]),
      X_diff  = abs(X_dist - lag(X_dist, default = X_dist[1])),
      
      # a break if ANY threshold is exceeded
      gap_flag = (qry_gap  >  maxgap) |
        (ref_gap  >  maxgap) |
        (X_diff   >  X_dist_diff),
      
      # start contig 1, bump by 1 whenever gap_flag==TRUE
      contig_id = cumsum(gap_flag) + 1
    ) %>%
    ungroup() %>%
    
    # 3) collapse each contig into one “merged alignment”
    group_by(rid, qid, strand, contig_id) %>%
    summarise(
      rs      = min(rs),
      re      = max(re),
      qs      = if (first(strand) == "+") min(qs) else max(qs),
      qe      = if (first(strand) == "+") max(qe) else min(qe),
      meanlen = sum(meanlen),
      X_dist  = mean(X_dist),
      slope   = if_else(re == rs, NA_real_, (qe - qs) / (re - rs)),
      "rlen" = unique(rlen),
      "qlen" = unique(qlen),
      .groups = "keep"
    ) %>%
    
    # 4) drop anything too small
    filter(meanlen > minlen)
    # 4) drop any tiny contigs
    
}
filter_delta(delta_table) %>% plot_delta

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
        XDF_diff2   = cumsum(abs(XDD) >= X_dist_diff),
        new_contigs = rleid(strand, XDF_diff2)
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
  function(delta_table) {
  # --- your existing filtering & rowwise min/max --------------------------
    fd <- delta_table %>%
      filter_delta() %>%             # your pipeline up through `fd`
      rowwise() %>%
      mutate(
        qs2 = min(qs, qe),
        qe2 = max(qs, qe)
      ) %>%
      ungroup()
    
    # --- annotate a simple row index so we can refer back by position ----------
    fd2 <- fd %>%
      arrange(rs) %>%
      mutate(idx = row_number())
    
    nested_idx <- fd2 %>%
      filter(strand == "+" &
               lag(strand, default = NA) == "-" &
               lead(strand, default = NA) == "-") %>%
      pull(idx)
    # should be “4” in your example
    
    # 2) bracketers: the minus‐strand just before/after that nested piece
    outer_above <- max(fd2$idx[fd2$strand == "-" & fd2$idx < nested_idx], na.rm = TRUE)
    outer_below <- min(fd2$idx[fd2$strand == "-" & fd2$idx > nested_idx], na.rm = TRUE)
    outer_idx   <- c(outer_above, outer_below)  # should be c(2,5)
    
    # 3) any other minus strand is a localized inversion
    localized_idx <- fd2$idx[fd2$strand == "-" & !(fd2$idx %in% outer_idx)]
    
    # 4) build one‐row summary for each event
    events <- 
      bind_rows(
      # localized
      fd2 %>% 
        filter(idx %in% localized_idx) %>%
        summarise(
          variant_specific      = "localized",
          rs         = min(rs),
          re         = max(re),
          qs         = min(qs2),
          qe         = max(qe2),
          contigs    = paste(new_contigs, collapse = ","),
          strand= unique(strand),
          X_dist = mean(X_dist),
          width = (re-rs+1),
          rid = unique(rid),
          qid = unique(qid)
        ),
      # outer
      fd2 %>% filter(idx %in% outer_idx) %>%
        summarise(
          variant_specific      = "outer",
          rs         = min(rs),
          re         = max(re),
          qs         = min(qs2),
          qe         = max(qe2),
          contigs    = paste(new_contigs, collapse = ","),
          strand= unique(strand),
          X_dist = mean(X_dist),
          width = (re-rs+1),
          rid = unique(rid),
          qid = unique(qid)
        ),
      # nested
      fd2 %>% filter(idx == nested_idx) %>%
        summarise(
          variant_specific      = "nested",
          rs         = rs,
          re         = re,
          qs         = qs2,
          qe         = qe2,
          contigs    = paste(new_contigs, collapse = ","),
          strand=unique(strand),
          X_dist = mean(X_dist),
          width = (re-rs+1),
          rid = unique(rid),
          qid = unique(qid)
        )
    )
    
    return(events)
}

events$contigs





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







