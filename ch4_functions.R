##Ch4 4 functions


load_arrange <- 
  function(
    assembly_path,
    sync_dir_the_rest
  ){
    
    fasta_path <- 
      sprintf(
        "%s/%s_genomic.fna.gz",
        assembly_path, 
        gsub(".*/", "", assembly_path)
      )
    
    gff_path <- 
      sprintf(
        "%s/%s_genomic.gff.gz",
        assembly_path,
        gsub(".*/", "", assembly_path)
      )
    
    fna_exists <- 
      file_exists_at_url(fasta_path)
    gff_exists <- 
      file_exists_at_url(gff_path)
    
    dnaA_table <- 
      dnaA_from_gff(
        gff_path
      )
    if (!fna_exists | !gff_exists){
      position_tib <- 
        tibble( 
          asm_name = 
            basename(assembly_path),
          note = "missing fasta/gff"
        )
    } else {
      position_tib <-   
        dnaA_sync(fasta_path, dnaA_table, sync_dir_the_rest)
    }
    return(position_tib)
  }

sync_dir_the_rest %>% list.files
load_arrange()

arrange_align <- 
  function(
    assembly_path,
    ref_path,
    alignment_dir_the_rest,
    sync_dir_the_rest
  ){
    
    fasta_path <- 
      sprintf(
        "%s/%s_genomic.fna.gz",
        assembly_path, 
        gsub(".*/", "", assembly_path)
      )
    
    gff_path <- 
      sprintf(
        "%s/%s_genomic.gff.gz",
        assembly_path,
        gsub(".*/", "", assembly_path)
      )
    
    fna_exists <- 
      file_exists_at_url(fasta_path)
    gff_exists <- 
      file_exists_at_url(gff_path)
    
    dnaA_table <- 
      dnaA_from_gff(
        gff_path
      )
    if (!fna_exists | !gff_exists){
      position_tib <- 
        tibble( 
          asm_name = 
            basename(assembly_path),
          note = "missing fasta/gff"
        )
    } else {
      position_tib <-   
        dnaA_sync(fasta_path, dnaA_table, sync_dir_the_rest)
    }
    
    out_str <- 
      sprintf(
        "%s_v_%s",
        sub(".txt", "", basename(ref_path)), 
        sub(".txt", "", basename(position_tib$out_file))
      )
    ## align to the ref
    simple_nucmer(
      ref = ref_path, 
      qry = position_tib$out_file,
      output = 
        out_str,
      alignment_dir = alignment_dir_the_rest
    ) %>% system()
    
    delta_filt_file <- 
      sprintf(
        "%s/%s_filtered.delta",
        alignment_dir_the_rest,
        out_str
        )
        
    all_delta_unfilt_files <- 
      sprintf(
        "%s/unfiltered/%s.delta",
        alignment_dir_the_rest,
        out_str
      ) 
    
    delta_SVs <- 
      delta_all_SVs(
        delta_file = filtered_delta, 
        unfiltered_delta_file = unfiltered_delta,
        species_name = species_name
      )
    SRs <- 
      delta_SVs %>% 
      filter(variant=="Structural rearrangement")
    
    if (delta_SVs$variant_specific[1]=="none"){
      file.remove(position_tib$out_file)
    }
    
    
    
    out_tib <- 
      position_tib 
    
    out_list <- 
      list(
        position_tib,
        delta_SVs
      )
    
    return(position_tib)
  }






filter_combined_delta <- 
  function(
    dd, 
    minimum_alignment_percentage, minimum_inversion_length
  ){
    
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



## in align coverage is not enforced yet
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
