# ==========================================================================
# SVMC: OriC localization & genome orientation
# dnaA detection, dnaA-box scoring, GC-skew disparity, OriC localization.
# ==========================================================================

dnaA_from_gff <- 
  function(gff_path, genome_acc){
    
    # Open connection and ensure it's closed even if read.gff fails
    gff_gzcon <- 
      gzcon(url(gff_path))
    on.exit(try(close(gff_gzcon), silent = TRUE), add = TRUE)
    
    # Defensive read — return empty data frame on network/parse failure
    gff_table <- 
      tryCatch(
        read.gff(gff_gzcon),
        error = function(e) {
          warning(
            "dnaA_from_gff: failed to read ", gff_path,
            " — ", conditionMessage(e),
            call. = FALSE
          )
          NULL
        }
      )
    
    if (is.null(gff_table)) {
      return(data.frame())
    }
    
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
    
    if (nrow(dnaA_vect) > 1) {
      
      if (
        any(
          grepl(
            "gene=DnaA",
            dnaA_vect$attributes,
            ignore.case = T
          )
        )
      ) {
        dnaA_vect <-
          dnaA_vect[
            grepl(
              "gene=DnaA",
              dnaA_vect$attributes,
              ignore.case = T
            ),
          ]
      }
      
      if (all(pull(dnaA_vect, strand) == "+")) {
        dnaA_vect <-
          dnaA_vect %>%
          filter(start == min(start)) %>%
          dplyr::slice(1)
      } else if (all(pull(dnaA_vect, strand) == "-")) {
        dnaA_vect <-
          dnaA_vect %>%
          filter(start == max(start)) %>% 
          dplyr::slice(1)
      } else {
        dnaA_vect <-
          dnaA_vect %>% 
          dplyr::slice(1)
      }
    }
    
    # Type-consistent return: emit a warning but keep as 0-row data frame
    # (do NOT replace with a character string — breaks downstream callers)
    if (nrow(dnaA_vect) == 0) {
      warning(
        "dnaA_from_gff: no dnaA annotation found in ",
        gff_path,
        call. = FALSE
      )
    }
    
    return(dnaA_vect)
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


id_dnaA_boxes <- 
  function(
    genome_seq,
    dnaA_box = 
      DNAString("TTATCCACA")
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
 
    
    return(out_tib)
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
    
    #dnaA_box_tib$b %>% boxplot
    
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
    
    
    
    ### agreement metric
    # ori_score_tib <- 
    #     OriC_score(
    #         GC_disp, dnaA_box_tib, dnaA_gene_pos= fasta_len, 
    #         genome_len = fasta_len
    #         )
    
    # FIXME(SVMC): oriC_scores_v2() is called below but is NOT defined anywhere in
    # the original functions.R. Recreate it or remove this call before the
    # package can pass R CMD check. See handoff 'Known gaps' #1.
    ori_scores <- 
      oriC_scores_v2(
        dbox = dbox,
        GC_tib = GC_tib,
        dgene = round(fasta_len / 2),
        genome_len = fasta_len
      ) 
    #
    # ori_score_plot <- 
    #     ori_scores %>%
    #     select(mid, box_z, gc_z, gene_z, ori_score) %>%
    #     pivot_longer(cols = -mid, names_to = "component", values_to = "value") %>%
    #     ggplot(aes(x = mid / 1000, y = value, color = component)) +
    #     geom_line(size = 1) +
    #     labs(x = "Genome position (Kb)", y = "Normalized score",
    #          title = "Component contributions to OriC agreement score") +
    #     theme_classic()
    
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
    top_regions
    plot_OriC_markers
    
    
    
    dnaA_gene_box_GC_disp <- 
      plot_OriC_markers(
        GC_disp,
        dnaA_box_tib,
        fasta_len,
        top_regions
      )
    
    oriC_plot_dir <- 
      sprintf(
        "%s/OriC",
        plot_dir
      )
    if (!dir.exists(oriC_plot_dir)){dir.create(oriC_plot_dir)}
    ggsave(
      sprintf(
        "%s/%s.png",
        oriC_plot_dir,
        st2$accession[j]
      ),
      plot = dnaA_gene_box_GC_disp, 
      width = 8, height = 6, dpi = 300
    )
  }


compute_disparities <- function(fasta_file, output_file, sample_n = 10000) {
  #' Compute cumulative G-C and A-T disparities along a genome sequence.
  #' Drop-in replacement for nucleotide_disparities.py — runs ~10-20x faster.
  #'
  #' @param fasta_file Path to FASTA file (.txt or .fna or .fna.gz)
  #' @param output_file Path to write CSV with columns G-C and A-T
  #' @param sample_n Sample every N bases (default 10000)
  #' @return Tibble with G-C and A-T columns (also written to output_file)
  
  
  seq <- readDNAStringSet(fasta_file)[[1]]
  seq_len <- nchar(seq)
  n_windows <- floor(seq_len / sample_n)
  
  if (n_windows == 0) {
    af <- alphabetFrequency(seq)  # named vector for DNAString
    out <- data.frame(`G-C` = af["G"] - af["C"],
                      `A-T` = af["A"] - af["T"],
                      check.names = FALSE)
    write.csv(out, output_file, row.names = FALSE)
    return(out)
  }
  
  # Vectorized: split into windows and count all at once
  starts <- seq(1, n_windows * sample_n, by = sample_n)
  ends   <- starts + sample_n - 1
  
  # Use Views for zero-copy subsetting
  v <- Views(seq, start = starts, end = pmin(ends, seq_len))
  af <- alphabetFrequency(v)
  
  # Cumulative sums (matching the Python behavior)
  gc <- cumsum(af[, "G"] - af[, "C"])
  at <- cumsum(af[, "A"] - af[, "T"])
  
  out <- data.frame(`G-C` = as.integer(gc),
                    `A-T` = as.integer(at),
                    check.names = FALSE)
  
  write.csv(out, output_file, row.names = FALSE)
  return(out)
}


GC_check_R <- function(fasta_file, disparity_file) {
  #' Compute GC disparity and locate the GC-skew minimum.
  #' Drop-in replacement for GC_check() — no Python dependency.
  #'
  #' @param fasta_file Path to FASTA file
  #' @param disparity_file Path to write/read disparity CSV
  #' @return Tibble with accession, GC_ori_loc, and optionally ori_dist
  
  
  # Compute disparities (replaces python3 nucleotide_disparities.py)
  compute_disparities(fasta_file, disparity_file)
  
  skew_tib <- fread(disparity_file)
  
  acc_i <- basename(sub(".txt$", "", fasta_file))
  disparities_dir <- dirname(disparity_file)
  
  gc_values <- skew_tib$`G-C`
  min_idx <- which.min(gc_values)
  
  # Check if minimum is at position 1 (already oriented correctly)
  if (min_idx == 1) {
    return(
      tibble(
        accession  = acc_i,
        GC_ori_loc = 1
      )
    )
  }
  
  # Zoom in for precise localization
  minima_coords <- c(
    1e4 * (min_idx - 1),
    1e4 * (min_idx + 1)
  )
  
  original_ff <- readDNAStringSet(fasta_file)
  
  if (minima_coords[2] > width(original_ff)) {
    minima_coords[2] <- width(original_ff)
  }
  if (minima_coords[1] < 1) {
    minima_coords[1] <- 1
  }
  
  disp_regions_dir <- sprintf("%s/regions", disparities_dir)
  if (!dir.exists(disp_regions_dir)) dir.create(disp_regions_dir, recursive = TRUE)
  
  ori_region_file <- sprintf("%s/%s_ori.txt", disp_regions_dir, acc_i)
  
  ori_region <- subseq(original_ff, minima_coords[1], minima_coords[2])
  names(ori_region) <- paste0(acc_i, "_ori_region")
  writeXStringSet(ori_region, ori_region_file)
  
  # Fine-grained disparity at single-base resolution
  ori_disparity_file <- sprintf("%s/%s_ori_region.csv", disp_regions_dir, acc_i)
  compute_disparities(ori_region_file, ori_disparity_file, sample_n = 1)
  
  ori_disp <- fread(ori_disparity_file)
  GC_min_region_idx <- which.min(ori_disp$`G-C`)[1]
  ori_location <- minima_coords[1] + GC_min_region_idx
  
  tibble(
    accession     = acc_i,
    GC_ori_loc    = ori_location,
    ori_dist      = min(abs(1 - ori_location), abs(ori_location - width(original_ff))),
    ori_dist_frac = round(100 * ori_dist / width(original_ff), 4),
    genome_len    = as.numeric(width(original_ff))
  )
}


# --------------------------------------------------------------------------
# Internal placeholder for the multi-marker OriC scorer.
#
# locate_OriC() calls oriC_scores_v2() to combine dnaA-box density and GC-skew
# into a single OriC score, but the original functions.R never defined it
# (see handoff 'Known gaps' #1). This stub keeps the package self-consistent
# (no undefined globals) and fails loudly rather than returning wrong results.
# Recreate the real implementation, or supply OriC coordinates directly,
# before relying on locate_OriC().
# Not exported.
# --------------------------------------------------------------------------
oriC_scores_v2 <- function(dbox, GC_tib, dgene, genome_len) {
  stop(
    "oriC_scores_v2() is not yet implemented. It is required by locate_OriC() ",
    "to merge dnaA-box density and GC-skew into multi-marker OriC scores. ",
    "Recreate it (or provide OriC coordinates directly) before using ",
    "locate_OriC(). See the package issue tracker.",
    call. = FALSE
  )
}
