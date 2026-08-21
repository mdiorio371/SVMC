# ==========================================================================
# SVMC: whole-genome alignment (MUMmer/nucmer)
# nucmer command generation, delta parsing/filtering/plotting, benchmarking.
# ==========================================================================

simple_nucmer <- function(
    ref,
    qry,
    alignment_dir,
    output,
    minmatch   = 32,
    nuc_remove = FALSE,
    check_nuc  = TRUE
){
  unfiltered_dir <- file.path(alignment_dir, "unfiltered")
  dir.create(unfiltered_dir, showWarnings = FALSE, recursive = TRUE)
  
  prefix <- file.path(alignment_dir, output)
  delta  <- sprintf("%s.delta", prefix)
  filt   <- file.path(alignment_dir, sprintf("%s_filtered.delta", output))
  unf    <- file.path(unfiltered_dir, sprintf("%s.delta", output))
  
  # back to --mumreference (no --maxmatch)
  nuc_params <- sprintf("--maxgap=500 --mincluster=100 --minmatch=%d", minmatch)
  
  mkdir_cmd  <- sprintf("mkdir -p %s %s", alignment_dir, unfiltered_dir)
  nucmer_cmd <- sprintf("nucmer %s --prefix=%s %s %s",
                        nuc_params, prefix, ref, qry)
  filter_cmd <- sprintf("delta-filter -r -q %s > %s", delta, filt)
  post_cmd   <- if (isTRUE(nuc_remove)) {
    sprintf("rm -f %s", delta)
  } else {
    sprintf("mv %s %s", delta, unf)
  }
  
  core <- paste(mkdir_cmd, nucmer_cmd, filter_cmd, post_cmd, sep = " && ")
  
  if (isTRUE(check_nuc)) {
    sprintf("if [ ! -e %s ]; then %s ; fi", filt, core)
  } else {
    core
  }
}

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
          # Calculate raw X_dist
          X_dist_raw = 
            round(max(c(min(c(fwd_dist, rev_dist)), 0))),
          # Adjust for genome size difference
          avg_genome_size = mean(c(rlen, qlen)),
          genome_size_diff = abs(rlen - qlen),
          # Normalize: add expected offset from size difference to threshold tolerance
          X_dist = 
            round(X_dist_raw / (1 + genome_size_diff / avg_genome_size))
        ) %>%
        ungroup %>%
        arrange(rs) %>%
        dplyr::select(-c(
          refmid, qrymid, expected_offset_ref,
          expected_offset_qry, mean_offset,
          X_dist_raw, avg_genome_size, genome_size_diff
        ))
      
      return(delta_alignments)
    }
  }

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
        # Calculate raw X_dist
        X_dist_raw = 
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
          ),
        # Adjust for genome size difference
        avg_genome_size = mean(c(rlen, qlen)),
        genome_size_diff = abs(rlen - qlen),
        X_dist = round(X_dist_raw / (1 + genome_size_diff / avg_genome_size))
      ) %>% 
      ungroup() %>%
      dplyr::select(-X_dist_raw, -avg_genome_size, -genome_size_diff)
    
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
