# 04-alignment.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- read_delta  [from SVCM.R:695] ---
read_delta <- function(delta_path) {
  delta_lines <- (strsplit(readLines(delta_path), " "))[-c(1, 2)]
  if (length(delta_lines) == 0L) return(NULL)

  id_lines <- delta_lines[lengths(delta_lines) == 4][[1]]
  delta_alignments <- delta_lines[lengths(delta_lines) == 7] %>%
    unlist() %>%
    matrix(ncol = 7, byrow = TRUE)

  delta_alignments %>%
    apply(2, as.numeric) %>%
    matrix(ncol = 7) %>%
    `colnames<-`(c("rs", "re", "qs", "qe", "error", "e2", "zero")) %>%
    as_tibble() %>%
    dplyr::select(1:5) %>%
    mutate(
      strand = ifelse(qe - qs > 0, '+', '-'),
      rid  = sub(">", "", strsplit(id_lines, " ")[[1]]),
      qid  = strsplit(id_lines, " ")[[2]],
      rlen = as.numeric(strsplit(id_lines, " ")[[3]]),
      qlen = as.numeric(strsplit(id_lines, " ")[[4]]),
      rcov = abs(re - rs + 1),
      qcov = abs(qe - qs + 1),
      perc_error = round(100 * error / pmax(rcov, qcov), 2)
    ) %>%
    rowwise() %>%
    mutate(
      meanlen = ceiling(mean(c(rcov, qcov))),
      refmid  = (rs + re) / 2,
      qrymid  = (qs + qe) / 2,
      expected_offset_ref = (refmid * ((qlen / rlen) - 1)) / sqrt(2),
      expected_offset_qry = (qrymid * ((qlen / rlen) - 1)) / sqrt(2),
      mean_offset = mean(c(expected_offset_ref, expected_offset_qry)),
      fwd_dist = abs(abs(refmid - qrymid) / sqrt(2) - mean_offset),
      rev_dist = abs(((refmid + qrymid) - mean(c(rlen, qlen))) / sqrt(2)) - mean_offset,
      X_dist_raw = round(max(c(min(c(fwd_dist, rev_dist)), 0))),
      avg_genome_size = mean(c(rlen, qlen)),
      genome_size_diff = abs(rlen - qlen),
      X_dist = round(X_dist_raw / (1 + genome_size_diff / avg_genome_size))
    ) %>%
    ungroup() %>%
    arrange(rs) %>%
    dplyr::select(-c(refmid, qrymid, expected_offset_ref, expected_offset_qry,
                     mean_offset, X_dist_raw, avg_genome_size, genome_size_diff))
}

# --- filter_delta  [from SVCM.R:751] ---
filter_delta <- function(delta_table, maxgap = 1e4, minlen = 1e4, X_dist_diff = 5e4) {
  if (is.character(delta_table)) delta_table <- read_delta(delta_table)

  out_tibble <- delta_table %>%
    group_by(strand) %>%
    mutate(
      qry_gapsize = qs - lag(qe, default = qs[1]),
      ref_gapsize = rs - lag(re, default = rs[1]),
      X_diff      = X_dist - lag(X_dist, default = 0)
    ) %>%
    ungroup() %>%
    mutate(
      qry_gapsize  = ifelse(strand == "+", qry_gapsize, qry_gapsize * -1),
      qry_gaps_up  = cumsum(ifelse(qry_gapsize < maxgap, 0, 1)),
      qry_gaps_down = cumsum(ifelse(qry_gapsize > -maxgap, 0, 1)),
      ref_gaps     = cumsum(ifelse(abs(ref_gapsize) < maxgap, 0, 1)),
      XDD          = X_dist - lag(X_dist, default = 0),
      XDF          = cumsum(ifelse(abs(XDD) < X_dist_diff, 0, 1)),
      new_contigs  = data.table::rleid(strand, qry_gaps_up, qry_gaps_down, ref_gaps)
    ) %>%
    group_by(new_contigs, strand) %>%
    dplyr::summarise(
      X_dist  = mean(X_dist),
      meanlen = sum(meanlen),
      rs = min(rs), re = max(re),
      qs = unique(ifelse(strand == "+", min(qs), max(qs))),
      qe = unique(ifelse(strand == "+", max(qe), min(qe))),
      rid  = unique(rid), qid = unique(qid),
      slope = (qe - qs) / (re - rs),
      rlen = unique(rlen), qlen = unique(qlen),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    filter((meanlen > minlen) | (new_contigs == 1) | (new_contigs == max(new_contigs))) %>%
    mutate(
      qry_gapsize = qs - lag(qe, default = qs[1]),
      ref_gapsize = rs - lag(re, default = rs[1]),
      XDD         = X_dist - lag(X_dist, default = 0),
      XDF_diff2   = cumsum(abs(XDD) >= X_dist_diff),
      new_contigs = data.table::rleid(strand, XDF_diff2)
    ) %>%
    ungroup() %>%
    group_by(new_contigs, strand) %>%
    dplyr::summarise(
      rs = min(rs), re = max(re),
      qs = unique(ifelse(strand == "+", min(qs), max(qs))),
      qe = unique(ifelse(strand == "+", max(qe), min(qe))),
      rid = unique(rid), qid = unique(qid),
      slope  = (qe - qs) / (re - rs),
      rlen   = unique(rlen), qlen = unique(qlen),
      meanlen = mean(c(abs(qe - qs), abs(re - rs))),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(
      X_dist_raw = ifelse(
        strand == "+",
        mean(c(abs((rs - qs) / sqrt(2)), abs((re - qe) / sqrt(2)))),
        mean(c(abs((qs + rs - qlen) / sqrt(2)), abs((qe + re - qlen) / sqrt(2))))
      ),
      avg_genome_size = mean(c(rlen, qlen)),
      genome_size_diff = abs(rlen - qlen),
      X_dist = round(X_dist_raw / (1 + genome_size_diff / avg_genome_size))
    ) %>%
    ungroup() %>%
    dplyr::select(-X_dist_raw, -avg_genome_size, -genome_size_diff)

  out_tibble
}

# --- plot_delta  [from SVCM.R:823] ---
plot_delta <- function(delta_table, gtitle = "", xlb = NULL, ylb = NULL) {
  if (is.character(delta_table)) delta_table <- read_delta(delta_table)

  myColors <- brewer.pal(5, "Set1")[c(1, 2)]
  names(myColors) <- c("-", "+")
  colScale <- scale_colour_manual(name = "grp", values = myColors)

  p <- ggplot(
    delta_table %>% mutate(rid = sub("NZ_", "", rid), qid = sub("NZ_", "", qid)),
    aes(x = rs, xend = re, y = qs, yend = qe,
        colour = factor(strand, levels = c("-", "+")))
  ) +
    geom_segment(alpha = 1, linewidth = 1.5) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.35),
      axis.title = element_text(size = 12)
    ) +
    geom_abline(intercept = mean(c(delta_table$rlen[1], delta_table$qlen[1])),
                slope = -1, color = "darkgrey", linetype = "dashed", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1,
                color = "darkgrey", linetype = "dashed", linewidth = 0.5) +
    colScale +
    scale_x_continuous(limits = c(1, delta_table$rlen[1]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(1, delta_table$qlen[1]), expand = c(0, 0))

  if (is.null(xlb) & is.null(ylb)) {
    p <- p +
      xlab(sub("NZ_", "", unique(pull(delta_table, rid)))) +
      ylab(sub("NZ_", "", unique(pull(delta_table, qid))))
  } else {
    p <- p + xlab(xlb) + ylab(ylb)
  }
  p
}

# --- benchmark_alignments  [from SVCM.R:2734] ---
benchmark_alignments <- function(genome_matrix, delta_dir, sync_dir,
                                 strategy_name, mc.cores = 6L) {
  message("Starting benchmarking: ", strategy_name)
  cmds <- generate_nucmer_commands(genome_matrix, delta_dir, sync_dir,
                                   nuc_remove = FALSE, check_nuc = FALSE)
  results <- pbmcapply::pbmclapply(seq_len(nrow(genome_matrix)), function(i) {
    ref <- genome_matrix[i, 1]; qry <- genome_matrix[i, 2]
    pair_id <- paste(ref, qry, sep = "_vs_")
    start_time <- Sys.time()
    res <- tryCatch({
      exit_code <- system(cmds[i])
      end_time <- Sys.time()
      tibble::tibble(strategy = strategy_name, pair = pair_id, ref = ref, qry = qry,
                     runtime_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
                     exit_code = exit_code)
    }, error = function(e) {
      tibble::tibble(strategy = strategy_name, pair = pair_id, ref = ref, qry = qry,
                     runtime_sec = NA_real_, exit_code = NA_integer_)
    })
    res
  }, mc.cores = mc.cores)

  results <- Filter(Negate(is.null), results)
  dplyr::bind_rows(results)
}

# --- filter_combined_delta  [from SVCM.R:2678] ---
filter_combined_delta <- function(delta_dir) {
  filt_files <- list.files(delta_dir, pattern = "_filtered\\.delta$", full.names = TRUE)
  if (!length(filt_files)) return(tibble::tibble())

  res <- lapply(filt_files, function(f) {
    dt <- tryCatch(read_delta(f), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    fd <- tryCatch(filter_delta(dt), error = function(e) NULL)
    if (is.null(fd) || !nrow(fd)) return(NULL)
    fd
  })

  dplyr::bind_rows(Filter(Negate(is.null), res))
}

# --- svcm_delta_paths  [from SVCM.R:2912] ---
svcm_delta_paths <- function(alignment_dir, unfiltered_dir, pair_id) {
  list(
    filt = file.path(alignment_dir, sprintf("%s_filtered.delta", pair_id)),
    unf  = file.path(unfiltered_dir, sprintf("%s.delta", pair_id))
  )
}

# --- resolve_delta_paths  [from SVCM6.R:292] ---
resolve_delta_paths <- function(alignment_dir, unfiltered_dir, pair_id) {
  # primary expected contract
  filt1 <- file.path(alignment_dir, sprintf("%s_filtered.delta", pair_id))
  unf1  <- file.path(unfiltered_dir, sprintf("%s.delta", pair_id))
  
  # if missing, attempt find
  if (!file.exists(filt1)) {
    cand <- list.files(alignment_dir, pattern = paste0("^", gsub("\\|", "\\\\|", pair_id), ".*\\.delta$"),
                       full.names = TRUE)
    cand <- cand[!grepl("unfiltered", cand)]
    filt1 <- if (length(cand) > 0) cand[1] else filt1
  }
  if (!file.exists(unf1)) {
    cand <- list.files(unfiltered_dir, pattern = paste0("^", gsub("\\|", "\\\\|", pair_id), ".*\\.delta$"),
                       full.names = TRUE)
    unf1 <- if (length(cand) > 0) cand[1] else unf1
  }
  
  list(delta_filt = filt1, delta_unfilt = unf1)
}

# --- svcm_step2_combine_per_qid  [from SVCM6.R:835] ---
svcm_step2_combine_per_qid <- function(
    r_vars_dir,
    sv_by_qid_dir,
    drop_non_ref_rid = TRUE,
    quarantine_dirname = "_quarantine_bad_cache",
    rid_singleton = TRUE
) {
  expected_rid <- NA_character_
  ref_rid_path <- file.path(r_vars_dir, "ref_rid.rds")
  if (file.exists(ref_rid_path)) expected_rid <- readRDS(ref_rid_path)
  expected_rid <- if (is.na(expected_rid)) NA_character_ else sub(" .*", "", expected_rid)
  
  files <- list.files(sv_by_qid_dir, pattern="\\.rds$", full.names=TRUE, recursive=FALSE)
  if (!length(files)) stop("No per-qid SV files found in: ", sv_by_qid_dir)
  
  qdir <- file.path(sv_by_qid_dir, quarantine_dirname)
  dir.create(qdir, recursive=TRUE, showWarnings=FALSE)
  
  required_cols <- c("rid","qid","start","end","width","variant")
  
  read_one <- function(fp) {
    stem <- sub("\\.rds$", "", basename(fp))
    
    obj <- tryCatch(readRDS(fp), error=function(e) e)
    if (inherits(obj, "error") || !is.data.frame(obj)) {
      return(list(ok=FALSE, why="read error or not data.frame", fp=fp, stem=stem))
    }
    
    miss <- setdiff(required_cols, names(obj))
    if (length(miss)) {
      return(list(ok=FALSE, why=paste0("missing cols: ", paste(miss, collapse=", ")), fp=fp, stem=stem))
    }
    
    # enforce pipeline invariants
    obj$qid <- stem
    obj$rid <- sub(" .*", "", as.character(obj$rid))
    if (!"variant_specific" %in% names(obj)) obj$variant_specific <- NA_character_
    
    # optional: drop non-chromosome rid rows
    if (isTRUE(drop_non_ref_rid) && !is.na(expected_rid)) {
      obj <- obj[obj$rid == expected_rid, , drop=FALSE]
    }
    
    list(ok=TRUE, sv=tibble::as_tibble(obj), fp=fp, stem=stem)
  }
  
  parsed <- lapply(files, read_one)
  bad <- vapply(parsed, function(x) !isTRUE(x$ok), logical(1))
  
  if (any(bad)) {
    for (x in parsed[bad]) {
      to <- file.path(qdir, basename(x$fp))
      if (!file.exists(to)) file.rename(x$fp, to)
    }
    warning("Quarantined ", sum(bad), " bad cache files to: ", qdir)
  }
  
  svs <- dplyr::bind_rows(lapply(parsed[!bad], `[[`, "sv"))
  
  if (nrow(svs)) {
    svs <- dplyr::distinct(svs, qid, rid, variant, variant_specific, start, end, .keep_all=TRUE)
  }
  
  if (isTRUE(rid_singleton) && nrow(svs)) {
    u <- unique(svs$rid)
    if (!is.na(expected_rid)) {
      stopifnot(length(u) == 1L, u == expected_rid)
    } else {
      stopifnot(length(u) <= 1L)
    }
  }
  
  out_path <- file.path(r_vars_dir, "all_svs_raw.rds")
  saveRDS(svs, out_path)
  svs
}

# --- step2_align_call_resumable  [from SVCM6.R:640] ---
step2_align_call_resumable <- function(sync_dir,
                                       alignment_dir,
                                       unfiltered_dir,
                                       r_vars_dir,
                                       sv_by_qid_dir,
                                       failures_dir,
                                       species_name,
                                       ref_genome,
                                       ref_path,
                                       n_cores = 8L,
                                       keep_delta_bp = 50000L,
                                       svcm_version = "1.0.0",
                                       params_fields = list(nucmer_minmatch = 32L,
                                                            delta_filter_1to1 = TRUE),
                                       force_recompute = FALSE) {
  # prerequisites
  if (!exists("simple_nucmer", mode = "function")) stop("simple_nucmer() not found. Did you source('SVCM.R')?")
  if (!exists("delta_all_SVs", mode = "function")) stop("delta_all_SVs() not found. Did you source('SVCM.R')?")
  rid_expected <- sub(" .*", "", names(Biostrings::readDNAStringSet(ref_path))[1])
  
  dir.create(sv_by_qid_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(failures_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(unfiltered_dir, recursive = TRUE, showWarnings = FALSE)
  
  # staleness key
  params_path <- file.path(r_vars_dir, "step2_params.rds")
  step2_params <- c(list(svcm_version = svcm_version, ref_genome = ref_genome), params_fields)
  cached_params <- load_cache(params_path)
  params_stale <- force_recompute
  
  if (!is.null(cached_params) && !force_recompute) {
    # compare only keys that affect SV output
    for (nm in names(step2_params)) {
      if (!identical(cached_params[[nm]], step2_params[[nm]])) {
        message(sprintf("  WARNING: step2_params$%s changed (%s -> %s) => cache stale",
                        nm, as.character(cached_params[[nm]]), as.character(step2_params[[nm]])))
        params_stale <- TRUE
      }
    }
  }
  save_cache(step2_params, params_path)
  
  inv <- inventory_synced(sync_dir)
  qids <- setdiff(inv$sync_accs, ref_genome)
  gz_files <- file.path(sync_dir, sprintf("%s.txt.gz", qids))
  n_queries <- length(qids)
  message(sprintf("Step 2: %d queries vs ref %s", n_queries, ref_genome))
  
  # run summary (latest-per-qid)
  run_summary_path <- file.path(r_vars_dir, "run_summary.rds")
  run_summary <- load_cache(run_summary_path)
  if (is.null(run_summary)) run_summary <- tibble::tibble()
  
  pb <- utils::txtProgressBar(min = 1, max = n_queries, style = 3)
  
  rid_expected_ref <- sub(" .*", "", names(Biostrings::readDNAStringSet(ref_path))[1])
  saveRDS(rid_expected_ref, file.path(r_vars_dir, "ref_rid.rds"))
  
  for (i in seq_len(n_queries)) {
    qid    <- qids[i]
    qid_gz <- gz_files[i]
    sv_out <- file.path(sv_by_qid_dir, sprintf("%s.rds", qid))
    
    # Tier A: cached SVs exist + params not stale
    if (!params_stale && file.exists(sv_out) && file.size(sv_out) > 0) {
      utils::setTxtProgressBar(pb, i)
      next
    }
    
    t_start <- Sys.time()
    status <- "success"
    tier   <- "C"
    n_svs  <- 0L
    has_large <- FALSE
    err_msg <- NA_character_
    tmp_plain <- NULL
    
    tryCatch({
      if (!file.exists(qid_gz)) stop("Missing synced .txt.gz for qid: ", qid, " at ", qid_gz)
      
      pair_id <- sprintf("%s_v_%s", ref_genome, qid)
      dp <- resolve_delta_paths(alignment_dir, unfiltered_dir, pair_id)
      delta_filt   <- dp$delta_filt
      delta_unfilt <- dp$delta_unfilt
      
      # Tier B: reuse deltas only if filtered delta exists and non-empty AND (preferably) unfiltered exists
      if (file.exists(delta_filt) && file.size(delta_filt) > 0 &&
          file.exists(delta_unfilt) && file.size(delta_unfilt) > 0) {
        tier <- "B"
      } else {
        tier <- "C"
        
        # materialize qid plain fasta (temporary)
        tmp_plain <- file.path(sync_dir, sprintf("%s.txt", qid))
        if (!file.exists(tmp_plain) || file.size(tmp_plain) == 0) {
          suppressWarnings(system2("gunzip", args = c("-k", "-f", shQuote(qid_gz))))
        }
        if (!file.exists(tmp_plain) || file.size(tmp_plain) == 0) stop("Failed to materialize query .txt for: ", qid)
        
        # run nucmer + delta-filter (via your simple_nucmer)
        nuc_cmd <- simple_nucmer(ref = ref_path,
                                 qry = tmp_plain,
                                 output = pair_id,
                                 alignment_dir = alignment_dir,
                                 check_nuc = FALSE)
        system(nuc_cmd)
        
        # re-resolve in case naming differs
        dp <- resolve_delta_paths(alignment_dir, unfiltered_dir, pair_id)
        delta_filt   <- dp$delta_filt
        delta_unfilt <- dp$delta_unfilt
        
        if (!file.exists(delta_filt) || file.size(delta_filt) == 0) {
          stop("nucmer produced no filtered delta for qid=", qid, " expected: ", delta_filt)
        }
        if (!file.exists(delta_unfilt) || file.size(delta_unfilt) == 0) {
          # not fatal for all variants, but delta_duplications usually needs it
          message("  NOTE: unfiltered delta missing/empty for qid=", qid, "; duplication calls may be incomplete.")
        }
      }
      
      svs <- delta_all_SVs(
        delta_file            = delta_filt,
        unfiltered_delta_file = if (file.exists(delta_unfilt)) delta_unfilt else NA_character_,
        species_name          = species_name
      )
      
      if (is.data.frame(svs) && nrow(svs) > 0) {
        svs$rid <- sub(" .*", "", as.character(svs$rid))
        svs <- svs[svs$rid == rid_expected_ref, , drop = FALSE]  # FILTER
      }
      
      svs <- sv_standardize_schema(svs, qid = qid, rid_expected = NULL)  # qid overwritten inside function
      
      
      n_svs <- if (is.null(svs)) 0L else nrow(svs)
      has_large <- if (n_svs == 0) FALSE else any(svs$width >= keep_delta_bp, na.rm = TRUE)
      
      saveRDS(svs, sv_out)
      
      # cleanup deltas unless large SVs found
      if (!has_large) {
        if (file.exists(delta_filt))   file.remove(delta_filt)
        if (file.exists(delta_unfilt)) file.remove(delta_unfilt)
      }
      
      # cleanup temp query .txt
      if (!is.null(tmp_plain) && file.exists(tmp_plain)) file.remove(tmp_plain)
      
    }, error = function(e) {
      status <<- "FAILED"
      err_msg <<- conditionMessage(e)
      # write failure artifact
      err_path <- file.path(failures_dir, sprintf("%s.rds", qid))
      saveRDS(list(qid = qid, error = e, time = Sys.time()), err_path)
      if (!is.null(tmp_plain) && file.exists(tmp_plain)) file.remove(tmp_plain)
    })
    
    t_end <- Sys.time()
    row <- tibble::tibble(
      qid = qid,
      status = status,
      tier = tier,
      n_svs = n_svs,
      has_large = has_large,
      runtime_sec = as.numeric(difftime(t_end, t_start, units = "secs")),
      sv_out_path = sv_out,
      error_message = err_msg
    )
    
    # keep latest-per-qid
    run_summary <- run_summary |>
      dplyr::filter(.data$qid != !!qid) |>
      dplyr::bind_rows(row)
    
    # periodic save
    if (i %% 50 == 0 || i == n_queries) saveRDS(run_summary, run_summary_path)
    
    utils::setTxtProgressBar(pb, i)
  }
  
  close(pb)
  saveRDS(run_summary, run_summary_path)
  
  # report
  n_ok   <- sum(run_summary$status == "success", na.rm = TRUE)
  n_fail <- sum(run_summary$status == "FAILED", na.rm = TRUE)
  n_B    <- sum(run_summary$tier == "B" & run_summary$status == "success", na.rm = TRUE)
  n_C    <- sum(run_summary$tier == "C" & run_summary$status == "success", na.rm = TRUE)
  message(sprintf("Step 2 DONE: %d succeeded, %d failed | Tier B: %d | Tier C: %d", n_ok, n_fail, n_B, n_C))
  
  invisible(run_summary)
}
