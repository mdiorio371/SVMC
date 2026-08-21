# 12-orchestration.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- svcm_step2_post_summary  [from SVCM.R:3363] ---
svcm_step2_post_summary <- function(
    sync_dir,
    r_vars_dir,
    sv_by_qid_dir,
    species_name,
    ref_path,
    include = c("", "1v_the_rest"),
    exclude = c("mash_id", "ss20"),
    length_bins = c(0, 1e3, 1e4, 5e4, 1e5, 5e5, 1e6, Inf),
    recurrence_min_bp = 10e3,
    ro_min = 0.95,
    out_subdir = "step2_summary"
) {
  # Requires: svcm_sync_inventory(), svcm_token_from_fasta_header(), svcm_load_cache()/saveRDS
  if (!exists("svcm_sync_inventory", mode = "function")) stop("svcm_sync_inventory() not found.")
  if (!exists("svcm_token_from_fasta_header", mode = "function")) stop("svcm_token_from_fasta_header() not found.")
  if (!exists("sv_rec2", mode = "function")) stop("sv_rec2() not found (for SR recurrence).")
  if (!exists("group_svs", mode = "function")) stop("group_svs() not found (for non-SR recurrence).")
  
  sync_dir     <- normalizePath(sync_dir, winslash="/", mustWork=TRUE)
  r_vars_dir   <- normalizePath(r_vars_dir, winslash="/", mustWork=TRUE)
  sv_by_qid_dir<- normalizePath(sv_by_qid_dir, winslash="/", mustWork=TRUE)
  ref_path     <- normalizePath(ref_path, winslash="/", mustWork=TRUE)
  
  out_dir <- file.path(r_vars_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # inventory = "started with"
  inv <- svcm_sync_inventory(sync_dir, include = include, exclude = exclude, recursive = FALSE)
  if (nrow(inv) == 0) stop("No synced genomes found under: ", sync_dir)
  inv <- inv %>% dplyr::mutate(path = normalizePath(.data$path, winslash="/", mustWork=TRUE))
  
  # drop non-FASTA artifacts pre-emptively
  looks_fa <- vapply(inv$path, svcm_file_looks_fasta, logical(1))
  inv <- inv[looks_fa, , drop = FALSE]
  
  # run_summary = what "have"
  run_summary_path <- file.path(r_vars_dir, "run_summary.rds")
  if (!file.exists(run_summary_path)) stop("Missing run_summary.rds: ", run_summary_path)
  rs <- readRDS(run_summary_path)
  if (!("qid" %in% names(rs))) stop("run_summary.rds missing qid column")
  
  qid_started <- sort(unique(inv$qid))
  qid_success <- rs %>%
    dplyr::filter(.data$status == "success") %>%
    dplyr::pull(.data$qid) %>%
    unique() %>%
    intersect(qid_started) %>%
    sort()
  
  # cohort sv files (non-recursive)
  sv_files <- list.files(sv_by_qid_dir, pattern = "\\.rds$", full.names = TRUE, recursive = FALSE)
  sv_files <- sv_files[basename(sv_files) %in% paste0(qid_success, ".rds")]
  
  # read + combine (only cohort)
  sv_list <- lapply(sv_files, function(fp) {
    x <- readRDS(fp)
    if (is.null(x) || !is.data.frame(x)) return(NULL)
    x
  })
  sv_all <- dplyr::bind_rows(Filter(Negate(is.null), sv_list))
  
  # save combined
  all_svs_fp <- file.path(r_vars_dir, "all_svs_raw.rds")
  saveRDS(sv_all, all_svs_fp)
  
  # ref genome length + rid token for recurrence
  ref_seq <- Biostrings::readDNAStringSet(ref_path)
  if (length(ref_seq) != 1L) stop("ref_path must be single-record FASTA for summary.")
  L <- as.integer(Biostrings::width(ref_seq))
  ref_rid_token <- svcm_token_from_fasta_header(names(ref_seq)[1])
  
  # ---------- SV length distribution ----------
  len_tbl <- sv_all %>%
    dplyr::mutate(
      width = as.numeric(.data$width),
      len_bin = cut(.data$width, breaks = length_bins, right = FALSE, include.lowest = TRUE)
    ) %>%
    dplyr::count(.data$variant, .data$len_bin, name = "n_svs") %>%
    dplyr::group_by(.data$variant) %>%
    dplyr::mutate(prop = .data$n_svs / sum(.data$n_svs)) %>%
    dplyr::ungroup()
  
  # ---------- Recurrence (hybrid, to keep it fast) ----------
  # Filter to a size floor for recurrence reporting (you can set recurrence_min_bp = 0 if you want everything)
  sv_rec_in <- sv_all %>% dplyr::filter(as.numeric(.data$width) >= recurrence_min_bp)
  
  sr <- sv_rec_in %>% dplyr::filter(.data$variant == "Structural_rearrangement")
  non_sr <- sv_rec_in %>% dplyr::filter(.data$variant != "Structural_rearrangement")
  
  # SR recurrence: canonicalize + sv_rec2
  sr_rec <- tibble::tibble()
  if (nrow(sr) > 0) {
    sr2 <- canon_sr_short_arc(sr, genome_len = L) %>% dplyr::mutate(refmid = .data$refmid_circ)
    sr_rec <- sv_rec2(
      sr2,
      ref_acc = ref_rid_token,
      by_rid = TRUE,
      size_min_bp = 0L,
      match_variant_specific = TRUE,
      ro_thresh = ro_min,
      genome_len = L
    ) %>%
      dplyr::mutate(motif_uid = paste0("MOTIF|SR|", .data$recurrence_group))
  }
  
  # non-SR recurrence: group_svs per variant
  non_sr_rec <- tibble::tibble()
  if (nrow(non_sr) > 0) {
    non_sr_rec <- non_sr %>%
      dplyr::group_by(.data$variant) %>%
      dplyr::group_modify(~{
        group_svs(.x, min_overlap_bp = 1L, ro_min = ro_min, genome_len = L)
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        motif_uid = paste0("MOTIF|", .data$variant, "|", .data$sv_group),
        group_size = as.integer(.data$group_size)
      )
  }
  
  motifs <- dplyr::bind_rows(sr_rec, non_sr_rec) %>%
    dplyr::mutate(width = as.numeric(.data$width))
  
  motif_tbl <- motifs %>%
    dplyr::group_by(.data$motif_uid, .data$variant) %>%
    dplyr::summarise(
      n_present = dplyr::n_distinct(.data$qid),
      width_median = stats::median(.data$width, na.rm = TRUE),
      refmid_median = stats::median(as.numeric(.data$refmid), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_present), dplyr::desc(.data$width_median))
  
  recurrence_by_len <- motif_tbl %>%
    dplyr::mutate(
      len_bin = cut(.data$width_median, breaks = length_bins, right = FALSE, include.lowest = TRUE)
    ) %>%
    dplyr::count(.data$variant, .data$len_bin, name = "n_motifs") %>%
    dplyr::arrange(.data$variant, .data$len_bin)
  
  # ---------- cohort + pipeline summary ----------
  summary_tbl <- tibble::tibble(
    species = species_name,
    n_started = length(qid_started),
    n_success = length(qid_success),
    n_failed  = length(setdiff(qid_started, qid_success)),
    success_rate = if (length(qid_started) > 0) 100 * length(qid_success) / length(qid_started) else NA_real_,
    n_svs_total = nrow(sv_all),
    recurrence_min_bp = recurrence_min_bp,
    ro_min = ro_min,
    genome_len = L,
    ref_rid = ref_rid_token
  )
  
  cohort_fp <- file.path(out_dir, "qid_success.txt")
  writeLines(qid_success, cohort_fp)
  
  saveRDS(list(
    summary = summary_tbl,
    qid_started = qid_started,
    qid_success = qid_success,
    all_svs_path = all_svs_fp,
    length_distribution = len_tbl,
    motif_table = motif_tbl,
    recurrence_by_len = recurrence_by_len
  ), file.path(out_dir, "step2_summary_bundle.rds"))
  
  utils::write.csv(summary_tbl, file.path(out_dir, "step2_summary_overview.csv"), row.names = FALSE)
  utils::write.csv(len_tbl,     file.path(out_dir, "step2_sv_length_distribution.csv"), row.names = FALSE)
  utils::write.csv(motif_tbl,   file.path(out_dir, "step2_motif_recurrence_table.csv"), row.names = FALSE)
  utils::write.csv(recurrence_by_len, file.path(out_dir, "step2_recurrence_by_length.csv"), row.names = FALSE)
  
  message("Saved Step2 summary to: ", out_dir)
  invisible(list(out_dir = out_dir, summary = summary_tbl, qid_success = qid_success))
}

# --- svcm_step3_split  [from SVCM6.R:1004] ---
svcm_step3_split <- function(all_svs_raw,
                             kde_result,
                             min_large_bp = 50000L) {
  svcm_stop_if_missing(c("dplyr", "tibble"))
  
  svs <- tibble::as_tibble(all_svs_raw)
  if (!"width" %in% names(svs)) stop("all_svs_raw missing 'width' column", call. = FALSE)
  
  gene_pile <- svs %>%
    dplyr::filter(.data$width >= kde_result$lower_bp[1] & .data$width <= kde_result$upper_bp[1])
  
  large_svs <- svs %>%
    dplyr::filter(.data$width >= as.integer(min_large_bp))
  
  list(gene_pile = gene_pile, large_svs = large_svs)
}

# --- svcm_step3_kde_split_annotate_cds  [from SVCM6.R:1167] ---
svcm_step3_kde_split_annotate_cds <- function(r_vars_dir,
                                              prokka_dir,
                                              ref_genome,
                                              min_large_bp = 50000,
                                              search_bp = c(500, 5000),
                                              cap_bounds_bp = c(300, 8000),
                                              gene_thresh = 80,
                                              bp_pad = 0L) {
  svcm_stop_if_missing(c("dplyr", "tibble"))
  
  all_svs_path <- file.path(r_vars_dir, "all_svs_raw.rds")
  if (!file.exists(all_svs_path)) stop("Missing all_svs_raw.rds at: ", all_svs_path, call. = FALSE)
  all_svs_raw <- readRDS(all_svs_path) %>% tibble::as_tibble()
  
  # 3.1 KDE (cache)
  kde_path <- file.path(r_vars_dir, "kde_result.rds")
  kde <- svcm_load_cache(kde_path)
  if (is.null(kde)) {
    kde <- detect_gene_length_peak_kde(
      widths_bp = all_svs_raw$width,
      search_bp = search_bp,
      cap_bounds_bp = cap_bounds_bp
    )
    svcm_save_cache(kde, kde_path)
  }
  
  # 3.2 split (cache)
  gene_path  <- file.path(r_vars_dir, "gene_pile.rds")
  large_path <- file.path(r_vars_dir, "large_svs.rds")
  
  gene_pile <- svcm_load_cache(gene_path)
  large_svs <- svcm_load_cache(large_path)
  if (is.null(gene_pile) || is.null(large_svs)) {
    sp <- svcm_step3_split(all_svs_raw, kde_result = kde, min_large_bp = min_large_bp)
    gene_pile <- sp$gene_pile
    large_svs <- sp$large_svs
    svcm_save_cache(gene_pile, gene_path)
    svcm_save_cache(large_svs, large_path)
  }
  
  # 3.3 Prokka CDS GRanges (cache)
  cds_path <- file.path(r_vars_dir, "cds_gr.rds")
  cds_gr <- svcm_load_cache(cds_path)
  if (is.null(cds_gr)) {
    gff_path <- file.path(prokka_dir, paste0(ref_genome, ".gff"))
    if (!file.exists(gff_path)) {
      # Try to give a helpful pointer if the ref FASTA path was saved elsewhere
      ref_plain_guess <- file.path(dirname(prokka_dir), "..", "sync", "UNKNOWN", "UNKNOWN", "_ref_plain", paste0(ref_genome, ".txt"))
      msg <- paste0(
        "Missing Prokka GFF: ", gff_path, "\n\n",
        "You need to run Prokka on the reference FASTA and put the .gff here.\n",
        "Expected filename: ", paste0(ref_genome, ".gff"), "\n\n",
        "Example command (edit paths):\n  ",
        svcm_prokka_cmd(ref_fasta = "<PATH_TO_REF_FASTA.txt>", outdir = prokka_dir, prefix = ref_genome, cpus = 8L),
        "\n"
      )
      stop(msg, call. = FALSE)
    }
    cds_gr <- svcm_read_prokka_gff_cds(gff_path)
    svcm_save_cache(cds_gr, cds_path)
  }
  
  # 3.3 annotate gene pile (cache)
  ann_path <- file.path(r_vars_dir, "gene_annotations.rds")
  gene_annotations <- svcm_load_cache(ann_path)
  if (is.null(gene_annotations)) {
    gene_annotations <- svcm_annotate_gene_pile_cds(gene_pile, cds_gr, gene_thresh = gene_thresh)
    svcm_save_cache(gene_annotations, ann_path)
  }
  
  list(
    kde_result = kde,
    gene_pile = gene_pile,
    large_svs = large_svs,
    cds_gr = cds_gr,
    gene_annotations = gene_annotations
  )
}

# --- step1_setup_and_check2  [from SVCM.R:4861] ---
step1_setup_and_check2 <- function(ncbi_table, species_name, n_orient = 200L, n_cores = 4L) {
  cat("\nSTEP 1:", species_name, "\n", sep = " ")
  cat(strrep("=", 80), "\n")
  
  sa_entries <- ncbi_table %>%
    dplyr::filter(.data$species == sanitize_species(species_name))
  
  if (nrow(sa_entries) == 0) stop("No entries in ncbi_table for: ", sanitize_species(species_name))
  
  # DO NOT use slice(); it may be masked. Use base indexing.
  phylum_val <- first_scalar_chr(sa_entries$phylum_clade[1])
  if (is.na(phylum_val) || phylum_val == "") stop("phylum_clade missing/empty for: ", species_name)
  clade_name <- gsub(" ", "_", phylum_val)
  
  dirs <- list(
    sync      = file.path("processing/sync",      clade_name, species_name),
    alignment = file.path("processing/alignment", clade_name, species_name),
    identity  = file.path("processing/identity",  clade_name, species_name),
    r_vars    = file.path("processing/r_vars",    clade_name, species_name),
    results   = file.path("processing/results",   clade_name, species_name),
    temp      = file.path("processing/temp",      clade_name, species_name)
  )
  lapply(dirs, ensure_dir)
  
  step1_ck <- file.path(dirs$temp, "step1_checkpoint.rds")
  # if (rds_exists(step1_ck)) {
  #   cat("  - Step1 checkpoint exists; loading:", step1_ck, "\n")
  #   x <- readRDS(step1_ck)
  #   x$species_name <- species_name
  #   return(x)
  # }
  
  n_orient <- min(as.integer(n_orient), nrow(sa_entries))
  cat("  - Entries:", nrow(sa_entries), " | Orienting first:", n_orient, "\n")
  cat("  - Clade:", clade_name, "\n")
  cat("  - results:", dirs$results, "\n")
  
  sync_dir <- dirs$sync
  id_dir   <- dirs$identity
  
  assembly_table <- pbmcapply::pbmclapply(
    X = seq_len(n_orient),
    FUN = function(i) {
      assembly_path <- as.character(sa_entries$ftp_path[i])
      
      fasta_path <- sprintf("%s/%s_genomic.fna.gz", assembly_path, basename(assembly_path))
      gff_path   <- sprintf("%s/%s_genomic.gff.gz", assembly_path, basename(assembly_path))
      
      fna_exists <- file_exists_at_url(fasta_path)
      gff_exists <- file_exists_at_url(gff_path)
      
      if (!isTRUE(fna_exists) || !isTRUE(gff_exists)) {
        return(tibble(asm_name = basename(assembly_path), note = "missing fasta/gff"))
      }
      
      dnaA_table <- dnaA_from_gff(gff_path)
      pos_tib    <- dnaA_sync(fasta_path, dnaA_table, sync_dir)
      pos_tib %>% dplyr::mutate(asm_name = basename(assembly_path), note = "ok")
    },
    mc.cores = n_cores
  ) %>% dplyr::bind_rows()
  
  cat("  - Computing pairwise identities for reference selection...\n")
  species_identities <- 
    get_pairwise_identities(sync_dir = sync_dir, id_dir = id_dir, species_name = species_name)
  
  most_related <- species_identities$id_dists %>%
    dplyr::arrange(.data$Sequence_mean_distance) %>%
    dplyr::slice(1) %>%               # explicitly dplyr
    dplyr::pull(.data$Sequence)
  
  done_asm <- assembly_table$asm_name
  the_rest <- sa_entries %>% dplyr::filter(!(basename(.data$ftp_path) %in% done_asm))
  
  summary_table <- tibble(
    species = species_name,
    clade   = clade_name,
    n_total = nrow(sa_entries),
    n_oriented_seed = nrow(assembly_table),
    n_remaining = nrow(the_rest),
    reference = most_related
  )
  
  x <- list(
    species_name = species_name,
    dirs = dirs,
    assembly_table = assembly_table,
    species_table = sa_entries,
    the_rest = the_rest,
    most_related = most_related,
    ids = species_identities$mash_tib, 
    summary = summary_table
  )
  
  saveRDS(x, step1_ck)
  cat("  ✓ Saved:", step1_ck, "\n")
  print(summary_table)
  x
}

# --- step3_full_processing2  [from SVCM.R:4961] ---
step3_full_processing2 <- function(step1_data, batch_size = 50L, n_cores = 4L) {
  cat("\nSTEP 3:", step1_data$species_name, "\n", sep = " ")
  cat(strrep("=", 80), "\n")
  
  dirs <- step1_data$dirs
  out_all_svs <- file.path(dirs$results, "all_svs_raw.rds")
  out_ori     <- file.path(dirs$results, "ori_info_all.rds")
  out_batch   <- file.path(dirs$results, "batch_processing_summary.rds")
  out_stats   <- file.path(dirs$results, "overall_stats.rds")
  
  if (rds_exists(out_all_svs) && rds_exists(out_ori) && rds_exists(out_batch) && rds_exists(out_stats)) {
    cat("  - Step3 outputs exist; skipping.\n")
    return(invisible(TRUE))
  }
  
  alignment_dir <- file.path(dirs$alignment, "1v_the_rest")
  sync_dir      <- file.path(dirs$sync,      "1v_the_rest")
  ensure_dir(alignment_dir); ensure_dir(sync_dir); ensure_dir(dirs$temp); ensure_dir(dirs$results)
  
  ref_path <- file.path(dirs$sync, paste0(step1_data$most_related, ".txt"))
  if (!file.exists(ref_path)) stop("Reference sync file missing: ", ref_path)
  
  n_genomes <- nrow(step1_data$the_rest)
  n_batches <- ceiling(n_genomes / batch_size)
  
  cat("  - Genomes:", n_genomes, " | batch_size:", batch_size, " | batches:", n_batches, " | cores:", n_cores, "\n")
  
  all_batch_svs <- vector("list", n_batches)
  all_batch_ori <- vector("list", n_batches)
  batch_summary <- vector("list", n_batches)
  
  for (b in seq_len(n_batches)) {
    i1 <- (b - 1L) * batch_size + 1L
    i2 <- min(b * batch_size, n_genomes)
    batch_genomes <- step1_data$the_rest[i1:i2, , drop = FALSE]
    
    cat("  - Batch", b, "/", n_batches, " (", i1, "-", i2, ")\n", sep = "")
    
    t0 <- Sys.time()
    batch_res <- pbmcapply::pbmclapply(
      X = seq_len(nrow(batch_genomes)),
      FUN = function(i) {
        tryCatch({
          arrange_align(
            assembly_path          = batch_genomes$ftp_path[i],
            ref_path               = ref_path,
            alignment_dir_the_rest = alignment_dir,
            sync_dir_the_rest      = sync_dir,
            call_SVs               = TRUE,
            species_name           = step1_data$species_name
          )
        }, error = function(e) NULL)
      },
      mc.cores = n_cores
    )
    
    have <- !vapply(batch_res, is.null, logical(1))
    all_batch_svs[[b]] <- dplyr::bind_rows(lapply(which(have), \(ii) batch_res[[ii]][[2]]))
    all_batch_ori[[b]] <- dplyr::bind_rows(lapply(which(have), \(ii) batch_res[[ii]][[1]]))
    
    dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    batch_summary[[b]] <- tibble(
      batch = b,
      n_genomes = nrow(batch_genomes),
      n_successful = sum(have),
      n_failed = sum(!have),
      n_svs = nrow(all_batch_svs[[b]]),
      n_ori_rows = nrow(all_batch_ori[[b]]),
      time_mins = dt
    )
  }
  
  all_svs <- dplyr::bind_rows(all_batch_svs)
  ori_info_all <- dplyr::bind_rows(all_batch_ori)
  batch_summary_df <- dplyr::bind_rows(batch_summary)
  
  overall_stats <- batch_summary_df %>%
    dplyr::summarise(
      total_genomes = sum(.data$n_genomes),
      total_successful = sum(.data$n_successful),
      success_rate = 100 * total_successful / total_genomes,
      total_svs = sum(.data$n_svs),
      total_time_mins = sum(.data$time_mins),
      .groups = "drop"
    )
  
  saveRDS(all_svs, out_all_svs)
  saveRDS(ori_info_all, out_ori)
  saveRDS(batch_summary_df, out_batch)
  saveRDS(overall_stats, out_stats)
  
  cat("  ✓ Saved Step3 outputs to:", dirs$results, "\n")
  print(overall_stats)
  invisible(TRUE)
}

# --- run_result1_salmonella_specific  [from SVCM.R:4409] ---
run_result1_salmonella_specific <- function(
    stepC_fp = file.path(
      "processing", "results", "Gamma", "Salmonella_enterica",
      "sv_structure_stepBC", "stepC_sv_pcoa_bundle.rds"
    ),
    identity_dir = file.path(
      "processing", "identity", "Gamma", "Salmonella_enterica"
    ),
    out_dir = file.path(
      "processing", "results", "Gamma", "Salmonella_enterica",
      "sv_structure_stepBC", "subtype_discovery"
    ),
    mash_bundle_fp = NULL,
    prefer_full_matrix = TRUE,
    module_jaccard_thresh = 0.90,
    min_module_carriers = 10L,
    min_major_bp = 2e5,
    min_major_carriers = 20L,
    target_sizes_bp = c(0.86e6, 0.55e6, 1.53e6),
    max_target_rel_diff = 0.35,
    jaccard_exclusive_max = 0.25,
    dbscan_minPts = 5L,
    dbscan_max_clusters = 9L,
    nperm = 199L
) {
  species <- "Salmonella_enterica"
  clade <- "Gamma"
  species_lab <- gsub("_", " ", species, fixed = TRUE)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(mash_bundle_fp)) {
    mash_bundle_fp <- find_mash_bundle_fp(identity_dir)
  }
  
  b <- readRDS(stepC_fp)
  mb <- readRDS(mash_bundle_fp)
  
  M0 <- choose_M_matrix(b, prefer_full_matrix = prefer_full_matrix)
  motif_meta0 <- make_motif_meta(b, M0)
  D_mash <- extract_mash_dist(mb)
  
  h <- harmonize_M_and_D(M0, D_mash)
  M_int <- h$M
  Dmat_int <- h$Dmat
  D_int <- h$D
  
  cat("\nShared genomes:", nrow(M_int), "\n")
  print(flatness_report(D_int))
  
  ord <- make_ord_df(D_int)
  ord_df <- ord$df
  ord_df$qid <- rownames(Dmat_int)
  
  db_out <- summarize_dbscan_grid(ord_df[, c("MDS1", "MDS2")], minPts = dbscan_minPts)
  db_pick <- choose_dbscan_eps_sane(db_out, max_clusters = dbscan_max_clusters)
  eps_use <- db_pick$eps[1]
  
  fit_db <- dbscan::dbscan(ord_df[, c("MDS1", "MDS2")], eps = eps_use, minPts = dbscan_minPts)
  ord_df$cluster_raw <- as.character(fit_db$cluster)
  ord_df$cluster_plot <- collapse_dbscan_clusters(ord_df$cluster_raw, min_cluster_n = 5L)
  
  cluster_levels <- c(
    sort(unique(ord_df$cluster_plot[!ord_df$cluster_plot %in% c("Noise", "Other")])),
    intersect(c("Noise", "Other"), unique(ord_df$cluster_plot))
  )
  ord_df$cluster_plot <- factor(ord_df$cluster_plot, levels = cluster_levels)
  
  motif_meta_int <- motif_meta0 |>
    dplyr::filter(motif_uid %in% colnames(M_int)) |>
    dplyr::mutate(
      carriers_overlap = as.integer(colSums(M_int != 0)[match(motif_uid, colnames(M_int))])
    ) |>
    dplyr::filter(variant_collapsed == "Inversion") |>
    dplyr::arrange(dplyr::desc(carriers_overlap), motif_uid)
  
  mods <- build_modules_from_jaccard(
    M_mat = M_int[, motif_meta_int$motif_uid, drop = FALSE],
    motif_tbl = motif_meta_int,
    jaccard_thresh = module_jaccard_thresh
  )
  
  M_mod <- mods$M_module
  membership_tbl <- mods$membership_tbl
  feature_tbl <- mods$module_tbl |>
    dplyr::rename(feature_id = module_id)
  
  scan_tbl <- lapply(colnames(M_mod), function(fid) {
    scan_one_feature(
      feature_id = fid,
      M_bin = M_mod,
      Dmat = Dmat_int,
      cluster_raw = ord_df$cluster_raw,
      nperm = nperm
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::left_join(feature_tbl, by = "feature_id") |>
    dplyr::arrange(dplyr::desc(delta_z), dplyr::desc(carriers_overlap), feature_id)
  
  sel <- select_salmonella_major_modules(
    scan_tbl = scan_tbl,
    M_bin = M_mod,
    target_sizes_bp = target_sizes_bp,
    min_major_bp = min_major_bp,
    min_major_carriers = min_major_carriers,
    max_target_rel_diff = max_target_rel_diff,
    jaccard_exclusive_max = jaccard_exclusive_max
  )
  
  major_ids <- sel$selected_ids
  major_tbl <- sel$selected_tbl |>
    dplyr::mutate(
      label = target_label
    )
  
  label_map <- stats::setNames(major_tbl$label, major_tbl$feature_id)
  
  state_df <- assign_states(M_mod, major_ids, label_map)
  plot_df_all <- ord_df |>
    dplyr::left_join(state_df, by = "qid")
  
  ## add module presence columns for the selected states
  for (fid in major_ids) {
    plot_df_all[[fid]] <- as.integer(M_mod[plot_df_all$qid, fid] != 0)
  }
  
  lead_id <- major_ids[1]
  
  p_cluster <- ggplot2::ggplot(plot_df_all, ggplot2::aes(MDS1, MDS2)) +
    ggplot2::geom_point(color = "grey80", alpha = 0.35, size = 1.5) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df_all, !is.na(cluster_plot)),
      ggplot2::aes(fill = cluster_plot),
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.90,
      size = 2.4
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = species_lab,
      x = paste0("MDS1 (", sprintf("%.1f", ord$ord_pct["MDS1"]), "%)"),
      y = paste0("MDS2 (", sprintf("%.1f", ord$ord_pct["MDS2"]), "%)"),
      fill = "DBSCAN"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "italic", hjust = 0.5)
    )
  
  p_lead <- ggplot2::ggplot(plot_df_all, ggplot2::aes(MDS1, MDS2)) +
    ggplot2::geom_point(color = "grey80", alpha = 0.35, size = 1.5) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df_all, .data[[lead_id]] == 1L),
      shape = 21,
      fill = "black",
      color = "black",
      stroke = 0.25,
      alpha = 0.95,
      size = 2.4
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = species_lab,
      x = paste0("MDS1 (", sprintf("%.1f", ord$ord_pct["MDS1"]), "%)"),
      y = paste0("MDS2 (", sprintf("%.1f", ord$ord_pct["MDS2"]), "%)")
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "italic", hjust = 0.5)
    )
  
  p_state <- ggplot2::ggplot(plot_df_all, ggplot2::aes(MDS1, MDS2)) +
    ggplot2::geom_point(color = "grey80", alpha = 0.35, size = 1.5) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df_all, state != "None"),
      ggplot2::aes(fill = state),
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.95,
      size = 2.5
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = species_lab,
      x = paste0("MDS1 (", sprintf("%.1f", ord$ord_pct["MDS1"]), "%)"),
      y = paste0("MDS2 (", sprintf("%.1f", ord$ord_pct["MDS2"]), "%)"),
      fill = "State"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "italic", hjust = 0.5)
    )
  
  combo_plot <- p_cluster / p_lead / p_state
  
  plot_obj <- list(
    species = species,
    clade = clade,
    baseline = basename(mash_bundle_fp),
    ord_pct = ord$ord_pct,
    flatness = flatness_report(D_int),
    dbscan_summary = db_out,
    dbscan_pick = db_pick,
    dbscan_eps = eps_use,
    plot_df_all = plot_df_all,
    major_tbl = major_tbl,
    scan_tbl = scan_tbl,
    feature_tbl = feature_tbl,
    membership_tbl = membership_tbl,
    params = list(
      module_jaccard_thresh = module_jaccard_thresh,
      min_module_carriers = min_module_carriers,
      min_major_bp = min_major_bp,
      min_major_carriers = min_major_carriers,
      target_sizes_bp = target_sizes_bp,
      max_target_rel_diff = max_target_rel_diff,
      jaccard_exclusive_max = jaccard_exclusive_max,
      dbscan_minPts = dbscan_minPts,
      dbscan_max_clusters = dbscan_max_clusters,
      nperm = nperm
    )
  )
  
  saveRDS(plot_obj, file.path(out_dir, "plot_data_major_states.rds"))
  
  ggplot2::ggsave(
    filename = file.path(out_dir, "result1_salmonella_specific.png"),
    plot = combo_plot,
    width = 8,
    height = 13,
    dpi = 300
  )
  
  cat("\nChosen DBSCAN row:\n")
  print(db_pick, row.names = FALSE)
  
  cat("\nSelected Salmonella major states:\n")
  print(
    major_tbl |>
      dplyr::select(
        feature_id,
        label,
        target_size_bp,
        carriers_overlap,
        representative_motif,
        representative_width_median,
        delta_z,
        pct_dom_nonnoise,
        outside_dom_carriers_nonnoise
      )
  )
  
  invisible(
    list(
      plot = combo_plot,
      plot_obj = plot_obj
    )
  )
}

# --- select_salmonella_major_modules  [from SVCM.R:4323] ---
select_salmonella_major_modules <- function(
    scan_tbl,
    M_bin,
    target_sizes_bp = c(0.86e6, 0.55e6, 1.53e6),
    min_major_bp = 2e5,
    min_major_carriers = 20L,
    max_target_rel_diff = 0.35,
    jaccard_exclusive_max = 0.25
) {
  stopifnot(all(c("feature_id", "carriers_overlap", "delta_z", "representative_width_median") %in% names(scan_tbl)))
  
  feature_sets <- lapply(colnames(M_bin), function(fid) {
    rownames(M_bin)[as.vector(M_bin[, fid] != 0)]
  })
  names(feature_sets) <- colnames(M_bin)
  
  cand_tbl <- scan_tbl |>
    dplyr::filter(
      is.finite(delta_z),
      delta_z > 0,
      carriers_overlap >= min_major_carriers,
      is.finite(representative_width_median),
      representative_width_median >= min_major_bp
    ) |>
    dplyr::mutate(
      width_mbp = representative_width_median / 1e6
    )
  
  if (nrow(cand_tbl) == 0) {
    stop("No large inversion modules passed the Salmonella-specific candidate filter.")
  }
  
  selected_ids <- character(0)
  selected_targets <- numeric(0)
  
  for (target_bp in target_sizes_bp) {
    target_tbl <- cand_tbl |>
      dplyr::mutate(
        rel_diff_target = abs(representative_width_median - target_bp) / target_bp
      ) |>
      dplyr::filter(rel_diff_target <= max_target_rel_diff)
    
    if (length(selected_ids) > 0) {
      keep_idx <- vapply(target_tbl$feature_id, function(fid) {
        all(vapply(selected_ids, function(prev) {
          jaccard_sets(feature_sets[[fid]], feature_sets[[prev]]) <= jaccard_exclusive_max
        }, logical(1)))
      }, logical(1))
      target_tbl <- target_tbl[keep_idx, , drop = FALSE]
    }
    
    if (nrow(target_tbl) == 0) {
      next
    }
    
    pick <- target_tbl |>
      dplyr::arrange(
        rel_diff_target,
        dplyr::desc(carriers_overlap),
        dplyr::desc(delta_z),
        feature_id
      ) |>
      dplyr::slice(1)
    
    selected_ids <- c(selected_ids, pick$feature_id)
    selected_targets <- c(selected_targets, target_bp)
  }
  
  if (length(selected_ids) == 0) {
    stop("No Salmonella major modules matched the target inversion sizes.")
  }
  
  selected_tbl <- scan_tbl |>
    dplyr::filter(feature_id %in% selected_ids) |>
    dplyr::mutate(
      target_size_bp = selected_targets[match(feature_id, selected_ids)],
      target_label = sprintf("%.2f Mbp Inversion", target_size_bp / 1e6)
    ) |>
    dplyr::arrange(match(feature_id, selected_ids))
  
  list(
    selected_ids = selected_ids,
    selected_tbl = selected_tbl
  )
}

# --- svcm_self_check  [from SVCM6.R:486] ---
svcm_self_check <- function(sync_dir,
                            r_vars_dir,
                            id_dir,
                            sv_by_qid_dir,
                            alignment_dir,
                            unfiltered_dir,
                            expected_ref = NULL) {
  res <- list()
  
  inv <- inventory_synced(sync_dir)
  res$sync_has_gz <- length(inv$sync_gz_files) > 0
  
  res$has_position_table <- file.exists(file.path(r_vars_dir, "position_table.rds"))
  res$has_ref_rds        <- file.exists(file.path(r_vars_dir, "ref_genome.rds"))
  
  if (!is.null(expected_ref)) {
    res$ref_in_sync <- expected_ref %in% inv$sync_accs
    res$ref_plain_exists <- file.exists(file.path(sync_dir, sprintf("%s.txt", expected_ref))) ||
      file.exists(file.path(sync_dir, sprintf("%s.txt.gz", expected_ref)))
  } else {
    res$ref_in_sync <- NA
    res$ref_plain_exists <- NA
  }
  
  # seed-medoid artifacts
  res$has_seed_medoid <- dir.exists(file.path(id_dir, "seed_medoid")) ||
    file.exists(file.path(r_vars_dir, "mash_result_seed.rds")) ||
    file.exists(file.path(r_vars_dir, "mash_seed_result.rds"))
  
  # step2 cache
  res$has_sv_by_qid_dir <- dir.exists(sv_by_qid_dir)
  res$n_sv_by_qid <- if (dir.exists(sv_by_qid_dir)) length(list.files(sv_by_qid_dir, pattern="\\.rds$")) else 0L
  
  # delta dirs exist
  res$alignment_dir_exists <- dir.exists(alignment_dir)
  res$unfiltered_dir_exists <- dir.exists(unfiltered_dir)
  
  # report
  ok <- function(x) if (isTRUE(x)) "PASS" else if (identical(x, NA)) "NA" else "FAIL"
  message("SVCM self-check:")
  message("  Sync inventory exists:        ", ok(res$sync_has_gz))
  message("  position_table.rds present:   ", ok(res$has_position_table))
  message("  ref_genome.rds present:       ", ok(res$has_ref_rds))
  message("  seed-medoid artifacts exist:  ", ok(res$has_seed_medoid))
  message("  sv_by_qid dir exists:         ", ok(res$has_sv_by_qid_dir),
          " (n=", res$n_sv_by_qid, ")")
  message("  alignment_dir exists:         ", ok(res$alignment_dir_exists))
  message("  unfiltered_dir exists:        ", ok(res$unfiltered_dir_exists))
  if (!is.null(expected_ref)) {
    message("  ref in synced filenames:      ", ok(res$ref_in_sync))
    message("  ref materializable:           ", ok(res$ref_plain_exists))
  }
  invisible(res)
}
