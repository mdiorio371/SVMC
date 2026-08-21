# R/98-from-notebooks.R
# Functions recovered from Ch2.Rmd / SVCM.Rmd analysis notebooks.

# --- assign_gmm_submodes  [from Ch2.Rmd] ---
assign_gmm_submodes <- function(df) {
  n <- nrow(df)
  if(n < 10) return(rep(NA_character_, n))
  gmm <- recursive_gmm(df$log_length)
  # Handle empty or failed GMM fit gracefully
  if(nrow(gmm) == 0 || !"mean_log10" %in% colnames(gmm)) return(rep(NA_character_, n))
  map_chr(df$log_length, function(x) {
    idx <- which.min(abs(x - gmm$mean_log10))
    if(length(idx)==0) return(NA_character_)
    gmm$submode[idx]
  })
}

# --- file_acc  [from SVCM.Rmd] ---
file_acc <- function(x) {
  sub("\\.fna$", "", basename(x), ignore.case = TRUE)
}

# --- find_sync_fasta  [from SVCM.Rmd] ---
find_sync_fasta <- function(qid, sync_dir) {
  # Expected: <qid>.txt, but also allow .fna/.fa/.fasta and gz variants.
  pats <- c(
    paste0("^", qid, "\\.txt$"),
    paste0("^", qid, "\\.(fa|fna|fasta)$"),
    paste0("^", qid, "\\.(txt|fa|fna|fasta)\\.gz$")
  )
  files <- list.files(sync_dir, full.names = TRUE)
  hits <- character(0)
  for (p in pats) {
    hits <- c(hits, files[grepl(p, basename(files), ignore.case = TRUE)])
  }
  hits <- unique(hits)
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

# --- get_present_qids  [from SVCM.Rmd] ---
get_present_qids <- function(M_mod, fid) {
  if (!exists("M_mod", inherits = TRUE)) {
    return(character(0))
  }
  if (is.null(fid) || length(fid) == 0 || is.na(fid) || !(fid %in% colnames(M_mod))) {
    return(character(0))
  }
  rownames(M_mod)[as.vector(M_mod[, fid] != 0)]
}

# --- make_feature_labels  [from SVCM.Rmd] ---
make_feature_labels <- function(major_tbl) {
  base_lab <- sprintf("%.2f Mbp Inversion", major_tbl$representative_width_median / 1e6)
  base_lab[!is.finite(major_tbl$representative_width_median)] <- paste0("Inversion ", seq_len(sum(!is.finite(major_tbl$representative_width_median))))

  dup_idx <- ave(seq_along(base_lab), base_lab, FUN = seq_along)
  label <- ifelse(
    duplicated(base_lab) | duplicated(base_lab, fromLast = TRUE),
    paste0(base_lab, " [", dup_idx, "]"),
    base_lab
  )

  stats::setNames(label, major_tbl$feature_id)
}

# --- normalize_label  [from SVCM.Rmd] ---
normalize_label <- function(x) {
  x <- basename(x)
  x <- sub("\\.(fa|fna|fasta|txt)$", "", x, ignore.case = TRUE)
  x <- sub("\\.aln$", "", x, ignore.case = TRUE)
  x <- sub("\\.ref$", "", x, ignore.case = TRUE)
  x <- sub("\\s+.*$", "", x)
  x <- sub("\\|.*$", "", x)
  x
}

# --- pick_compact_by_xy  [from SVCM.Rmd] ---
pick_compact_by_xy <- function(df, n_pick, x_col = "MDS1", y_col = "MDS2") {
  if (nrow(df) == 0) {
    return(character(0))
  }
  if (nrow(df) <= n_pick) {
    return(df$qid)
  }

  xy <- as.matrix(df[, c(x_col, y_col), drop = FALSE])
  dmat <- as.matrix(stats::dist(xy))

  ## medoid = point with smallest total distance to all others
  medoid_idx <- which.min(rowSums(dmat))

  dist_to_medoid <- dmat[, medoid_idx]
  ord <- order(dist_to_medoid, decreasing = FALSE)

  df$qid[ord[seq_len(n_pick)]]
}

# --- pick_diverse_by_xy  [from SVCM.Rmd] ---
pick_diverse_by_xy <- function(df, n_pick, x_col = "MDS1", y_col = "MDS2") {
  if (nrow(df) == 0) {
    return(character(0))
  }
  if (nrow(df) <= n_pick) {
    return(df$qid)
  }

  xy <- as.matrix(df[, c(x_col, y_col), drop = FALSE])
  d0 <- as.matrix(stats::dist(xy))

  ## start at medoid
  start_idx <- which.min(rowSums(d0))
  picked <- start_idx

  while (length(picked) < n_pick) {
    remain <- setdiff(seq_len(nrow(df)), picked)
    d_to_picked <- vapply(remain, function(i) min(d0[i, picked]), numeric(1))
    next_idx <- remain[which.max(d_to_picked)]
    picked <- c(picked, next_idx)
  }

  df$qid[picked]
}

# --- pick_state_compact  [from SVCM.Rmd] ---
pick_state_compact <- function(df,
                               n_pick,
                               require_single_cluster = TRUE,
                               x_col = "MDS1",
                               y_col = "MDS2") {
  if (nrow(df) == 0) {
    return(character(0))
  }

  df1 <- df

  if ("cluster_plot" %in% names(df1)) {
    cl_tab <- df1 |>
      dplyr::count(cluster_plot, sort = TRUE)

    if (nrow(cl_tab) > 0) {
      dom_cluster <- cl_tab$cluster_plot[1]
      df_dom <- df1 |>
        dplyr::filter(cluster_plot == dom_cluster)

      ## for the figure, it is usually better to stay inside one compact cluster
      if (require_single_cluster && nrow(df_dom) >= min(3, n_pick)) {
        df1 <- df_dom
      }
    }
  }

  pick_compact_by_xy(df1, n_pick = n_pick, x_col = x_col, y_col = y_col)
}

# --- pick_state_reps  [from SVCM.Rmd] ---
pick_state_reps <- function(df,
                            n_pick,
                            prefer_cluster = TRUE,
                            x_col = "MDS1",
                            y_col = "MDS2") {
  if (nrow(df) == 0) {
    return(character(0))
  }

  df1 <- df

  if (prefer_cluster && "cluster_plot" %in% names(df1)) {
    cl_tab <- df1 |>
      dplyr::count(cluster_plot, sort = TRUE)

    if (nrow(cl_tab) > 0) {
      dom_cluster <- cl_tab$cluster_plot[1]
      df_dom <- df1 |>
        dplyr::filter(cluster_plot == dom_cluster)

      ## only restrict to dominant cluster if it is not too tiny
      if (nrow(df_dom) >= max(4, ceiling(n_pick * 0.7))) {
        df1 <- df_dom
      }
    }
  }

  ## choose diverse representatives within the retained state context
  pick_diverse_by_xy(df1, n_pick = n_pick, x_col = x_col, y_col = y_col)
}

# --- pick_target_module  [from SVCM.Rmd] ---
pick_target_module <- function(module_meta_tbl,
                               target_mbp,
                               tol_frac = 0.30,
                               min_carriers = 5L,
                               exclude_ids = character(0)) {
  cand <- module_meta_tbl |>
    dplyr::filter(
      !feature_id %in% exclude_ids,
      is.finite(width_mbp),
      carriers_module >= min_carriers,
      dplyr::between(width_mbp, target_mbp * (1 - tol_frac), target_mbp * (1 + tol_frac))
    ) |>
    dplyr::mutate(
      abs_diff = abs(width_mbp - target_mbp)
    ) |>
    dplyr::arrange(abs_diff, dplyr::desc(carriers_module), feature_id)

  if (nrow(cand) == 0) {
    return(NA_character_)
  }
  cand$feature_id[1]
}

# --- plot_candidate_space  [from SVCM.Rmd] ---
plot_candidate_space <- function(motif_tbl, variant_name = NULL) {
  df <- motif_tbl |>
    dplyr::filter(
      is.finite(width_median),
      carriers > 0
    )

  if (!is.null(variant_name)) {
    df <- df |>
      dplyr::filter(variant_collapsed == variant_name)
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = width_kb, y = carriers)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = if (is.null(variant_name)) gsub("_", " ", species, fixed = TRUE) else paste(gsub("_", " ", species, fixed = TRUE), "-", variant_name),
      x = "Median motif width (kb, log10)",
      y = "Carrier count (log10)"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "italic", hjust = 0.5)
    )
}

# --- plot_presence_panel  [from SVCM.Rmd] ---
plot_presence_panel <- function(plot_df, feature_id, title_text) {
  carriers <- plot_df$qid[plot_df[[feature_id]] == 1L]

  ggplot2::ggplot(plot_df, ggplot2::aes(MDS1, MDS2)) +
    ggplot2::geom_point(color = "grey80", alpha = 0.35, size = 1.5) +
    ggplot2::geom_point(
      data = dplyr::filter(plot_df, qid %in% carriers),
      shape = 21,
      fill = "black",
      color = "black",
      alpha = 0.95,
      stroke = 0.25,
      size = 2.4
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "italic", hjust = 0.5)
    )
}

# --- reciprocal_overlap  [from Ch2.Rmd] ---
reciprocal_overlap <- function(i, j) {
  ov_start <- max(r1v2$start[i], r1v2$start[j])
  ov_end   <- min(r1v2$end[i],   r1v2$end[j])
  ov_len   <- max(0, ov_end - ov_start + 1)
  len_i    <- r1v2$end[i] - r1v2$start[i] + 1
  len_j    <- r1v2$end[j] - r1v2$start[j] + 1
  c(ov_len / len_i, ov_len / len_j)
}

# --- run_result1_species  [from SVCM.Rmd] ---
run_result1_species <- function(
    species,
    clade,
    stepC_fp,
    out_dir,
    mash_bundle_fp = NULL,
    identity_dir = NULL,
    prefer_full_matrix = TRUE,
    collapse_modules = TRUE,
    module_jaccard_thresh = 0.90,
    max_states = 4L,
    jaccard_exclusive_max = 0.25,
    min_module_carriers = 10L,
    dbscan_minPts = 5L,
    dbscan_min_cluster_n = 5L,
    nperm = 199L
) {
  stopifnot(file.exists(stepC_fp))

  if (is.null(mash_bundle_fp)) {
    if (is.null(identity_dir)) {
      stop("Provide mash_bundle_fp or identity_dir.")
    }
    mash_bundle_fp <- find_mash_bundle_fp(identity_dir)
  }
  stopifnot(file.exists(mash_bundle_fp))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat_rule(paste("Result 1:", species))
  cat("Step C bundle:\n", stepC_fp, "\n", sep = "")
  cat("Mash bundle:\n", mash_bundle_fp, "\n", sep = "")

  b <- readRDS(stepC_fp)
  mb <- readRDS(mash_bundle_fp)

  M0 <- choose_M_matrix(b, prefer_full_matrix = prefer_full_matrix)
  motif_meta0 <- make_motif_meta(b, M0)
  D_mash <- extract_mash_dist(mb)

  h <- harmonize_M_and_D(M0, D_mash)
  M_int <- h$M
  Dmat_int <- h$Dmat
  D_int <- h$D

  ord <- make_ord_df(D_int)
  ord_df <- ord$df
  ord_df$qid <- rownames(Dmat_int)

  db_out <- summarize_dbscan_grid(ord_df[, c("MDS1", "MDS2")], minPts = dbscan_minPts)
  mid_i <- ceiling(nrow(db_out) / 2)
  eps_mid <- db_out$eps[mid_i]

  fit_mid <- dbscan::dbscan(ord_df[, c("MDS1", "MDS2")], eps = eps_mid, minPts = dbscan_minPts)
  ord_df$cluster_raw <- as.character(fit_mid$cluster)
  ord_df$cluster_plot <- collapse_dbscan_clusters(ord_df$cluster_raw, min_cluster_n = dbscan_min_cluster_n)

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

  if (nrow(motif_meta_int) == 0) {
    stop("No inversion-class motifs found after overlap harmonization.")
  }

  if (collapse_modules) {
    mods <- build_modules_from_jaccard(
      M_mat = M_int,
      motif_tbl = motif_meta_int,
      jaccard_thresh = module_jaccard_thresh
    )

    M_feat <- mods$M_module
    membership_tbl <- mods$membership_tbl
    feature_tbl <- mods$module_tbl |>
      dplyr::rename(feature_id = module_id)
  } else {
    M_feat <- M_int[, motif_meta_int$motif_uid, drop = FALSE]
    membership_tbl <- motif_meta_int |>
      dplyr::transmute(
        motif_uid = motif_uid,
        module_id = motif_uid
      )
    feature_tbl <- motif_meta_int |>
      dplyr::transmute(
        feature_id = motif_uid,
        n_motifs = 1L,
        carriers_module = carriers_overlap,
        representative_motif = motif_uid,
        representative_variant = variant,
        representative_width_median = width_median
      )
  }

  scan_tbl <- lapply(colnames(M_feat), function(fid) {
    scan_one_feature(
      feature_id = fid,
      M_bin = M_feat,
      Dmat = Dmat_int,
      cluster_raw = ord_df$cluster_raw,
      nperm = nperm
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::left_join(feature_tbl, by = "feature_id") |>
    dplyr::arrange(dplyr::desc(delta_z), dplyr::desc(carriers_overlap), feature_id)

  major_ids <- select_major_features(
    scan_tbl = scan_tbl,
    M_bin = M_feat,
    max_states = max_states,
    jaccard_exclusive_max = jaccard_exclusive_max,
    min_carriers = min_module_carriers
  )

  if (length(major_ids) == 0) {
    major_ids <- scan_tbl |>
      dplyr::filter(carriers_overlap >= min_module_carriers) |>
      dplyr::slice(1:min(3, dplyr::n())) |>
      dplyr::pull(feature_id)
  }

  major_tbl <- scan_tbl |>
    dplyr::filter(feature_id %in% major_ids) |>
    dplyr::arrange(match(feature_id, major_ids))

  label_map <- make_feature_labels(major_tbl)

  major_tbl <- major_tbl |>
    dplyr::mutate(
      label = unname(label_map[feature_id])
    )

  state_df <- assign_states(M_feat, major_ids, label_map)
  plot_df_all <- ord_df |>
    dplyr::left_join(state_df, by = "qid")

  lead_id <- major_tbl$feature_id[1]
  plot_df_feat <- plot_df_all
  plot_df_feat[[lead_id]] <- as.integer(M_feat[plot_df_feat$qid, lead_id] != 0)

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
      title = species_label(species),
      x = paste0("MDS1 (", sprintf("%.1f", ord$ord_pct["MDS1"]), "%)"),
      y = paste0("MDS2 (", sprintf("%.1f", ord$ord_pct["MDS2"]), "%)"),
      fill = "DBSCAN"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "italic", hjust = 0.5)
    )

  p_lead <- plot_presence_panel(
    plot_df = plot_df_feat,
    feature_id = lead_id,
    title_text = species_label(species)
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
      title = species_label(species),
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
    dbscan_eps = eps_mid,
    plot_df_all = plot_df_all,
    major_tbl = major_tbl,
    scan_tbl = scan_tbl,
    feature_tbl = feature_tbl,
    membership_tbl = membership_tbl,
    params = list(
      prefer_full_matrix = prefer_full_matrix,
      collapse_modules = collapse_modules,
      module_jaccard_thresh = module_jaccard_thresh,
      max_states = max_states,
      jaccard_exclusive_max = jaccard_exclusive_max,
      min_module_carriers = min_module_carriers,
      dbscan_minPts = dbscan_minPts,
      dbscan_min_cluster_n = dbscan_min_cluster_n,
      nperm = nperm
    )
  )

  saveRDS(plot_obj, file.path(out_dir, "plot_data_major_states.rds"))

  ggplot2::ggsave(
    filename = file.path(out_dir, "result1_panels.png"),
    plot = combo_plot,
    width = 8,
    height = 13,
    dpi = 300
  )

  cat_rule("Selected major inversion states")
  print(
    major_tbl |>
      dplyr::select(
        feature_id,
        label,
        carriers_overlap,
        representative_motif,
        representative_width_median,
        delta_z,
        pct_dom_nonnoise,
        outside_dom_carriers_nonnoise
      )
  )

  cat_rule("DBSCAN summary")
  print(db_out, row.names = FALSE)

  invisible(
    list(
      plot = combo_plot,
      plot_obj = plot_obj
    )
  )
}

# --- select_major_features  [from SVCM.Rmd] ---
select_major_features <- function(scan_tbl, M_bin,
                                  max_states = 4L,
                                  jaccard_exclusive_max = 0.25,
                                  min_carriers = 10L) {
  scan_tbl_filt <- scan_tbl |>
    dplyr::filter(
      carriers_overlap >= min_carriers,
      is.finite(delta_z)
    ) |>
    dplyr::arrange(dplyr::desc(delta_z), dplyr::desc(carriers_overlap), feature_id)

  if (nrow(scan_tbl_filt) == 0) {
    return(character(0))
  }

  feature_sets <- lapply(colnames(M_bin), function(fid) {
    rownames(M_bin)[as.vector(M_bin[, fid] != 0)]
  })
  names(feature_sets) <- colnames(M_bin)

  chosen <- character(0)

  for (fid in scan_tbl_filt$feature_id) {
    if (length(chosen) >= max_states) {
      break
    }

    keep <- TRUE
    if (length(chosen) > 0) {
      for (prev in chosen) {
        jac <- jaccard_sets(feature_sets[[fid]], feature_sets[[prev]])
        if (is.finite(jac) && jac > jaccard_exclusive_max) {
          keep <- FALSE
          break
        }
      }
    }

    if (keep) {
      chosen <- c(chosen, fid)
    }
  }

  chosen
}

# --- stage_one_fasta  [from SVCM.Rmd] ---
stage_one_fasta <- function(qid, sync_dir, input_dir) {
  src <- find_sync_fasta(qid, sync_dir)
  if (is.na(src) || !file.exists(src)) return(list(qid = qid, staged = FALSE, src = NA_character_))

  dest <- file.path(input_dir, paste0(qid, ".fna"))

  # If gz, decompress to .fna (parsnp is often picky)
  if (grepl("\\.gz$", src, ignore.case = TRUE)) {
    ok <- tryCatch({
      con_in <- gzfile(src, open = "rt")
      on.exit(close(con_in), add = TRUE)
      txt <- readLines(con_in)
      writeLines(txt, dest)
      TRUE
    }, error = function(e) FALSE)
    return(list(qid = qid, staged = ok, src = src))
  }

  ok <- file.copy(src, dest, overwrite = TRUE)
  list(qid = qid, staged = ok, src = src)
}
