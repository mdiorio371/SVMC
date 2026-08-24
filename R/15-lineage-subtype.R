# 15-lineage-subtype.R
# Lineage-association scoring via MASH-subtype purity (Chapter 4 method).
# Extracted from Ch4.Rmd. Produces subtype_purity = the lineage-association score.
# subtype_purity = (carriers in dominant MASH subtype) / (total carriers).

# --- upper_tri_values ---
upper_tri_values <- function(mat) {
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(numeric(0))
  mat[upper.tri(mat, diag = FALSE)]
}

# --- extract_dist_object  [from Ch4.Rmd:3788]  (brace +0) ---
extract_dist_object <- function(x) {
  if (inherits(x, "dist")) return(x)

  if (is.list(x)) {
    for (nm in c("dist", "D", "mash_dist", "mash_dist_obj")) {
      if (nm %in% names(x) && inherits(x[[nm]], "dist")) {
        return(x[[nm]])
      }
    }

    if ("dist_matrix" %in% names(x) && is.matrix(x$dist_matrix)) {
      return(stats::as.dist(x$dist_matrix))
    }
  }

  stop("Could not find a dist object in the supplied mash object.", call. = FALSE)
}

# --- compute_set_metrics_once  [from Ch4.Rmd:3811]  (brace +0) ---
compute_set_metrics_once <- function(Dmat, subtype_vec, carrier_idx) {
  all_idx <- seq_len(nrow(Dmat))
  noncarrier_idx <- setdiff(all_idx, carrier_idx)

  within_vals <- upper_tri_values(Dmat[carrier_idx, carrier_idx, drop = FALSE])
  between_vals <- if (length(noncarrier_idx) > 0L) {
    as.vector(Dmat[carrier_idx, noncarrier_idx, drop = FALSE])
  } else {
    numeric(0)
  }

  carrier_subtypes <- subtype_vec[carrier_idx]
  subtype_tab <- sort(table(carrier_subtypes), decreasing = TRUE)

  dominant_subtype <- names(subtype_tab)[1]
  dominant_n <- as.integer(subtype_tab[1])
  purity <- dominant_n / length(carrier_idx)

  within_median <- if (length(within_vals) > 0L) stats::median(within_vals, na.rm = TRUE) else NA_real_
  between_median <- if (length(between_vals) > 0L) stats::median(between_vals, na.rm = TRUE) else NA_real_
  separation <- if (is.finite(within_median) && is.finite(between_median)) between_median - within_median else NA_real_

  list(
    dominant_subtype = dominant_subtype,
    dominant_n = dominant_n,
    purity = purity,
    within_median = within_median,
    between_median = between_median,
    separation = separation
  )
}

# --- permute_same_size_metrics  [from Ch4.Rmd:3843]  (brace +0) ---
permute_same_size_metrics <- function(Dmat, subtype_vec, n_carriers, n_perm = 1000L, seed = 1L) {
  set.seed(seed)
  n_total <- nrow(Dmat)

  out <- replicate(n_perm, {
    idx <- sort(sample.int(n_total, size = n_carriers, replace = FALSE))
    met <- compute_set_metrics_once(Dmat, subtype_vec, idx)
    c(
      purity = met$purity,
      within_median = met$within_median,
      separation = met$separation
    )
  })

  as.data.frame(t(out))
}

# --- derive_mash_subtypes  [from Ch4.Rmd:3860]  (brace +0) ---
derive_mash_subtypes <- function(mash_obj, k = NULL, k_min = 2L, k_max = 8L, hc_method = "ward.D2") {
  D <- extract_dist_object(mash_obj)
  labels <- attr(D, "Labels")
  stopifnot(!is.null(labels), length(labels) >= 3L)

  hc <- stats::hclust(D, method = hc_method)

  if (is.null(k)) {
    candidate_k <- seq.int(k_min, min(k_max, length(labels) - 1L))
    if (length(candidate_k) == 0L) stop("No valid k values for subtype derivation.", call. = FALSE)

    sil_tbl <- lapply(candidate_k, function(kk) {
      grp <- stats::cutree(hc, k = kk)
      if (length(unique(grp)) < 2L) {
        return(tibble::tibble(k = kk, avg_sil = NA_real_))
      }
      sil <- cluster::silhouette(grp, D)
      tibble::tibble(k = kk, avg_sil = mean(sil[, "sil_width"], na.rm = TRUE))
    }) %>%
      dplyr::bind_rows()

    k <- sil_tbl %>%
      dplyr::filter(!is.na(avg_sil)) %>%
      dplyr::arrange(dplyr::desc(avg_sil), k) %>%
      dplyr::slice(1) %>%
      dplyr::pull(k)

    if (length(k) == 0L || is.na(k)) stop("Failed to choose subtype k from silhouette.", call. = FALSE)
  } else {
    sil_tbl <- tibble::tibble(k = k, avg_sil = NA_real_)
  }

  grp <- stats::cutree(hc, k = k)

  subtype_tbl <- tibble::tibble(
    genome_id = labels,
    mash_subtype = paste0("S", grp)
  )

  list(
    k = k,
    subtype_tbl = subtype_tbl,
    silhouette_tbl = sil_tbl,
    hc = hc
  )
}

# --- run_exclusive_motif_mash_subtype_test  [from Ch4.Rmd:3911]  (brace +0) ---
run_exclusive_motif_mash_subtype_test <- function(
  M_all,
  mash_obj,
  exclusivity_tbl,
  exclusive_col = "candidate_exclusive_primary_state_all",
  motif_col = "motif_uid",
  variant_col = "variant",
  carriers_col = "carriers",
  min_carriers = 5L,
  subtype_tbl = NULL,
  k = NULL,
  k_min = 2L,
  k_max = 8L,
  n_perm = 1000L,
  seed = 1L,
  out_dir = NULL
) {
  stopifnot(inherits(M_all, "dgCMatrix"))
  stopifnot(is.data.frame(exclusivity_tbl))
  stopifnot(all(c(motif_col, variant_col, carriers_col, exclusive_col) %in% names(exclusivity_tbl)))

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  D <- extract_dist_object(mash_obj)
  mash_labels <- attr(D, "Labels")
  stopifnot(!is.null(mash_labels))

  subtype_info <- NULL
  if (is.null(subtype_tbl)) {
    subtype_info <- derive_mash_subtypes(
      mash_obj = D,
      k = k,
      k_min = k_min,
      k_max = k_max
    )
    subtype_tbl <- subtype_info$subtype_tbl
  }

  common_genomes <- Reduce(intersect, list(
    rownames(M_all),
    mash_labels,
    subtype_tbl$genome_id
  ))

  if (length(common_genomes) < 3L) {
    stop("Too few common genomes among M_all, Mash labels, and subtype table.", call. = FALSE)
  }

  M_use <- M_all[common_genomes, , drop = FALSE]
  Dmat <- as.matrix(D)[common_genomes, common_genomes, drop = FALSE]
  subtype_vec <- subtype_tbl$mash_subtype[match(common_genomes, subtype_tbl$genome_id)]
  names(subtype_vec) <- common_genomes

  cand_tbl <- exclusivity_tbl %>%
    dplyr::transmute(
      motif_uid = as.character(.data[[motif_col]]),
      variant = as.character(.data[[variant_col]]),
      carriers = as.integer(.data[[carriers_col]]),
      is_exclusive = as.logical(.data[[exclusive_col]])
    ) %>%
    dplyr::filter(is_exclusive) %>%
    dplyr::filter(!is.na(carriers), carriers >= min_carriers) %>%
    dplyr::filter(motif_uid %in% colnames(M_use)) %>%
    dplyr::distinct()

  if (nrow(cand_tbl) == 0L) {
    stop("No exclusive motifs passed the carrier threshold in the aligned evaluation cohort.", call. = FALSE)
  }

  result_list <- vector("list", nrow(cand_tbl))

  for (i in seq_len(nrow(cand_tbl))) {
    motif_id <- cand_tbl$motif_uid[i]
    carrier_idx <- which(as.vector(M_use[, motif_id]) > 0)

    if (length(carrier_idx) < min_carriers) next

    obs <- compute_set_metrics_once(Dmat, subtype_vec, carrier_idx)

    noncarrier_idx <- setdiff(seq_len(nrow(Dmat)), carrier_idx)
    dom_noncarrier_n <- sum(subtype_vec[noncarrier_idx] == obs$dominant_subtype)

    fisher_mat <- matrix(
      c(
        obs$dominant_n,
        length(carrier_idx) - obs$dominant_n,
        dom_noncarrier_n,
        length(noncarrier_idx) - dom_noncarrier_n
      ),
      nrow = 2,
      byrow = TRUE
    )

    fisher_p <- suppressWarnings(stats::fisher.test(fisher_mat, alternative = "greater")$p.value)

    null_tbl <- permute_same_size_metrics(
      Dmat = Dmat,
      subtype_vec = subtype_vec,
      n_carriers = length(carrier_idx),
      n_perm = n_perm,
      seed = seed + i
    )

    p_perm_purity <- mean(null_tbl$purity >= obs$purity, na.rm = TRUE)
    p_perm_within <- mean(null_tbl$within_median <= obs$within_median, na.rm = TRUE)   # lower is better
    p_perm_separation <- mean(null_tbl$separation >= obs$separation, na.rm = TRUE)     # higher is better

    result_list[[i]] <- tibble::tibble(
      motif_uid = motif_id,
      variant = cand_tbl$variant[i],
      carriers = length(carrier_idx),
      dominant_subtype = obs$dominant_subtype,
      dominant_subtype_n = obs$dominant_n,
      subtype_purity = obs$purity,
      fisher_p = fisher_p,
      within_median = obs$within_median,
      between_median = obs$between_median,
      separation = obs$separation,
      p_perm_purity = p_perm_purity,
      p_perm_within = p_perm_within,
      p_perm_separation = p_perm_separation
    )
  }

  motif_results_tbl <- dplyr::bind_rows(result_list) %>%
    dplyr::mutate(
      q_fisher = stats::p.adjust(fisher_p, method = "BH"),
      q_perm_purity = stats::p.adjust(p_perm_purity, method = "BH"),
      q_perm_within = stats::p.adjust(p_perm_within, method = "BH"),
      q_perm_separation = stats::p.adjust(p_perm_separation, method = "BH"),
      supported_by_subtype = q_fisher < 0.05 & q_perm_purity < 0.05,
      supported_by_distance = q_perm_within < 0.05 & q_perm_separation < 0.05,
      strong_subtype_marker = supported_by_subtype & supported_by_distance
    ) %>%
    dplyr::arrange(variant, dplyr::desc(subtype_purity), dplyr::desc(separation), motif_uid)

  variant_summary_tbl <- motif_results_tbl %>%
    dplyr::group_by(variant) %>%
    dplyr::summarise(
      n_tested = dplyr::n(),
      n_subtype_supported = sum(supported_by_subtype, na.rm = TRUE),
      n_distance_supported = sum(supported_by_distance, na.rm = TRUE),
      n_strong_markers = sum(strong_subtype_marker, na.rm = TRUE),
      median_purity = stats::median(subtype_purity, na.rm = TRUE),
      median_within = stats::median(within_median, na.rm = TRUE),
      median_separation = stats::median(separation, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_strong_markers), dplyr::desc(median_purity), variant)

  subtype_membership_tbl <- subtype_tbl %>%
    dplyr::filter(genome_id %in% common_genomes) %>%
    dplyr::count(mash_subtype, name = "n_genomes") %>%
    dplyr::arrange(dplyr::desc(n_genomes), mash_subtype)

  out <- list(
    common_genomes = common_genomes,
    subtype_tbl = subtype_tbl,
    subtype_info = subtype_info,
    subtype_membership_tbl = subtype_membership_tbl,
    motif_results_tbl = motif_results_tbl,
    variant_summary_tbl = variant_summary_tbl
  )

  if (!is.null(out_dir)) {
    saveRDS(out, file.path(out_dir, "exclusive_motif_mash_subtype_test.rds"))
    readr::write_tsv(subtype_membership_tbl, file.path(out_dir, "subtype_membership_tbl.tsv"))
    readr::write_tsv(motif_results_tbl, file.path(out_dir, "exclusive_motif_mash_subtype_results.tsv"))
    readr::write_tsv(variant_summary_tbl, file.path(out_dir, "exclusive_motif_mash_subtype_variant_summary.tsv"))
    if (!is.null(subtype_info) && !is.null(subtype_info$silhouette_tbl)) {
      readr::write_tsv(subtype_info$silhouette_tbl, file.path(out_dir, "mash_subtype_silhouette_tbl.tsv"))
    }
  }

  out
}

# --- make_motif_lineage_role_tbl  [from Ch4.Rmd:518]  (brace +0) ---
make_motif_lineage_role_tbl <- function(species_res) {
  excl <- species_res$exclusivity_refined$motif_summary_refined_tbl %>%
    dplyr::rename(
      carriers_excl = carriers
    )

  st <- species_res$subtype_test$motif_results_tbl %>%
    dplyr::rename(
      carriers_test = carriers
    )

  out <- excl %>%
    dplyr::left_join(
      st,
      by = c("motif_uid", "variant")
    ) %>%
    dplyr::mutate(
      # Basic support flags
      context_concentrated = dplyr::coalesce(subtype_purity >= 0.90, FALSE),
      context_moderate     = dplyr::coalesce(subtype_purity >= 0.70, FALSE),

      subtype_supported = dplyr::coalesce(
        q_fisher <= 0.05 | q_perm_purity <= 0.05,
        FALSE
      ),

      cohesion_supported = dplyr::coalesce(q_perm_within <= 0.05, FALSE),

      separation_supported = dplyr::coalesce(
        q_perm_separation <= 0.05 & separation > 0,
        FALSE
      ),

      derivative_like = dplyr::coalesce(
        exclusivity_label_all %in% c("nested_or_derivative"),
        FALSE
      ),

      background_like = dplyr::coalesce(
        exclusivity_label_all %in% c(
          "cooccurring_or_background",
          "cooccurring_with_other_same_class_motifs"
        ),
        FALSE
      ),

      exclusive_primary = dplyr::coalesce(
        candidate_exclusive_primary_state_all,
        FALSE
      ),

      # Ordinal lineage-role classification
      lineage_role = dplyr::case_when(
        # Strong global / subtype-partitioning state
        dplyr::coalesce(strong_subtype_marker, FALSE) ~ "major_lineage_defining",

        exclusive_primary &
          context_concentrated &
          subtype_supported &
          (cohesion_supported | separation_supported) ~ "major_lineage_defining",

        # Strongly localized and cohesive, but not globally subtype-discriminative
        context_concentrated &
          !subtype_supported &
          (cohesion_supported | separation_supported) &
          !derivative_like ~ "nested_lineage_associated",

        # Recurrent and somewhat localized, but spread across contexts
        context_moderate &
          !derivative_like &
          !background_like ~ "multi_context_recurrent",

        # Default weak end
        TRUE ~ "derivative_or_background"
      ),

      lineage_role_ord = dplyr::case_when(
        lineage_role == "derivative_or_background"   ~ 0L,
        lineage_role == "multi_context_recurrent"    ~ 1L,
        lineage_role == "nested_lineage_associated"  ~ 2L,
        lineage_role == "major_lineage_defining"     ~ 3L,
        TRUE ~ NA_integer_
      ),

      lineage_role = factor(
        lineage_role,
        levels = c(
          "derivative_or_background",
          "multi_context_recurrent",
          "nested_lineage_associated",
          "major_lineage_defining"
        ),
        ordered = TRUE
      ),

      lineage_role_reason = dplyr::case_when(
        lineage_role == "major_lineage_defining" &
          dplyr::coalesce(strong_subtype_marker, FALSE) ~
          "exclusive and strongly supported by subtype purity/cohesion/separation",

        lineage_role == "major_lineage_defining" ~
          "exclusive primary state with concentrated subtype context and distance support",

        lineage_role == "nested_lineage_associated" ~
          "highly concentrated and cohesive within one background, but not globally subtype-discriminative",

        lineage_role == "multi_context_recurrent" ~
          "recurrent motif with some context concentration, but not a strong lineage partition",

        TRUE ~
          "derivative, background, or weakly localized recurrent motif"
      )
    )

  out
}
