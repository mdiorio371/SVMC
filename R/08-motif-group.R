# 08-motif-group.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- assign_reciprocal_overlap_groups  [from SVCM.R:2395] ---
assign_reciprocal_overlap_groups <- function(df, min_ro = 0.5) {
  if (nrow(df) == 0L) { df$ro_group <- integer(0); return(df) }
  if (nrow(df) == 1L) { df$ro_group <- 1L; return(df) }

  n <- nrow(df)
  rng <- IRanges::IRanges(start = as.integer(df$start), end = as.integer(df$end))
  hits <- IRanges::findOverlaps(rng, drop.self = TRUE, drop.redundant = TRUE)

  if (!length(hits)) { df$ro_group <- seq_len(n); return(df) }

  qi <- S4Vectors::queryHits(hits); si <- S4Vectors::subjectHits(hits)
  w <- as.integer(df$width)
  ovl <- pmax.int(0L, pmin(IRanges::end(rng)[qi], IRanges::end(rng)[si]) -
                    pmax(IRanges::start(rng)[qi], IRanges::start(rng)[si]) + 1L)
  ro <- pmin(ovl / w[qi], ovl / w[si])
  keep <- ro >= min_ro

  if (!any(keep)) { df$ro_group <- seq_len(n); return(df) }

  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  g <- igraph::add_edges(g, as.vector(rbind(qi[keep], si[keep])))
  df$ro_group <- as.integer(igraph::components(g)$membership)
  df
}

# --- sv_rec2  [from SVCM.R:2462] ---
sv_rec2 <- function(
    sv_tbl,
    ref_acc                = NULL,
    by_rid                 = TRUE,
    size_min_bp            = 2000L,
    ro_thresh_noninsert    = 0.5,
    ro_thresh_inv          = 0.2,
    min_overlap_bp         = 50L,
    match_variant_specific = FALSE,
    debug                  = FALSE
) {
  require(dplyr)
  require(data.table)

  # SV class classifier
  .classify_sv <- function(variant, variant_specific) {
    var_lower <- tolower(variant)
    spec_lower <- tolower(variant_specific)
    if (grepl("translocation", var_lower)) return("trans")
    if (grepl("insertion", spec_lower)) return("ins")
    if (grepl("deletion", spec_lower)) return("del")
    if (grepl("duplication", var_lower)) return("dup")
    if (grepl("inversion", var_lower)) return("inv")
    return("other")
  }

  # Adaptive tolerance functions (scale per SV class)
  .pos_tol_vec <- function(len, cls) {
    sapply(seq_along(len), function(i) {
      l <- len[i]; c <- cls[i]
      if (c == "ins") return(as.integer(pmax(100L, l * 0.2)))
      if (c == "trans") return(as.integer(pmax(500L, l * 0.15)))
      return(as.integer(pmax(50L, l * 0.1)))
    })
  }
  .len_tol_vec <- function(len, cls) as.integer(len * 0.2)
  .bd_tol_vec  <- function(len, cls) as.integer(pmax(1000L, len * 0.05))

  # DSU functions
  dsu_init <- function(n) list(parent = seq_len(n), rank = rep(1L, n))

  dsu_find <- function(dsu, x) {
    if (dsu$parent[x] != x) dsu$parent[x] <- dsu_find(dsu, dsu$parent[x])
    dsu$parent[x]
  }

  dsu_union <- function(dsu, x, y) {
    px <- dsu_find(dsu, x); py <- dsu_find(dsu, y)
    if (px != py) {
      if (dsu$rank[px] < dsu$rank[py]) dsu$parent[px] <- py
      else if (dsu$rank[px] > dsu$rank[py]) dsu$parent[py] <- px
      else { dsu$parent[py] <- px; dsu$rank[px] <- dsu$rank[px] + 1L }
    }
    dsu
  }

  # Main processing
  tbl <- sv_tbl
  if (!is.null(ref_acc)) tbl <- tbl[tbl$rid == ref_acc, , drop = FALSE]
  rid_groups <- if (by_rid) split(tbl, tbl$rid) else list(tbl)

  res <- list()
  gid_offset <- 0L

  for (rg in seq_along(rid_groups)) {
    dat <- rid_groups[[rg]]
    if (nrow(dat) == 0) next

    x <- dat %>%
      mutate(
        len   = as.integer(width),
        mid   = as.integer(refmid),
        sv_id = sprintf("%s|%s|%d-%d|%s", rid, qid, start, end, variant),
        cls   = mapply(.classify_sv, variant, variant_specific, USE.NAMES = FALSE),
        posw  = .pos_tol_vec(len, cls),
        lent  = .len_tol_vec(len, cls),
        bdw   = .bd_tol_vec(len, cls)
      ) %>%
      filter(!is.na(start), !is.na(end), start <= end, len >= size_min_bp)

    if (nrow(x) < 2) {
      x$recurrence_group <- sprintf("rg_%06d", gid_offset + 1L)
      x$group_size <- 1L
      gid_offset <- gid_offset + 1L
      res[[rg]] <- x
      next
    }

    blk_cols <- if (match_variant_specific) c("variant", "variant_specific") else "variant"
    x <- x %>% mutate(blk = do.call(paste, c(across(all_of(blk_cols)), sep = "|")))

    dsu <- dsu_init(nrow(x))
    XDT <- as.data.table(x)

    for (b in unique(XDT$blk)) {
      DT <- XDT[blk == b]
      if (nrow(DT) < 2) next

      D_ins <- copy(DT[cls == "ins"])
      D_oth <- copy(DT[cls != "ins"])

      # Insertions — fast midpoint-proximity + length-similarity path
      if (nrow(D_ins) > 1) {
        D_ins[, idx := .I]
        D_ins[, gidx := match(sv_id, x$sv_id)]
        D_ins[, `:=`(lo = mid - posw, hi = mid + posw)]
        setkey(D_ins, mid)

        C <- D_ins[D_ins, on = .(mid >= lo, mid <= hi), nomatch = 0L, allow.cartesian = TRUE]
        C <- C[i.idx < idx]

        if (nrow(C)) {
          C[, mid_ok := abs(i.mid - mid) <= pmax(i.posw, posw)]
          C[, len_ok := abs(i.len - len) <= pmax(i.lent, lent)]
          E <- C[mid_ok & len_ok, .(a = i.gidx, b = gidx)]
          if (nrow(E)) {
            for (k in seq_len(nrow(E))) dsu <- dsu_union(dsu, E$a[k], E$b[k])
          }
        }
      }

      # Non-insertions — RO + boundary-OR + large-SV special case
      if (nrow(D_oth) > 1) {
        D_oth[, idx := .I]
        D_oth[, gidx := match(sv_id, x$sv_id)]
        D_oth[, `:=`(lo = mid - posw, hi = mid + posw)]
        setkey(D_oth, mid)

        C <- D_oth[D_oth, on = .(mid >= lo, mid <= hi), nomatch = 0L, allow.cartesian = TRUE]
        C <- C[i.idx < idx]

        if (nrow(C)) {
          # Reciprocal overlap
          C[, ovl := pmax(0L, pmin(i.end, end) - pmax(i.start, start) + 1L)]
          C[, ro1 := ovl / i.len]
          C[, ro2 := ovl / len]

          # Adaptive RO threshold for translocations
          C[, ro_thr := ifelse(cls == "trans" & i.cls == "trans", 0.3, ro_thresh_noninsert)]
          C[, ro_ok := (ovl >= min_overlap_bp) & (ro1 >= ro_thr) & (ro2 >= ro_thr)]

          # Tolerances
          C[, posw_pair := pmax(i.posw, posw)]
          C[, mid_ok := abs(i.mid - mid) <= posw_pair]
          C[, lent_pair := pmax(i.lent, lent)]
          C[, len_ok := abs(i.len - len) <= lent_pair]

          # Boundary check for translocations (OR logic: either start or end matches)
          C[, bdw_pair := pmax(i.bdw, bdw)]
          C[, boundary_ok := (cls == "trans" & i.cls == "trans") &
              ((abs(i.start - start) <= bdw_pair) | (abs(i.end - end) <= bdw_pair))]

          # Large SV special case
          C[, large_sv_ok := (i.len > 500000L | len > 500000L) &
              abs(i.mid - mid) <= pmax(50000L, as.integer(0.1 * pmax(i.len, len)))]

          # Match decision: (ro_ok | boundary_ok | large_sv_ok) AND mid_ok AND len_ok
          E <- C[(ro_ok | boundary_ok | large_sv_ok) & mid_ok & len_ok, .(a = i.gidx, b = gidx)]

          if (debug && any(C$cls == "trans")) {
            cat("\nTranslocation comparisons:\n")
            print(C[cls == "trans", .(i.variant_specific, variant_specific,
                                       ro_ok, boundary_ok, mid_ok, len_ok)])
          }

          if (nrow(E)) {
            for (k in seq_len(nrow(E))) dsu <- dsu_union(dsu, E$a[k], E$b[k])
          }
        }
      }
    }

    comp <- vapply(seq_len(nrow(x)), function(i) dsu_find(dsu, i), integer(1))
    cf <- as.integer(factor(comp))
    x$recurrence_group <- sprintf("rg_%06d", gid_offset + cf)
    gid_offset <- gid_offset + length(unique(cf))

    res[[rg]] <- as.data.frame(x)
  }

  # Combine and reorder by group size
  result_df <- bind_rows(res)

  if (nrow(result_df) > 0) {
    group_order <- result_df %>%
      group_by(recurrence_group) %>%
      summarise(size = n(), .groups = 'drop') %>%
      arrange(desc(size)) %>%
      mutate(new_group = sprintf("rg_%06d", seq_len(n())))

    result_df <- result_df %>%
      left_join(group_order[, c("recurrence_group", "new_group", "size")],
                by = "recurrence_group") %>%
      mutate(
        recurrence_group = new_group,
        group_size = size
      ) %>%
      select(rid, qid, species, variant, variant_specific,
             start, end, refmid, width, len, mid, sv_id,
             recurrence_group, group_size)
  }

  return(result_df)
}

# --- jaccard_sets  [from SVCM_result1_helpers.R:199] ---
jaccard_sets <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(0)
  u <- length(union(a, b))
  if (u == 0) return(0)
  length(intersect(a, b)) / u
}

# --- build_modules_from_jaccard  [from SVCM_result1_helpers.R:225] ---
build_modules_from_jaccard <- function(M_mat, motif_tbl,
                                       jaccard_thresh = 0.90) {
  cand_uids <- colnames(M_mat)
  if (length(cand_uids) == 0) {
    stop("M_mat has no columns (no candidate motifs).")
  }

  # Carrier sets per motif
  carrier_sets <- lapply(seq_along(cand_uids), function(j) {
    rownames(M_mat)[which(M_mat[, j] != 0)]
  })
  names(carrier_sets) <- cand_uids

  # Pairwise Jaccard — keep only edges >= threshold
  edges <- list()
  n <- length(cand_uids)
  if (n > 1) {
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        jac <- jaccard_sets(carrier_sets[[i]], carrier_sets[[j]])
        if (jac >= jaccard_thresh) {
          edges[[length(edges) + 1L]] <- c(cand_uids[i], cand_uids[j], jac)
        }
      }
    }
  }

  if (length(edges) > 0L) {
    edge_df <- do.call(rbind, edges) |> as.data.frame(stringsAsFactors = FALSE)
    names(edge_df) <- c("from", "to", "jaccard")
    edge_df$jaccard <- as.numeric(edge_df$jaccard)
    g <- igraph::graph_from_data_frame(edge_df[, c("from", "to")],
                                       vertices = cand_uids,
                                       directed = FALSE)
  } else {
    edge_df <- data.frame(from = character(0), to = character(0),
                          jaccard = numeric(0), stringsAsFactors = FALSE)
    g <- igraph::make_empty_graph(directed = FALSE)
    g <- igraph::add_vertices(g, length(cand_uids), name = cand_uids)
  }

  comp <- igraph::components(g)$membership
  comp_id <- as.integer(factor(comp))                  # contiguous 1..K
  module_id <- sprintf("INV_MOD_%03d", comp_id)

  membership_tbl <- tibble::tibble(
    motif_uid = names(comp),
    module_id = module_id[match(names(comp), cand_uids)]
  )

  # M_module: row-wise OR of motifs sharing a module
  unique_mods <- sort(unique(membership_tbl$module_id))
  M_module <- matrix(0L,
                     nrow = nrow(M_mat),
                     ncol = length(unique_mods),
                     dimnames = list(rownames(M_mat), unique_mods))
  for (mod in unique_mods) {
    motifs_in_mod <- membership_tbl$motif_uid[membership_tbl$module_id == mod]
    cols <- which(colnames(M_mat) %in% motifs_in_mod)
    if (length(cols) == 1L) {
      M_module[, mod] <- as.integer(M_mat[, cols] != 0)
    } else {
      M_module[, mod] <- as.integer(rowSums(M_mat[, cols, drop = FALSE] != 0) > 0)
    }
  }

  # Module summary table
  mod_carrier_counts <- colSums(M_module != 0)

  # Pick representative motif per module (most-carried)
  motif_tbl_loc <- motif_tbl
  if (!"carriers_overlap" %in% names(motif_tbl_loc)) {
    motif_tbl_loc$carriers_overlap <- NA_integer_
  }
  if (!"representative_width_median" %in% names(motif_tbl_loc) &&
      "width_median" %in% names(motif_tbl_loc)) {
    motif_tbl_loc$representative_width_median <- motif_tbl_loc$width_median
  }
  if (!"variant_collapsed" %in% names(motif_tbl_loc)) {
    motif_tbl_loc$variant_collapsed <- NA_character_
  }

  module_tbl <- membership_tbl %>%
    dplyr::left_join(motif_tbl_loc, by = "motif_uid") %>%
    dplyr::group_by(.data$module_id) %>%
    dplyr::arrange(dplyr::desc(.data$carriers_overlap),
                   .data$motif_uid, .by_group = TRUE) %>%
    dplyr::summarise(
      n_motifs_in_module          = dplyr::n(),
      representative_motif        = dplyr::first(.data$motif_uid),
      representative_width_median = dplyr::first(.data$representative_width_median),
      variant_collapsed           = dplyr::first(.data$variant_collapsed),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      carriers_module = as.integer(mod_carrier_counts[.data$module_id])
    )

  list(
    M_module       = M_module,
    membership_tbl = membership_tbl,
    module_tbl     = module_tbl,
    edge_df        = edge_df
  )
}

# --- make_motif_meta  [from SVCM_result1_helpers.R:117] ---
make_motif_meta <- function(b, M) {
  motif_ids <- colnames(M)

  if (!is.null(b$motif_meta) && "motif_uid" %in% names(b$motif_meta)) {
    meta <- b$motif_meta %>%
      dplyr::filter(.data$motif_uid %in% motif_ids)
  } else {
    meta <- tibble::tibble(motif_uid = motif_ids)
  }

  # variant_collapsed
  if (!"variant_collapsed" %in% names(meta)) {
    pc <- if ("plot_class" %in% names(meta)) meta$plot_class else NA_character_
    meta$variant_collapsed <- dplyr::case_when(
      grepl("Inversion|Structural_rearrangement|Substructural_inversion", pc) ~ "Inversion",
      grepl("Translocation", pc) ~ "Translocation",
      grepl("Duplication",    pc) ~ "Duplication",
      grepl("Insertion|Deletion|Indel", pc) ~ "Indel",
      TRUE ~ pc
    )
  }

  # representative_width_median
  if (!"representative_width_median" %in% names(meta) &&
      "width_median" %in% names(meta)) {
    meta$representative_width_median <- meta$width_median
  }

  # representative_motif
  if (!"representative_motif" %in% names(meta)) {
    meta$representative_motif <- meta$motif_uid
  }

  # carriers_overlap (used by later code)
  if (!"carriers_overlap" %in% names(meta) &&
      "n_carriers" %in% names(meta)) {
    meta$carriers_overlap <- meta$n_carriers
  }
  meta
}

# --- choose_M_matrix  [from SVCM_result1_helpers.R:107] ---
choose_M_matrix <- function(b, prefer_full_matrix = TRUE) {
  if (prefer_full_matrix && !is.null(b$M_all)) return(b$M_all)
  if (!is.null(b$M_sel)) return(b$M_sel)
  if (!is.null(b$M_all)) return(b$M_all)
  stop("Step C bundle has neither $M_all nor $M_sel.")
}

# --- harmonize_M_and_D  [from SVCM_result1_helpers.R:161] ---
harmonize_M_and_D <- function(M, D, normalize_extensions = TRUE) {
  normalize_id <- function(x) {
    if (!normalize_extensions) return(as.character(x))
    x <- basename(as.character(x))
    x <- sub("\\.gz$",    "", x)
    x <- sub("\\.fasta$", "", x)
    x <- sub("\\.fna$",   "", x)
    x <- sub("\\.fa$",    "", x)
    x <- sub("\\.txt$",   "", x)
    x
  }

  Dmat <- if (inherits(D, "dist")) as.matrix(D) else as.matrix(D)

  m_ids <- normalize_id(rownames(M))
  d_ids <- normalize_id(rownames(Dmat))

  rownames(M)    <- m_ids
  rownames(Dmat) <- d_ids
  colnames(Dmat) <- d_ids

  shared <- intersect(m_ids, d_ids)
  if (length(shared) < 5) {
    stop(sprintf(
      "Too few shared IDs between motif matrix (n=%d) and distance matrix (n=%d): shared=%d. Check ID normalization.",
      length(m_ids), length(d_ids), length(shared)
    ))
  }

  M_h    <- M[shared, , drop = FALSE]
  Dmat_h <- Dmat[shared, shared, drop = FALSE]

  list(M = M_h, Dmat = Dmat_h, D = stats::as.dist(Dmat_h))
}

# --- mean_within_dist  [from SVCM_result1_helpers.R:209] ---
mean_within_dist <- function(idx, Dmat) {
  if (length(idx) < 2) return(NA_real_)
  sub <- Dmat[idx, idx, drop = FALSE]
  mean(sub[upper.tri(sub)], na.rm = TRUE)
}

# --- mean_between_dist  [from SVCM_result1_helpers.R:215] ---
mean_between_dist <- function(idx1, idx2, Dmat) {
  if (length(idx1) == 0 || length(idx2) == 0) return(NA_real_)
  mean(Dmat[idx1, idx2], na.rm = TRUE)
}

# --- .compute_width  [from ch4_functions.R:1324] ---
.compute_width <- function(df) {
  if ("width" %in% names(df) && is.numeric(df$width)) {
    as.numeric(df$width)
  } else if (all(c("start","end") %in% names(df))) {
    as.numeric(abs(df$end - df$start) + 1)
  } else {
    stop("sv_rec must contain either numeric 'width' or ('start','end').")
  }
}

# --- .make_presence_matrix  [from ch4_functions.R:1309] ---
.make_presence_matrix <- function(df, sample_col = "qid", event_col = "event_id") {
  stopifnot(all(c(sample_col, event_col) %in% names(df)))
  key <- df %>% distinct(.data[[sample_col]], .data[[event_col]])
  samp <- factor(key[[sample_col]])
  ev   <- factor(key[[event_col]])
  i <- as.integer(samp)
  j <- as.integer(ev)
  x <- rep(1L, nrow(key))
  P <- sparseMatrix(i = i, j = j, x = x,
                    dims = c(nlevels(samp), nlevels(ev)),
                    dimnames = list(levels(samp), levels(ev)))
  P
}

# --- svcm_ro_bin  [from SVCM6.R:1256] ---
svcm_ro_bin <- function(ro) {
  dplyr::case_when(
    is.na(ro)            ~ NA_character_,
    ro >= 0.90           ~ ">=90%",
    ro >= 0.50 & ro < 0.90~ "50–90%",
    ro < 0.50            ~ "<50%"
  )
}

# --- svcm_best_cds_ro  [from SVCM6.R:1267] ---
svcm_best_cds_ro <- function(sv_tbl, cds_gr,
                             sv_id_cols = c("rid","qid","start","end","width","variant","variant_specific","species")) {
  
  sv <- tibble::as_tibble(sv_tbl)
  
  # ensure required cols exist
  need <- c("rid","start","end","width")
  miss <- setdiff(need, names(sv))
  if (length(miss)) stop("SV table missing cols: ", paste(miss, collapse=", "))
  
  # row ids for stable mapping
  sv$.sv_i <- seq_len(nrow(sv))
  
  # build SV GRanges
  sv_gr <- GenomicRanges::GRanges(
    seqnames = sv$rid,
    ranges   = IRanges::IRanges(start = as.integer(sv$start), end = as.integer(sv$end))
  )
  S4Vectors::mcols(sv_gr)$.sv_i <- sv$.sv_i
  S4Vectors::mcols(sv_gr)$sv_width <- as.integer(sv$width)
  
  # overlaps
  hits <- GenomicRanges::findOverlaps(sv_gr, cds_gr, ignore.strand = TRUE)
  
  # initialize outputs
  best_ro    <- rep(NA_real_, length(sv_gr))
  best_gene  <- rep(NA_character_, length(sv_gr))
  best_prod  <- rep(NA_character_, length(sv_gr))
  best_locus <- rep(NA_character_, length(sv_gr))
  best_cds_w <- rep(NA_integer_, length(sv_gr))
  
  if (length(hits) > 0) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    
    inter <- GenomicRanges::pintersect(
      GenomicRanges::ranges(sv_gr)[qh],
      GenomicRanges::ranges(cds_gr)[sh]
    )
    inter_w <- IRanges::width(inter)
    
    cds_w <- IRanges::width(GenomicRanges::ranges(cds_gr))[sh]
    sv_w  <- S4Vectors::mcols(sv_gr)$sv_width[qh]
    
    ro <- inter_w / pmax(1L, pmin(sv_w, cds_w))
    
    mg <- GenomicRanges::mcols(cds_gr)
    gene    <- if ("gene" %in% names(mg)) as.character(mg$gene[sh]) else NA_character_
    product <- if ("product" %in% names(mg)) as.character(mg$product[sh]) else NA_character_
    locus   <- if ("locus_tag" %in% names(mg)) as.character(mg$locus_tag[sh]) else NA_character_
    
    # keep max RO per SV (ties: first)
    ord <- order(qh, -ro)
    qh2 <- qh[ord]; ro2 <- ro[ord]
    gene2 <- gene[ord]; prod2 <- product[ord]; locus2 <- locus[ord]; cdsw2 <- cds_w[ord]
    
    first_idx <- !duplicated(qh2)
    q_keep <- qh2[first_idx]
    
    best_ro[q_keep]    <- ro2[first_idx]
    best_gene[q_keep]  <- gene2[first_idx]
    best_prod[q_keep]  <- prod2[first_idx]
    best_locus[q_keep] <- locus2[first_idx]
    best_cds_w[q_keep] <- as.integer(cdsw2[first_idx])
  }
  
  out <- sv %>%
    mutate(
      best_ro = best_ro,
      ro_bin  = svcm_ro_bin(best_ro),
      best_gene = best_gene,
      best_product = best_prod,
      best_locus_tag = best_locus,
      best_cds_width = best_cds_w
    )
  
  # keep only expected columns + new annotations (retain originals too)
  out
}
