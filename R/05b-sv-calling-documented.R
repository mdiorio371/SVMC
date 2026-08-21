# ==========================================================================
# SVMC: structural variation calling
# Replichore domains and the five SV classes; SV aggregation and grouping.
# ==========================================================================

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
) {
  if (is.character(delta_table)) delta_table <- read_delta(delta_table)
  
  dt <- delta_table %>% 
    dplyr::mutate(qs2 = pmin(qs, qe), qe2 = pmax(qs, qe))
  qlen <- dt$qlen[1]
  rlen <- dt$rlen[1]
  rid1 <- as.character(dt$rid[1])
  qid1 <- as.character(dt$qid[1])
  
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
      rid = rid1, 
      qid = qid1,
      start = as.integer(IRanges::start(ins_gr)),
      end   = as.integer(IRanges::end(ins_gr)),
      width = as.integer(IRanges::width(ins_gr)),
      ref_insertion_point = as.numeric(ref_pt),
      variant = "Insertion",
      variant_specific = "Insertion"
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
      rid = rid1, 
      qid = qid1,
      start = as.integer(IRanges::start(ref_gaps)),
      end   = as.integer(IRanges::end(ref_gaps)),
      width = as.integer(IRanges::width(ref_gaps)),
      ref_insertion_point = NA_real_,
      variant = "Deletion",
      variant_specific = "Deletion"
    )
  } else tibble::tibble()
  
  # -------- Combine, guard empty schema ----
  out <- dplyr::bind_rows(insertions_df, deletions_df)
  
  if (nrow(out) == 0L || ncol(out) == 0L) {
    out <- tibble::tibble(
      rid = character(), 
      qid = character(),
      start = integer(),  
      end = integer(),  
      width = integer(),
      ref_insertion_point = numeric(),
      variant = character(),
      variant_specific = character()
    )
    if (isTRUE(add_size_class)) {
      out <- tibble::add_column(out, size_class = character())
    }
    return(out)
  }
  
  if (!is.null(minlen)) out <- dplyr::filter(out, width >= minlen)
  
  if (add_size_class) {
    labs <- c(
      paste0("<", size_bands[1]),
      paste0(size_bands[1], "–", size_bands[2] - 1),
      paste0("≥", size_bands[2])
    )
    out <- out %>%
      dplyr::mutate(size_class = cut(
        width, 
        breaks = c(-Inf, size_bands, Inf), 
        labels = labs, 
        right = FALSE
      ))
  }
  
  dplyr::arrange(out, start)
}

delta_duplications <- function(
    delta_table,
    unfiltered_delta_table,
    minoverlap          = 200,
    max_tandem_gap      = 10000,
    min_prop            = 0.9,
    drop_within_primary = FALSE,   # conservative: drop trivial self/contained inside same primary block
    primary_cover_tol   = 0.95,
    include_query_side  = FALSE    # also detect two-query→one-ref (query-side duplications)
) {
  # ----------- IO / parsing -----------
  if (is.character(delta_table)) {
    delta_table <- read_delta(delta_table)
  }
  if (is.character(unfiltered_delta_table)) {
    unfiltered_delta_table <- read_delta(unfiltered_delta_table)
  }
  
  # Early shrink: normalize qs/qe, drop duplicate alignments, drop tiny
  ufd <- unfiltered_delta_table %>%
    dplyr::mutate(qs2 = pmin(qs, qe), qe2 = pmax(qs, qe)) %>%
    dplyr::distinct(qid, qs2, qe2, rid, rs, re, strand, .keep_all = TRUE) %>%
    dplyr::filter((qe2 - qs2 + 1L) >= minoverlap) %>%
    dplyr::ungroup()
  
  if (!nrow(ufd)) {
    return(
      dplyr::tibble(
        rid = delta_table$rid[1],
        qid = delta_table$qid[1],
        start = NA_integer_, end = NA_integer_, width = NA_integer_,
        orig_ref_start = NA_integer_, orig_ref_end = NA_integer_,
        dup_ref_start  = NA_integer_, dup_ref_end  = NA_integer_,
        dup_len = NA_integer_, dup_gap = NA_integer_,
        domains = NA_character_,
        variant = "Duplication", variant_specific = NA_character_,
        orientation = NA_character_,
        collapsed = 0L, duplication_found = FALSE,
        refmid = NA_integer_
      )
    )
  }
  
  stopifnot(all(c("qs", "qe", "rs", "re") %in% names(ufd)))
  ref_len <- unique(ufd$rlen)
  stopifnot(length(ref_len) == 1)
  
  # Only build primary blocks if we might use them
  primary_gr <- NULL
  if (drop_within_primary) {
    primary_gr <- GenomicRanges::GRanges(
      seqnames = delta_table$rid,
      ranges   = IRanges::IRanges(start = delta_table$rs, end = delta_table$re)
    )
    if (length(primary_gr) > 0) {
      S4Vectors::mcols(primary_gr)$block_id <- seq_along(primary_gr)
    }
  }
  
  # GRanges in query space (for reference-side detection)
  ufd_qry_gr <- GenomicRanges::GRanges(
    seqnames  = ufd$qid,
    ranges    = IRanges::IRanges(start = ufd$qs2, end = ufd$qe2),
    strand    = ufd$strand,
    rs        = ufd$rs,      re        = ufd$re,   # reference coords of each aln
    qry_start = ufd$qs2,     qry_end   = ufd$qe2,
    qry_name  = ufd$qid,     ref_name  = ufd$rid
  )
  
  # GRanges in reference space (for query-side detection)
  ufd_ref_gr <- GenomicRanges::GRanges(
    seqnames  = ufd$rid,
    ranges    = IRanges::IRanges(start = ufd$rs, end = ufd$re),
    strand    = ufd$strand,
    rs        = ufd$rs,      re        = ufd$re,
    qry_start = ufd$qs2,     qry_end   = ufd$qe2,  # query coords of each aln
    qry_name  = ufd$qid,     ref_name  = ufd$rid
  )
  
  # cached domain assigner
  have_assign_domains <- exists("assign_domains", mode = "function")
  safe_assign_domains <- function(pos, ref_len) {
    if (have_assign_domains) {
      assign_domains(pos, ref_len)
    } else {
      cut(pos / ref_len, breaks = c(-Inf, .25, .5, .75, Inf),
          labels = c("Q1","Q2","Q3","Q4"), right = TRUE)
    }
  }
  
  # ---------- pair finding: reference-side (two ref loci hit by same query locus)
  find_dups_refside <- function(gr) {
    h <- GenomicRanges::findOverlaps(gr, gr, minoverlap = minoverlap, ignore.strand = TRUE)
    if (length(h) == 0) return(dplyr::tibble())
    
    row1 <- S4Vectors::queryHits(h); row2 <- S4Vectors::subjectHits(h)
    keep <- row1 < row2
    row1 <- row1[keep]; row2 <- row2[keep]
    if (!length(row1)) return(dplyr::tibble())
    
    starts  <- IRanges::start(gr); ends <- IRanges::end(gr)
    meta    <- S4Vectors::mcols(gr)
    strands <- as.character(GenomicRanges::strand(gr))
    
    # overlap length in QUERY space
    dup_len <- pmax(0L, pmin(ends[row1], ends[row2]) - pmax(starts[row1], starts[row2]) + 1L)
    
    # project overlapped sub-interval (of row1) back to REFERENCE space
    denom <- pmax(1L, (ends[row1] - starts[row1]))
    frac  <- (pmax(starts[row1], starts[row2]) - starts[row1]) / denom
    orig_ref_dup_start <- round(meta$rs[row1] + frac * (meta$re[row1] - meta$rs[row1]))
    orig_ref_dup_end   <- orig_ref_dup_start + dup_len - 1L
    
    # gap between the two REF hits → tandem vs dispersed (relative to reference)
    orig_ref_start <- meta$rs[row1]; orig_ref_end <- meta$re[row1]
    dup_ref_start  <- meta$rs[row2]; dup_ref_end  <- meta$re[row2]
    dup_gap <- dplyr::case_when(
      dup_ref_start > orig_ref_end  ~ dup_ref_start - orig_ref_end,
      orig_ref_start > dup_ref_end  ~ orig_ref_start - dup_ref_end,
      TRUE ~ 0L
    )
    
    dplyr::tibble(
      rid = meta$ref_name[row1], qid = meta$qry_name[row1],
      orig_ref_dup_start, orig_ref_dup_end,
      orig_ref_start, orig_ref_end,
      dup_ref_start,  dup_ref_end,
      dup_len, dup_gap,
      duplication_type = dplyr::if_else(dup_gap <= max_tandem_gap, "tandem", "dispersed"),
      orientation      = dplyr::if_else(strands[row1] == strands[row2], "same", "inverted")
    )
  }
  
  # ---------- pair finding: query-side (two query loci hit the same ref locus)
  find_dups_queryside <- function(gr) {
    h <- GenomicRanges::findOverlaps(gr, gr, minoverlap = minoverlap, ignore.strand = TRUE)
    if (length(h) == 0) return(dplyr::tibble())
    
    row1 <- S4Vectors::queryHits(h); row2 <- S4Vectors::subjectHits(h)
    keep <- row1 < row2
    row1 <- row1[keep]; row2 <- row2[keep]
    if (!length(row1)) return(dplyr::tibble())
    
    starts  <- IRanges::start(gr); ends <- IRanges::end(gr)
    meta    <- S4Vectors::mcols(gr)
    strands <- as.character(GenomicRanges::strand(gr))
    
    # overlap length in REFERENCE space (the shared ref locus)
    orig_ref_dup_start <- pmax(starts[row1], starts[row2])
    orig_ref_dup_end   <- pmin(ends[row1],   ends[row2])
    dup_len <- pmax(0L, orig_ref_dup_end - orig_ref_dup_start + 1L)
    
    # For query-side, tandem/dispersed is judged on QUERY separation
    q1s <- meta$qry_start[row1]; q1e <- meta$qry_end[row1]
    q2s <- meta$qry_start[row2]; q2e <- meta$qry_end[row2]
    dup_gap_q <- dplyr::case_when(
      q2s > q1e ~ q2s - q1e,
      q1s > q2e ~ q1s - q2e,
      TRUE ~ 0L
    )
    
    dplyr::tibble(
      rid = meta$ref_name[row1], qid = meta$qry_name[row1],
      orig_ref_dup_start, orig_ref_dup_end,
      # keep these for symmetry/context (ref spans of both alignments)
      orig_ref_start = meta$rs[row1], orig_ref_end = meta$re[row1],
      dup_ref_start  = meta$rs[row2], dup_ref_end  = meta$re[row2],
      dup_len, dup_gap = dup_gap_q,
      duplication_type = dplyr::if_else(dup_gap_q <= max_tandem_gap, "tandem", "dispersed"),
      orientation      = dplyr::if_else(strands[row1] == strands[row2], "same", "inverted")
    )
  }
  
  # Optional conservative primary-block filter
  primary_filter <- function(tbl) {
    if (!drop_within_primary || is.null(primary_gr) || nrow(tbl) == 0) return(tbl)
    
    cover_stats <- function(starts, ends, seqnames_vec) {
      gr <- GenomicRanges::GRanges(
        seqnames = seqnames_vec,
        ranges   = IRanges::IRanges(starts, ends)
      )
      hits <- GenomicRanges::findOverlaps(gr, primary_gr)
      if (length(hits) == 0) {
        return(list(block = rep(NA_integer_, length(gr)), prop = rep(0, length(gr))))
      }
      ovw <- IRanges::width(
        GenomicRanges::pintersect(
          gr[S4Vectors::queryHits(hits)],
          primary_gr[S4Vectors::subjectHits(hits)]
        )
      )
      idx <- tapply(
        seq_along(ovw),
        S4Vectors::queryHits(hits),
        function(i) i[which.max(ovw[i])]
      )
      block <- rep(NA_integer_, length(gr))
      prop  <- rep(0, length(gr))
      qi <- as.integer(names(idx))
      si <- S4Vectors::subjectHits(hits)[unlist(idx, use.names = FALSE)]
      block[qi] <- S4Vectors::mcols(primary_gr)$block_id[si]
      prop[qi]  <- ovw[unlist(idx, use.names = FALSE)] / IRanges::width(gr[qi])
      list(block = block, prop = prop)
    }
    
    cs_orig <- cover_stats(tbl$orig_ref_dup_start, tbl$orig_ref_dup_end, seqnames_vec = tbl$rid)
    cs_dup  <- cover_stats(tbl$dup_ref_start,      tbl$dup_ref_end,      seqnames_vec = tbl$rid)
    
    inter_w <- pmax(0L, pmin(tbl$orig_ref_dup_end, tbl$dup_ref_end) -
                      pmax(tbl$orig_ref_dup_start, tbl$dup_ref_start) + 1L)
    small_len <- pmin(tbl$orig_ref_dup_end - tbl$orig_ref_dup_start + 1L,
                      tbl$dup_ref_end      - tbl$dup_ref_start      + 1L)
    overlap_prop_small <- ifelse(small_len > 0, inter_w / small_len, 0)
    
    keep <- !(
      !is.na(cs_orig$block) & !is.na(cs_dup$block) &
        (cs_orig$block == cs_dup$block) &
        (cs_orig$prop  >= primary_cover_tol) &
        (cs_dup$prop   >= primary_cover_tol) &
        (overlap_prop_small >= primary_cover_tol)
    )
    tbl[keep, , drop = FALSE]
  }
  
  # Consolidate → score → thin to high-confidence, then cluster
  consolidate_and_cluster <- function(pair_tbl) {
    if (nrow(pair_tbl) == 0) return(pair_tbl)
    
    all_dups <- pair_tbl %>%
      dplyr::mutate(
        orig_mid      = round((orig_ref_start + orig_ref_end) / 2),
        dup_mid       = round((dup_ref_start  + dup_ref_end)  / 2),
        mid_domain    = safe_assign_domains(orig_mid, ref_len),
        refmid_domain = safe_assign_domains(dup_mid,  ref_len),
        domains       = paste(mid_domain, refmid_domain, sep = "_"),
        orig_ref_len  = orig_ref_end - orig_ref_start + 1L,
        dup_ref_len   = dup_ref_end  - dup_ref_start  + 1L,
        prop_orig     = dup_len / pmax(1L, orig_ref_len),
        prop_dup      = dup_len / pmax(1L, dup_ref_len),
        prop_min      = dup_len / pmax(1L, pmin(orig_ref_len, dup_ref_len))
      ) %>%
      dplyr::filter(prop_min > min_prop) %>%
      dplyr::select(
        rid, qid,
        orig_ref_dup_start, orig_ref_dup_end,
        orig_ref_start, orig_ref_end,
        dup_ref_start, dup_ref_end,
        dup_len, dup_gap,
        domains, duplication_type, orientation
      ) %>%
      dplyr::distinct(
        rid, qid,
        orig_ref_dup_start, orig_ref_dup_end,
        dup_ref_start, dup_ref_end,
        duplication_type, orientation, domains,
        .keep_all = TRUE
      )
    
    if (nrow(all_dups) == 0) return(all_dups)
    
    # FAST EDGES + CLUSTER (no igraph)
    N <- nrow(all_dups)
    proj_st <- all_dups$orig_ref_dup_start; proj_en <- all_dups$orig_ref_dup_end
    targ_st <- all_dups$dup_ref_start;     targ_en <- all_dups$dup_ref_end
    
    # candidate pairs in projected/source space (reference coords)
    h1 <- IRanges::findOverlaps(
      IRanges::IRanges(proj_st, proj_en),
      IRanges::IRanges(proj_st, proj_en),
      minoverlap = 1L
    )
    q1 <- S4Vectors::queryHits(h1); s1 <- S4Vectors::subjectHits(h1)
    keep1 <- q1 < s1
    q1 <- q1[keep1]; s1 <- s1[keep1]
    
    w1    <- pmax(0L, pmin(proj_en[q1], proj_en[s1]) - pmax(proj_st[q1], proj_st[s1]) + 1L)
    minw1 <- pmin(proj_en[q1]-proj_st[q1]+1L, proj_en[s1]-proj_st[s1]+1L)
    ok1   <- (w1 >= (min_prop * minw1))
    
    # candidate pairs in target space (reference coords of partner copy)
    h2 <- IRanges::findOverlaps(
      IRanges::IRanges(targ_st, targ_en),
      IRanges::IRanges(targ_st, targ_en),
      minoverlap = 1L
    )
    q2 <- S4Vectors::queryHits(h2); s2 <- S4Vectors::subjectHits(h2)
    keep2 <- q2 < s2
    q2 <- q2[keep2]; s2 <- s2[keep2]
    
    w2    <- pmax(0L, pmin(targ_en[q2], targ_en[s2]) - pmax(targ_st[q2], targ_st[s2]) + 1L)
    minw2 <- pmin(targ_en[q2]-targ_st[q2]+1L, targ_en[s2]-targ_st[s2]+1L)
    ok2   <- (w2 >= (min_prop * minw2))
    
    keyN <- N + 1L
    key1 <- q1*keyN + s1
    key2 <- q2*keyN + s2
    if (requireNamespace("fastmatch", quietly = TRUE)) {
      in_both <- ok1 & (fastmatch::fmatch(key1, key2, nomatch = 0L) > 0L)
    } else {
      in_both <- ok1 & !is.na(match(key1, key2))
    }
    from <- q1[in_both]; to <- s1[in_both]
    
    if (!length(from)) {
      out0 <- all_dups %>%
        dplyr::mutate(
          start = orig_ref_dup_start,
          end   = orig_ref_dup_end,
          width = end - start + 1L,
          variant = "Duplication",
          variant_specific = duplication_type,
          collapsed = 1L,
          duplication_found = TRUE,
          refmid = floor((start + end)/2)
        ) %>%
        dplyr::select(
          rid, qid, start, end, width,
          orig_ref_start, orig_ref_end,
          dup_ref_start, dup_ref_end,
          dup_len, dup_gap,
          domains, variant, variant_specific, orientation,
          collapsed, duplication_found, refmid
        ) %>%
        dplyr::distinct(rid, qid, start, end, width, variant, variant_specific, .keep_all = TRUE)
      return(out0)
    }
    
    # union–find components
    parent <- seq_len(N); rank <- integer(N)
    uf_find <- function(x){
      while (parent[x] != x) {
        parent[x] <<- parent[parent[x]]
        x <- parent[x]
      }
      x
    }
    uf_union <- function(a,b){
      ra <- uf_find(a); rb <- uf_find(b)
      if (ra == rb) return(invisible(NULL))
      if (rank[ra] < rank[rb]) parent[ra] <<- rb
      else if (rank[ra] > rank[rb]) parent[rb] <<- ra
      else { parent[rb] <<- ra; rank[ra] <<- rank[ra] + 1L }
      invisible(NULL)
    }
    for (i in seq_along(from)) uf_union(from[i], to[i])
    for (i in seq_len(N)) parent[i] <- uf_find(i)
    cluster <- match(parent, unique(parent))
    
    all_dups %>%
      dplyr::mutate(cluster = cluster) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        rid = dplyr::first(rid),
        qid = dplyr::first(qid),
        start = min(orig_ref_dup_start),
        end   = max(orig_ref_dup_end),
        orig_ref_start = min(orig_ref_start),
        orig_ref_end   = max(orig_ref_end),
        dup_ref_start  = min(dup_ref_start),
        dup_ref_end    = max(dup_ref_end),
        width   = end - start + 1L,
        dup_len = width,
        dup_gap = min(dup_gap),
        domains = paste(sort(unique(domains)), collapse = ";"),
        variant = "Duplication",
        variant_specific = if (dplyr::n_distinct(duplication_type) == 1L)
          dplyr::first(duplication_type) else "mixed",
        orientation = if (dplyr::n_distinct(orientation) == 1L)
          dplyr::first(orientation) else "mixed",
        collapsed = dplyr::n(),
        duplication_found = TRUE,
        .groups = "drop"
      ) %>%
      dplyr::mutate(refmid = floor((start + end)/2)) %>%
      dplyr::distinct(rid, qid, start, end, width, variant, variant_specific, .keep_all = TRUE)
  }
  
  # ---- Run both passes ----
  ref_pairs  <- find_dups_refside(ufd_qry_gr)   # duplication on reference
  ref_pairs  <- primary_filter(ref_pairs)
  ref_calls  <- consolidate_and_cluster(ref_pairs)
  
  if (include_query_side) {
    qry_pairs  <- find_dups_queryside(ufd_ref_gr)  # duplication on query
    qry_pairs  <- primary_filter(qry_pairs)
    qry_calls  <- consolidate_and_cluster(qry_pairs)
  } else {
    qry_calls <- dplyr::tibble()
  }
  
  # If neither side produced calls, return the empty scaffold
  if ((nrow(ref_calls) == 0) && (nrow(qry_calls) == 0)) {
    return(
      dplyr::tibble(
        rid = delta_table$rid[1],
        qid = delta_table$qid[1],
        start = NA_integer_, end = NA_integer_, width = NA_integer_,
        orig_ref_start = NA_integer_, orig_ref_end = NA_integer_,
        dup_ref_start  = NA_integer_, dup_ref_end  = NA_integer_,
        dup_len = NA_integer_, dup_gap = NA_integer_,
        domains = NA_character_,
        variant = "Duplication", variant_specific = NA_character_,
        orientation = NA_character_,
        collapsed = 0L, duplication_found = FALSE,
        refmid = NA_integer_
      )
    )
  }
  
  # Combine both sides; deduplicate on the canonical output subset
  out <- dplyr::bind_rows(ref_calls, qry_calls) %>%
    dplyr::distinct(rid, qid, start, end, width, variant, variant_specific, .keep_all = TRUE)
  
  out
}

delta_translocations <- function(
    delta_table,
    X_dist_threshold = 0.95
) {
  
  if (is.character(delta_table)) {
    delta_table <- read_delta(delta_table)
  }
  
  d2 <- delta_table %>%
    rowwise() %>%
    mutate(
      qs2 = min(c(qs, qe)),
      qe2 = max(c(qs, qe))
    ) %>% 
    ungroup()
  
  tlocs <- delta_table %>%
    filter(
      !is.na(X_dist),
      X_dist > quantile(X_dist, X_dist_threshold)
    ) %>%
    mutate(
      refmid = round((rs + re) / 2),
      qrymid = round((qs + qe) / 2),
      ref_domain = assign_domains(refmid, rlen),
      qry_domain = assign_domains(qrymid, qlen),
      replichores = paste(ref_domain, qry_domain, sep = "-")
    ) %>%
    rowwise() %>%
    mutate(
      start = rs,
      end = re,
      width = abs(start - end + 1),
      variant = "Translocation",
      variant_specific = replichores
    ) %>%
    ungroup() %>%
    dplyr::select(c(
      rid, qid,
      start, end, width, 
      refmid, qrymid,
      strand, 
      variant, variant_specific,
      X_dist
    ))
  
  return(tlocs)
}

delta_structural_rearrangements <- function(
    delta_table,
    symmetric_window   = 0.10,
    noise_bp           = 10000,
    xdist_tol          = 0.25,
    min_run_bp         = 1000,
    min_inv_bp         = 10000,
    max_denoise_passes = 3L,
    emit_flanks        = FALSE,
    emit_megacompound  = TRUE
){
  # ---------- Utilities ----------
  empty_out <- function() {
    tibble::tibble(
      rid = character(), qid = character(),
      start = numeric(), end = numeric(), width = numeric(),
      strand = character(), variant = character(),
      variant_specific = character(), X_dist = numeric()
    )
  }
  safe_arrange <- function(df) if (nrow(df)) dplyr::arrange(df, .data$start) else df
  
  # Fuse - (+tiny) - and + (-tiny) +
  denoise_once <- function(runs, noise_bp, xdist_tol, min_run_bp) {
    n <- nrow(runs); if (n < 3) return(list(runs = runs, changed = FALSE))
    s   <- runs$strand
    len <- runs$re - runs$rs + 1L
    mid <- 2:(n - 1L); L <- mid - 1L; R <- mid + 1L
    
    caseA <- s[L] == "-" & s[mid] == "+" & s[R] == "-" &
      len[mid] <= noise_bp &
      pmax(runs$X_dist[L], runs$X_dist[R], na.rm = TRUE) > 0 &
      abs(runs$X_dist[L] - runs$X_dist[R]) <=
      xdist_tol * pmax(runs$X_dist[L], runs$X_dist[R], na.rm = TRUE)
    
    caseB <- s[L] == "+" & s[mid] == "-" & s[R] == "+" &
      len[mid] <= noise_bp
    
    idx <- which(caseA | caseB)
    if (!length(idx)) return(list(runs = runs, changed = FALSE))
    
    new <- runs[L[idx], ]
    new$rs     <- pmin(runs$rs[L[idx]],  runs$rs[R[idx]])
    new$re     <- pmax(runs$re[L[idx]],  runs$re[R[idx]])
    new$qs2    <- pmin(runs$qs2[L[idx]], runs$qs2[R[idx]])
    new$qe2    <- pmax(runs$qe2[L[idx]], runs$qe2[R[idx]])
    new$X_dist <- rowMeans(cbind(runs$X_dist[L[idx]], runs$X_dist[R[idx]]), na.rm = TRUE)
    new$strand <- runs$strand[L[idx]]
    
    keep <- rep(TRUE, n)
    keep[c(rbind(L[idx], mid[idx], R[idx]))] <- FALSE
    
    runs2 <- dplyr::bind_rows(runs[keep, ], new) %>%
      dplyr::arrange(.data$rs)
    
    r2id <- c(1L, 1L + cumsum(runs2$strand[-1L] != runs2$strand[-nrow(runs2)]))
    runs2 <- runs2 %>%
      dplyr::mutate(r2id = r2id) %>%
      dplyr::group_by(.data$r2id, .data$strand) %>%
      dplyr::summarise(
        rs     = min(.data$rs),  re = max(.data$re),
        qs2    = min(.data$qs2), qe2 = max(.data$qe2),
        X_dist = mean(.data$X_dist, na.rm = TRUE),
        rid    = dplyr::first(.data$rid),
        qid    = dplyr::first(.data$qid),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$rs) %>%
      dplyr::filter((.data$re - .data$rs + 1L) >= min_run_bp) %>%
      dplyr::select(-.data$r2id)
    
    list(runs = runs2, changed = TRUE)
  }
  
  # Find maximal alternating families
  find_alt_families <- function(runs, loc_thresh) {
    n <- nrow(runs); if (n < 3) return(list())
    s <- runs$strand
    used <- rep(FALSE, n)
    fams <- list()
    i <- 2L
    while (i <= n - 2L) {
      if (used[i]) { i <- i + 1L; next }
      j <- i
      while (j < n && runs$strand[j + 1L] != runs$strand[j]) j <- j + 1L
      if ((j - i + 1L) >= 3L) {
        cand <- which(s[i:j] == s[i]) + i - 1L
        cand <- cand[cand >= (i + 2L)]
        if (length(cand)) {
          cand_ok <- cand[runs$X_dist[i] <= loc_thresh & runs$X_dist[cand] <= loc_thresh]
          if (length(cand_ok)) {
            R <- max(cand_ok)
            fam_idx <- i:R
            used[fam_idx] <- TRUE
            fams[[length(fams) + 1L]] <- fam_idx
            i <- R + 1L
            next
          }
        }
      }
      i <- j + 1L
    }
    fams
  }
  
  # Island compounds inside a family (emits B-like spans)
  emit_island_compounds <- function(runs, fam_idx, loc_thresh, min_inv_bp) {
    if (length(fam_idx) < 3) return(NULL)
    L <- fam_idx[1L]; R <- fam_idx[length(fam_idx)]
    outer_strand <- runs$strand[L]
    opp_strand   <- if (outer_strand == "+") "-" else "+"
    
    inner_idx <- fam_idx[2:(length(fam_idx)-1L)]
    if (!length(inner_idx)) return(NULL)
    
    flanks <- inner_idx[runs$strand[inner_idx] == outer_strand]
    if (length(flanks) < 2L) return(NULL)
    
    out <- list(); a_pos <- 1L
    while (a_pos < length(flanks)) {
      left <- flanks[a_pos]
      cand <- flanks[(a_pos+1L):length(flanks)]
      if (!length(cand)) break
      ok <- sapply(cand, function(rr){
        between <- inner_idx[inner_idx > left & inner_idx < rr]
        if (!length(between)) return(FALSE)
        opp_between <- between[runs$strand[between] == opp_strand]
        if (!length(opp_between)) return(FALSE)
        j <- min(opp_between); k <- max(opp_between)
        (runs$X_dist[j] <= loc_thresh) && (runs$X_dist[k] <= loc_thresh)
      })
      cand_ok <- cand[ok]
      if (length(cand_ok)) {
        right <- max(cand_ok)
        between <- inner_idx[inner_idx > left & inner_idx < right]
        opp_between <- between[runs$strand[between] == opp_strand]
        j <- min(opp_between); k <- max(opp_between)
        
        span_start <- runs$rs[j]; span_end <- runs$re[k]
        span_len   <- as.numeric(span_end - span_start + 1)
        if (span_len >= min_inv_bp) {
          out[[length(out)+1L]] <- tibble::tibble(
            rid = runs$rid[j], qid = runs$qid[j],
            start = span_start, end = span_end, width = span_len,
            strand = outer_strand,                            # keep consistent with “inversion = -”
            variant = "Structural_rearrangement",
            variant_specific = "inner_compound",
            X_dist = mean(c(runs$X_dist[left], runs$X_dist[right]), na.rm = TRUE)
          )
        }
        a_pos <- which(flanks == right)
      }
      a_pos <- a_pos + 1L
    }
    if (length(out)) dplyr::bind_rows(out) else NULL
  }
  
  # Nested symmetric calls between interior same-strand flanks (emits D)
  # Returns list(calls=tibble, used_flanks=integer idx of flank runs)
  emit_nested_symmetric <- function(runs, fam_idx, loc_thresh, min_inv_bp) {
    if (length(fam_idx) < 3) return(list(calls=NULL, used_flanks=integer()))
    L <- fam_idx[1L]; R <- fam_idx[length(fam_idx)]
    outer_strand <- runs$strand[L]
    inner_idx <- fam_idx[2:(length(fam_idx)-1L)]
    if (!length(inner_idx)) return(list(calls=NULL, used_flanks=integer()))
    
    flanks <- inner_idx[runs$strand[inner_idx] == outer_strand]
    if (length(flanks) < 2L) return(list(calls=NULL, used_flanks=integer()))
    
    out <- list(); used <- integer()
    a_pos <- 1L
    while (a_pos < length(flanks)) {
      left <- flanks[a_pos]
      cand <- flanks[(a_pos+1L):length(flanks)]
      if (!length(cand)) break
      ok <- sapply(cand, function(right){
        between <- inner_idx[inner_idx > left & inner_idx < right]
        if (!length(between)) return(FALSE)
        have_opp <- any(runs$strand[between] != outer_strand)
        x_ok <- (runs$X_dist[left] <= loc_thresh) && (runs$X_dist[right] <= loc_thresh)
        have_opp && x_ok
      })
      cand_ok <- cand[ok]
      if (length(cand_ok)) {
        right <- max(cand_ok)
        start <- min(runs$rs[left], runs$rs[right])
        end   <- max(runs$re[left], runs$re[right])
        span  <- as.numeric(end - start + 1)
        if (span >= min_inv_bp) {
          out[[length(out)+1L]] <- tibble::tibble(
            rid = runs$rid[left], qid = runs$qid[left],
            start = start, end = end, width = span,
            strand = outer_strand,
            variant = "Structural_rearrangement",
            variant_specific = "inner_symmetric",
            X_dist = mean(c(runs$X_dist[left], runs$X_dist[right]), na.rm = TRUE)
          )
          used <- c(used, left, right)       # mark flank runs so we can suppress them from inner_nested
        }
        a_pos <- which(flanks == right)
      }
      a_pos <- a_pos + 1L
    }
    list(calls = (if (length(out)) dplyr::bind_rows(out) else NULL),
         used_flanks = unique(used))
  }
  
  # ---------- Load & normalize ----------
  delta_table <- if (is.character(delta_table)) read_delta(delta_table) else delta_table
  if (is.null(delta_table) || !nrow(delta_table)) return(empty_out())
  
  needed <- c("strand","rs","re","qs","qe","rid","qid","rlen","X_dist")
  missing <- setdiff(needed, colnames(delta_table))
  if (length(missing)) {
    warning("delta_structural_rearrangements: missing columns: ", paste(missing, collapse = ", "))
    return(empty_out())
  }
  
  fd <- delta_table %>%
    filter_delta() %>%
    dplyr::mutate(
      rs = as.numeric(.data$rs), re = as.numeric(.data$re),
      qs = as.numeric(.data$qs), qe = as.numeric(.data$qe),
      X_dist = as.numeric(.data$X_dist),
      strand = as.character(.data$strand),
      rid = as.character(.data$rid), qid = as.character(.data$qid),
      qs2 = pmin(.data$qs, .data$qe),
      qe2 = pmax(.data$qs, .data$qe)
    ) %>%
    dplyr::arrange(.data$rs)
  
  if (!nrow(fd)) return(empty_out())
  if (length(unique(fd$rlen)) != 1L) warning("Multiple rlen values detected; using the first.")
  loc_thresh <- round(symmetric_window * fd$rlen[1])
  
  # Collapse same-strand runs
  run_id <- c(1L, 1L + cumsum(fd$strand[-1L] != fd$strand[-nrow(fd)]))
  runs <- fd %>%
    dplyr::mutate(run_id = run_id) %>%
    dplyr::group_by(.data$run_id, .data$strand) %>%
    dplyr::summarise(
      rs     = min(.data$rs),  re = max(.data$re),
      qs2    = min(.data$qs2), qe2 = max(.data$qe2),
      X_dist = mean(.data$X_dist, na.rm = TRUE),
      rid    = dplyr::first(.data$rid),
      qid    = dplyr::first(.data$qid),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$rs) %>%
    dplyr::filter((.data$re - .data$rs + 1L) >= min_run_bp)
  
  if (!nrow(runs)) return(empty_out())
  
  # Denoise fractured flanks
  passes <- 0L; changed <- TRUE
  while (changed && nrow(runs) >= 3L && passes < max_denoise_passes) {
    dd <- denoise_once(runs, noise_bp, xdist_tol, min_run_bp)
    runs <- dd$runs; changed <- dd$changed; passes <- passes + 1L
  }
  
  # ---------- Build families & call events ----------
  fams <- find_alt_families(runs, loc_thresh)
  out <- list()
  
  if (length(fams)) {
    fam_calls <- lapply(fams, function(idx) {
      L <- idx[1L]; R <- idx[length(idx)]
      outer_strand <- runs$strand[L]
      opp_strand   <- if (outer_strand == "+") "-" else "+"
      
      # OUTER
      outer <- tibble::tibble(
        rid = runs$rid[L], qid = runs$qid[L],
        start = min(runs$rs[L], runs$rs[R]),
        end   = max(runs$re[L], runs$re[R]),
        width = as.numeric(max(runs$re[c(L,R)]) - min(runs$rs[c(L,R)]) + 1),
        strand = outer_strand,
        variant = "Structural_rearrangement",
        variant_specific = "outer",
        X_dist = mean(c(runs$X_dist[L], runs$X_dist[R]), na.rm = TRUE)
      )
      
      inner_idx <- idx[2:(length(idx)-1L)]
      
      # Mega interior (union of all interior opposite-strand runs)
      inner_megacompound <- NULL
      if (isTRUE(emit_megacompound)) {
        opp_runs <- inner_idx[runs$strand[inner_idx] == opp_strand]
        if (length(opp_runs)) {
          mc_start <- min(runs$rs[opp_runs]); mc_end <- max(runs$re[opp_runs])
          mc_w     <- as.numeric(mc_end - mc_start + 1)
          if (mc_w >= min_inv_bp) {
            inner_megacompound <- tibble::tibble(
              rid = runs$rid[L], qid = runs$qid[L],
              start = mc_start, end = mc_end, width = mc_w,
              strand = outer_strand,
              variant = "Structural_rearrangement",
              variant_specific = "inner_compound_all",
              X_dist = mean(c(runs$X_dist[L], runs$X_dist[R]), na.rm = TRUE)
            )
          }
        }
      }
      
      # Island compounds (e.g., B)
      inner_islands <- emit_island_compounds(runs, idx, loc_thresh, min_inv_bp)
      
      # Nested symmetric (e.g., D) + which flanks were used
      ns <- emit_nested_symmetric(runs, idx, loc_thresh, min_inv_bp)
      nested_sym <- ns$calls
      used_flanks <- ns$used_flanks
      
      # Inner nested (minus runs) excluding those used as nested-symmetric flanks
      inner_minus <- inner_idx[runs$strand[inner_idx] == "-"]
      if (length(used_flanks)) inner_minus <- setdiff(inner_minus, used_flanks)
      inner_nested <- NULL
      if (length(inner_minus)) {
        inner_nested <- tibble::tibble(
          rid   = runs$rid[inner_minus],
          qid   = runs$qid[inner_minus],
          start = runs$rs[inner_minus],
          end   = runs$re[inner_minus],
          width = as.numeric(runs$re[inner_minus] - runs$rs[inner_minus] + 1),
          strand = "-",
          variant = "Structural_rearrangement",
          variant_specific = "inner_nested",
          X_dist = runs$X_dist[inner_minus]
        )
      }
      
      flanks <- NULL
      if (isTRUE(emit_flanks)) {
        flanks <- tibble::tibble(
          rid   = c(runs$rid[L], runs$rid[R]),
          qid   = c(runs$qid[L], runs$qid[R]),
          start = c(runs$rs[L],  runs$rs[R]),
          end   = c(runs$re[L],  runs$re[R]),
          width = as.numeric(c(runs$re[L]-runs$rs[L]+1, runs$re[R]-runs$rs[R]+1)),
          strand = c(runs$strand[L], runs$strand[R]),
          variant = "Structural_rearrangement",
          variant_specific = "flank",
          X_dist = c(runs$X_dist[L], runs$X_dist[R])
        )
      }
      
      dplyr::bind_rows(
        outer,
        if (!is.null(inner_megacompound)) inner_megacompound,
        if (!is.null(inner_islands))      inner_islands,
        if (!is.null(nested_sym))         nested_sym,
        if (!is.null(inner_nested))       inner_nested,
        if (!is.null(flanks))             flanks
      )
    })
    out$symmetric <- dplyr::bind_rows(fam_calls)
  }
  
  # Singles: minus runs not in any family
  in_family <- rep(FALSE, nrow(runs))
  if (length(fams)) in_family[unlist(fams)] <- TRUE
  singles <- runs %>%
    dplyr::mutate(run_idx = dplyr::row_number()) %>%
    dplyr::filter(.data$strand == "-", !in_family) %>%
    dplyr::transmute(
      rid, qid,
      start = .data$rs, end = .data$re,
      width = as.numeric(.data$re - .data$rs + 1),
      strand,
      variant = "Structural_rearrangement",
      variant_specific = dplyr::if_else(.data$X_dist > loc_thresh, "localized", "symmetric_single"),
      X_dist
    )
  if (nrow(singles)) out$single <- singles
  
  # ---------- Combine & filter ----------
  combined <- if (length(out)) suppressWarnings(dplyr::bind_rows(out)) else NULL
  if (is.null(combined) || !nrow(combined)) return(empty_out())
  
  res <- combined %>%
    dplyr::mutate(width = as.numeric(.data$width)) %>%
    dplyr::filter(.data$width >= min_inv_bp) %>%
    safe_arrange()
  
  if (!nrow(res)) return(empty_out())
  res
}

delta_substructural_inversions <- function(
    delta_table,
    delta_SR,
    symmetric_window = 0.10,
    include_outside = TRUE,
    fully_contained = TRUE,
    min_inv_bp = 1000
) {
  # Empty output template
  empty_out <- function() {
    tibble::tibble(
      rid = character(),
      qid = character(),
      start = double(),
      end = double(),
      width = double(),
      meanlen = double(),
      variant = character(),
      variant_specific = character(),
      context = character()
    )
  }
  
  # ---- Load & normalize ----
  dt <- if (is.character(delta_table)) read_delta(delta_table) else delta_table
  if (is.null(dt) || !nrow(dt)) return(empty_out())
  
  # Normalize coordinates
  dt <- dt %>%
    dplyr::mutate(
      qs2 = pmin(qs, qe),
      qe2 = pmax(qs, qe),
      rs = as.numeric(rs),
      re = as.numeric(re),
      qs = as.numeric(qs),
      qe = as.numeric(qe),
      width = re - rs + 1
    ) %>%
    dplyr::filter(re > rs)
  
  # Ensure meanlen exists
  if (!"meanlen" %in% names(dt)) {
    dt <- dt %>% dplyr::mutate(meanlen = (width + abs(qe2 - qs2)) / 2)
  }
  
  # Ensure X_dist exists
  if (!"X_dist" %in% names(dt)) {
    if ("fwd_dist" %in% names(dt) || "rev_dist" %in% names(dt)) {
      dt <- dt %>% 
        dplyr::mutate(X_dist = dplyr::coalesce(fwd_dist, rev_dist, 0))
    } else {
      warning("X_dist not found - setting to 0")
      dt <- dt %>% dplyr::mutate(X_dist = 0)
    }
  }
  
  # Infer rlen if missing
  if (!"rlen" %in% names(dt) || is.na(dt$rlen[1])) {
    dt <- dt %>%
      dplyr::group_by(rid) %>%
      dplyr::mutate(rlen = max(re, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }
  
  localized_threshold <- round(symmetric_window * dt$rlen[1])
  
  # ---- Early exit if no SR ----
  if (is.null(delta_SR) || !nrow(delta_SR)) {
    invs <- dt %>%
      dplyr::filter(strand == "-", width >= min_inv_bp) %>%
      dplyr::mutate(
        variant = dplyr::if_else(
          X_dist < localized_threshold, 
          "Substructural_inversion",
          "Substructural_inversion"
        ),
        variant_specific = dplyr::if_else(
          X_dist < localized_threshold,
          "symmetric",
          "localized"
        ),
        context = "outside_all_SR"
      ) %>%
      dplyr::select(rid, qid, start = rs, end = re, width, meanlen, variant, variant_specific, context)
    return(invs)
  }
  
  # ---- Prepare SR windows ----
  sr <- delta_SR %>%
    dplyr::filter(variant_specific %in% c("symmetric_outer", "localized", "symmetric_single")) %>%
    dplyr::transmute(
      sr_id = dplyr::row_number(),
      sr_rs = start,
      sr_re = end,
      sr_strand = strand,
      sr_variant = variant_specific
    ) %>%
    dplyr::arrange(sr_rs)
  
  if (!nrow(sr)) {
    invs <- dt %>%
      dplyr::filter(strand == "-", width >= min_inv_bp) %>%
      dplyr::mutate(
        variant = "Substructural_inversion",
        variant_specific = dplyr::if_else(
          X_dist < localized_threshold,
          "symmetric",
          "localized"
        ),
        context = "outside_all_SR"
      ) %>%
      dplyr::select(rid, qid, start = rs, end = re, width, meanlen, variant, variant_specific, context)
    return(invs)
  }
  
  # ---- 1) Inside SR (substructural events) ----
  dt_sr <- dt %>%
    dplyr::mutate(seg_id = dplyr::row_number()) %>%
    tidyr::crossing(sr) %>%
    dplyr::filter(
      if (fully_contained) {
        rs >= sr_rs & re <= sr_re
      } else {
        pmax(rs, sr_rs) <= pmin(re, sr_re)
      }
    ) %>%
    dplyr::filter(strand != sr_strand) %>%
    dplyr::group_by(seg_id) %>%
    dplyr::slice_min(order_by = (sr_re - sr_rs), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  inside_calls <- dt_sr %>%
    dplyr::distinct(seg_id, .keep_all = TRUE) %>%
    dplyr::filter(width >= min_inv_bp) %>%
    dplyr::mutate(
      variant = "Substructural_inversion",
      variant_specific = dplyr::case_when(
        strand == "+" & sr_strand == "-" & X_dist < localized_threshold ~ "restored_orientation",
        strand == "+" & sr_strand == "-" ~ "nested_inversion_localized",
        strand == "-" & sr_strand == "+" & X_dist < localized_threshold ~ "symmetric",
        strand == "-" & sr_strand == "+" ~ "localized",
        TRUE ~ "complex_rearrangement"
      ),
      context = "inside_SR"
    ) %>%
    dplyr::select(rid, qid, start = rs, end = re, width, meanlen, variant, variant_specific, context)
  
  # ---- 2) Outside all SR (optional) ----
  outside_calls <- tibble::tibble()
  if (include_outside) {
    ovlp <- dt %>%
      dplyr::filter(strand == "-") %>%
      dplyr::mutate(seg_id = dplyr::row_number()) %>%
      tidyr::crossing(sr) %>%
      dplyr::mutate(overlap = pmax(rs, sr_rs) <= pmin(re, sr_re)) %>%
      dplyr::group_by(seg_id) %>%
      dplyr::summarise(any_overlap = any(overlap), .groups = "drop")
    
    outside_ids <- ovlp %>% dplyr::filter(!any_overlap) %>% dplyr::pull(seg_id)
    
    outside_calls <- dt %>%
      dplyr::filter(strand == "-") %>%
      dplyr::mutate(seg_id = dplyr::row_number()) %>%
      dplyr::filter(seg_id %in% outside_ids, width >= min_inv_bp) %>%
      dplyr::mutate(
        variant = "Substructural_inversion",
        variant_specific = dplyr::if_else(
          X_dist < localized_threshold,
          "symmetric",
          "localized"
        ),
        context = "outside_all_SR"
      ) %>%
      dplyr::select(rid, qid, start = rs, end = re, width, meanlen, variant, variant_specific, context)
  }
  
  # ---- Combine & deduplicate ----
  dplyr::bind_rows(inside_calls, outside_calls) %>%
    dplyr::distinct(start, end, rid, qid, variant, context, .keep_all = TRUE) %>%
    dplyr::arrange(start)
}

delta_all_SVs <- function(
    delta_file, 
    unfiltered_delta_file, 
    species_name
) {
  
  ### Indels
  delta_indels_i <- delta_indels(delta_file) %>%
    transmute(
      rid, qid, 
      start, end, width,
      refmid = round((start + end) / 2),
      variant,
      variant_specific,
      species = species_name
    )
  
  if (nrow(delta_indels_i) == 0) {
    delta_indels_i <- tibble(
      rid = character(),
      qid = character(),
      start = numeric(),
      end = numeric(),
      width = numeric(),
      refmid = numeric(),
      variant = character(),
      variant_specific = character(),
      species = character()
    )
  }
  
  ### Duplications
  delta_dups_i <- delta_duplications(
    delta_table = delta_file, 
    unfiltered_delta_table = unfiltered_delta_file
  ) %>% 
    transmute(
      rid, qid,
      start, end, width,
      refmid = round((start + end) / 2),
      variant,
      variant_specific,
      species = species_name
    )
  
  if (nrow(delta_dups_i) == 0) {
    delta_dups_i <- delta_indels_i[0, ]
  }
  
  ### Translocations
  delta_tlocs_i <- delta_translocations(delta_file) %>%
    transmute(
      rid, qid, 
      start, end, width,
      refmid = round((start + end) / 2),
      variant,
      variant_specific,
      species = species_name
    )
  
  if (nrow(delta_tlocs_i) == 0) {
    delta_tlocs_i <- delta_indels_i[0, ]
  }
  
  ### Structural rearrangements
  delta_SR_temp <- delta_structural_rearrangements(delta_file)
  
  delta_subSR_i <- delta_substructural_inversions(delta_file, delta_SR_temp) %>%
    ungroup() %>%
    transmute(
      rid, qid,
      start, end, width = meanlen,
      refmid = round((start + end) / 2),
      variant,
      variant_specific,
      species = species_name
    )
  
  if (nrow(delta_subSR_i) == 0) {
    delta_subSR_i <- delta_indels_i[0, ]
  }
  
  delta_SR_i <- delta_SR_temp %>%
    transmute(
      rid, qid,
      start, end, width,
      refmid = round((start + end) / 2),
      variant,
      variant_specific,
      species = species_name
    )
  
  if (nrow(delta_SR_i) == 0) {
    delta_SR_i <- delta_indels_i[0, ]
  }
  
  ### Combine all SVs
  all_SVs <- bind_rows(
    delta_indels_i,
    delta_dups_i,
    delta_tlocs_i,
    delta_subSR_i,
    delta_SR_i
  ) %>% 
    arrange(desc(width))
  
  return(all_SVs)
}

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

group_svs <- function(sv_tib, min_overlap_bp = 1L, ro_min = 0.95) {
  stopifnot(all(c("start","end","width") %in% names(sv_tib)))
  n <- nrow(sv_tib); if (n == 0L) return(sv_tib)
  
  start <- as.integer(sv_tib$start)
  end   <- as.integer(sv_tib$end)
  w     <- as.integer(sv_tib$width)
  
  ok <- !is.na(start) & !is.na(end) & !is.na(w) & (start <= end) & (w > 0L)
  sv_tib <- sv_tib[ok, , drop = FALSE]; start <- start[ok]; end <- end[ok]; w <- w[ok]
  n <- nrow(sv_tib); if (n == 0L) return(sv_tib)
  
  rng  <- IRanges::IRanges(start = start, end = end)
  hits <- IRanges::findOverlaps(rng, drop.self = TRUE, drop.redundant = TRUE)
  
  if (length(hits) == 0L) {
    sv_tib$sv_group <- seq_len(n); sv_tib$group_size <- 1L; return(sv_tib)
  }
  
  i <- as.integer(S4Vectors::queryHits(hits))
  j <- as.integer(S4Vectors::subjectHits(hits))
  
  # overlap (bp) on original ranges
  ovl <- pmax.int(0L, pmin(end[i], end[j]) - pmax(start[i], start[j]) + 1L)
  ro  <- pmin(ovl / w[i], ovl / w[j])
  
  keep <- (ovl >= min_overlap_bp) & (ro >= ro_min)
  if (!any(keep)) {
    sv_tib$sv_group <- seq_len(n); sv_tib$group_size <- 1L; return(sv_tib)
  }
  
  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  igraph::V(g)$name <- as.character(seq_len(n))
  g <- igraph::add_edges(g, as.vector(rbind(i[keep], j[keep])))
  
  memb <- igraph::components(g)$membership
  grp  <- as.integer(memb[seq_len(n)])
  
  sv_tib$sv_group <- grp
  sizes <- tabulate(grp, nbins = max(grp, na.rm = TRUE))
  sv_tib$group_size <- sizes[grp]
  sv_tib
}
