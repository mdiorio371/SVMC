# 07-annotation.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- annotate_SVs  [from SVCM.R:2075] ---
annotate_SVs <- function(
    SV_granges, gene_granges,
    mode = c("both","sv","gene","either","jaccard"),
    sv_thresh = 80, gene_thresh = 80, jaccard_thresh = 0.5,
    keep = c("best","all"),
    bp_pad = 100L, breakpoint_min_width = 50000L,
    genome_len = NULL
) {
  mode <- match.arg(mode); keep <- match.arg(keep)
  stopifnot(inherits(SV_granges, "GRanges"), inherits(gene_granges, "GRanges"))
  if (length(SV_granges) == 0) return(tibble::tibble())

  common <- intersect(GenomeInfoDb::seqlevels(SV_granges), GenomeInfoDb::seqlevels(gene_granges))
  if (!length(common)) return(tibble::tibble())
  SVg   <- GenomeInfoDb::keepSeqlevels(SV_granges, common, pruning.mode = "coarse")
  Geneg <- GenomeInfoDb::keepSeqlevels(gene_granges, common, pruning.mode = "coarse")

  .mcol <- function(gr, nm) if (nm %in% names(S4Vectors::mcols(gr))) S4Vectors::mcols(gr)[[nm]] else rep(NA_character_, length(gr))
  sv_species <- .mcol(SVg, "species")
  sv_variant <- { v <- .mcol(SVg, "variant"); if (all(is.na(v))) .mcol(SVg, "variant_specific") else v }
  sv_refmid  <- .mcol(SVg, "refmid")
  gene_gene  <- .mcol(Geneg, "gene")
  gene_prod  <- .mcol(Geneg, "product")
  sv_w <- IRanges::width(SVg); gene_w <- IRanges::width(Geneg)

  # ---- SPAN overlaps ----
  hits <- GenomicRanges::findOverlaps(SVg, Geneg, ignore.strand = TRUE)
  if (length(hits)) {
    qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
    ov <- IRanges::pintersect(SVg[qh], Geneg[sh])
    ov_bp <- IRanges::width(ov)
    pct_sv <- 100 * ov_bp / sv_w[qh]; pct_gene <- 100 * ov_bp / gene_w[sh]
    jacc <- ov_bp / (sv_w[qh] + gene_w[sh] - ov_bp)

    fully <- switch(mode,
      "both"    = (pct_sv >= sv_thresh) & (pct_gene >= gene_thresh),
      "sv"      = (pct_sv >= sv_thresh),
      "gene"    = (pct_gene >= gene_thresh),
      "either"  = (pct_sv >= sv_thresh) | (pct_gene >= gene_thresh),
      "jaccard" = (jacc >= jaccard_thresh)
    )

    span_df <- tibble::tibble(
      sv_idx = qh,
      gene_accession = as.character(GenomeInfoDb::seqnames(Geneg))[sh],
      gene = ifelse(fully, gene_gene[sh], NA_character_),
      product = ifelse(fully, gene_prod[sh], NA_character_),
      overlap_bp = as.integer(ov_bp),
      pct_overlap_sv = as.numeric(pct_sv),
      pct_overlap_gene = as.numeric(pct_gene),
      jaccard = as.numeric(jacc),
      span_class = ifelse(fully, "genic", "partially_genic")
    )
    if (keep == "best") {
      span_df <- span_df %>%
        dplyr::group_by(.data$sv_idx) %>%
        dplyr::slice_max(order_by = .data$overlap_bp, with_ties = FALSE) %>%
        dplyr::ungroup()
    }
  } else {
    span_df <- tibble::tibble()
  }

  base <- tibble::tibble(
    sv_idx = seq_along(SVg),
    sv_accession = as.character(GenomeInfoDb::seqnames(SVg)),
    species = sv_species, variant = sv_variant,
    sv_start = BiocGenerics::start(SVg), sv_end = BiocGenerics::end(SVg),
    sv_width = sv_w, refmid = sv_refmid
  )

  out <- base %>%
    dplyr::left_join(span_df, by = "sv_idx") %>%
    dplyr::mutate(
      span_class = dplyr::coalesce(.data$span_class, "intergenic"),
      overlap_bp = dplyr::coalesce(.data$overlap_bp, 0L)
    )

  # ---- BREAKPOINT genic/intergenic for large SVs ----
  bp_pad <- as.integer(bp_pad)
  if (isTRUE(bp_pad > 0L)) {
    big <- which(is.finite(out$sv_width) & out$sv_width >= as.numeric(breakpoint_min_width))
    out$bp1_is_genic <- FALSE; out$bp2_is_genic <- FALSE
    out$bp1_gene <- NA_character_; out$bp2_gene <- NA_character_
    out$bp1_product <- NA_character_; out$bp2_product <- NA_character_
    out$bp_any_genic <- FALSE; out$bp_both_genic <- FALSE

    if (length(big)) {
      b_start <- out$sv_start[big]; b_end <- out$sv_end[big]
      gr1 <- GenomicRanges::GRanges(
        seqnames = out$sv_accession[big],
        ranges = IRanges::IRanges(start = pmax(1L, b_start - bp_pad), end = b_start + bp_pad)
      )
      gr2 <- GenomicRanges::GRanges(
        seqnames = out$sv_accession[big],
        ranges = IRanges::IRanges(start = pmax(1L, b_end - bp_pad), end = b_end + bp_pad)
      )
      if (!is.null(genome_len)) {
        genome_len <- as.integer(genome_len)
        BiocGenerics::end(gr1) <- pmin(BiocGenerics::end(gr1), genome_len)
        BiocGenerics::end(gr2) <- pmin(BiocGenerics::end(gr2), genome_len)
      }

      .best_bp_hit <- function(bp_gr) {
        h <- GenomicRanges::findOverlaps(bp_gr, Geneg, ignore.strand = TRUE)
        if (!length(h)) return(list(hit = rep(FALSE, length(bp_gr)),
                                    gene = rep(NA_character_, length(bp_gr)),
                                    product = rep(NA_character_, length(bp_gr))))
        qh <- S4Vectors::queryHits(h); sh <- S4Vectors::subjectHits(h)
        ov <- IRanges::width(IRanges::pintersect(bp_gr[qh], Geneg[sh]))
        df <- tibble::tibble(q = qh, s = sh, ov = ov) %>%
          dplyr::group_by(.data$q) %>%
          dplyr::slice_max(.data$ov, with_ties = FALSE) %>%
          dplyr::ungroup()
        hit <- rep(FALSE, length(bp_gr)); hit[df$q] <- TRUE
        gene <- rep(NA_character_, length(bp_gr)); gene[df$q] <- gene_gene[df$s]
        prod <- rep(NA_character_, length(bp_gr)); prod[df$q] <- gene_prod[df$s]
        list(hit = hit, gene = gene, product = prod)
      }

      h1 <- .best_bp_hit(gr1); h2 <- .best_bp_hit(gr2)
      out$bp1_is_genic[big] <- h1$hit; out$bp2_is_genic[big] <- h2$hit
      out$bp1_gene[big] <- h1$gene; out$bp2_gene[big] <- h2$gene
      out$bp1_product[big] <- h1$product; out$bp2_product[big] <- h2$product
      out$bp_any_genic[big]  <- out$bp1_is_genic[big] | out$bp2_is_genic[big]
      out$bp_both_genic[big] <- out$bp1_is_genic[big] & out$bp2_is_genic[big]
    }
  }
  out
}

# --- classify_svs_one  [from SVCM.R:2220] ---
classify_svs_one <- function(
    sv_tbl, cds_gr, bp_pad = 100L, min_large_bp = 50000L,
    genome_len = NULL,
    type_col = c("variant_specific","variant","sv_type","type","svtype")
) {
  stopifnot(is.data.frame(sv_tbl), all(c("start","end","width") %in% names(sv_tbl)))
  stopifnot(inherits(cds_gr, "GRanges"))

  type_col <- intersect(type_col, names(sv_tbl))
  type_col <- if (length(type_col)) type_col[1] else NULL

  x <- sv_tbl %>%
    dplyr::mutate(
      start = as.integer(.data$start), end = as.integer(.data$end),
      width = as.numeric(.data$width),
      sv_type = if (!is.null(type_col)) as.character(.data[[type_col]]) else "Unknown",
      is_large = is.finite(.data$width) & .data$width >= as.numeric(min_large_bp)
    )

  SVg <- GenomicRanges::GRanges(
    seqnames = "chr",
    ranges = IRanges::IRanges(start = pmin(x$start, x$end), end = pmax(x$start, x$end))
  )
  S4Vectors::mcols(SVg)$variant <- x$sv_type

  ann <- annotate_SVs(
    SV_granges = SVg, gene_granges = cds_gr,
    keep = "best", bp_pad = bp_pad,
    breakpoint_min_width = min_large_bp, genome_len = genome_len
  )

  x$bp1_is_genic  <- ann$bp1_is_genic  %||% FALSE
  x$bp2_is_genic  <- ann$bp2_is_genic  %||% FALSE
  x$bp_any_genic  <- ann$bp_any_genic  %||% (x$bp1_is_genic | x$bp2_is_genic)
  x$bp_both_genic <- ann$bp_both_genic %||% (x$bp1_is_genic & x$bp2_is_genic)
  x
}

# --- enrichment_from_categories  [from SVCM.R:2268] ---
enrichment_from_categories <- function(
    classified_sv, cds_gr, genome_len, bp_pad = 100L,
    group_cols = c("sv_type"), filter_large_only = TRUE
) {
  stopifnot(is.data.frame(classified_sv), inherits(cds_gr, "GRanges"))
  genome_len <- as.integer(genome_len)
  stopifnot(genome_len > 0L)
  pad <- as.integer(bp_pad)

  # Expected probability from CDS coverage (expanded by ±pad)
  cds_exp <- GenomicRanges::resize(cds_gr, width = IRanges::width(cds_gr) + 2L*pad, fix = "center")
  GenomicRanges::start(cds_exp) <- pmax(1L, GenomicRanges::start(cds_exp))
  GenomicRanges::end(cds_exp)   <- pmin(genome_len, GenomicRanges::end(cds_exp))
  cov_bp <- sum(IRanges::width(IRanges::reduce(GenomicRanges::ranges(cds_exp))))
  p_exp <- cov_bp / genome_len

  df <- classified_sv
  if (filter_large_only && "is_large" %in% names(df)) df <- dplyr::filter(df, .data$is_large)

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_svs = dplyr::n(),
      n_any  = sum(.data$bp_any_genic, na.rm = TRUE),
      n_both = sum(.data$bp_both_genic, na.rm = TRUE),
      frac_any  = n_any / n_svs,
      frac_both = n_both / n_svs,
      p_exp  = p_exp,
      p_any  = stats::binom.test(n_any, n_svs, p = p_exp)$p.value,
      p_both = stats::binom.test(n_both, n_svs, p = p_exp^2)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      q_any  = stats::p.adjust(.data$p_any, method = "BH"),
      q_both = stats::p.adjust(.data$p_both, method = "BH")
    )
}

# --- make_cds_granges  [from SVCM6.R:403] ---
make_cds_granges <- function(gff_path, sv_rid = NULL) {
  if (exists("read_prokka_gff2", mode = "function")) {
    prokka <- read_prokka_gff2(gff_path)
    cds <- prokka$cds
  } else {
    cds <- read_prokka_cds_fallback(gff_path)
  }
  
  # Optional remap of seqid convention to sv_rid (only if needed)
  if (!is.null(sv_rid)) {
    seqs <- unique(cds$seqid)
    if (!sv_rid %in% seqs) {
      message(sprintf("  NOTE: Prokka seqid(s) do not include SV rid '%s'. Remapping CDS seqid -> sv_rid.", sv_rid))
      cds$seqid <- sv_rid
    }
  }
  
  GenomicRanges::GRanges(
    seqnames = cds$seqid,
    ranges   = IRanges::IRanges(start = cds$start, end = cds$end),
    strand   = cds$strand,
    gene     = cds$gene,
    product  = cds$product
  )
}

# --- annotate_gene_pile_cds  [from SVCM6.R:430] ---
annotate_gene_pile_cds <- function(gene_pile_df, cds_gr, gene_thresh = 80) {
  # gene_pile_df must have rid, start, end, qid, variant, width (at least)
  need <- c("rid","start","end","qid","variant","width")
  miss <- setdiff(need, names(gene_pile_df))
  if (length(miss) > 0) stop("gene_pile missing columns: ", paste(miss, collapse = ", "))
  
  sv_gr <- GenomicRanges::GRanges(
    seqnames = gene_pile_df$rid,
    ranges   = IRanges::IRanges(start = gene_pile_df$start, end = gene_pile_df$end)
  )
  
  hits <- GenomicRanges::findOverlaps(sv_gr, cds_gr, ignore.strand = TRUE)
  
  # percent CDS coverage of each SV interval (union CDS overlaps per SV)
  sv_len <- IRanges::width(GenomicRanges::ranges(sv_gr))
  covered <- numeric(length(sv_gr))
  best_gene <- rep(NA_character_, length(sv_gr))
  best_prod <- rep(NA_character_, length(sv_gr))
  best_olap <- numeric(length(sv_gr))
  
  if (length(hits) > 0) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    
    ol <- IRanges::pintersect(GenomicRanges::ranges(sv_gr)[qh],
                              GenomicRanges::ranges(cds_gr)[sh])
    ol_w <- IRanges::width(ol)
    
    # accumulate coverage per SV via max non-overlapping union approximation:
    # For CDS-only classification, a conservative sum of overlaps is okay
    # (Prokka CDS features rarely heavily overlap).
    for (k in seq_along(ol_w)) {
      i <- qh[k]
      covered[i] <- covered[i] + ol_w[k]
      if (ol_w[k] > best_olap[i]) {
        best_olap[i] <- ol_w[k]
        best_gene[i] <- mcols(cds_gr)$gene[sh[k]]
        best_prod[i] <- mcols(cds_gr)$product[sh[k]]
      }
    }
  }
  
  pct <- 100 * covered / pmax(sv_len, 1)
  span_class <- ifelse(pct >= gene_thresh, "genic",
                       ifelse(pct <= 0, "intergenic", "partial"))
  
  out <- gene_pile_df
  out$cds_cov_bp   <- covered
  out$cds_cov_pct  <- pct
  out$span_class   <- span_class
  out$best_gene    <- best_gene
  out$best_product <- best_prod
  out
}

# --- read_prokka_cds_fallback  [from SVCM6.R:376] ---
read_prokka_cds_fallback <- function(gff_path) {
  if (!file.exists(gff_path)) stop("Prokka GFF not found: ", gff_path)
  
  # read GFF3, keep CDS
  x <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    x <- data.table::fread(gff_path, sep = "\t", header = FALSE, data.table = FALSE, fill = TRUE, quote = "")
  } else {
    x <- utils::read.delim(gff_path, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  }
  if (ncol(x) < 9) stop("GFF parse failed (expected 9 columns): ", gff_path)
  colnames(x)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  x <- x[x$type == "CDS", , drop = FALSE]
  if (nrow(x) == 0) stop("No CDS rows found in Prokka GFF: ", gff_path)
  
  # parse attributes: ID, gene, product
  get_attr <- function(attr, key) {
    m <- regmatches(attr, regexpr(paste0(key, "=[^;]+"), attr))
    m <- sub(paste0("^", key, "="), "", m)
    ifelse(nchar(m) == 0, NA_character_, m)
  }
  x$gene    <- get_attr(x$attributes, "gene")
  x$product <- get_attr(x$attributes, "product")
  
  x
}

# --- svcm_read_prokka_gff_cds  [from SVCM6.R:1024] ---
svcm_read_prokka_gff_cds <- function(gff_path) {
  svcm_stop_if_missing(c("tibble"))
  
  if (!file.exists(gff_path)) stop("Missing Prokka GFF: ", gff_path, call. = FALSE)
  
  # Prefer rtracklayer if available (most robust)
  if (requireNamespace("rtracklayer", quietly = TRUE) &&
      requireNamespace("GenomicRanges", quietly = TRUE) &&
      requireNamespace("IRanges", quietly = TRUE) &&
      requireNamespace("S4Vectors", quietly = TRUE)) {
    
    gr <- rtracklayer::import(gff_path)
    cds <- gr[GenomicRanges::mcols(gr)$type %in% c("CDS")]
    
    if (length(cds) == 0L) stop("No CDS features found in Prokka GFF: ", gff_path, call. = FALSE)
    
    # Prokka attributes typically include: ID, locus_tag, gene, product
    m <- GenomicRanges::mcols(cds)
    gene    <- if ("gene" %in% names(m)) as.character(m$gene) else NA_character_
    product <- if ("product" %in% names(m)) as.character(m$product) else NA_character_
    locus   <- if ("locus_tag" %in% names(m)) as.character(m$locus_tag) else NA_character_
    
    GenomicRanges::mcols(cds)$gene <- gene
    GenomicRanges::mcols(cds)$product <- product
    GenomicRanges::mcols(cds)$locus_tag <- locus
    
    return(cds)
  }
  
  # Fallback: tab parse (works, but attributes parsing is simplistic)
  svcm_stop_if_missing(c("data.table", "GenomicRanges", "IRanges", "S4Vectors"))
  x <- data.table::fread(gff_path, sep = "\t", header = FALSE, data.table = FALSE,
                         fill = TRUE, quote = "", comment.char = "#")
  if (ncol(x) < 9) stop("GFF parse failed; expected >=9 columns: ", gff_path, call. = FALSE)
  colnames(x)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  x <- x[x$type == "CDS", , drop = FALSE]
  if (!nrow(x)) stop("No CDS rows in Prokka GFF: ", gff_path, call. = FALSE)
  
  get_attr <- function(attr, key) {
    m <- regmatches(attr, regexpr(paste0(key, "=[^;]+"), attr))
    m <- sub(paste0("^", key, "="), "", m)
    ifelse(nchar(m) == 0, NA_character_, m)
  }
  
  gene    <- get_attr(x$attributes, "gene")
  product <- get_attr(x$attributes, "product")
  locus   <- get_attr(x$attributes, "locus_tag")
  
  GenomicRanges::GRanges(
    seqnames = x$seqid,
    ranges   = IRanges::IRanges(start = x$start, end = x$end),
    strand   = x$strand,
    gene     = gene,
    product  = product,
    locus_tag= locus
  )
}

# --- svcm_annotate_gene_pile_cds  [from SVCM6.R:1085] ---
svcm_annotate_gene_pile_cds <- function(gene_pile_df,
                                        cds_gr,
                                        gene_thresh = 80) {
  svcm_stop_if_missing(c("GenomicRanges", "IRanges", "S4Vectors", "tibble"))
  
  gp <- tibble::as_tibble(gene_pile_df)
  need <- c("rid","start","end","qid","variant","width")
  miss <- setdiff(need, names(gp))
  if (length(miss)) stop("gene_pile missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  
  # SV ranges
  sv_gr <- GenomicRanges::GRanges(
    seqnames = gp$rid,
    ranges   = IRanges::IRanges(start = gp$start, end = gp$end)
  )
  sv_len <- IRanges::width(GenomicRanges::ranges(sv_gr))
  
  hits <- GenomicRanges::findOverlaps(sv_gr, cds_gr, ignore.strand = TRUE)
  
  covered_bp <- numeric(length(sv_gr))
  best_gene  <- rep(NA_character_, length(sv_gr))
  best_prod  <- rep(NA_character_, length(sv_gr))
  best_ov    <- numeric(length(sv_gr))
  
  if (length(hits) > 0) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    
    # intersection widths for each hit
    inter <- GenomicRanges::pintersect(
      GenomicRanges::ranges(sv_gr)[qh],
      GenomicRanges::ranges(cds_gr)[sh]
    )
    inter_w <- IRanges::width(inter)
    
    # union coverage per SV via reduce() on intersections
    inter_ir <- IRanges::IRanges(IRanges::start(inter), IRanges::end(inter))
    split_ir <- split(inter_ir, qh)
    
    qs <- as.integer(names(split_ir))
    covered_bp[qs] <- vapply(split_ir, function(ir) sum(IRanges::width(IRanges::reduce(ir))), numeric(1))
    
    # best gene/product by max overlap hit
    for (k in seq_along(inter_w)) {
      i <- qh[k]
      if (inter_w[k] > best_ov[i]) {
        best_ov[i] <- inter_w[k]
        mg <- GenomicRanges::mcols(cds_gr)
        best_gene[i] <- if ("gene" %in% names(mg)) as.character(mg$gene[sh[k]]) else NA_character_
        best_prod[i] <- if ("product" %in% names(mg)) as.character(mg$product[sh[k]]) else NA_character_
      }
    }
  }
  
  pct <- 100 * covered_bp / pmax(sv_len, 1)
  span_class <- dplyr::case_when(
    pct >= gene_thresh ~ "genic",
    covered_bp <= 0    ~ "intergenic",
    TRUE               ~ "partial"
  )
  
  gp$cds_cov_bp   <- covered_bp
  gp$cds_cov_pct  <- pct
  gp$span_class   <- span_class
  gp$best_gene    <- best_gene
  gp$best_product <- best_prod
  gp
}

# --- build_prokka_cmd  [from i9.R:12] ---
build_prokka_cmd <- function(
    ref_fasta,
    outdir,
    prefix,
    force = TRUE,
    genus = NULL,
    species = NULL,
    strain = NULL,
    kingdom = "Bacteria",
    cpus = 1L
) {
  parts <- c(
    "prokka",
    "--outdir", shQuote(outdir),
    "--prefix", shQuote(prefix),
    "--kingdom", shQuote(kingdom),
    "--cpus", as.integer(cpus)
  )
  
  if (isTRUE(force)) {
    parts <- c(parts, "--force")
  }
  if (!is.null(genus)) {
    parts <- c(parts, "--genus", shQuote(genus))
  }
  if (!is.null(species)) {
    parts <- c(parts, "--species", shQuote(species))
  }
  if (!is.null(strain)) {
    parts <- c(parts, "--strain", shQuote(strain))
  }
  
  parts <- c(parts, shQuote(ref_fasta))
  paste(parts, collapse = " ")
}

# --- svcm_prokka_cmd  [from SVCM6.R:1157] ---
svcm_prokka_cmd <- function(ref_fasta, outdir, prefix, cpus = 8L) {
  sprintf(
    'prokka --outdir %s --prefix %s --force --cpus %d %s',
    shQuote(outdir), shQuote(prefix), as.integer(cpus), shQuote(ref_fasta)
  )
}

# --- annotate_breakpoints_simple  [from i9.R:659] ---
annotate_breakpoints_simple <- function(x, cds_tbl) {
  x <- tibble::as_tibble(x)
  
  if (!"seqid" %in% names(x)) {
    x$seqid <- unique(cds_tbl$seqid)[1]
  }
  
  gene_gr <- GenomicRanges::GRanges(
    seqnames = cds_tbl$seqid,
    ranges   = IRanges::IRanges(start = cds_tbl$start, end = cds_tbl$end),
    strand   = cds_tbl$strand
  )
  S4Vectors::mcols(gene_gr)$locus_tag <- cds_tbl$locus_tag
  S4Vectors::mcols(gene_gr)$gene      <- cds_tbl$gene
  S4Vectors::mcols(gene_gr)$product   <- cds_tbl$product
  
  bp1_gr <- GenomicRanges::GRanges(
    seqnames = x$seqid,
    ranges   = IRanges::IRanges(start = x$start, end = x$start)
  )
  bp2_gr <- GenomicRanges::GRanges(
    seqnames = x$seqid,
    ranges   = IRanges::IRanges(start = x$end, end = x$end)
  )
  
  hit1 <- GenomicRanges::findOverlaps(bp1_gr, gene_gr, select = "first")
  hit2 <- GenomicRanges::findOverlaps(bp2_gr, gene_gr, select = "first")
  
  x %>%
    mutate(
      bp1_in_gene = !is.na(hit1),
      bp2_in_gene = !is.na(hit2),
      bp1_locus_tag = ifelse(is.na(hit1), NA_character_, S4Vectors::mcols(gene_gr)$locus_tag[hit1]),
      bp2_locus_tag = ifelse(is.na(hit2), NA_character_, S4Vectors::mcols(gene_gr)$locus_tag[hit2]),
      bp1_gene = ifelse(is.na(hit1), NA_character_, S4Vectors::mcols(gene_gr)$gene[hit1]),
      bp2_gene = ifelse(is.na(hit2), NA_character_, S4Vectors::mcols(gene_gr)$gene[hit2]),
      bp1_product = ifelse(is.na(hit1), NA_character_, S4Vectors::mcols(gene_gr)$product[hit1]),
      bp2_product = ifelse(is.na(hit2), NA_character_, S4Vectors::mcols(gene_gr)$product[hit2]),
      breakpoint_context = case_when(
        bp1_in_gene & bp2_in_gene ~ "both_genic",
        bp1_in_gene | bp2_in_gene ~ "one_genic",
        TRUE ~ "intergenic"
      )
    )
}

# --- summarize_breakpoint_context_by_role  [from i9.R:117] ---
summarize_breakpoint_context_by_role <- function(x) {
  x %>%
    group_by(role) %>%
    summarise(
      n = n(),
      prop_both_genic = mean(breakpoint_context == "both_genic"),
      prop_one_genic  = mean(breakpoint_context == "one_genic"),
      prop_intergenic = mean(breakpoint_context == "intergenic"),
      .groups = "drop"
    )
}

# --- summarize_breakpoint_context_by_motif  [from i9.R:129] ---
summarize_breakpoint_context_by_motif <- function(x) {
  x %>%
    group_by(motif_uid, role) %>%
    summarise(
      n_events = n(),
      prop_both_genic = mean(breakpoint_context == "both_genic"),
      prop_one_genic  = mean(breakpoint_context == "one_genic"),
      prop_intergenic = mean(breakpoint_context == "intergenic"),
      .groups = "drop"
    )
}
