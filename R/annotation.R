# ==========================================================================
# SVMC: gene annotation of SVs
# Prokka GFF import and SV-to-gene annotation (canonical annotate_SVs).
# ==========================================================================

# Note: annotate_SVs() here is the former annotate_SVs4 (renamed).
# annotate_SVs2() and annotate_SVs3() were removed as superseded.

read_prokka_gff2 <- function(path,
                             keep_types = c("gene","CDS","rRNA","tRNA","tmRNA")) {
  gr <- rtracklayer::import(path)                 # auto-detects GFF3
  if (!is.null(keep_types)) {
    has_type <- "type" %in% names(mcols(gr))
    if (has_type) gr <- gr[mcols(gr)$type %in% keep_types]
  }
  
  df <- as_tibble(as.data.frame(gr))
  # Ensure expected columns exist
  for (nm in c("ID","Parent","gene","Name","locus_tag","product","phase","type","source")) {
    if (!nm %in% names(df)) df[[nm]] <- NA
  }
  
  features <- df %>%
    transmute(
      seqid     = as.character(seqnames),
      source    = as.character(source),
      type      = as.character(type),
      start     = as.integer(start),
      end       = as.integer(end),
      strand    = as.character(strand),
      phase     = suppressWarnings(as.integer(phase)),
      ID        = as.character(ID),
      Parent    = as.character(Parent),
      # prefer gene, else Name, else locus_tag/product as a last resort
      gene      = coalesce_chr(as.character(gene), as.character(Name)),
      locus_tag = as.character(locus_tag),
      product   = as.character(product)
    )
  
  # --- collapse CDS fragments to one row per locus_tag (gene-level) ---
  cds_collapsed <- features %>%
    filter(type == "CDS") %>%
    group_by(seqid, locus_tag) %>%
    summarise(
      start     = min(start, na.rm = TRUE),
      end       = max(end, na.rm = TRUE),
      strand    = na.omit(strand)[1],
      n_parts   = dplyr::n(),
      gene      = na.omit(gene)[1],
      product   = na.omit(product)[1],
      .groups   = "drop"
    ) %>%
    mutate(length = end - start + 1L)
  
  list(features = features, cds = cds_collapsed)
}


annotate_SVs <- 
  function(
    SV_granges, 
    gene_granges,
    sv_thresh = 80, 
    gene_thresh = 80,
    mode = c("both","sv","gene","either","jaccard"),
    jaccard_thresh = 0.5,
    keep = c("best","all"),
    bp_pad = 0L,
    inv_like = c("inversion","local inversion","rearrangement","translocation","RAI"),
    verbose = FALSE) {
  
  mode <- match.arg(mode)
  keep <- match.arg(keep)
  
  # Input validation
  if (!is(SV_granges, "GRanges")) stop("SV_granges must be a GRanges object")
  if (!is(gene_granges, "GRanges")) stop("gene_granges must be a GRanges object")
  if (length(SV_granges) == 0) return(tibble())
  if (length(gene_granges) == 0) {
    warning("No genes provided - all SVs will be marked as intergenic")
    return(tibble(
      hit_type = "none",
      sv_accession = as.character(seqnames(SV_granges)),
      sv_start = start(SV_granges),
      sv_end = end(SV_granges),
      sv_width = width(SV_granges),
      func_group = "intergenic"
    ))
  }
  
  # --- 0) Align seqlevels safely
  common_seq <- intersect(seqlevels(SV_granges), seqlevels(gene_granges))
  
  if (length(common_seq) == 0) {
    warning("No common sequence names between SVs and genes")
    return(tibble())
  }
  
  SVg <- keepSeqlevels(SV_granges, common_seq, pruning.mode = "coarse")
  Geneg <- keepSeqlevels(gene_granges, common_seq, pruning.mode = "coarse")
  
  if (verbose) {
    cat("Processing", length(SVg), "SVs against", length(Geneg), "genes\n")
    cat("Common sequences:", length(common_seq), "\n")
  }
  
  # Helper function with better error handling
  get_col <- function(gr, nm) {
    if (nm %in% names(mcols(gr))) {
      mcols(gr)[[nm]]
    } else {
      rep(NA_character_, length(gr))
    }
  }
  
  # Extract metadata
  sv_species <- get_col(SVg, "species")
  sv_variant <- {
    v <- get_col(SVg, "variant")
    if (all(is.na(v))) get_col(SVg, "variant_specific") else v
  }
  sv_refmid <- get_col(SVg, "refmid")
  gene_gene <- get_col(Geneg, "gene")
  gene_prod <- get_col(Geneg, "product")
  
  # Precompute widths
  sv_w <- width(SVg)
  gene_w <- width(Geneg)
  
  # --- 1A) SPAN overlaps
  span_hits <- findOverlaps(SVg, Geneg, ignore.strand = TRUE)
  
  if (length(span_hits) > 0) {
    qh <- queryHits(span_hits)
    sh <- subjectHits(span_hits)
    span_ov <- pintersect(SVg[qh], Geneg[sh])
    span_bp <- width(span_ov)
    pct_sv <- 100 * span_bp / sv_w[qh]
    pct_gene <- 100 * span_bp / gene_w[sh]
    jacc <- span_bp / (sv_w[qh] + gene_w[sh] - span_bp)
    
    fully_genic <- switch(
      mode,
      "both"    = (pct_sv >= sv_thresh) & (pct_gene >= gene_thresh),
      "sv"      = (pct_sv >= sv_thresh),
      "gene"    = (pct_gene >= gene_thresh),
      "either"  = (pct_sv >= sv_thresh) | (pct_gene >= gene_thresh),
      "jaccard" = (jacc >= jaccard_thresh)
    )
    
    span_df <- tibble(
      hit_type         = "span",
      gene_accession   = as.character(seqnames(Geneg))[sh],
      gene             = ifelse(fully_genic, gene_gene[sh], NA_character_),
      product          = ifelse(fully_genic, gene_prod[sh], NA_character_),
      sv_accession     = as.character(seqnames(SVg))[qh],
      species          = sv_species[qh],
      variant          = sv_variant[qh],
      sv_start         = start(SVg)[qh],
      sv_end           = end(SVg)[qh],
      sv_width         = sv_w[qh],
      refmid           = sv_refmid[qh],
      overlap_bp       = as.integer(span_bp),
      pct_overlap_sv   = as.numeric(pct_sv),
      pct_overlap_gene = as.numeric(pct_gene),
      jaccard          = as.numeric(jacc),
      fully_genic      = fully_genic
    )
  } else {
    span_df <- tibble()
  }
  
  # --- 1B) BREAKPOINT overlaps
  if (bp_pad < 0L) bp_pad <- 0L
  
  # Improved breakpoint range creation
  if (bp_pad > 0 || length(grep(paste(tolower(inv_like), collapse="|"), 
                                tolower(sv_variant), ignore.case=TRUE)) > 0) {
    
    # Create breakpoint ranges
    bp_starts <- GRanges(
      seqnames = seqnames(SVg),
      ranges = IRanges(start = pmax(1, start(SVg) - bp_pad),
                       end = start(SVg) + bp_pad)
    )
    bp_ends <- GRanges(
      seqnames = seqnames(SVg),
      ranges = IRanges(start = pmax(1, end(SVg) - bp_pad),
                       end = end(SVg) + bp_pad)
    )
    
    # Combine and add SV index
    BPg <- c(bp_starts, bp_ends)
    mcols(BPg)$sv_idx <- c(seq_along(SVg), seq_along(SVg))
    mcols(BPg)$bp_side <- c(rep("start", length(SVg)), rep("end", length(SVg)))
    
    bp_hits <- findOverlaps(BPg, Geneg, ignore.strand = TRUE)
    
    if (length(bp_hits) > 0) {
      bq_idx <- mcols(BPg)$sv_idx[queryHits(bp_hits)]
      bs <- subjectHits(bp_hits)
      bp_side <- mcols(BPg)$bp_side[queryHits(bp_hits)]
      
      bp_ov <- pintersect(BPg[queryHits(bp_hits)], Geneg[bs])
      bp_bp <- width(bp_ov)
      
      bp_df <- tibble(
        hit_type         = "breakpoint",
        gene_accession   = as.character(seqnames(Geneg))[bs],
        gene             = gene_gene[bs],
        product          = gene_prod[bs],
        sv_accession     = as.character(seqnames(SVg))[bq_idx],
        species          = sv_species[bq_idx],
        variant          = sv_variant[bq_idx],
        sv_start         = start(SVg)[bq_idx],
        sv_end           = end(SVg)[bq_idx],
        sv_width         = sv_w[bq_idx],
        refmid           = sv_refmid[bq_idx],
        overlap_bp       = as.integer(bp_bp),
        pct_overlap_sv   = NA_real_,
        pct_overlap_gene = 100 * as.numeric(bp_bp) / gene_w[bs],
        jaccard          = NA_real_,
        fully_genic      = TRUE,
        bp_side          = bp_side  # Added: which breakpoint
      )
    } else {
      bp_df <- tibble()
    }
  } else {
    bp_df <- tibble()
  }
  
  # --- 2) Enhanced functional classifier
  classify_func <- function(g, p, fully) {
    if (!isTRUE(fully)) return("partially genic")
    if (is.na(p) && is.na(g)) return(NA_character_)
    
    # Convert to lowercase for matching
    p_lower <- tolower(as.character(p))
    g_lower <- tolower(as.character(g))
    
    # Transposable elements (highest priority)
    if (grepl("transpos|\\bis\\d+|\\btn\\d+|insertion sequence", p_lower)) return("Transposase")
    
    # Phage-related
    if (grepl("phage|prophage|capsid|tail|terminase|integrase|portal", p_lower)) return("Phage-associated")
    
    # Hypothetical/Unknown
    if (grepl("hypothetical|uncharacterized|unknown function", p_lower)) return("Hypothetical protein")
    if (grepl("\\bputative\\b|probable|predicted", p_lower)) return("Putative protein")
    
    # RNA genes
    if (grepl("ribosomal\\s+rna|\\brrna\\b|23s|16s|5s", p_lower)) return("rRNA")
    if (grepl("\\btrna\\b|transfer rna", p_lower)) return("tRNA")
    
    # Specific systems
    if (grepl("esx|esat-6|type vii secretion", p_lower)) return("ESX secretion family")
    if (grepl("toxin|antitoxin|\\bta\\b", p_lower)) return("Toxin-antitoxin system")
    if (grepl("restriction|modification|methylase", p_lower)) return("Restriction-modification")
    
    # Transporters
    if (grepl("abc transporter|atp-binding cassette", p_lower)) return("ABC transporter")
    if (grepl("transporter|permease|efflux|uptake|channel", p_lower)) return("Transporter")
    
    # Enzymes (more specific)
    if (grepl("polymerase|helicase|recombinase|topoisomerase", p_lower)) return("DNA metabolism")
    if (grepl("dehydrogenase|reductase|oxidase|catalase|peroxidase", p_lower)) return("Oxidoreductase")
    if (grepl("synthase|synthetase|ligase", p_lower)) return("Synthase/Ligase")
    if (grepl("hydrolase|peptidase|protease|nuclease", p_lower)) return("Hydrolase")
    if (grepl("kinase|phosphatase", p_lower)) return("Signaling")
    if (grepl("transferase|aminotransferase", p_lower)) return("Transferase")
    
    # Structural/Functional
    if (grepl("ribosomal protein|\\brps|\\brpl", p_lower)) return("Ribosomal protein")
    if (grepl("chaperone|heat shock|\\bhsp|\\bdnak", p_lower)) return("Chaperone")
    if (grepl("flagell|pili|fimbr|adhesin", p_lower)) return("Surface structure")
    
    # Regulatory
    if (grepl("regulator|repressor|activator|\\bfur\\b|\\blexr", p_lower)) return("Regulatory")
    
    # Gene-specific classifications
    if (!is.na(g)) {
      if (g_lower %in% c("rho", "nusa", "nusg")) return("Transcription")
      if (g_lower %in% c("rpfa", "rpfb", "rpfc", "rpfd", "rpfe")) return("Resuscitation")
      if (grepl("^pe\\d+|^ppe\\d+", g_lower)) return("PE/PPE family")
    }
    
    # Catch-all for other enzymes
    if (grepl("\\base$|\\base\\b", p_lower)) return("Enzyme")
    
    "Other"
  }
  
  # Combine span and breakpoint results
  both_df <- bind_rows(span_df, bp_df)
  
  if (nrow(both_df) > 0) {
    # Variant-aware prioritization
    inv_mask <- !is.na(both_df$variant) & 
      tolower(both_df$variant) %in% tolower(inv_like)
    
    both_df$hit_priority <- case_when(
      inv_mask & both_df$hit_type == "breakpoint" ~ 3L,
      both_df$hit_type == "span" & both_df$fully_genic ~ 2L,
      TRUE ~ 1L
    )
    
    # Add functional classification
    both_df$func_group <- mapply(classify_func, 
                                 both_df$gene, 
                                 both_df$product, 
                                 both_df$fully_genic, 
                                 USE.NAMES = FALSE)
  }
  
  # --- 3) Add intergenic rows for SVs with no hits
  if (nrow(both_df) > 0) {
    sv_key <- paste(as.character(seqnames(SVg)), start(SVg), end(SVg), sep = "_")
    hit_key <- paste(both_df$sv_accession, both_df$sv_start, both_df$sv_end, sep = "_")
    nohit_idx <- which(!sv_key %in% hit_key)
  } else {
    nohit_idx <- seq_along(SVg)
  }
  
  if (length(nohit_idx) > 0) {
    intergenic_df <- tibble(
      hit_type         = "none",
      gene_accession   = NA_character_,
      gene             = NA_character_,
      product          = NA_character_,
      sv_accession     = as.character(seqnames(SVg))[nohit_idx],
      species          = sv_species[nohit_idx],
      variant          = sv_variant[nohit_idx],
      sv_start         = start(SVg)[nohit_idx],
      sv_end           = end(SVg)[nohit_idx],
      sv_width         = sv_w[nohit_idx],
      refmid           = sv_refmid[nohit_idx],
      overlap_bp       = 0L,
      pct_overlap_sv   = 0,
      pct_overlap_gene = 0,
      jaccard          = 0,
      fully_genic      = FALSE,
      hit_priority     = 0L,
      func_group       = "intergenic"
    )
    
    out <- bind_rows(both_df, intergenic_df)
  } else {
    out <- both_df
  }
  
  # --- 4) Keep best or all hits
  if (keep == "best" && nrow(out) > 0) {
    out <- out %>%
      arrange(desc(hit_priority), desc(overlap_bp)) %>%
      group_by(sv_accession, sv_start, sv_end) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      arrange(sv_accession, sv_start)
  } else if (nrow(out) > 0) {
    out <- arrange(out, sv_accession, desc(hit_priority), desc(overlap_bp))
  }
  
  if (verbose && nrow(out) > 0) {
    cat("\nAnnotation summary:\n")
    print(table(out$hit_type))
    cat("\nFunctional groups:\n")
    print(table(out$func_group))
  }
  
  return(out)
}
