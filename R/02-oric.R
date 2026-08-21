# 02-oric.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- dnaA_from_gff  [from SVCM.R:100] ---
dnaA_from_gff <- function(gff_path) {
  gff_gzcon <- gzcon(url(gff_path))
  gff_table <- ape::read.gff(gff_gzcon)
  close(gff_gzcon)

  product_string   <- "product="
  dnaAstring       <- "DnaA"
  replication_string <- " replication initiator "

  dnaA_strings <- c(
    sprintf("%schromosome%s%s", product_string, replication_string, dnaAstring),
    sprintf("%schromosome%sprotein %s", product_string, replication_string, dnaAstring)
  )

  dnaA_vect <- gff_table %>%
    as_tibble() %>%
    filter(
      (grepl("gene=dnaA", attributes, ignore.case = TRUE) |
         grepl(dnaA_strings[1], attributes, ignore.case = TRUE) |
         grepl(dnaA_strings[2], attributes, ignore.case = TRUE) |
         grepl("product=DnaA", attributes, ignore.case = TRUE)) &
        type == "CDS"
    ) %>%
    as.data.frame()

  if (nrow(dnaA_vect) > 1L) {
    if (any(grepl("gene=DnaA", dnaA_vect$attributes, ignore.case = TRUE))) {
      dnaA_vect <- dnaA_vect[grepl("gene=DnaA", dnaA_vect$attributes, ignore.case = TRUE), ]
    }
    if (all(pull(dnaA_vect, strand) == "+")) {
      dnaA_vect <- dnaA_vect %>% filter(start == min(start)) %>% dplyr::slice(1)
    } else if (all(pull(dnaA_vect, strand) == "-")) {
      dnaA_vect <- dnaA_vect %>% filter(start == max(start)) %>% dplyr::slice(1)
    } else {
      dnaA_vect <- dnaA_vect %>% dplyr::slice(1)
    }
  }

  if (nrow(dnaA_vect) == 0L) return("no dnaA annotation")
  dnaA_vect
}

# --- dnaA_sync  [from SVCM.R:148] ---
dnaA_sync <- function(fasta_path, dnaA_table, out_dir) {
  genome_seq <- Biostrings::readDNAStringSet(fasta_path)
  genome_seq <- genome_seq[order(Biostrings::width(genome_seq), decreasing = TRUE), ][1]
  genome_length <- Biostrings::width(genome_seq)
  genome_acc <- (strsplit(names(genome_seq), " "))[[1]][1]

  gff_ori_strand <- pull(dnaA_table, strand)

  if (gff_ori_strand == "+") {
    ori_start <- pull(dnaA_table, start)
    out_seq <- Biostrings::xscat(
      Biostrings::subseq(genome_seq, start = ori_start, end = genome_length),
      Biostrings::subseq(genome_seq, start = 1, end = (ori_start - 1))
    ) %>% `names<-`(genome_acc)
  } else {
    gene_end <- pull(dnaA_table, end)
    if (gene_end > Biostrings::width(genome_seq)) gene_end <- Biostrings::width(genome_seq)
    ori_start <- genome_length - gene_end + 1
    out_seq <- Biostrings::xscat(
      Biostrings::subseq(Biostrings::reverseComplement(genome_seq), start = ori_start, end = genome_length),
      Biostrings::subseq(Biostrings::reverseComplement(genome_seq), start = 1, end = (ori_start - 1))
    ) %>% `names<-`(genome_acc)
  }

  ensure_dir(out_dir)
  out_file <- sprintf("%s/%s.txt", out_dir, genome_acc)
  Biostrings::writeXStringSet(out_seq, filepath = out_file)

  tibble(
    asm_name   = basename(dirname(fasta_path)),
    accession  = genome_acc,
    ori_pos    = dnaA_table$start[1],
    ori_strand = dnaA_table$strand[1],
    len        = genome_length,
    out_file   = out_file,
    note       = NA
  )
}

# --- load_arrange  [from SVCM.R:195] ---
load_arrange <- function(assembly_path, sync_dir_the_rest) {
  fasta_path <- sprintf("%s/%s_genomic.fna.gz",
                        assembly_path, gsub(".*/", "", assembly_path))
  gff_path   <- sprintf("%s/%s_genomic.gff.gz",
                        assembly_path, gsub(".*/", "", assembly_path))

  fna_exists <- file_exists_at_url(fasta_path)
  gff_exists <- file_exists_at_url(gff_path)

  dnaA_table <- dnaA_from_gff(gff_path)

  if (!fna_exists | !gff_exists) {
    position_tib <- tibble(
      asm_name = basename(assembly_path),
      note     = "missing fasta/gff"
    )
  } else {
    position_tib <- dnaA_sync(fasta_path, dnaA_table, sync_dir_the_rest)
  }
  return(position_tib)
}

# --- arrange_align  [from SVCM.R:229] ---
arrange_align <- function(
    assembly_path,
    ref_path,
    alignment_dir_the_rest,
    sync_dir_the_rest,
    call_SVs     = FALSE,
    species_name = NULL
) {
  fasta_path <- sprintf("%s/%s_genomic.fna.gz",
                        assembly_path, gsub(".*/", "", assembly_path))
  gff_path   <- sprintf("%s/%s_genomic.gff.gz",
                        assembly_path, gsub(".*/", "", assembly_path))

  fna_exists <- file_exists_at_url(fasta_path)
  gff_exists <- file_exists_at_url(gff_path)

  if (!fna_exists | !gff_exists) {
    position_tib <- tibble(
      asm_name = basename(assembly_path),
      note     = "missing fasta/gff",
      out_file = NA_character_
    )
    return(list(position_tib, NULL))
  }

  dnaA_table   <- dnaA_from_gff(gff_path)
  position_tib <- dnaA_sync(fasta_path, dnaA_table, sync_dir_the_rest)

  if (!isTRUE(call_SVs)) return(list(position_tib, NULL))

  if (is.null(species_name)) stop("call_SVs = TRUE requires species_name.")

  out_str <- sprintf(
    "%s_v_%s",
    sub("\\.txt$", "", basename(ref_path)),
    sub("\\.txt$", "", basename(position_tib$out_file))
  )

  # Align to the reference
  simple_nucmer(
    ref           = ref_path,
    qry           = position_tib$out_file,
    output        = out_str,
    alignment_dir = alignment_dir_the_rest,
    check_nuc     = FALSE
  ) %>% system()

  delta_filt_file  <- sprintf("%s/%s_filtered.delta", alignment_dir_the_rest, out_str)
  unfiltered_delta <- sprintf("%s/unfiltered/%s.delta", alignment_dir_the_rest, out_str)

  delta_SVs <- delta_all_SVs(
    delta_file            = delta_filt_file,
    unfiltered_delta_file = unfiltered_delta,
    species_name          = species_name
  )

  list(position_tib, delta_SVs)
}
