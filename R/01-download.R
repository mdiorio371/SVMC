# 01-download.R
# Rebuilt cleanly during merge.

# --- svcm_sync_inventory  [from SVCM.R] ---
svcm_sync_inventory <- function(
    sync_dir,
    include = c("", "1v_the_rest"),
    exclude = c("mash_id", "ss20"),
    recursive = FALSE
) {
  sync_dir <- normalizePath(sync_dir, winslash = "/", mustWork = TRUE)
  
  include <- unique(as.character(include))
  exclude <- unique(as.character(exclude))
  
  # Directories to scan
  scan_dirs <- unique(file.path(sync_dir, include))
  scan_dirs[include == ""] <- sync_dir
  scan_dirs <- scan_dirs[file.exists(scan_dirs)]
  
  if (!length(scan_dirs)) {
    return(tibble::tibble(qid = character(), bucket = character(), path = character()))
  }
  
  # Filter out excluded buckets by basename
  if (length(exclude)) {
    scan_dirs <- scan_dirs[!(basename(scan_dirs) %in% exclude)]
  }
  if (!length(scan_dirs)) {
    return(tibble::tibble(qid = character(), bucket = character(), path = character()))
  }
  
  inv_list <- lapply(scan_dirs, function(d) {
    ff <- list.files(
      d,
      pattern = "\\.txt(\\.gz)?$",
      full.names = TRUE,
      recursive = recursive
    )
    if (!length(ff)) return(NULL)
    
    # bucket = relative path component immediately under sync_dir
    rel_dir <- sub(paste0("^", sync_dir, "/?"), "", dirname(ff))
    bucket <- ifelse(rel_dir == "" | is.na(rel_dir), "", sub("/.*$", "", rel_dir))
    
    qid <- basename(ff)
    qid <- sub("\\.txt(\\.gz)?$", "", qid, ignore.case = TRUE)
    
    tibble::tibble(
      qid = as.character(qid),
      bucket = as.character(bucket),
      path = normalizePath(ff, winslash = "/", mustWork = TRUE)
    )
  })
  
  inv <- dplyr::bind_rows(inv_list)
  inv <- inv %>% dplyr::distinct(.data$qid, .keep_all = TRUE)
  inv
}

# --- svcm_list_txt_like  [from SVCM.R] ---
svcm_list_txt_like <- function(dir, recursive = FALSE) {
  if (!dir.exists(dir)) return(character())
  list.files(dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE, recursive = recursive, ignore.case = TRUE)
}

# --- svcm_token_from_fasta_header  [from SVCM.R] ---
svcm_token_from_fasta_header <- function(h) sub(" .*", "", as.character(h))

# --- svcm_materialize_plain  [from SVCM.R] ---
svcm_materialize_plain <- function(src_path, tmp_dir, keep = FALSE) {
  stopifnot(file.exists(src_path))
  svcm_ensure_dir(tmp_dir)
  
  if (grepl("\\.txt$", src_path, ignore.case = TRUE) && !grepl("\\.gz$", src_path, ignore.case = TRUE)) {
    # already plain
    return(normalizePath(src_path, winslash = "/", mustWork = TRUE))
  }
  
  svcm_stop_if_missing_tools(c("gunzip"))
  out_path <- file.path(tmp_dir, paste0(svcm_id_from_file(src_path), ".txt"))
  # overwrite existing
  if (file.exists(out_path)) file.remove(out_path)
  code <- suppressWarnings(system2("gunzip", args = c("-c", shQuote(src_path)), stdout = out_path, stderr = TRUE))
  if (!file.exists(out_path) || file.size(out_path) == 0) {
    stop("Failed to materialize plain fasta from: ", src_path, "\n", paste(code, collapse = "\n"))
  }
  if (!isTRUE(keep)) {
    # caller will clean up
  }
  normalizePath(out_path, winslash = "/", mustWork = TRUE)
}

# --- materialize_ref_plain  [from SVCM6.R] ---
materialize_ref_plain <- function(sync_dir, ref_genome) {
  ref_gz   <- file.path(sync_dir, sprintf("%s.txt.gz", ref_genome))
  ref_path <- file.path(sync_dir, sprintf("%s.txt",    ref_genome))
  if (!file.exists(ref_path)) {
    if (!file.exists(ref_gz)) stop("Reference .txt.gz not found: ", ref_gz)
    # Keep a persistent .txt for reference
    suppressWarnings(system2("gunzip", args = c("-k", "-f", shQuote(ref_gz))))
  }
  if (!file.exists(ref_path) || file.size(ref_path) == 0) stop("Reference .txt not materialized: ", ref_path)
  ref_path
}

# --- svcm_ensure_ref_plain  [from SVCM6.R] ---
svcm_ensure_ref_plain <- function(sync_dir, ref_genome,
                                  include = c(""), exclude = character(),
                                  out_subdir = "_ref_plain",
                                  keep = c("largest")) {
  stopifnot(length(ref_genome) == 1L, nzchar(ref_genome))
  
  # Find candidate file in sync_dir or included buckets
  buckets <- unique(include)
  buckets <- buckets[!buckets %in% exclude]
  
  candidates <- unlist(lapply(buckets, function(b) {
    d <- file.path(sync_dir, b)
    c(file.path(d, paste0(ref_genome, ".txt")),
      file.path(d, paste0(ref_genome, ".txt.gz")))
  }), use.names = FALSE)
  
  in_path <- candidates[file.exists(candidates)][1]
  if (is.na(in_path)) {
    stop("svcm_ensure_ref_plain: cannot find ref genome file for ", ref_genome,
         " under sync_dir (include=", paste(include, collapse=","), ")")
  }
  
  # Read (gz ok)
  ref_seqs <- Biostrings::readDNAStringSet(in_path)
  if (length(ref_seqs) == 0L) stop("svcm_ensure_ref_plain: reference FASTA is empty: ", in_path)
  
  # Keep largest contig only (recommended)
  if ("largest" %in% keep && length(ref_seqs) > 1L) {
    ref_seqs <- ref_seqs[order(Biostrings::width(ref_seqs), decreasing = TRUE)][1]
  }
  names(ref_seqs) <- sub(" .*", "", names(ref_seqs))
  # Write a persistent single-record plain FASTA
  out_dir  <- file.path(sync_dir, out_subdir)
  ensure_dir(out_dir)
  out_path <- file.path(out_dir, paste0(ref_genome, ".txt"))
  
  Biostrings::writeXStringSet(ref_seqs, out_path)
  out_path
}

# --- gunzip_to_file  [from SVCM6.R] ---
gunzip_to_file <- function(src_gz, out_path) {
  stop_if_missing_tools(c("gunzip"))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # overwrite target if exists
  if (file.exists(out_path)) file.remove(out_path)
  
  # use gunzip -c and redirect stdout to file
  # system2(stdout=filename) is reliable; capture stderr to temp for debugging
  errf <- tempfile("gunzip_stderr_", fileext = ".txt")
  on.exit(unlink(errf), add = TRUE)
  
  code <- suppressWarnings(system2("gunzip", args = c("-c", shQuote(src_gz)),
                                   stdout = out_path, stderr = errf))
  if (!identical(code, 0L)) {
    err <- paste(readLines(errf, warn = FALSE), collapse = "\n")
    stop("gunzip failed for: ", src_gz, "\nExit code: ", code, "\n", err)
  }
  if (!file.exists(out_path) || file.size(out_path) == 0) {
    err <- paste(readLines(errf, warn = FALSE), collapse = "\n")
    stop("Failed to materialize plain file: ", out_path, "\n", err)
  }
  invisible(out_path)
}

# --- svcm_file_looks_fasta  [from SVCM.R] ---
svcm_file_looks_fasta <- function(path, n_check = 5L) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  x <- tryCatch(readLines(con, n = n_check, warn = FALSE), error = function(e) character())
  x <- trimws(x); x <- x[nzchar(x)]
  if (!length(x)) return(FALSE)
  startsWith(x[1], ">")
}

# --- svcm_id_from_file  [from SVCM.R] ---
svcm_id_from_file <- function(x) sub("\\.txt(\\.gz)?$", "", basename(x), ignore.case = TRUE)

# --- step1_sync  [from SVCM6.R] ---
step1_sync <- function(sync_dir, r_vars_dir, assembly_paths_rds = NULL, n_cores = 8L) {
  inv <- inventory_synced(sync_dir)
  pos_path <- file.path(r_vars_dir, "position_table.rds")
  
  cached_pos <- load_cache(pos_path)
  
  # If synced gz exist, skip sync regardless of assembly_paths presence.
  if (length(inv$sync_files) > 0) {
    if (!is.null(cached_pos)) {
      message(sprintf("Step 1.1 SKIP: %d synced genomes on disk; using cached position_table.rds", length(inv$sync_files)))
      return(list(position_table = cached_pos, sync_accs = inv$sync_accs, sync_files = inv$sync_files))
    }
    # create a minimal placeholder position table (so later steps can proceed)
    message(sprintf("Step 1.1 SKIP: %d synced genomes on disk; position_table.rds missing -> creating minimal placeholder", length(inv$sync_files)))
    position_table <- tibble::tibble(out_file = basename(inv$sync_files),
                                     len = NA_real_,
                                     note = NA_character_)
    save_cache(position_table, pos_path)
    return(list(position_table = position_table, sync_accs = inv$sync_accs, sync_files = inv$sync_files))
  }
  
  # Otherwise, need assembly paths + load_arrange()
  if (is.null(assembly_paths_rds)) assembly_paths_rds <- file.path(r_vars_dir, "assembly_paths.rds")
  if (!file.exists(assembly_paths_rds)) stop("No synced genomes found and missing assembly_paths.rds at: ", assembly_paths_rds)
  if (!exists("load_arrange", mode = "function")) stop("load_arrange() not found. Did you source('SVCM.R')?")
  
  assembly_paths <- readRDS(assembly_paths_rds)
  message(sprintf("Step 1.1: syncing %d assemblies (cores=%d) ...", length(assembly_paths), n_cores))
  
  if (!requireNamespace("pbmcapply", quietly = TRUE)) stop("Missing package: pbmcapply")
  sync_results <- pbmcapply::pbmclapply(assembly_paths, function(ap) {
    load_arrange(assembly_path = ap, sync_dir_the_rest = sync_dir)
  }, mc.cores = n_cores)
  
  position_table <- dplyr::bind_rows(sync_results)
  save_cache(position_table, pos_path)
  
  inv2 <- inventory_synced(sync_dir)
  n_ok <- sum(is.na(position_table$note))
  message(sprintf("Step 1.1 DONE: Synced %d / %d | On disk now: %d .txt.gz",
                  n_ok, nrow(position_table), length(inv2$sync_files)))
  
  list(position_table = position_table, sync_accs = inv2$sync_accs, sync_files = inv2$sync_files)
}

# --- inventory_synced  [from SVCM6.R] ---
inventory_synced <- function(sync_dir,
                             include = c(""),
                             exclude = character()) {
  buckets <- unique(include)
  buckets <- buckets[!buckets %in% exclude]
  
  files <- unlist(lapply(buckets, function(b) {
    d <- file.path(sync_dir, b)
    if (!dir.exists(d)) return(character())
    list.files(d, pattern = "\\.txt(\\.gz)?$", full.names = TRUE, recursive = FALSE)
  }), use.names = FALSE)
  
  accs <- sub("\\.txt(\\.gz)?$", "", basename(files), ignore.case = TRUE)
  
  list(
    sync_files = files,
    sync_accs  = accs,
    sync_txt_files = files[grepl("\\.txt$", files, ignore.case = TRUE)],
    sync_gz_files  = files[grepl("\\.txt\\.gz$", files, ignore.case = TRUE)]
  )
}
