# 00-utils.R
# Core utilities, cache helpers, small predicates. Deduplicated during merge.

# --- ensure_dir ---
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# --- rds_exists ---
rds_exists <- function(path) file.exists(path) && !isTRUE(file.info(path)$isdir)

# --- sanitize_species ---
sanitize_species <- function(species_name) sub("_", " ", species_name)

# --- matrix_density ---
matrix_density <- function(M) Matrix::nnzero(M) / (nrow(M) * ncol(M))

# --- first_existing ---
first_existing <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm) == 0L) return(NA_character_)
  nm[1]
}

# --- first_scalar_chr ---
first_scalar_chr <- function(x) {
  if (length(x) == 0L) return(NA_character_)
  x <- x[[1]]
  if (length(x) == 0L) NA_character_ else as.character(x)[1]
}

# --- id_from_file ---
id_from_file <- function(f) sub("\\.(fna|fa|fasta|txt)(\\.gz)?$", "", basename(f), ignore.case = TRUE)

# --- svcm_ensure_dir ---
svcm_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# --- svcm_rds_exists ---
svcm_rds_exists <- function(path) file.exists(path) && !isTRUE(file.info(path)$isdir)

# --- svcm_load_cache ---
svcm_load_cache <- function(path) if (svcm_rds_exists(path)) readRDS(path) else NULL

# --- svcm_save_cache ---
svcm_save_cache <- function(obj, path) {
  svcm_ensure_dir(dirname(path))
  saveRDS(obj, path)
  invisible(path)
}

# --- svcm_stop_if_missing_tools ---
svcm_stop_if_missing_tools <- function(tools) {
  miss <- tools[Sys.which(tools) == ""]
  if (length(miss)) stop("Missing required tools: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

# --- svcm_stop_if_missing ---
svcm_stop_if_missing <- function(x, what = "object") {
  if (is.null(x) || (is.data.frame(x) && nrow(x) == 0L))
    stop("Missing required ", what, call. = FALSE)
  invisible(TRUE)
}

# --- additional cache/schema helpers (from SVCM6.R) ---

# --- load_cache ---
load_cache <- function(path) {
  if (file.exists(path)) readRDS(path) else NULL
}

# --- save_cache ---
save_cache <- function(obj, path, label = basename(path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path)
  fs <- tryCatch(file.size(path), error = function(e) NA_real_)
  msg <- if (is.na(fs)) "" else sprintf(" (%s bytes)", format(fs, big.mark = ","))
  message(sprintf("  [cached] %s%s", label, msg))
  invisible(path)
}

# --- stop_if_missing_tools ---
stop_if_missing_tools <- function(tools) {
  missing <- tools[Sys.which(tools) == ""]
  if (length(missing) > 0) {
    stop("Missing required command-line tools: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# --- sv_standardize_schema ---
sv_standardize_schema <- function(sv_tbl, qid, rid_expected = NULL) {
  if (is.null(sv_tbl)) return(NULL)
  if (!is.data.frame(sv_tbl)) stop("SV output is not a data.frame for qid=", qid)
  
  # ensure essential columns exist
  must <- c("rid", "start", "end", "width", "variant")
  missing <- setdiff(must, names(sv_tbl))
  if (length(missing) > 0) {
    stop("SV table missing columns for qid=", qid, ": ", paste(missing, collapse = ", "))
  }
  
  sv_tbl$qid <- qid
  sv_tbl$rid <- sub(" .*", "", as.character(sv_tbl$rid))
  
  if (!"variant_specific" %in% names(sv_tbl)) sv_tbl$variant_specific <- NA_character_
  
  # enforce rid invariant if asked
  if (!is.null(rid_expected) && nrow(sv_tbl) > 0) {
    if (any(!is.na(sv_tbl$rid)) && any(sv_tbl$rid != rid_expected, na.rm = TRUE)) {
      stop("rid mismatch for qid=", qid, " expected=", rid_expected,
           " observed=", paste(unique(sv_tbl$rid), collapse = ","))
    }
  }
  
  sv_tbl
}
