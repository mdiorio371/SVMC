# 03-identity.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- get_dists  [from SVCM.R:293] ---
get_dists <- function(distance_matrix) {
  diag(distance_matrix) <- NA
  sequence_averages <- rowMeans(distance_matrix, na.rm = TRUE)
  overall_mean_distance <- mean(as.vector(distance_matrix), na.rm = TRUE)
  difference_from_mean <- sequence_averages - overall_mean_distance

  tibble(
    Sequence = names(sequence_averages),
    Sequence_mean_distance = sequence_averages,
    Mean_distance = overall_mean_distance,
    Difference_from_Mean = difference_from_mean
  )
}

# --- get_pairwise_identities_p2  [from SVCM.R:425] ---
get_pairwise_identities_p2 <- function(sync_dir, id_dir, species_name,
                                    clust_method = "complete",
                                    files = NULL,
                                    threads = 8L,
                                    sketch_size = 10000L) {
  
  mash_file <- file.path(id_dir, sprintf("%s_mash.txt", species_name))
  
  if (!file.exists(mash_file)) {
    mash_bash <- mash_commands(species_name, id_dir, sync_dir,
                               files = files, threads = threads,
                               sketch_size = sketch_size)
    cat(mash_bash[1], "\n", mash_bash[2], "\n")
    system(mash_bash[1])  # do NOT ignore stderr
    system(mash_bash[2])
  }
  
  id_table <- data.table::fread(mash_file, select = 1:3) %>%
    mutate(
      V1 = sub("\\.txt$", "", basename(V1)),
      V2 = sub("\\.txt$", "", basename(V2)),
      identity = (1 - as.numeric(V3)) * 100
    ) %>%
    dplyr::select(ref = V1, qry = V2, identity) %>%
    distinct(ref, qry, .keep_all = TRUE) %>%   # safety
    as_tibble()
  
  id_wide <- id_table %>%
    tidyr::pivot_wider(names_from = qry, values_from = identity) %>%
    tibble::column_to_rownames("ref") %>%
    as.matrix()
  
  # if anything is missing, fill from transpose (should usually be complete already)
  if (anyNA(id_wide)) {
    id_wide2 <- t(id_wide)
    id_wide[is.na(id_wide)] <- id_wide2[is.na(id_wide)]
  }
  
  id_wide_dist <- 100 - id_wide
  id_dist_lt <- as.dist(id_wide_dist)
  
  list(
    mash_tib = id_table,
    mash_D = id_dist_lt
  )
}

# --- mash_seed_medoid  [from SVCM6.R:149] ---
mash_seed_medoid <- function(sync_dir,
                             id_dir_seed,
                             seed_accs,
                             n_threads = 8L,
                             mash_k = 21L,
                             mash_s = 1000L) {
  stop_if_missing_tools(c("mash", "gunzip"))
  
  mash_in_dir <- file.path(id_dir_seed, "mash_in")
  dir.create(mash_in_dir, recursive = TRUE, showWarnings = FALSE)
  
  # materialize seed FASTA files into mash_in_dir with controlled basenames
  # output extension doesn't matter; mash reads FASTA.
  seed_out <- file.path(mash_in_dir, sprintf("%s.fa", seed_accs))
  
  for (i in seq_along(seed_accs)) {
    acc <- seed_accs[i]
    src <- file.path(sync_dir, sprintf("%s.txt.gz", acc))
    if (!file.exists(src)) stop("Missing synced genome: ", src)
    gunzip_to_file(src_gz = src, out_path = seed_out[i])
  }
  
  # write input list for mash
  list_file <- file.path(id_dir_seed, "mash_inputs.txt")
  writeLines(seed_out, con = list_file)
  
  prefix <- file.path(id_dir_seed, "seed")
  msh    <- sprintf("%s.msh", prefix)
  distf  <- file.path(id_dir_seed, "seed_mash_dist.txt")
  
  # sketch
  sketch_cmd <- c("sketch",
                  "-p", as.integer(n_threads),
                  "-k", as.integer(mash_k),
                  "-s", as.integer(mash_s),
                  "-o", shQuote(prefix),
                  "-l", shQuote(list_file))
  code1 <- suppressWarnings(system2("mash", args = sketch_cmd, stdout = TRUE, stderr = TRUE))
  if (!file.exists(msh) || file.size(msh) == 0) {
    stop("mash sketch failed (no .msh produced). Output:\n", paste(code1, collapse = "\n"))
  }
  
  # dist (all-vs-all within same sketch)
  code2 <- suppressWarnings(system2("mash",
                                    args = c("dist", shQuote(msh), shQuote(msh)),
                                    stdout = distf, stderr = TRUE))
  if (!file.exists(distf) || file.size(distf) == 0) {
    stop("mash dist failed (empty dist file).")
  }
  
  # parse dist: ref  query  dist  p  shared
  dt <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(distf, header = FALSE, select = 1:3, showProgress = FALSE)
    if (ncol(dt) < 3) stop("mash dist parse failed; dist file malformed: ", distf)
    names(dt) <- c("ref", "qry", "dist")
  } else {
    x <- read.table(distf, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) < 3) stop("mash dist parse failed; dist file malformed: ", distf)
    dt <- data.frame(ref = x[[1]], qry = x[[2]], dist = x[[3]], stringsAsFactors = FALSE)
  }
  
  # normalize IDs: strip paths + strip extensions
  norm_id <- function(z) {
    b <- basename(z)
    b <- sub("\\.fa(sta)?$|\\.fna$|\\.txt(\\.gz)?$", "", b, ignore.case = TRUE)
    b
  }
  dt$ref_id <- norm_id(dt$ref)
  dt$qry_id <- norm_id(dt$qry)
  
  # compute mean distance per genome without building full matrix
  ids <- sort(unique(c(dt$ref_id, dt$qry_id)))
  idx <- setNames(seq_along(ids), ids)
  
  iref <- idx[dt$ref_id]
  iqry <- idx[dt$qry_id]
  d    <- as.numeric(dt$dist)
  
  ok <- !is.na(iref) & !is.na(iqry) & (iref != iqry)
  iref <- iref[ok]; iqry <- iqry[ok]; d <- d[ok]
  
  sums   <- numeric(length(ids))
  counts <- numeric(length(ids))
  
  # add to both endpoints
  for (k in seq_along(d)) {
    a <- iref[k]; b <- iqry[k]; val <- d[k]
    sums[a]   <- sums[a] + val
    sums[b]   <- sums[b] + val
    counts[a] <- counts[a] + 1
    counts[b] <- counts[b] + 1
  }
  
  mean_dist <- sums / pmax(counts, 1)
  out <- tibble::tibble(
    id        = ids,
    mean_dist = mean_dist,
    n_pairs   = counts
  ) |>
    dplyr::arrange(mean_dist)
  
  ref <- out$id[1]
  
  list(
    ref_genome = ref,
    mean_dists = out,
    mash_in_dir = mash_in_dir,
    mash_sketch = msh,
    mash_dist_file = distf,
    seed_accs = seed_accs
  )
}

# --- step1_seed_medoid  [from SVCM6.R:588] ---
step1_seed_medoid <- function(sync_dir,
                              id_dir,
                              r_vars_dir,
                              species_name,
                              n_seed_subset = 200L,
                              n_threads = 8L,
                              force = FALSE) {
  ref_path_rds <- file.path(r_vars_dir, "ref_genome.rds")
  cached_ref <- load_cache(ref_path_rds)
  
  inv <- inventory_synced(sync_dir)
  if (length(inv$sync_accs) == 0) stop("No synced genomes found in: ", sync_dir)
  
  if (!force && !is.null(cached_ref) && cached_ref %in% inv$sync_accs) {
    message(sprintf("Step 1.2 SKIP: reusing cached reference %s", cached_ref))
    ref_genome <- cached_ref
  } else {
    set.seed(42)
    seed_accs <- if (length(inv$sync_accs) > n_seed_subset) sample(inv$sync_accs, n_seed_subset) else inv$sync_accs
    
    id_dir_seed <- file.path(id_dir, "seed_medoid")
    dir.create(id_dir_seed, recursive = TRUE, showWarnings = FALSE)
    
    message(sprintf("Step 1.2: seed-medoid MASH on %d genomes into: %s", length(seed_accs), id_dir_seed))
    mash_res <- mash_seed_medoid(sync_dir = sync_dir,
                                 id_dir_seed = id_dir_seed,
                                 seed_accs = seed_accs,
                                 n_threads = n_threads)
    
    ref_genome <- mash_res$ref_genome
    if (!ref_genome %in% inv$sync_accs) {
      stop("Seed medoid produced ref_genome not present in sync_dir inventory: ", ref_genome)
    }
    
    save_cache(ref_genome, ref_path_rds)
    save_cache(mash_res, file.path(r_vars_dir, "mash_seed_result.rds"))
    message(sprintf("Step 1.2 DONE: selected reference %s", ref_genome))
  }
  
  # materialize persistent ref .txt for nucmer
  ref_plain <- svcm_ensure_ref_plain(
    sync_dir   = sync_dir,
    ref_genome = ref_genome,
    include    = c(""),      # or your include buckets
    exclude    = character()
  )
  ref_rid <- sub(" .*", "", names(Biostrings::readDNAStringSet(ref_plain))[1])
  save_cache(ref_rid, file.path(r_vars_dir, "ref_rid.rds"), label = "ref_rid.rds")
  list(ref_genome = ref_genome, ref_path = ref_plain)
}

# --- svcm_run_mash_dist_to_file  [from SVCM.R:3266] ---
svcm_run_mash_dist_to_file <- function(sketch_msh, out_txt, threads = 8L) {
  mash <- Sys.which("mash")
  if (mash == "") stop("mash not found on PATH")
  
  sketch_msh <- normalizePath(sketch_msh, mustWork = TRUE)
  out_txt    <- normalizePath(out_txt,    mustWork = FALSE)
  dir.create(dirname(out_txt), recursive = TRUE, showWarnings = FALSE)
  
  err_txt <- sub("\\.txt$", ".stderr.txt", out_txt)
  
  args <- c("dist", "-p", as.integer(threads), "-t", sketch_msh, sketch_msh)
  status <- system2(mash, args, stdout = out_txt, stderr = err_txt)
  
  # salvage if Mash wrote table-ish output to stderr (rare, but you saw it)
  if ((!file.exists(out_txt) || is.na(file.info(out_txt)$size) || file.info(out_txt)$size == 0) &&
      file.exists(err_txt) && file.info(err_txt)$size > 0) {
    x <- readLines(err_txt, warn = FALSE)
    tab_ok <- any(vapply(strsplit(x, "\t", fixed = TRUE),
                         function(z) length(z) >= 3 && suppressWarnings(is.finite(as.numeric(z[3]))),
                         logical(1)))
    if (tab_ok) writeLines(x, out_txt)
  }
  
  if (!identical(status, 0L)) stop("mash dist failed (exit=", status, "). See: ", err_txt)
  if (!file.exists(out_txt) || is.na(file.info(out_txt)$size) || file.info(out_txt)$size == 0) {
    stop("mash dist produced empty output: ", out_txt, "\nSee: ", err_txt)
  }
  out_txt
}

# --- svcm_parse_mash_dist  [from SVCM.R:3296] ---
svcm_parse_mash_dist <- function(distf) {
  # Returns list: dt (long), D (square matrix), ids (vector)
  dt <- tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      x <- data.table::fread(distf, header = FALSE, select = 1:3, showProgress = FALSE)
      names(x) <- c("ref", "qry", "dist")
      as.data.frame(x)
    } else {
      x <- utils::read.table(distf, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      data.frame(ref = x[[1]], qry = x[[2]], dist = x[[3]], stringsAsFactors = FALSE)
    }
  }, error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) stop("mash dist parse failed/empty: ", distf)
  
  norm_id <- function(z) {
    b <- basename(z)
    sub("\\.fa(sta)?$|\\.fna$|\\.txt(\\.gz)?$", "", b, ignore.case = TRUE)
  }
  
  dt$ref_id <- norm_id(dt$ref)
  dt$qry_id <- norm_id(dt$qry)
  dt$dist   <- as.numeric(dt$dist)
  
  ids <- sort(unique(c(dt$ref_id, dt$qry_id)))
  idx <- setNames(seq_along(ids), ids)
  
  i <- idx[dt$ref_id]; j <- idx[dt$qry_id]
  ok <- !is.na(i) & !is.na(j) & is.finite(dt$dist)
  i <- i[ok]; j <- j[ok]; d <- dt$dist[ok]
  
  n <- length(ids)
  D <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  D[cbind(i,j)] <- d
  diag(D) <- 0
  # symmetrize if needed
  if (anyNA(D)) {
    Dt <- t(D)
    D[is.na(D)] <- Dt[is.na(D)]
  }
  
  list(dt = dt, D = D, ids = ids)
}

# --- svcm_step_mash_full_on_cohort  [from SVCM.R:3616] ---
svcm_step_mash_full_on_cohort <- function(
    sync_dir,
    id_dir,
    species_name,
    qids,
    include = c("", "1v_the_rest"),
    exclude = c("mash_id", "ss20"),
    threads = 8L,
    mash_k = 21L,
    mash_s = 10000L,
    make_plot = TRUE,
    plot_format = c("png","pdf")[1]
) {
  svcm_stop_if_missing_tools(c("mash", "gunzip"))
  if (!exists("svcm_sync_inventory", mode = "function")) stop("svcm_sync_inventory() not found.")
  
  sync_dir <- normalizePath(sync_dir, winslash="/", mustWork=TRUE)
  id_dir   <- normalizePath(id_dir,   winslash="/", mustWork=TRUE)
  
  mash_dir <- file.path(id_dir, "mash_full")
  mash_in  <- file.path(mash_dir, "mash_in")
  dir.create(mash_in, recursive = TRUE, showWarnings = FALSE)
  
  inv <- svcm_sync_inventory(sync_dir)
  inv <- inv %>% dplyr::mutate(path = normalizePath(.data$path, winslash="/", mustWork=TRUE))
  inv <- inv %>% dplyr::filter(.data$qid %in% qids)
  
  if (nrow(inv) == 0) stop("No qids matched inventory for mash_full.")
  
  # materialize .fa inputs
  out_files <- file.path(mash_in, paste0(inv$qid, ".fa"))
  
  for (i in seq_len(nrow(inv))) {
    src <- inv$path[i]
    dest <- out_files[i]
    if (grepl("\\.gz$", src, ignore.case = TRUE)) {
      if (file.exists(dest)) file.remove(dest)
      err <- tempfile(fileext = ".stderr.txt")
      status <- system2("gunzip", args = c("-c", src), stdout = dest, stderr = err)
      if (!identical(status, 0L) || !file.exists(dest) || is.na(file.size(dest)) || file.size(dest) == 0) {
        msg <- if (file.exists(err)) paste(readLines(err, warn=FALSE), collapse="\n") else ""
        stop("Failed to materialize mash input from gz: ", src, "\n", msg)
      }
    } else {
      ok <- file.copy(src, dest, overwrite = TRUE)
      if (!ok || !file.exists(dest) || is.na(file.size(dest)) || file.size(dest) == 0) {
        stop("Failed to copy mash input: ", src)
      }
    }
  }
  
  list_file <- file.path(mash_dir, "mash_inputs.txt")
  writeLines(normalizePath(out_files, mustWork = TRUE), con = list_file)
  
  prefix <- file.path(mash_dir, "mash_full")
  msh    <- paste0(prefix, ".msh")
  distf  <- file.path(mash_dir, "mash_full_dist.txt")
  err1   <- file.path(mash_dir, "mash_full_sketch.stderr.txt")
  
  # sketch (avoid shQuote in system2 args)
  sketch_args <- c(
    "sketch",
    "-p", as.integer(threads),
    "-k", as.integer(mash_k),
    "-s", as.integer(mash_s),
    "-o", prefix,
    "-l", list_file
  )
  status1 <- system2("mash", args = sketch_args, stdout = TRUE, stderr = err1)
  if (!file.exists(msh) || is.na(file.size(msh)) || file.size(msh) == 0) {
    msg <- if (file.exists(err1)) paste(readLines(err1, warn=FALSE), collapse="\n") else ""
    stop("mash sketch failed: ", msh, "\n", msg)
  }
  
  # dist
  svcm_run_mash_dist_to_file(msh, distf, threads = threads)
  
  parsed <- svcm_parse_mash_dist(distf)
  D <- parsed$D
  ids <- parsed$ids
  
  # mean distances (diversity quick read)
  mean_dist <- rowMeans(D, na.rm = TRUE)
  md_tbl <- tibble::tibble(qid = names(mean_dist), mean_dist = as.numeric(mean_dist)) %>%
    dplyr::arrange(.data$mean_dist)
  
  # optional heatmap artifact
  plot_path <- NA_character_
  if (isTRUE(make_plot) && requireNamespace("ggplot2", quietly = TRUE)) {
    ID <- (1 - D) * 100
    dlt <- stats::as.dist(100 - ID)
    hc <- stats::hclust(dlt, method = "ward.D")
    ord <- hc$labels[hc$order]
    
    long <- as.data.frame(as.table(ID), stringsAsFactors = FALSE)
    colnames(long) <- c("ref","qry","identity")
    
    p <- ggplot2::ggplot(long, ggplot2::aes(x=factor(.data$ref, ord), y=factor(.data$qry, ord), fill=.data$identity)) +
      ggplot2::geom_tile() +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill="white"),
        plot.title = ggplot2::element_text(size=16, face="bold", hjust=0.5),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
      ) +
      ggplot2::ggtitle(sprintf("%s\nPairwise MASH similarity (cohort n=%d)", gsub("_"," ", species_name), length(ids)))
    
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      p <- p + ggplot2::scale_fill_viridis_c(option="B", limits=c(95,100))
    } else {
      p <- p + ggplot2::scale_fill_gradient(limits=c(95,100))
    }
    
    plot_path <- file.path(mash_dir, paste0("mash_full_heatmap.", plot_format))
    ggplot2::ggsave(plot_path, p, width=7.5, height=6.5, dpi=300)
  }
  
  saveRDS(list(
    species = species_name,
    qids = ids,
    mash_k = mash_k,
    mash_s = mash_s,
    threads = threads,
    sketch_file = msh,
    dist_file = distf,
    dist_matrix = D,
    mean_dist = md_tbl,
    heatmap_path = plot_path
  ), file.path(mash_dir, "mash_full_bundle.rds"))
  
  utils::write.csv(md_tbl, file.path(mash_dir, "mash_full_mean_dist.csv"), row.names = FALSE)
  
  message("Saved mash_full artifacts to: ", mash_dir)
  invisible(list(mash_dir = mash_dir, mean_dist = md_tbl, heatmap_path = plot_path))
}

# --- find_mash_bundle_fp  [from SVCM_result1_helpers.R:36] ---
find_mash_bundle_fp <- function(identity_dir,
                                pattern = "mash.*bundle.*\\.rds$") {
  if (!dir.exists(identity_dir)) {
    stop("identity_dir does not exist: ", identity_dir)
  }
  candidates <- list.files(
    identity_dir, pattern = pattern,
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  if (!length(candidates)) {
    stop("No MASH bundle files found under: ", identity_dir,
         "\n  Pattern: ", pattern)
  }
  # Prefer "mash_full_bundle" filenames, deprioritize "subsplit",
  # then prefer shorter paths.
  score <- ifelse(grepl("mash_full_bundle", basename(candidates),
                        ignore.case = TRUE), 0L, 1L) +
           ifelse(grepl("subsplit", candidates, ignore.case = TRUE), 2L, 0L)
  candidates[order(score, nchar(candidates))[1]]
}

# --- extract_mash_dist  [from SVCM_result1_helpers.R:60] ---
extract_mash_dist <- function(mb) {
  if (!is.null(mb$mash_D)) {
    D <- mb$mash_D
    if (inherits(D, "dist")) return(D)
    if (is.matrix(D))       return(stats::as.dist(D))
  }
  if (!is.null(mb$D)) {
    D <- mb$D
    if (inherits(D, "dist")) return(D)
    if (is.matrix(D))       return(stats::as.dist(D))
  }
  if (is.null(mb$mash_tib)) {
    stop("MASH bundle has neither $mash_D, $D, nor $mash_tib.")
  }

  tib <- mb$mash_tib
  col_id1 <- intersect(c("ref", "id1", "ID1", "qid1", "Reference"), names(tib))[1]
  col_id2 <- intersect(c("query", "id2", "ID2", "qid2", "Query"),     names(tib))[1]
  col_d   <- intersect(c("dist", "distance", "mash_dist", "Distance"), names(tib))[1]
  if (any(is.na(c(col_id1, col_id2, col_d)))) {
    stop("mash_tib missing required columns. Have: ",
         paste(names(tib), collapse = ", "))
  }

  ids <- sort(unique(c(tib[[col_id1]], tib[[col_id2]])))
  n <- length(ids)
  M <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  diag(M) <- 0

  ii <- match(tib[[col_id1]], ids)
  jj <- match(tib[[col_id2]], ids)
  M[cbind(ii, jj)] <- tib[[col_d]]
  M[cbind(jj, ii)] <- tib[[col_d]]

  if (anyNA(M[upper.tri(M)])) {
    n_missing <- sum(is.na(M[upper.tri(M)]))
    stop(sprintf(
      "Reconstructed MASH matrix has %d missing distances in upper triangle. Source mash_tib is incomplete.",
      n_missing
    ))
  }
  stats::as.dist(M)
}

# --- canon_sr_short_arc  [from SVCM.R:3340] ---
canon_sr_short_arc <- function(sr_tbl, genome_len) {
  stopifnot(is.numeric(genome_len), length(genome_len) == 1, genome_len > 0)
  
  sr_tbl %>%
    dplyr::mutate(
      L = as.integer(genome_len),
      width = as.integer(width),
      start = as.integer(start),
      end   = as.integer(end),
      width_short = dplyr::if_else(width > (L/2), as.integer(L - width), width),
      start2 = dplyr::if_else(width > (L/2), as.integer(end + 1L), start),
      end2_raw = dplyr::if_else(width > (L/2), as.integer(start - 1L), end),
      end2 = dplyr::if_else(width > (L/2) & end2_raw < start2, as.integer(end2_raw + L), end2_raw),
      start = start2,
      end   = end2,
      refmid_lin  = as.integer(floor((start + end) / 2)),
      refmid_circ = as.integer(((refmid_lin - 1L) %% L) + 1L),
      width = width_short
    ) %>%
    dplyr::select(-L, -width_short, -start2, -end2_raw, -end2, -refmid_lin)
}
