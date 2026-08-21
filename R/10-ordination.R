# 10-ordination.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- select_across_pcoa  [from SVCM.R:3865] ---
select_across_pcoa <- function(scores_df, n=100L, nbins=10L, seed=1L) {
  set.seed(seed)
  df <- scores_df %>% dplyr::filter(is.finite(.data$Axis.1), is.finite(.data$Axis.2), nzchar(.data$qid))
  if (nrow(df) <= n) return(df$qid)
  
  bx <- unique(stats::quantile(df$Axis.1, probs=seq(0,1,length.out=nbins+1), na.rm=TRUE))
  by <- unique(stats::quantile(df$Axis.2, probs=seq(0,1,length.out=nbins+1), na.rm=TRUE))
  if (length(bx) < 3 || length(by) < 3) {
    nb2 <- max(3L, floor(nbins/2))
    bx <- unique(stats::quantile(df$Axis.1, probs=seq(0,1,length.out=nb2+1), na.rm=TRUE))
    by <- unique(stats::quantile(df$Axis.2, probs=seq(0,1,length.out=nb2+1), na.rm=TRUE))
  }
  
  df$binx <- cut(df$Axis.1, breaks=bx, include.lowest=TRUE, right=FALSE)
  df$biny <- cut(df$Axis.2, breaks=by, include.lowest=TRUE, right=FALSE)
  df$cell <- interaction(df$binx, df$biny, drop=TRUE)
  
  reps <- df %>% dplyr::group_by(.data$cell) %>% dplyr::slice_sample(n=1) %>% dplyr::ungroup()
  
  farthest_first <- function(all_df, seed_df, target_n) {
    all_xy <- as.matrix(all_df[,c("Axis.1","Axis.2")])
    seed_xy <- as.matrix(seed_df[,c("Axis.1","Axis.2")])
    dmin <- rep(Inf, nrow(all_df))
    for (k in seq_len(nrow(seed_xy))) {
      dx <- all_xy[,1]-seed_xy[k,1]; dy <- all_xy[,2]-seed_xy[k,2]
      dmin <- pmin(dmin, dx*dx + dy*dy)
    }
    chosen <- unique(seed_df$qid)
    while (length(chosen) < target_n) {
      already <- all_df$qid %in% chosen
      dmin2 <- dmin; dmin2[already] <- -Inf
      j <- which.max(dmin2)
      chosen <- c(chosen, all_df$qid[j])
      dx <- all_xy[,1]-all_xy[j,1]; dy <- all_xy[,2]-all_xy[j,2]
      dmin <- pmin(dmin, dx*dx + dy*dy)
    }
    chosen
  }
  
  if (nrow(reps) >= n) farthest_first(reps, reps[1,,drop=FALSE], n) else farthest_first(df, reps, n)
}

# --- choose_dbscan_eps_sane  [from SVCM.R:3950] ---
choose_dbscan_eps_sane <- function(db_out, max_clusters = 9L) {
  cand <- db_out |>
    dplyr::filter(n_clusters >= 2L, n_clusters <= max_clusters)
  
  if (nrow(cand) == 0) {
    cand <- db_out |>
      dplyr::filter(n_clusters >= 2L) |>
      dplyr::mutate(cluster_dev = abs(n_clusters - max_clusters)) |>
      dplyr::arrange(cluster_dev, dplyr::desc(eps), pct_noise)
    return(cand |> dplyr::slice(1))
  }
  
  cand |>
    dplyr::arrange(dplyr::desc(eps), pct_noise, n_clusters) |>
    dplyr::slice(1)
}

# --- summarize_dbscan_grid  [from SVCM_result1_helpers.R:439] ---
summarize_dbscan_grid <- function(coords, minPts = 5L, n_grid = 30L) {
  coords <- as.matrix(coords)
  d <- as.numeric(stats::dist(coords))
  qs <- stats::quantile(d, probs = seq(0.01, 0.30, length.out = n_grid),
                        na.rm = TRUE)
  eps_grid <- unname(as.numeric(qs))

  rows <- lapply(eps_grid, function(eps) {
    fit <- dbscan::dbscan(coords, eps = eps, minPts = minPts)
    n_clust <- length(unique(fit$cluster[fit$cluster > 0]))
    pct_noise <- mean(fit$cluster == 0)
    tibble::tibble(eps = eps, n_clusters = n_clust, pct_noise = pct_noise)
  })
  dplyr::bind_rows(rows)
}

# --- collapse_dbscan_clusters  [from SVCM_result1_helpers.R:458] ---
collapse_dbscan_clusters <- function(cluster_raw, min_cluster_n = 5L) {
  out <- as.character(cluster_raw)
  out[out == "0"] <- "Noise"
  tab <- table(out)
  small <- names(tab)[tab < min_cluster_n & !names(tab) %in% c("Noise")]
  out[out %in% small] <- "Other"
  out
}

# --- make_ord_df  [from SVCM_result1_helpers.R:470] ---
make_ord_df <- function(D) {
  if (!inherits(D, "dist")) D <- stats::as.dist(D)
  fit <- stats::cmdscale(D, k = 2, eig = TRUE)
  pts <- as.data.frame(fit$points)
  names(pts) <- c("MDS1", "MDS2")
  ev_pos <- pmax(fit$eig, 0)
  tot <- sum(ev_pos)
  pct <- ev_pos[1:2] / tot * 100
  list(
    df       = pts,
    ord_pct  = c(MDS1 = pct[1], MDS2 = pct[2]),
    fit      = fit
  )
}

# --- flatness_report  [from SVCM_result1_helpers.R:489] ---
flatness_report <- function(D, ord = NULL) {
  if (is.null(ord)) ord <- make_ord_df(D)
  ev_pos <- pmax(ord$fit$eig, 0)
  tot <- sum(ev_pos)
  if (tot <= 0) {
    return(c(pct_axis1 = NA, pct_axis2 = NA, pct_axis1_2 = NA,
             n_axes_for_80pct = NA_integer_))
  }
  pct1 <- ev_pos[1] / tot * 100
  pct2 <- ev_pos[2] / tot * 100
  cum_frac <- cumsum(ev_pos) / tot
  c(
    pct_axis1        = pct1,
    pct_axis2        = pct2,
    pct_axis1_2      = pct1 + pct2,
    n_axes_for_80pct = as.integer(which(cum_frac >= 0.80)[1])
  )
}

# --- get_stepC_inputs  [from SVCM.R:5057] ---
get_stepC_inputs <- function(species_name, clade_name) {
  step1_ck <- file.path(
    "processing/temp",
    clade_name,
    species_name,
    "step1_checkpoint.rds"
  )
  stopifnot(file.exists(step1_ck))
  
  step1_data <- readRDS(step1_ck)
  dirs <- step1_data$dirs
  ref_qid <- as.character(step1_data$most_related)
  ref_path <- file.path(dirs$sync, paste0(ref_qid, ".txt"))
  out_sv_dir <- file.path(dirs$results, "sv_structure_stepBC")
  sv_fp <- file.path(dirs$results, "all_svs_raw.rds")
  
  stopifnot(file.exists(ref_path))
  stopifnot(file.exists(sv_fp))
  dir.create(out_sv_dir, recursive = TRUE, showWarnings = FALSE)
  
  list(
    step1_data = step1_data,
    dirs = dirs,
    ref_qid = ref_qid,
    ref_path = ref_path,
    out_sv_dir = out_sv_dir,
    sv_fp = sv_fp
  )
}
