# ==========================================================================
# SVMC: SV size-mode decomposition
# KDE peak assignment and recursive Gaussian-mixture mode finding.
# ==========================================================================

assign_kde_peaks <- function(df, core_thresh = 0.1, shoulder_thresh = 0.3) {
  if (nrow(df) < 5 || !"width" %in% names(df)) return(NULL)
  
  df <- df %>%
    filter(!is.na(width)) %>%
    mutate(log_length = log10(width + 1))
  
  kde <- density(df$log_length, bw = "nrd0")
  peaks <- findpeaks(kde$y, minpeakdistance = 5)
  if (is.null(peaks)) return(NULL)
  
  peak_pos <- kde$x[peaks[, 2]]
  peak_ids <- paste0("Peak_", seq_along(peak_pos))
  
  df %>%
    rowwise() %>%
    mutate(
      nearest = which.min(abs(log_length - peak_pos)),
      mode_group = peak_ids[nearest],
      mode_mean_log = peak_pos[nearest],
      mode_mean_bp  = 10^mode_mean_log,
      log_distance  = abs(log_length - mode_mean_log),
      zone = case_when(
        log_distance <= core_thresh ~ "core",
        log_distance <= shoulder_thresh ~ "shoulder",
        TRUE ~ "outlier"
      ),
      n_modes = length(peak_pos)
    ) %>%
    ungroup()
}

recursive_gmm <- function(
    x,
    max_depth = 5,
    depth = 0,
    min_n = 100,
    G_range = 1:6,
    proportion_cutoff = 0.1,
    spread_ratio_cutoff = 1.0,
    modelNames = "V"
) {
  # ---- input hygiene ---------------------------------------------------------
  x <- x[is.finite(x)]
  if (!is.numeric(x)) stop("`x` must be a numeric vector of log10-lengths.", call. = FALSE)
  if (length(x) < min_n || depth >= max_depth) {
    return(tibble::tibble())
  }
  
  # ---- fit mixture -----------------------------------------------------------
  fit <- try(mclust::Mclust(x, G = G_range, modelNames = modelNames), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$parameters$mean)) {
    return(tibble::tibble())
  }
  
  means <- as.numeric(fit$parameters$mean)
  # for spherical models, 'sigmasq' is a scalar per component
  sigsq <- fit$parameters$variance$sigmasq
  sds   <- sqrt(as.numeric(sigsq))
  props <- as.numeric(fit$parameters$pro)
  
  # guard lengths
  k <- length(means)
  if (any(lengths(list(sds, props)) != k)) {
    return(tibble::tibble())
  }
  
  mean_bp <- 10^means
  sd_bp   <- 10^(means + sds) - mean_bp  # +1σ offset in bp (for banding)
  
  summary_table <- tibble::tibble(
    depth      = depth,
    submode    = paste0("D", depth, ".", seq_len(k)),
    mean_log10 = means,
    mean_bp    = mean_bp,
    sd_log10   = sds,
    sd_bp      = sd_bp,
    proportion = props,
    n          = length(x)
  )
  
  # ---- decide which components to recurse into -------------------------------
  nested_tables <- list(summary_table)
  
  for (i in seq_len(k)) {
    prop_i <- props[i]
    spread_ratio_i <- ifelse(mean_bp[i] > 0, sd_bp[i] / mean_bp[i], 0)
    
    if (is.finite(prop_i) && is.finite(spread_ratio_i) &&
        prop_i > proportion_cutoff && spread_ratio_i > spread_ratio_cutoff) {
      
      idx <- abs(x - means[i]) < 2 * sds[i]  # ±2 SD window in log10 space
      if (sum(idx) >= min_n) {
        nested_result <- recursive_gmm(
          x = x[idx],
          max_depth = max_depth,
          depth = depth + 1,
          min_n = min_n,
          G_range = if (length(G_range) > 0) pmax(1, pmin(3, G_range)) else 1:3,
          proportion_cutoff = proportion_cutoff,
          spread_ratio_cutoff = spread_ratio_cutoff,
          modelNames = modelNames
        )
        if (nrow(nested_result) > 0) nested_tables <- append(nested_tables, list(nested_result))
      }
    }
  }
  
  dplyr::bind_rows(nested_tables)
}
