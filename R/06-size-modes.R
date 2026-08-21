# 06-size-modes.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- .find_local_maxima  [from SVCM.R:1898] ---
.find_local_maxima <- function(y, min_peak_distance = 5L) {
  y <- as.numeric(y)
  if (length(y) < 3) return(integer())
  idx <- which(diff(sign(diff(y))) == -2) + 1L
  if (!length(idx)) return(integer())
  idx <- idx[order(y[idx], decreasing = TRUE)]
  keep <- integer()
  for (k in idx) {
    if (!length(keep) || all(abs(k - keep) >= min_peak_distance)) keep <- c(keep, k)
  }
  sort(keep)
}

# --- detect_gene_length_peak_kde  [from SVCM.R:1919] ---
detect_gene_length_peak_kde <- function(
    widths_bp,
    search_bp    = c(500, 5000),
    bw           = "nrd0",
    grid_n       = 2048L,
    cap_bounds_bp = c(300, 8000)
) {
  widths_bp <- as.numeric(widths_bp)
  widths_bp <- widths_bp[is.finite(widths_bp) & widths_bp > 0]
  if (length(widths_bp) < 30) {
    return(tibble::tibble(
      mode_bp = NA_real_, lower_bp = NA_real_, upper_bp = NA_real_,
      n = length(widths_bp), n_in_pile = NA_integer_, frac_in_pile = NA_real_
    ))
  }

  logw <- log10(widths_bp)
  wmin <- log10(search_bp[1]); wmax <- log10(search_bp[2])
  logw_win <- logw[logw >= wmin & logw <= wmax]
  if (length(logw_win) < 20) {
    return(tibble::tibble(
      mode_bp = NA_real_, lower_bp = NA_real_, upper_bp = NA_real_,
      n = length(widths_bp), n_in_pile = NA_integer_, frac_in_pile = NA_real_
    ))
  }

  kde <- stats::density(logw_win, bw = bw, n = grid_n)
  k <- which.max(kde$y)
  mode_log <- kde$x[k]
  half <- kde$y[k] / 2

  left <- k; while (left > 1 && kde$y[left] >= half) left <- left - 1L
  right <- k; while (right < length(kde$y) && kde$y[right] >= half) right <- right + 1L

  lower_bp <- max(10^kde$x[max(1L, left)], cap_bounds_bp[1])
  upper_bp <- min(10^kde$x[min(length(kde$x), right)], cap_bounds_bp[2])

  in_pile <- widths_bp >= lower_bp & widths_bp <= upper_bp

  tibble::tibble(
    mode_bp = 10^mode_log, lower_bp = lower_bp, upper_bp = upper_bp,
    n = length(widths_bp), n_in_pile = sum(in_pile), frac_in_pile = mean(in_pile)
  )
}

# --- detect_gene_length_peak_kde_safe  [from SVCM6.R:362] ---
detect_gene_length_peak_kde_safe <- function(widths_bp,
                                             search_bp = c(500, 5000),
                                             cap_bounds_bp = c(300, 8000)) {
  if (exists("detect_gene_length_peak_kde", mode = "function")) {
    return(detect_gene_length_peak_kde(widths_bp = widths_bp,
                                       search_bp = search_bp,
                                       cap_bounds_bp = cap_bounds_bp))
  }
  detect_gene_length_peak_kde_fallback(widths_bp = widths_bp,
                                       search_bp = search_bp,
                                       cap_bounds_bp = cap_bounds_bp)
}

# --- detect_gene_length_peak_kde_fallback  [from SVCM6.R:314] ---
detect_gene_length_peak_kde_fallback <- function(widths_bp,
                                                 search_bp = c(500, 5000),
                                                 cap_bounds_bp = c(300, 8000),
                                                 n = 4096) {
  w <- as.numeric(widths_bp)
  w <- w[is.finite(w) & w > 0]
  if (length(w) < 100) stop("Too few SV widths for KDE peak detection (n=", length(w), ")")
  
  # restrict to plausible gene-length window for mode search
  w_search <- w[w >= search_bp[1] & w <= search_bp[2]]
  if (length(w_search) < 50) stop("Too few widths inside search window for KDE (n=", length(w_search), ")")
  
  x <- log10(w_search)
  den <- stats::density(x, n = n)
  
  # mode in log space -> bp
  i_mode <- which.max(den$y)
  mode_log <- den$x[i_mode]
  mode_bp  <- 10^mode_log
  
  # FWHM bounds in log space: y >= half max around mode
  half <- den$y[i_mode] / 2
  # find left boundary
  left_idx <- max(which(den$x <= mode_log & den$y <= half), na.rm = TRUE)
  if (!is.finite(left_idx)) left_idx <- 1
  # find right boundary
  right_idx <- min(which(den$x >= mode_log & den$y <= half), na.rm = TRUE)
  if (!is.finite(right_idx)) right_idx <- length(den$x)
  
  lower_bp <- 10^den$x[left_idx]
  upper_bp <- 10^den$x[right_idx]
  
  # cap bounds
  lower_bp <- max(lower_bp, cap_bounds_bp[1])
  upper_bp <- min(upper_bp, cap_bounds_bp[2])
  
  # report fraction in pile using full widths
  in_pile <- w >= lower_bp & w <= upper_bp
  list(
    mode_bp = as.numeric(mode_bp),
    lower_bp = as.numeric(lower_bp),
    upper_bp = as.numeric(upper_bp),
    n = length(w),
    n_in_pile = sum(in_pile),
    frac_in_pile = sum(in_pile) / length(w)
  )
}
