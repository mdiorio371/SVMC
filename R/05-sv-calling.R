# 05-sv-calling.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- assign_domains  [from SVCM.R:869] ---
assign_domains <- function(midpoints, genome_length) {
  pct <- (midpoints / genome_length) * 100
  dplyr::case_when(
    pct > 87.5 | pct <= 12.5  ~ "O",
    pct > 12.5  & pct <= 37.5 ~ "R",
    pct > 37.5  & pct <= 62.5 ~ "T",
    pct > 62.5  & pct <= 87.5 ~ "L",
    TRUE ~ NA_character_
  )
}

# --- get_all_SVs  [from SVCM.R:1863] ---
get_all_SVs <- function(alignment_directory, species_name) {
  all_delta_filt_files <- list.files(alignment_directory, full.names = TRUE, pattern = "filtered.delta")
  all_delta_unfilt_files <- list.files(sprintf("%s/unfiltered", alignment_directory), full.names = TRUE)

  filt_core_names   <- sub("_filtered\\.delta$", "", basename(all_delta_filt_files))
  unfilt_core_names <- sub("\\.delta$", "", basename(all_delta_unfilt_files))
  unfilt_lookup     <- setNames(all_delta_unfilt_files, unfilt_core_names)

  n <- length(all_delta_filt_files)
  all_res <- vector("list", n)
  pb <- txtProgressBar(min = 1, max = n, style = 3)

  for (i in seq_len(n)) {
    core_name <- filt_core_names[i]
    if (!core_name %in% names(unfilt_lookup)) {
      warning(sprintf("No unfiltered file found for %s, skipping...", all_delta_filt_files[i]))
      next
    }
    all_res[[i]] <- delta_all_SVs(
      delta_file = all_delta_filt_files[i],
      unfiltered_delta_file = unfilt_lookup[[core_name]],
      species_name = species_name
    )
    setTxtProgressBar(pb, i)
  }
  close(pb)

  dplyr::bind_rows(Filter(Negate(is.null), all_res))
}
