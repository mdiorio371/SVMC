# ==========================================================================
# SVMC: motif cohesion scoring
# Lineage-defining structural-variation motif identification.
# ==========================================================================

#' Structural-variation motif cohesion
#'
#' Quantifies how genomically cohesive the carriers of a structural-variation
#' motif are relative to the cohort as a whole. Carriers of a lineage-defining
#' motif are more similar to one another (smaller pairwise MASH distances) than
#' two genomes drawn at random from the cohort, which yields a high cohesion
#' score.
#'
#' The score is defined as
#' \deqn{1 - \frac{\bar{d}_{within}}{\bar{d}_{overall}}}
#' where \eqn{\bar{d}_{within}} is the mean pairwise distance among the motif's
#' carriers and \eqn{\bar{d}_{overall}} is the mean pairwise distance across the
#' whole cohort. A score near `1` indicates a tightly clustered, likely
#' lineage-defining motif; a score near `0` indicates carriers no more similar
#' than the cohort average; negative scores indicate carriers that are *more*
#' divergent than the cohort average.
#'
#' Carriers are matched to the distance matrix **by name**, never by position,
#' so `M_all` and `dist_matrix` need not share row ordering.
#'
#' @param motif_uid Identifier of the motif in `M_all`: either a column name
#'   (character) or a column index (integer).
#' @param M_all Carrier matrix. Rows are genomes (`rownames` = genome IDs),
#'   columns are motifs (`colnames` = motif UIDs). A genome counts as a carrier
#'   of a motif when its entry is greater than zero.
#' @param dist_matrix Square, symmetric pairwise-distance matrix (e.g. MASH
#'   distances). Row and column names must be genome IDs matching
#'   `rownames(M_all)`.
#'
#' @return A single numeric cohesion score, or `NA_real_` when fewer than two
#'   carriers are present in `dist_matrix`, or when the cohort-wide mean
#'   distance is zero.
#'
#' @examples
#' ## toy cohort: 5 genomes, 2 motifs
#' ids <- paste0("g", 1:5)
#' M <- matrix(c(1, 1, 1, 0, 0,    # motif A carried by g1-g3
#'               1, 0, 0, 0, 1),   # motif B carried by g1 and g5
#'             nrow = 5, dimnames = list(ids, c("A", "B")))
#' ## g1-g3 are mutually close; g4 and g5 are outliers
#' D <- as.matrix(dist(c(g1 = 0, g2 = 0.01, g3 = 0.02, g4 = 0.5, g5 = 0.9)))
#' svmc_cohesion("A", M, D)   # high: tight cluster
#' svmc_cohesion("B", M, D)   # low / negative: carriers far apart
#'
#' @seealso [svmc_cohesion_screen()] to score every motif at once.
#' @export
svmc_cohesion <- function(motif_uid, M_all, dist_matrix) {
  carriers <- rownames(M_all)[M_all[, motif_uid] > 0]
  carriers <- carriers[!is.na(carriers)]
  if (length(carriers) < 2L) return(NA_real_)

  ## Match carriers to the distance matrix BY NAME, never positionally:
  ## M_all and dist_matrix may be ordered differently.
  carrier_idx <- match(carriers, rownames(dist_matrix))
  carrier_idx <- carrier_idx[!is.na(carrier_idx)]
  if (length(carrier_idx) < 2L) return(NA_real_)

  within <- dist_matrix[carrier_idx, carrier_idx, drop = FALSE]
  diag(within) <- NA
  mean_within  <- mean(within, na.rm = TRUE)
  mean_overall <- mean(dist_matrix[upper.tri(dist_matrix)], na.rm = TRUE)

  if (!is.finite(mean_overall) || mean_overall == 0) return(NA_real_)
  1 - (mean_within / mean_overall)
}

#' Score every motif in a carrier matrix
#'
#' Convenience wrapper that applies [svmc_cohesion()] to every motif (column)
#' of `M_all` and returns the results as a tibble. This is the function behind
#' the cross-species cohesion screen: motifs in the high-cohesion tail (by
#' default `>= 0.40`) are the lineage-defining candidates.
#'
#' @param M_all Carrier matrix; see [svmc_cohesion()].
#' @param dist_matrix Pairwise-distance matrix; see [svmc_cohesion()].
#' @param min_cohesion Optional numeric threshold. When supplied, only motifs
#'   with `cohesion >= min_cohesion` (and non-`NA`) are returned. Default
#'   `NULL` returns all motifs.
#'
#' @return A tibble with one row per motif and columns `motif_uid`,
#'   `n_carriers`, and `cohesion`, sorted by descending `cohesion`
#'   (`NA` scores sort last).
#'
#' @examples
#' ids <- paste0("g", 1:5)
#' M <- matrix(c(1, 1, 1, 0, 0,  1, 0, 0, 0, 1), nrow = 5,
#'             dimnames = list(ids, c("A", "B")))
#' D <- as.matrix(dist(c(g1 = 0, g2 = 0.01, g3 = 0.02, g4 = 0.5, g5 = 0.9)))
#' svmc_cohesion_screen(M, D)
#' svmc_cohesion_screen(M, D, min_cohesion = 0.40)
#'
#' @seealso [svmc_cohesion()]
#' @importFrom tibble tibble
#' @export
svmc_cohesion_screen <- function(M_all, dist_matrix, min_cohesion = NULL) {
  motif_uids <- colnames(M_all)
  if (is.null(motif_uids)) {
    stop("`M_all` must have column names (motif UIDs).", call. = FALSE)
  }

  cohesion <- vapply(
    motif_uids,
    function(u) svmc_cohesion(u, M_all, dist_matrix),
    numeric(1)
  )
  n_carriers <- vapply(
    motif_uids,
    function(u) sum(M_all[, u] > 0, na.rm = TRUE),
    integer(1)
  )

  out <- tibble::tibble(
    motif_uid  = motif_uids,
    n_carriers = n_carriers,
    cohesion   = unname(cohesion)
  )
  out <- out[order(out$cohesion, decreasing = TRUE), ]

  if (!is.null(min_cohesion)) {
    out <- out[!is.na(out$cohesion) & out$cohesion >= min_cohesion, ]
  }
  out
}
