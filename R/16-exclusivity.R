# 16-exclusivity.R
# Mutual-exclusivity classification of structural-variant motifs.
# The central SVMC concept: partition recurrent inversions into those that are
# mutually exclusive (alternative states, never co-occur) versus those that
# co-occur -- computed from genome structure alone, blind to lineage.

#' Classify structural variants by mutual exclusivity
#'
#' Partitions a set of structural-variant motifs into *exclusive* and
#' *co-occurring* classes based purely on their co-occurrence across genomes.
#' Two motifs are mutually exclusive when (almost) no genome carries both; a
#' motif is "co-occurring" when it shares carriers with another motif above a
#' tolerance. This partition uses only the presence/absence matrix -- no
#' lineage, serovar, or phylogenetic information -- which is what makes any
#' downstream association with lineage a genuine (non-circular) finding.
#'
#' The exclusivity score for each motif is the largest fraction of its carriers
#' that it shares with any *other* motif in the set:
#' \deqn{\mathrm{frac\_shared}_i = \max_{j \ne i} \frac{n_{ij}}{n_i}}
#' where \eqn{n_{ij}} is the number of genomes carrying both motifs \eqn{i} and
#' \eqn{j}, and \eqn{n_i} is the number carrying motif \eqn{i}. A motif with
#' `frac_shared` below `max_shared` is labelled `"exclusive"`.
#'
#' @param M A presence matrix (genomes x motifs). Rows are genomes, columns are
#'   motif IDs. A genome carries a motif when its entry is greater than zero.
#'   May be a base matrix or a sparse `Matrix::dgCMatrix`.
#' @param motif_ids Optional character vector of motif IDs (columns of `M`) to
#'   restrict the analysis to. Defaults to all columns of `M`.
#' @param max_shared Numeric tolerance in `[0, 1]`. A motif is called
#'   `"exclusive"` when its `frac_shared` is strictly below this value. The
#'   default `0.05` allows a small number of doubly-carried genomes (e.g. from
#'   miscalls) without breaking exclusivity.
#'
#' @return A tibble with one row per motif, sorted by descending `n_carriers`:
#'   \describe{
#'     \item{motif_uid}{the motif ID}
#'     \item{n_carriers}{number of genomes carrying the motif}
#'     \item{max_cooccur}{largest number of carriers shared with any other motif}
#'     \item{top_partner}{the motif ID it shares the most carriers with}
#'     \item{frac_shared}{`max_cooccur / n_carriers`}
#'     \item{class}{`"exclusive"` or `"co-occurring"`}
#'   }
#'
#' @examples
#' ## three exclusive states + one co-occurring passenger
#' ids <- paste0("g", 1:10)
#' M <- matrix(0, nrow = 10, ncol = 4,
#'             dimnames = list(ids, c("A", "B", "C", "P")))
#' M[1:4,  "A"] <- 1                 # state A
#' M[5:7,  "B"] <- 1                 # state B (exclusive to A)
#' M[8:10, "C"] <- 1                 # state C (exclusive to A, B)
#' M[c(1,2,5,6), "P"] <- 1           # passenger P: rides on A and B
#' svmc_inversion_exclusivity(M)
#'
#' @seealso [svmc_state_serovar_validation()] to test whether the exclusive
#'   states predict an external typing scheme.
#' @importFrom tibble tibble
#' @export
svmc_inversion_exclusivity <- function(M, motif_ids = NULL, max_shared = 0.05) {
  if (is.null(colnames(M))) {
    stop("`M` must have column names (motif IDs).", call. = FALSE)
  }
  ids <- if (is.null(motif_ids)) colnames(M) else intersect(motif_ids, colnames(M))
  if (length(ids) < 2L) {
    stop("Need at least two motifs present in `M` to assess exclusivity.",
         call. = FALSE)
  }

  X <- (M[, ids, drop = FALSE] > 0) * 1
  # motif x motif co-occurrence counts; diagonal = carrier counts
  cooc <- as.matrix(Matrix::crossprod(X))
  n_carriers <- diag(cooc)

  k <- length(ids)
  max_cooccur <- integer(k)
  top_partner <- rep(NA_character_, k)
  for (i in seq_len(k)) {
    others <- cooc[i, -i]
    if (length(others) && max(others) > 0) {
      max_cooccur[i] <- as.integer(max(others))
      top_partner[i] <- ids[-i][which.max(others)]
    } else {
      max_cooccur[i] <- 0L
    }
  }

  frac_shared <- ifelse(n_carriers > 0, max_cooccur / n_carriers, 0)

  out <- tibble::tibble(
    motif_uid   = ids,
    n_carriers  = as.integer(n_carriers),
    max_cooccur = max_cooccur,
    top_partner = top_partner,
    frac_shared = round(frac_shared, 3),
    class       = ifelse(frac_shared < max_shared, "exclusive", "co-occurring")
  )
  out[order(out$n_carriers, decreasing = TRUE), ]
}


#' Test whether a set of motifs is mutually exclusive
#'
#' Convenience check: returns `TRUE` when no genome carries more than one of the
#' supplied motifs (perfect mutual exclusivity), along with the number of
#' multiply-carrying genomes. Useful for confirming that a candidate set of
#' "primary" states genuinely partitions the cohort.
#'
#' @param M A presence matrix (genomes x motifs); see
#'   [svmc_inversion_exclusivity()].
#' @param motif_ids Character vector of motif IDs to test together.
#'
#' @return A list with `mutually_exclusive` (logical), `n_multiple` (number of
#'   genomes carrying more than one of the motifs), and `cooccurrence` (the
#'   motif-by-motif co-occurrence matrix; off-diagonal zeros indicate
#'   exclusivity).
#'
#' @examples
#' ids <- paste0("g", 1:9)
#' M <- matrix(0, 9, 3, dimnames = list(ids, c("A", "B", "C")))
#' M[1:3, "A"] <- 1; M[4:6, "B"] <- 1; M[7:9, "C"] <- 1
#' svmc_are_mutually_exclusive(M, c("A", "B", "C"))$mutually_exclusive  # TRUE
#'
#' @export
svmc_are_mutually_exclusive <- function(M, motif_ids) {
  ids <- intersect(motif_ids, colnames(M))
  if (length(ids) < 2L) {
    stop("Need at least two motifs present in `M`.", call. = FALSE)
  }
  X <- (M[, ids, drop = FALSE] > 0) * 1
  cooc <- as.matrix(Matrix::crossprod(X))
  n_multiple <- sum(Matrix::rowSums(X) > 1)
  off_diag <- cooc[upper.tri(cooc)]
  list(
    mutually_exclusive = all(off_diag == 0),
    n_multiple = as.integer(n_multiple),
    cooccurrence = cooc
  )
}
