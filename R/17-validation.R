# 17-validation.R
# Validate structural-variant states against an external typing scheme
# (e.g. SeqSero2 serovars). Precision/recall of state -> lineage.

#' Validate structural-variant states against an external typing scheme
#'
#' For each structural-variant *state* -- defined by a marker motif whose
#' carriers constitute the state -- computes precision and recall against a
#' reference lineage assignment (for example, SeqSero2 serovar calls). This is
#' the independent test of whether a structurally-defined state corresponds to
#' a real biological lineage.
#'
#' \strong{Precision} = fraction of the state's carriers that belong to the
#' target lineage (does carrying the inversion predict the lineage?).
#' \strong{Recall} = fraction of the target lineage's genomes that carry the
#' state (does the inversion capture the whole lineage, or a sublineage?).
#'
#' The exclusivity partition (from [svmc_inversion_exclusivity()]) is computed
#' from structure alone; the reference typing here is independent data, so a
#' state that scores high precision is a genuine structure-to-lineage
#' correspondence, not a circular restatement.
#'
#' @param M Presence matrix (genomes x motifs); carrier if entry > 0.
#' @param serovar Either a named character vector (names = genome IDs, values =
#'   reference call) or a data frame with columns `qid` and `serotype`.
#' @param state_map Named character vector: names are motif UIDs (columns of
#'   `M`), values are human-readable state labels. Each named motif defines a
#'   state consisting of its carriers.
#' @param lineage_map Named list mapping each lineage label to a definition:
#'   `list("Typhi" = list(state = "B", serovars = "Typhi"), ...)`, where `state`
#'   is a value used in `state_map` and `serovars` is a character vector of
#'   reference calls composing the lineage.
#'
#' @return A tibble with one row per lineage: `lineage`, `state`, precision
#'   (`precision_n`, `precision_d`, `precision_pct`) and recall (`recall_n`,
#'   `recall_d`, `recall_pct`).
#'
#' @examples
#' ids <- paste0("g", 1:8)
#' M <- matrix(0, 8, 2, dimnames = list(ids, c("mA", "mB")))
#' M[1:4, "mA"] <- 1        # state A carriers
#' M[5:7, "mB"] <- 1        # state B carriers
#' sero <- c(g1="X", g2="X", g3="X", g4="Y",   # A is mostly serovar X
#'           g5="Z", g6="Z", g7="Z", g8="X")   # B is serovar Z
#' state_map <- c(mA = "A", mB = "B")
#' lineage_map <- list(
#'   "X-lineage" = list(state = "A", serovars = "X"),
#'   "Z-lineage" = list(state = "B", serovars = "Z")
#' )
#' svmc_state_serovar_validation(M, sero, state_map, lineage_map)
#'
#' @seealso [svmc_inversion_exclusivity()]
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @export
svmc_state_serovar_validation <- function(M, serovar, state_map, lineage_map) {
  if (is.data.frame(serovar)) {
    sv <- stats::setNames(as.character(serovar$serotype), serovar$qid)
  } else {
    sv <- serovar
  }

  motif_ids <- names(state_map)
  missing_motifs <- setdiff(motif_ids, colnames(M))
  if (length(missing_motifs)) {
    stop("These state_map motifs are not columns of `M`: ",
         paste(missing_motifs, collapse = ", "), call. = FALSE)
  }

  # per-genome primary state (label from state_map, or NA)
  carrier_state <- vapply(rownames(M), function(g) {
    hits <- motif_ids[which(M[g, motif_ids] > 0)]
    if (length(hits) == 0L) NA_character_
    else if (length(hits) > 1L) "Multiple"
    else unname(state_map[hits])
  }, character(1))

  genome_sv <- sv[rownames(M)]

  dplyr::bind_rows(lapply(names(lineage_map), function(lin) {
    def         <- lineage_map[[lin]]
    state_label <- def$state
    serovars    <- def$serovars

    in_state        <- which(carrier_state == state_label)
    tp              <- sum(genome_sv[in_state] %in% serovars, na.rm = TRUE)
    n_state         <- length(in_state)
    n_lineage_total <- sum(genome_sv %in% serovars, na.rm = TRUE)

    tibble::tibble(
      lineage       = lin,
      state         = state_label,
      precision_n   = tp,
      precision_d   = n_state,
      precision_pct = if (n_state > 0) round(100 * tp / n_state, 1) else NA_real_,
      recall_n      = tp,
      recall_d      = n_lineage_total,
      recall_pct    = if (n_lineage_total > 0) round(100 * tp / n_lineage_total, 1) else NA_real_
    )
  }))
}
