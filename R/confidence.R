# Confidence bookkeeping.
#
# Two things can go wrong when you interpret a mutation on a predicted
# structure, and this file makes both of them visible rather than invisible:
#
#  1. The residue's own coordinates may be unreliable (low pLDDT).  Below 70 the
#     local geometry is a guess; below 50 the residue is very likely to be in an
#     intrinsically disordered region and its position in space is meaningless.
#  2. Two residues may each be confidently modelled, yet their *relative*
#     placement may not be (high inter-domain PAE).  Any statement of the form
#     "these two mutations cluster in 3D" is then unsupported.
#
# Everything here is reported, never silently filtered, with an explicit
# `interpretable` flag the user can choose to act on.

#' Segment a structure into confident and low-confidence stretches
#'
#' @param structure a `protein_structure`
#' @param threshold pLDDT below which a residue counts as low confidence
#' @param min_len ignore runs shorter than this when reporting
#' @return data.table: `start`, `end`, `length`, `band`, `mean_plddt`,
#'   `interpretable`
#' @export
confidence_segments <- function(structure, threshold = 70, min_len = 5L) {
  r <- structure$residues[order(pos)]
  lowr <- r$plddt < threshold
  run <- rle(lowr)
  ends <- cumsum(run$lengths); starts <- ends - run$lengths + 1L
  out <- data.table(
    start = r$pos[starts], end = r$pos[ends],
    length = run$lengths,
    low = run$values
  )
  out[, mean_plddt := vapply(seq_len(.N), function(i)
    mean(r$plddt[r$pos >= start[i] & r$pos <= end[i]], na.rm = TRUE), 0)]
  out[, band := as.character(plddt_band(mean_plddt))]
  out[, interpretable := !low]
  out[length >= min_len | !low][]
}

#' Attach confidence information to mapped mutations
#'
#' @param mapped output of [map_mutations_to_structure()]
#' @param structure a `protein_structure`
#' @param pae optional PAE matrix (from the `af_model`), used to report how
#'   confidently each mutated residue is placed relative to the rest of the chain
#' @param threshold pLDDT below which a residue is flagged as not interpretable
#' @return `mapped` with added columns `plddt`, `plddt_band`, `burial`,
#'   `n_neighbours`, `pae_global`, `interpretable`, `confidence_note`
#' @export
annotate_confidence <- function(mapped, structure, pae = NULL, threshold = 70) {
  m <- as.data.table(mapped)
  r <- structure$residues
  i <- match(m$struct_pos, r$pos)
  m[, `:=`(plddt = r$plddt[i], plddt_band = r$plddt_band[i],
           burial = r$burial[i], n_neighbours = r$n_neighbours[i])]
  if (!is.null(pae)) {
    pp <- pae_profile(pae)
    off <- if (!is.null(structure$model)) structure$model$uniprot_start - 1L else 0L
    m[, pae_global := pp$pae_global[match(struct_pos - off, pp$pos)]]
  } else {
    m[, pae_global := NA_real_]
  }
  m[, interpretable := !is.na(plddt) & plddt >= threshold]
  m[, confidence_note := data.table::fcase(
    is.na(struct_pos),  "not placed on the structure",
    is.na(plddt),       "residue absent from the model",
    plddt < 50,         "very low pLDDT - likely disordered; do not interpret 3D position",
    plddt < threshold,  "low pLDDT - local geometry uncertain",
    !is.na(pae_global) & pae_global > 15,
      "confident locally but poorly placed relative to the rest of the chain (high PAE)",
    default = "confident")]
  m[]
}

#' One-paragraph modelling-confidence summary for a gene
#'
#' @param structure a `protein_structure`
#' @param model the `af_model`
#' @param mapped output of [annotate_confidence()] (optional)
#' @return list with `text` (character) and the underlying numbers
#' @export
confidence_summary <- function(structure, model, mapped = NULL) {
  r <- structure$residues
  fr <- table(r$plddt_band) / nrow(r)
  getf <- function(nm) as.numeric(fr[nm] %||% 0)
  segs <- confidence_segments(structure)
  low_segs <- segs[low == TRUE & length >= 15]
  n_mut <- if (!is.null(mapped)) sum(!is.na(mapped$struct_pos)) else NA_integer_
  n_low <- if (!is.null(mapped)) sum(!is.na(mapped$struct_pos) & !mapped$interpretable) else NA_integer_

  txt <- sprintf(
    paste0("AlphaFold model %s (%s, v%s) for %s: mean pLDDT %.1f; %.0f%% of residues very high, ",
           "%.0f%% confident, %.0f%% low, %.0f%% very low."),
    model$entry_id, model$accession, as.character(model$version), model$organism %||% "",
    mean(r$plddt, na.rm = TRUE), 100 * getf("Very high (>90)"), 100 * getf("Confident (70-90)"),
    100 * getf("Low (50-70)"), 100 * getf("Very low (<50)"))
  if (nrow(low_segs)) {
    txt <- paste0(txt, sprintf(
      " %d low-confidence stretch(es) of >=15 residues (%s) should be treated as unmodelled, most likely intrinsically disordered.",
      nrow(low_segs), paste(sprintf("%d-%d", low_segs$start, low_segs$end), collapse = ", ")))
  }
  if (!is.na(n_mut)) {
    txt <- paste0(txt, sprintf(" %d of %d placed mutations (%.0f%%) fall in regions below pLDDT 70 and are flagged as not interpretable.",
                               n_low, n_mut, 100 * n_low / max(n_mut, 1)))
  }
  if (is.null(model$pae)) {
    txt <- paste0(txt, " PAE was not loaded, so inter-domain placement confidence was not assessed.")
  }
  list(text = txt, segments = segs, fractions = fr,
       n_mutations = n_mut, n_low_confidence = n_low)
}
