# Functional domains on the structure, with cross-species transfer.
#
# The problem this solves: a macaque, dog or naked mole-rat UniProt entry
# usually has no curated active sites or DNA-binding regions, and sometimes no
# reviewed entry at all. The human ortholog almost always does. So we take the
# human annotation and walk it across the ortholog protein alignment onto the
# residue numbering of whichever structure we are drawing.
#
# Transferred features are never silently presented as native: every row keeps
# `source`, `transferred` and `transfer_identity` so the figure and the report
# can say where the annotation came from.

#' Build the annotated domain table for a structure
#'
#' @param accession UniProt accession of the protein being drawn
#' @param human_accession UniProt accession of the human ortholog, used as the
#'   annotation donor when `transfer = TRUE`. `NULL` disables transfer.
#' @param map optional residue map *from* the human protein *to* the structure
#'   protein (see [residue_map()]). Required for transfer; if `NULL` and both
#'   accessions are given, one is computed from the two UniProt sequences.
#' @param transfer transfer human annotation across the alignment
#' @param interpro also fetch InterPro signatures for the structure's own entry
#' @param transfer_identity percent identity of the alignment used for transfer,
#'   recorded on every transferred row
#' @return data.table: `start`, `end`, `type`, `class`, `description`, `point`,
#'   `evidence`, `source`, `transferred`, `transfer_identity`, `label`
#' @export
domain_table <- function(accession, human_accession = NULL, map = NULL,
                         transfer = TRUE, interpro = TRUE,
                         transfer_identity = NA_real_) {
  own <- uniprot_features(accession)
  own[, `:=`(transferred = FALSE, transfer_identity = NA_real_)]

  if (isTRUE(interpro)) {
    ip <- interpro_domains(accession)
    if (nrow(ip)) {
      ip[, `:=`(transferred = FALSE, transfer_identity = NA_real_)]
      own <- rbindlist(list(own, ip), fill = TRUE)
    }
  }

  if (isTRUE(transfer) && !is.null(human_accession) &&
      !identical(human_accession, accession)) {
    hu <- uniprot_features(human_accession)
    if (nrow(hu)) {
      if (is.null(map)) {
        hs <- uniprot_entry(human_accession)$sequence
        ts <- uniprot_entry(accession)$sequence
        map <- residue_map(hs, ts)
        transfer_identity <- attr(map, "pid")
      }
      hu[, `:=`(start = .apply_map(start, map), end = .apply_map(end, map))]
      lost <- sum(is.na(hu$start) | is.na(hu$end))
      hu <- hu[!is.na(start) & !is.na(end) & end >= start]
      if (lost) .msg(sprintf("%d human feature(s) could not be transferred (aligned to gaps).", lost))
      hu[, `:=`(source = paste0(source, " (human ", human_accession, ")"),
                transferred = TRUE, transfer_identity = transfer_identity)]
      own <- rbindlist(list(own, hu), fill = TRUE)
    }
  }

  if (!nrow(own)) {
    return(data.table(start = integer(), end = integer(), type = character(),
                      class = character(), description = character(),
                      point = logical(), evidence = character(), source = character(),
                      transferred = logical(), transfer_identity = numeric(),
                      label = character()))
  }
  # prefer a native annotation over a transferred duplicate of the same span
  setorder(own, start, end, transferred)
  own <- own[!duplicated(own[, .(start, end, class, description)])]
  own[, label := .short_label(description, class)]
  own[]
}

.short_label <- function(description, class) {
  d <- as.character(description)
  d[is.na(d) | !nzchar(d)] <- as.character(class)[is.na(d) | !nzchar(d)]
  d <- sub("\\s*\\(IPR\\d+\\)$", "", d)
  d <- sub("^(.{26}).{3,}$", "\\1\u2026", d)
  d
}

#' Which domain(s) does each residue fall in?
#'
#' @param pos integer vector of residue positions
#' @param domains table from [domain_table()]
#' @param classes feature classes to consider (default: the structural ones;
#'   point features such as active sites are handled by [key_sites()])
#' @return character vector, one `"; "`-separated string per position
#' @export
annotate_positions <- function(pos, domains, classes = c("Domain", "Region", "Repeat",
                                                         "Zinc finger", "DNA binding",
                                                         "Coiled coil", "Motif",
                                                         "Transmembrane", "Signal peptide")) {
  d <- as.data.table(domains)[class %in% classes]
  vapply(pos, function(p) {
    if (is.na(p)) return(NA_character_)
    h <- d[start <= p & end >= p]
    if (!nrow(h)) return(NA_character_)
    paste(unique(h$label), collapse = "; ")
  }, "")
}

#' Key functional sites (single-residue features)
#'
#' @param domains table from [domain_table()]
#' @return data.table: `pos`, `class`, `description`, `source`, `transferred`
#' @export
key_sites <- function(domains) {
  d <- as.data.table(domains)[point == TRUE]
  if (!nrow(d)) return(data.table(pos = integer(), class = character(),
                                  description = character(), source = character(),
                                  transferred = logical()))
  out <- d[, .(pos = seq.int(start, end)), by = .(class, description, source, transferred)]
  unique(out[, .(pos, class, description, source, transferred)])[order(pos)]
}

#' Flag mutations that hit, or sit next to, a key functional site
#'
#' A mutation is reported as *functionally relevant* when it falls inside an
#' annotated domain, directly hits a key site, or lies within `radius`
#' Angstrom in 3D of a key site even though it is distant in sequence. The
#' third case is the one a lollipop plot cannot show you.
#'
#' @param pos residue positions of the mutations
#' @param structure a `protein_structure`
#' @param domains table from [domain_table()]
#' @param radius Angstrom for the 3D-proximity test
#' @return data.table: `pos`, `domain`, `site_hit`, `site_near`, `site_near_dist`
#' @export
flag_functional <- function(pos, structure, domains, radius = 8) {
  ks <- key_sites(domains)
  d <- structure$dist
  dom <- annotate_positions(pos, domains)
  hit <- vapply(pos, function(p) {
    h <- ks[pos == p]
    if (!nrow(h)) NA_character_ else paste(unique(sprintf("%s: %s", h$class, h$description)), collapse = "; ")
  }, "")
  near <- rep(NA_character_, length(pos)); near_d <- rep(NA_real_, length(pos))
  if (nrow(ks)) {
    sites <- intersect(as.character(ks$pos), colnames(d))
    for (i in seq_along(pos)) {
      p <- as.character(pos[i])
      if (is.na(pos[i]) || !p %in% colnames(d) || !length(sites)) next
      dd <- d[p, sites, drop = TRUE]
      dd <- dd[dd > 0 & dd <= radius]
      if (!length(dd)) next
      j <- as.integer(names(sort(dd))[1])
      h <- ks[pos == j]
      near[i]   <- paste(unique(sprintf("%s: %s (res %d)", h$class, h$description, j)), collapse = "; ")
      near_d[i] <- min(dd)
    }
  }
  data.table(pos = pos, domain = dom, site_hit = hit,
             site_near = near, site_near_dist = round(near_d, 2))
}
