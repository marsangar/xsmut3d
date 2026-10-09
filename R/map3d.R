# Placing mutations onto structure residue numbering.
#
# The anchor is always the sequence of the AlphaFold model (a UniProt sequence).
# Every protein whose coordinates we need to place -- the species transcript the
# calls were annotated against, the human ortholog, the COSMIC transcript -- is
# aligned to that one anchor, so there is a single, auditable coordinate system
# and a recorded percent identity for each hop.
#
#   dNdScv aa_pos (species protein ENSxxxP)
#        | pairwise alignment to the model sequence       <- identity recorded
#        v
#   structure residue number (== UniProt position)
#        | coordinates, pLDDT, PAE, domains, key sites
#        v
#   report + 3D view

ens_host <- function(human_assembly = c("GRCh38", "GRCh37")) {
  human_assembly <- match.arg(human_assembly)
  if (human_assembly == "GRCh37") "https://grch37.rest.ensembl.org" else "https://rest.ensembl.org"
}

#' Fetch a protein sequence from Ensembl by protein ID
#' @param protein_id Ensembl protein ID
#' @param host Ensembl REST host
#' @return character(1) amino-acid sequence
#' @export
ens_protein_seq <- function(protein_id, host = ens_host()) {
  z <- http_json(paste0(host, "/sequence/id/", protein_id),
                 query = list(type = "protein", `content-type` = "application/json"))
  z$seq
}

.ungap <- function(x) gsub("-", "", x, fixed = TRUE)

#' Resolve which protein structure to draw for a gene in a species
#'
#' @param gene gene symbol (human nomenclature; the ortholog is looked up)
#' @param species Ensembl production name, or `NULL` for human
#' @param structure_species `"auto"` (use the species' own AlphaFold model when
#'   one exists, otherwise fall back to the human ortholog's model), `"query"`
#'   (require the species' own model) or `"human"` (always use the human model)
#' @param human_assembly `"GRCh38"` or `"GRCh37"`
#' @return list: `accession`, `human_accession`, `species_accession`,
#'   `anchor` (`"query"` or `"human"`), `ortholog` (from [xsmut::get_ortholog()]
#'   or `NULL` for human input), `entry`
#' @export
resolve_structure_target <- function(gene, species = NULL,
                                     structure_species = c("auto", "query", "human"),
                                     human_assembly = c("GRCh38", "GRCh37")) {
  structure_species <- match.arg(structure_species)
  host <- ens_host(match.arg(human_assembly))
  human_like <- is.null(species) || species %in% c("homo_sapiens", "human")

  orth <- NULL
  sp_pid <- NA_character_; hu_pid <- NA_character_
  if (!human_like) {
    orth <- xsmut::get_ortholog(gene, species)
    sp_pid <- orth$target$protein_id
    hu_pid <- orth$human$protein_id
  } else {
    g <- http_json(paste0(host, "/lookup/symbol/homo_sapiens/", gene),
                   query = list(expand = 1, `content-type` = "application/json"))
    canon <- Filter(function(x) isTRUE(as.logical(x$is_canonical)), g$Transcript)
    tx <- if (length(canon)) canon[[1]] else g$Transcript[[1]]
    hu_pid <- tx$Translation$id
  }

  hu_acc <- uniprot_accession(gene, "homo_sapiens", protein_id = hu_pid)
  sp_acc <- if (!human_like) uniprot_accession(gene, species, protein_id = sp_pid) else hu_acc

  pick <- switch(structure_species,
    human = hu_acc,
    query = sp_acc,
    auto  = {
      if (!is.null(sp_acc) && !is.null(af_query(sp_acc$accession))) sp_acc else {
        if (!human_like) .msg(sprintf(
          "No AlphaFold model for the %s protein; falling back to the human ortholog's model (%s).",
          species, hu_acc$accession))
        hu_acc
      }
    })
  if (is.null(pick)) stop("Could not resolve a UniProt accession for ", gene,
                          if (!human_like) paste0(" in ", species) else "")
  if (structure_species == "query" && is.null(af_query(pick$accession)))
    stop("structure_species = 'query' but AlphaFold has no model for ", pick$accession)

  anchor <- if (!is.null(sp_acc) && identical(pick$accession, sp_acc$accession) && !human_like) "query" else "human"
  list(accession = pick$accession,
       human_accession = if (!is.null(hu_acc)) hu_acc$accession else NA_character_,
       species_accession = if (!is.null(sp_acc)) sp_acc$accession else NA_character_,
       anchor = anchor, ortholog = orth, entry = pick,
       human_protein_id = hu_pid, species_protein_id = sp_pid)
}

#' Map mutation calls onto the residue numbering of a structure
#'
#' Accepts a dNdScv `annotmuts` data frame (auto-detected and converted with
#' [xsmut::as_xsmut()]) or any table with `gene` plus `aa_change` or `aa_pos`.
#' Each distinct source protein (`species_protein_id`, i.e. the dNdScv `pid`) is
#' aligned to the structure's sequence once; the identity of that alignment is
#' recorded per row so that weak mappings are visible.
#'
#' @param mutations data.frame of calls
#' @param gene gene symbol to subset to
#' @param model an `af_model`
#' @param species Ensembl production name of the calls (`NULL` = human)
#' @param fallback_protein_id protein ID to use for rows lacking one
#' @param min_identity warn when a source protein aligns to the structure below
#'   this percent identity
#' @return data.table with the input columns plus `struct_pos`, `struct_aa`,
#'   `aa_ref`, `aa_alt`, `ref_match`, `align_identity`, `map_note`
#' @export
map_mutations_to_structure <- function(mutations, gene, model, species = NULL,
                                       fallback_protein_id = NULL,
                                       min_identity = 60) {
  stopifnot(inherits(model, "af_model"))
  mt <- as.data.table(mutations)
  if (xsmut::is_dndscv(mt)) mt <- xsmut::as_xsmut(mt)
  if (!"gene" %in% names(mt)) stop("`mutations` needs a `gene` column")
  mt <- mt[toupper(gene) == toupper(mt$gene)]
  if (!nrow(mt)) {
    warning("No rows for gene ", gene, " in the mutation table")
  }
  if (!"aa_pos" %in% names(mt)) mt[, aa_pos := xsmut::parse_aa_pos(aa_change)]
  mt[, aa_pos := as.integer(aa_pos)]
  sub <- parse_aa_sub(mt$aa_change)
  mt[, `:=`(aa_ref = sub$ref, aa_alt = sub$alt)]

  if (!"species_protein_id" %in% names(mt)) mt[, species_protein_id := NA_character_]
  if (!is.null(fallback_protein_id)) {
    mt[is.na(species_protein_id) | !nzchar(species_protein_id),
       species_protein_id := fallback_protein_id]
  }

  anchor_seq <- model$sequence
  offset <- model$uniprot_start - 1L

  pids <- unique(stats::na.omit(mt$species_protein_id))
  maps <- list(); ids <- numeric()
  for (pid in pids) {
    s <- tryCatch(ens_protein_seq(pid), error = function(e) NULL)
    if (is.null(s)) { maps[[pid]] <- NULL; next }
    m <- residue_map(s, anchor_seq)
    maps[[pid]] <- m
    ids[pid] <- attr(m, "pid")
    if (!is.na(ids[pid]) && ids[pid] < min_identity) {
      warning(sprintf("%s aligns to %s at only %.1f%% identity; residue mapping is unreliable.",
                      pid, model$entry_id, ids[pid]))
    }
  }
  # rows without a protein ID: align the species ortholog protein once
  if (any(is.na(mt$species_protein_id)) && !is.null(species)) {
    orth <- tryCatch(xsmut::get_ortholog(gene, species), error = function(e) NULL)
    if (!is.null(orth)) {
      s <- .ungap(orth$target$seq)
      maps[["__ortholog__"]] <- residue_map(s, anchor_seq)
      ids["__ortholog__"] <- attr(maps[["__ortholog__"]], "pid")
      mt[is.na(species_protein_id), species_protein_id := "__ortholog__"]
    }
  }
  # human input with no protein ID: the anchor is already the right numbering
  mt[is.na(species_protein_id), species_protein_id := "__identity__"]
  if ("__identity__" %in% mt$species_protein_id) {
    maps[["__identity__"]] <- seq_len(nchar(anchor_seq))
    ids["__identity__"] <- 100
  }

  mt[, `:=`(struct_pos = NA_integer_, align_identity = NA_real_, map_note = NA_character_)]
  for (pid in unique(mt$species_protein_id)) {
    i <- which(mt$species_protein_id == pid)
    m <- maps[[pid]]
    if (is.null(m)) { mt[i, map_note := paste("protein sequence unavailable:", pid)]; next }
    mt[i, struct_pos := .apply_map(aa_pos, m) + offset]
    mt[i, align_identity := unname(ids[pid])]
    mt[i, map_note := data.table::fcase(
      !is.na(struct_pos), sprintf("aligned to %s", model$entry_id),
      is.na(aa_pos),      "no protein position (non-coding / splice)",
      aa_pos > length(m), "position beyond the source protein",
      default = "residue aligns to a gap in the model sequence")]
  }

  res <- model$sequence
  mt[, struct_aa := NA_character_]
  ok <- !is.na(mt$struct_pos) & mt$struct_pos >= model$uniprot_start & mt$struct_pos <= model$uniprot_end
  mt[ok, struct_aa := substring(res, struct_pos - offset, struct_pos - offset)]
  mt[, ref_match := !is.na(aa_ref) & !is.na(struct_aa) & aa_ref == struct_aa]
  mt[!is.na(aa_ref) & !is.na(struct_aa) & !ref_match,
     map_note := paste0(map_note, "; reference residue differs (", aa_ref, "->", struct_aa, ")")]
  mt[]
}

#' Collapse mutations to one row per residue
#'
#' @param mapped output of [map_mutations_to_structure()]
#' @param by extra grouping columns to keep (default none)
#' @return data.table: `struct_pos`, `n_mutations`, `n_samples`, `aa_changes`,
#'   `consequences`, `top_consequence`
#' @export
collapse_by_residue <- function(mapped, by = character()) {
  m <- as.data.table(mapped)[!is.na(struct_pos)]
  if (!nrow(m)) return(data.table(struct_pos = integer(), n_mutations = integer(),
                                  n_samples = integer(), aa_changes = character(),
                                  consequences = character(), top_consequence = character()))
  if (!"sample_id" %in% names(m)) m[, sample_id := NA_character_]
  if (!"consequence" %in% names(m)) m[, consequence := NA_character_]
  m[, .(n_mutations = .N,
        n_samples = data.table::uniqueN(sample_id),
        aa_changes = paste(sort(unique(stats::na.omit(aa_change))), collapse = ", "),
        consequences = paste(sort(unique(stats::na.omit(consequence))), collapse = ", "),
        top_consequence = names(sort(table(consequence), decreasing = TRUE))[1] %||% NA_character_),
    by = c("struct_pos", by)][order(struct_pos)]
}

#' Place human COSMIC mutations on the structure
#'
#' @param cosmic normalised COSMIC table (see [xsmut::cosmic_gene_mutations()])
#' @param model an `af_model`
#' @param human_accession UniProt accession of the human protein COSMIC
#'   positions refer to; when it differs from the model's accession the two
#'   sequences are aligned
#' @return data.table: `struct_pos`, `n_samples`, `aa_change`, `consequence`
#' @export
cosmic_on_structure <- function(cosmic, model, human_accession = NULL) {
  cs <- as.data.table(cosmic)
  if (!nrow(cs)) return(data.table(struct_pos = integer(), n_samples = integer(),
                                   aa_change = character(), consequence = character()))
  if (!"aa_pos" %in% names(cs)) cs[, aa_pos := xsmut::parse_aa_pos(aa_change)]
  cs <- cs[!is.na(aa_pos)]
  map <- if (!is.null(human_accession) && !identical(human_accession, model$accession)) {
    residue_map(uniprot_entry(human_accession)$sequence, model$sequence)
  } else seq_len(nchar(model$sequence))
  cs[, struct_pos := .apply_map(aa_pos, map) + (model$uniprot_start - 1L)]
  cs <- cs[!is.na(struct_pos)]
  cs[, .(n_samples = .N), by = .(struct_pos, aa_change, consequence)][order(struct_pos)]
}
