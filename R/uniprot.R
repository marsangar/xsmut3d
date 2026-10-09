# Resolving a gene (in any species) to a UniProt accession, and fetching the
# functional annotation that makes a structure interpretable.
#
# Two independent routes, tried in order:
#   1. Ensembl xrefs for the canonical translation (exact, isoform-aware)
#   2. UniProt REST search on gene symbol + organism, reviewed entries first
#
# Both are cheap and cached for the session.

.uniprot_base  <- "https://rest.uniprot.org"
.ebi_proteins  <- "https://www.ebi.ac.uk/proteins/api"
.interpro_base <- "https://www.ebi.ac.uk/interpro/api"

# --- 1. via Ensembl ---------------------------------------------------------

.uniprot_from_ensembl <- function(protein_id, host = "https://rest.ensembl.org") {
  dbs <- c("UniProt/SWISSPROT", "UniProt/SPTREMBL")
  for (db in dbs) {
    x <- tryCatch(
      http_json(paste0(host, "/xrefs/id/", protein_id),
                query = list(external_db = db, `content-type` = "application/json")),
      error = function(e) NULL
    )
    if (length(x)) {
      acc <- vapply(x, function(z) as.character(z$primary_id %||% NA), "")
      acc <- acc[nzchar(acc) & !is.na(acc)]
      if (length(acc)) return(list(accession = acc[1], reviewed = db == "UniProt/SWISSPROT",
                                   source = paste0("Ensembl xref (", db, ")")))
    }
  }
  NULL
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- 2. via UniProt search --------------------------------------------------

.uniprot_search <- function(gene, taxon_id = NULL, organism = NULL) {
  q <- sprintf("(gene:%s)", gene)
  if (!is.null(taxon_id)) q <- paste0(q, sprintf(" AND (taxonomy_id:%s)", taxon_id))
  else if (!is.null(organism)) q <- paste0(q, sprintf(" AND (organism_name:\"%s\")", organism))
  res <- tryCatch(
    http_json(paste0(.uniprot_base, "/uniprotkb/search"),
              query = list(query = q, format = "json", size = 25,
                           fields = "accession,id,reviewed,protein_name,gene_names,length,organism_id")),
    error = function(e) NULL
  )
  if (is.null(res) || !length(res$results)) return(NULL)
  r <- res$results
  rev <- vapply(r, function(z) identical(z$entryType, "UniProtKB reviewed (Swiss-Prot)"), TRUE)
  len <- vapply(r, function(z) as.integer(z$sequence$length %||% z$length %||% NA), 1L)
  ord <- order(!rev, -len)
  z <- r[[ord[1]]]
  list(accession = z$primaryAccession, reviewed = rev[ord[1]], source = "UniProt search")
}

#' Resolve a gene to a UniProt accession for a given species
#'
#' @param gene gene symbol
#' @param species Ensembl production name (e.g. `"macaca_mulatta"`), or `NULL`
#'   to search UniProt without a species restriction
#' @param protein_id optional Ensembl protein ID (`ENSP...`, `ENSMMUP...`). When
#'   given this is tried first and gives the isoform-exact accession.
#' @param taxon_id optional NCBI taxonomy id, used for the UniProt search route
#' @return list with `accession`, `reviewed`, `source`, or `NULL` if nothing found
#' @export
uniprot_accession <- function(gene, species = NULL, protein_id = NULL, taxon_id = NULL) {
  if (!is.null(protein_id) && !is.na(protein_id) && nzchar(protein_id)) {
    hit <- .uniprot_from_ensembl(protein_id)
    if (!is.null(hit)) return(hit)
  }
  organism <- if (!is.null(species)) gsub("_", " ", species) else NULL
  .uniprot_search(gene, taxon_id = taxon_id, organism = organism)
}

#' Fetch the UniProt sequence and metadata for an accession
#'
#' @param accession UniProt accession
#' @return list with `accession`, `id`, `name`, `gene`, `organism`, `taxon_id`,
#'   `length`, `sequence`, `reviewed`
#' @export
uniprot_entry <- function(accession) {
  z <- http_json(sprintf("%s/uniprotkb/%s", .uniprot_base, accession),
                 query = list(format = "json"))
  list(
    accession = z$primaryAccession,
    id        = z$uniProtkbId,
    name      = z$proteinDescription$recommendedName$fullName$value %||%
                 z$proteinDescription$submissionNames[[1]]$fullName$value %||% NA_character_,
    gene      = z$genes[[1]]$geneName$value %||% NA_character_,
    organism  = z$organism$scientificName %||% NA_character_,
    taxon_id  = z$organism$taxonId %||% NA_integer_,
    length    = z$sequence$length,
    sequence  = z$sequence$value,
    reviewed  = identical(z$entryType, "UniProtKB reviewed (Swiss-Prot)")
  )
}

# ---------------------------------------------------------------------------
# Functional features
# ---------------------------------------------------------------------------

# Which UniProt feature types are worth drawing / flagging, and how to group them.
.feature_classes <- c(
  DOMAIN           = "Domain",
  REGION           = "Region",
  REPEAT           = "Repeat",
  ZN_FING          = "Zinc finger",
  DNA_BIND         = "DNA binding",
  COILED           = "Coiled coil",
  MOTIF            = "Motif",
  ACT_SITE         = "Active site",
  BINDING          = "Binding site",
  SITE             = "Site",
  METAL            = "Metal binding",
  DISULFID         = "Disulfide",
  MOD_RES          = "Modified residue",
  CARBOHYD         = "Glycosylation",
  CROSSLNK         = "Cross-link",
  LIPID            = "Lipidation",
  TRANSMEM         = "Transmembrane",
  SIGNAL           = "Signal peptide",
  PROPEP           = "Propeptide",
  MUTAGEN          = "Mutagenesis",
  VARIANT          = "Variant"
)

# Point features are single residues; the rest are ranges.
.point_features <- c("Active site", "Binding site", "Site", "Metal binding",
                     "Modified residue", "Glycosylation", "Cross-link",
                     "Lipidation", "Mutagenesis", "Variant")

#' Fetch UniProt functional features for an accession
#'
#' Uses the EBI Proteins API (`/features`), which exposes the same curated
#' feature table as the UniProt entry page, including active sites, binding
#' sites, DNA-binding regions, zinc fingers, disulfides and curated domains.
#'
#' @param accession UniProt accession
#' @param types optional character vector of UniProt feature type codes to keep
#'   (default: the interpretable subset listed in `xsmut3d:::.feature_classes`)
#' @return data.table: `start`, `end`, `type`, `class`, `description`, `point`,
#'   `evidence`, `source`
#' @export
uniprot_features <- function(accession, types = names(.feature_classes)) {
  z <- tryCatch(http_json(sprintf("%s/features/%s", .ebi_proteins, accession)),
                error = function(e) NULL)
  if (is.null(z) || !length(z$features)) {
    return(data.table(start = integer(), end = integer(), type = character(),
                      class = character(), description = character(),
                      point = logical(), evidence = character(), source = character()))
  }
  ft <- rbindlist(lapply(z$features, function(f) {
    data.table(
      start = suppressWarnings(as.integer(f$begin)),
      end   = suppressWarnings(as.integer(f$end)),
      type  = as.character(f$type),
      description = as.character(f$description %||% f$ftId %||% ""),
      evidence = if (length(f$evidences)) {
        paste(unique(vapply(f$evidences, function(e) as.character(e$code %||% ""), "")), collapse = ",")
      } else NA_character_
    )
  }), fill = TRUE)
  ft <- ft[type %in% types & !is.na(start) & !is.na(end)]
  ft[, class := unname(.feature_classes[type])]
  ft[, point := class %in% .point_features | start == end]
  ft[, source := "UniProt"]
  ft[order(start, end)]
}

#' Fetch InterPro domain assignments for an accession
#'
#' A useful cross-check / fallback when a species entry has little curated
#' UniProt annotation: InterPro signatures are computed, so they exist for
#' essentially every proteome.
#'
#' @param accession UniProt accession
#' @param databases signature databases to keep
#' @return data.table in the same schema as [uniprot_features()]
#' @export
interpro_domains <- function(accession, databases = c("pfam", "smart", "profile", "cdd")) {
  out <- lapply(databases, function(db) {
    z <- tryCatch(http_json(sprintf("%s/entry/%s/protein/uniprot/%s/",
                                    .interpro_base, db, accession),
                            query = list(page_size = 200)),
                  error = function(e) NULL)
    if (is.null(z) || !length(z$results)) return(NULL)
    rbindlist(lapply(z$results, function(r) {
      acc  <- r$metadata$accession
      nm   <- r$metadata$name
      locs <- r$proteins[[1]]$entry_protein_locations
      rbindlist(lapply(locs, function(l) rbindlist(lapply(l$fragments, function(fr)
        data.table(start = as.integer(fr$start), end = as.integer(fr$end),
                   type = toupper(db), description = sprintf("%s (%s)", nm, acc))))))
    }), fill = TRUE)
  })
  ft <- rbindlist(Filter(Negate(is.null), out), fill = TRUE)
  if (!nrow(ft)) {
    return(data.table(start = integer(), end = integer(), type = character(),
                      class = character(), description = character(),
                      point = logical(), evidence = character(), source = character()))
  }
  ft[, `:=`(class = "Domain", point = FALSE, evidence = NA_character_, source = "InterPro")]
  unique(ft)[order(start, end)]
}
