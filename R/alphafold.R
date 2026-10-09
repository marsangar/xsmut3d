# AlphaFold Protein Structure Database (https://alphafold.ebi.ac.uk).
#
# One REST call gives the model URLs plus the global confidence summary; the
# mmCIF (or PDB) file carries per-residue pLDDT in the B-factor column, and the
# companion JSON files carry the predicted aligned error (PAE) matrix and the
# AlphaMissense substitution scores.
#
# Everything is cached on disk (see xsmut3d_cache_dir()) so repeat runs and
# multi-gene loops are offline after the first fetch.

.af_api <- "https://alphafold.ebi.ac.uk/api/prediction"

#' Query the AlphaFold DB API for a UniProt accession
#'
#' @param accession UniProt accession (e.g. `"P04637"`)
#' @return list of prediction entries (one per fragment), or `NULL` if the
#'   accession has no AlphaFold model
#' @export
af_query <- function(accession) {
  z <- tryCatch(http_json(paste0(.af_api, "/", accession)), error = function(e) NULL)
  if (is.null(z) || !length(z)) return(NULL)
  z
}

#' Fetch an AlphaFold model (structure + confidence) for a UniProt accession
#'
#' Downloads the mmCIF coordinate file, the predicted aligned error matrix and,
#' when `alphamissense = TRUE` and the entry has one, the AlphaMissense
#' substitution table. All files are cached under [xsmut3d_cache_dir()].
#'
#' @param accession UniProt accession
#' @param fragment which model fragment to use when the protein is longer than
#'   the AlphaFold fragment limit (~2700 aa). Default `1`.
#' @param format `"cif"` (default) or `"pdb"`
#' @param pae download and parse the PAE matrix (used for domain-level
#'   confidence). Can be slow for very long proteins.
#' @param alphamissense download AlphaMissense pathogenicity scores if available
#' @param cache_dir directory for cached downloads
#' @return list of class `af_model`:
#'   `entry_id`, `accession`, `gene`, `organism`, `version`, `sequence`,
#'   `file` (path to the coordinate file), `format`, `n_fragments`,
#'   `offset` (residue number of the first residue of the fragment minus one),
#'   `global_plddt`, `plddt_fractions`, `pae` (matrix or NULL),
#'   `alphamissense` (data.table or NULL), `url`
#' @export
af_model <- function(accession, fragment = 1L, format = c("cif", "pdb"),
                     pae = TRUE, alphamissense = TRUE,
                     cache_dir = xsmut3d_cache_dir()) {
  format <- match.arg(format)
  z <- af_query(accession)
  if (is.null(z)) stop("No AlphaFold model for UniProt accession ", accession,
                       ". Try structure_species = 'human' to use the ortholog's model.")
  n_frag <- length(z)
  if (fragment > n_frag) stop("Only ", n_frag, " fragment(s) for ", accession)
  e <- z[[fragment]]
  if (n_frag > 1L) {
    .msg(sprintf("AlphaFold has %d fragments for %s; using %s. Residue numbering is offset.",
                 n_frag, accession, e$entryId))
  }

  url  <- if (format == "cif") e$cifUrl else e$pdbUrl
  file <- .download_cached(url, file.path(cache_dir, "models", basename(url)))

  pae_mat <- NULL
  if (isTRUE(pae) && !is.null(e$paeDocUrl)) {
    pf <- tryCatch(.download_cached(e$paeDocUrl, file.path(cache_dir, "pae", basename(e$paeDocUrl))),
                   error = function(err) NULL)
    if (!is.null(pf)) pae_mat <- .parse_pae(pf)
  }

  am <- NULL
  if (isTRUE(alphamissense) && !is.null(e$amAnnotationsUrl)) {
    af <- tryCatch(.download_cached(e$amAnnotationsUrl,
                                    file.path(cache_dir, "am", basename(e$amAnnotationsUrl))),
                   error = function(err) NULL)
    if (!is.null(af)) am <- .parse_alphamissense(af)
  }

  structure(list(
    entry_id     = e$entryId,
    accession    = e$uniprotAccession %||% accession,
    gene         = e$gene %||% NA_character_,
    organism     = e$organismScientificName %||% NA_character_,
    version      = e$latestVersion %||% NA,
    sequence     = e$uniprotSequence %||% NA_character_,
    uniprot_start = as.integer(e$uniprotStart %||% 1L),
    uniprot_end   = as.integer(e$uniprotEnd %||% nchar(e$uniprotSequence %||% "")),
    file         = file,
    format       = format,
    n_fragments  = n_frag,
    fragment     = fragment,
    global_plddt = as.numeric(e$globalMetricValue %||% NA),
    plddt_fractions = c(
      very_high = as.numeric(e$fractionPlddtVeryHigh %||% NA),
      confident = as.numeric(e$fractionPlddtConfident %||% NA),
      low       = as.numeric(e$fractionPlddtLow %||% NA),
      very_low  = as.numeric(e$fractionPlddtVeryLow %||% NA)
    ),
    pae           = pae_mat,
    alphamissense = am,
    url           = url
  ), class = "af_model")
}

#' @export
print.af_model <- function(x, ...) {
  cat(sprintf("<af_model> %s  %s (%s)  %s\n", x$entry_id, x$accession,
              x$gene, x$organism))
  cat(sprintf("  %d residues (UniProt %d-%d), model v%s, fragment %d/%d\n",
              nchar(x$sequence), x$uniprot_start, x$uniprot_end,
              as.character(x$version), x$fragment, x$n_fragments))
  cat(sprintf("  mean pLDDT %.1f  [very high %.0f%% | confident %.0f%% | low %.0f%% | very low %.0f%%]\n",
              x$global_plddt, 100 * x$plddt_fractions[["very_high"]],
              100 * x$plddt_fractions[["confident"]], 100 * x$plddt_fractions[["low"]],
              100 * x$plddt_fractions[["very_low"]]))
  cat(sprintf("  PAE: %s   AlphaMissense: %s\n",
              if (is.null(x$pae)) "not loaded" else sprintf("%d x %d", nrow(x$pae), ncol(x$pae)),
              if (is.null(x$alphamissense)) "not available" else sprintf("%d substitutions", nrow(x$alphamissense))))
  invisible(x)
}

# --- companion file parsers -------------------------------------------------

.parse_pae <- function(path) {
  z <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  if (is.data.frame(z)) z <- as.list(z)
  if (!is.null(z$predicted_aligned_error)) {
    m <- z$predicted_aligned_error
    if (is.list(m)) m <- m[[1]]
    return(as.matrix(if (is.list(m)) do.call(rbind, m) else m))
  }
  # legacy "residue1 / residue2 / distance" long form
  if (!is.null(z$residue1)) {
    r1 <- unlist(z$residue1); r2 <- unlist(z$residue2); d <- unlist(z$distance)
    n <- max(r1)
    m <- matrix(NA_real_, n, n)
    m[cbind(r1, r2)] <- d
    return(m)
  }
  NULL
}

.parse_alphamissense <- function(path) {
  # AlphaFold ships AlphaMissense as a small CSV:
  #   protein_variant,am_pathogenicity,am_class
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  nm <- tolower(names(dt))
  v  <- names(dt)[grep("variant", nm)][1]
  p  <- names(dt)[grep("patho", nm)][1]
  cl <- names(dt)[grep("class", nm)][1]
  if (is.na(v) || is.na(p)) return(NULL)
  out <- data.table(variant = as.character(dt[[v]]),
                    am_pathogenicity = as.numeric(dt[[p]]),
                    am_class = if (!is.na(cl)) as.character(dt[[cl]]) else NA_character_)
  sub <- parse_aa_sub(out$variant)
  out[, `:=`(ref = sub$ref, uniprot_pos = sub$pos, alt = sub$alt)]
  out[!is.na(uniprot_pos)]
}

# ---------------------------------------------------------------------------
# pLDDT banding
# ---------------------------------------------------------------------------

#' AlphaFold pLDDT confidence bands
#'
#' The standard AlphaFold DB bands. Regions below 70 should not be interpreted
#' as reliable local structure; regions below 50 are frequently intrinsically
#' disordered and their coordinates are close to meaningless.
#'
#' @param plddt numeric vector
#' @return ordered factor with levels
#'   `"Very low (<50)"`, `"Low (50-70)"`, `"Confident (70-90)"`, `"Very high (>90)"`
#' @export
plddt_band <- function(plddt) {
  lv <- c("Very low (<50)", "Low (50-70)", "Confident (70-90)", "Very high (>90)")
  factor(data.table::fcase(
    is.na(plddt),   NA_character_,
    plddt >= 90,    lv[4],
    plddt >= 70,    lv[3],
    plddt >= 50,    lv[2],
    default = lv[1]
  ), levels = lv, ordered = TRUE)
}

#' Colours for the pLDDT bands (the AlphaFold DB palette)
#' @export
plddt_palette <- c(
  "Very high (>90)"   = "#0053D6",
  "Confident (70-90)" = "#65CBF3",
  "Low (50-70)"       = "#FFDB13",
  "Very low (<50)"    = "#FF7D45"
)

#' Collapse a PAE matrix into per-residue domain confidence
#'
#' For each residue, the median PAE to all other residues within `window`
#' sequence separation is a cheap proxy for how well-packed the local unit is,
#' while the median PAE to the rest of the chain tells you whether that unit is
#' confidently placed *relative to* the rest of the protein. Pairs of residues
#' with high inter-PAE should not be compared by their 3D distance.
#'
#' @param pae PAE matrix from [af_model()]
#' @param window sequence separation defining "local"
#' @return data.table: `pos`, `pae_local`, `pae_global`
#' @export
pae_profile <- function(pae, window = 12L) {
  if (is.null(pae)) return(NULL)
  n <- nrow(pae)
  i <- seq_len(n)
  loc <- vapply(i, function(k) {
    j <- max(1, k - window):min(n, k + window)
    stats::median(pae[k, j], na.rm = TRUE)
  }, 0)
  glo <- vapply(i, function(k) stats::median(pae[k, ], na.rm = TRUE), 0)
  data.table(pos = i, pae_local = loc, pae_global = glo)
}
