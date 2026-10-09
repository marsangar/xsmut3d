#' @importFrom data.table data.table as.data.table rbindlist setnames fread setorder :=
#' @importFrom stats setNames p.adjust median quantile
NULL

# ---------------------------------------------------------------------------
# Tiny HTTP helper (same manners as xsmut: retry, polite UA, session memoise)
# ---------------------------------------------------------------------------

.ua <- "xsmut3d (https://github.com/marsangar/xsmut3d)"

.http_json_raw <- function(url, query = list(), accept = "application/json") {
  req <- httr2::request(url)
  if (length(query)) req <- httr2::req_url_query(req, !!!query)
  req <- httr2::req_headers(req, Accept = accept)
  req <- httr2::req_user_agent(req, .ua)
  req <- httr2::req_retry(req, max_tries = 4, backoff = ~ 2)
  resp <- httr2::req_perform(req)
  jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
}

#' @keywords internal
http_json <- memoise::memoise(.http_json_raw)

# Download a file to the cache directory, re-using it if already present.
.download_cached <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  req <- httr2::request(url)
  req <- httr2::req_user_agent(req, .ua)
  req <- httr2::req_retry(req, max_tries = 4, backoff = ~ 2)
  tmp <- paste0(dest, ".part")
  httr2::req_perform(req, path = tmp)
  file.rename(tmp, dest)
  dest
}

#' Directory used to cache AlphaFold models and API responses
#'
#' Defaults to `tools::R_user_dir("xsmut3d", "cache")`; override with the
#' `XSMUT3D_CACHE` environment variable or the `cache_dir` argument of
#' [af_model()].
#'
#' @return character(1) path (created if missing)
#' @export
xsmut3d_cache_dir <- function() {
  d <- Sys.getenv("XSMUT3D_CACHE", unset = "")
  if (!nzchar(d)) d <- tools::R_user_dir("xsmut3d", "cache")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# ---------------------------------------------------------------------------
# Amino-acid helpers
# ---------------------------------------------------------------------------

.aa3to1 <- c(
  ALA = "A", ARG = "R", ASN = "N", ASP = "D", CYS = "C", GLN = "Q", GLU = "E",
  GLY = "G", HIS = "H", ILE = "I", LEU = "L", LYS = "K", MET = "M", PHE = "F",
  PRO = "P", SER = "S", THR = "T", TRP = "W", TYR = "Y", VAL = "V",
  SEC = "U", PYL = "O", MSE = "M", UNK = "X"
)

aa3to1 <- function(x) {
  out <- unname(.aa3to1[toupper(x)])
  out[is.na(out)] <- "X"
  out
}

# "p.R175H" -> list(ref = "R", pos = 175, alt = "H"); "p.Q331*" -> alt "*"
#' Split an HGVS protein substitution into reference residue, position and alternate
#'
#' Non-substitutions (frameshifts, indels, splice) return `NA` for `ref`/`alt`
#' but still carry the position parsed by [xsmut::parse_aa_pos()].
#'
#' @param x character vector of HGVS protein strings
#' @return data.table with `ref`, `pos`, `alt`
#' @export
parse_aa_sub <- function(x) {
  s <- sub("^p\\.\\(?", "", as.character(x))
  s[is.na(s)] <- ""
  m <- regmatches(s, regexec("^([A-Z])(\\d+)([A-Z*=])$", s))
  ref <- vapply(m, function(z) if (length(z)) z[2] else NA_character_, "")
  alt <- vapply(m, function(z) if (length(z)) z[4] else NA_character_, "")
  syn <- !is.na(alt) & alt == "="
  alt[syn] <- ref[syn]
  data.table(ref = ref, pos = xsmut::parse_aa_pos(x), alt = alt)
}

# ---------------------------------------------------------------------------
# Pairwise protein alignment and residue maps (same approach as xsmut/crossmap)
# ---------------------------------------------------------------------------

.align_proteins <- function(seq_a, seq_b) {
  # pwalign (Bioc >= 3.19) or the older Biostrings home of pairwiseAlignment
  pa <- if (requireNamespace("pwalign", quietly = TRUE)) asNamespace("pwalign") else asNamespace("Biostrings")
  aln <- pa$pairwiseAlignment(
    Biostrings::AAString(gsub("[^A-Z]", "", toupper(seq_b))),
    Biostrings::AAString(gsub("[^A-Z]", "", toupper(seq_a))),
    substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 0.5,
    type = "global"
  )
  list(a = as.character(pa$alignedSubject(aln)),
       b = as.character(pa$alignedPattern(aln)),
       pid = tryCatch(pa$pid(aln), error = function(e) NA_real_))
}

# position in ungapped B -> position in ungapped A
.aln_map <- function(gapped_a, gapped_b) {
  a <- strsplit(gapped_a, "")[[1]]
  b <- strsplit(gapped_b, "")[[1]]
  stopifnot(length(a) == length(b))
  ia <- cumsum(a != "-"); ib <- cumsum(b != "-")
  keep <- b != "-"
  map <- rep(NA_integer_, sum(keep))
  map[ib[keep]] <- ifelse(a[keep] != "-", ia[keep], NA_integer_)
  map
}

#' Build a residue map between two protein sequences
#'
#' @param from character(1) sequence whose positions you have
#' @param to character(1) sequence whose positions you want
#' @return integer vector of length `nchar(from)`; element i is the position in
#'   `to` aligned to position i of `from`, or `NA` if it aligns to a gap.
#'   The attribute `"pid"` holds the percent identity of the alignment.
#' @export
residue_map <- function(from, to) {
  if (identical(gsub("[^A-Z]", "", toupper(from)), gsub("[^A-Z]", "", toupper(to)))) {
    m <- seq_len(nchar(from)); attr(m, "pid") <- 100; return(m)
  }
  aln <- .align_proteins(to, from)
  m <- .aln_map(aln$a, aln$b)
  attr(m, "pid") <- aln$pid
  m
}

# apply a residue map to a vector of positions
.apply_map <- function(pos, map) {
  out <- rep(NA_integer_, length(pos))
  ok <- !is.na(pos) & pos >= 1L & pos <= length(map)
  out[ok] <- map[pos[ok]]
  out
}

.msg <- function(...) if (!isTRUE(getOption("xsmut3d.quiet", FALSE))) message(...)
