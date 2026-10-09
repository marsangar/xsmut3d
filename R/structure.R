# Minimal, dependency-free readers for the two file formats AlphaFold serves,
# plus the geometry needed to talk about residues in 3D: representative atom
# coordinates, burial (neighbour count), and contact neighbourhoods.
#
# AlphaFold models are single-chain (A), with residue numbers equal to UniProt
# positions (offset by uniprot_start for multi-fragment entries) and pLDDT in
# the B-factor column.

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

.read_pdb_atoms <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[substr(ln, 1, 4) == "ATOM"]
  data.table(
    atom     = trimws(substr(ln, 13, 16)),
    res_name = trimws(substr(ln, 18, 20)),
    chain    = trimws(substr(ln, 22, 22)),
    res_num  = as.integer(substr(ln, 23, 26)),
    x        = as.numeric(substr(ln, 31, 38)),
    y        = as.numeric(substr(ln, 39, 46)),
    z        = as.numeric(substr(ln, 47, 54)),
    bfactor  = as.numeric(substr(ln, 61, 66))
  )
}

.read_cif_atoms <- function(path) {
  ln <- readLines(path, warn = FALSE)
  # locate the _atom_site loop and its column order
  hdr_idx <- grep("^_atom_site\\.", ln)
  if (!length(hdr_idx)) stop("No _atom_site loop found in ", path)
  cols <- sub("^_atom_site\\.", "", ln[hdr_idx])
  first <- max(hdr_idx) + 1L
  last  <- first
  while (last <= length(ln) && !grepl("^(#|loop_|_)", ln[last])) last <- last + 1L
  body <- ln[first:(last - 1L)]
  body <- body[nzchar(trimws(body))]
  m <- do.call(rbind, strsplit(trimws(body), "[[:space:]]+"))
  colnames(m) <- cols
  if ("group_PDB" %in% cols) m <- m[m[, "group_PDB"] == "ATOM", , drop = FALSE]
  get <- function(nm, alt = NULL) {
    if (nm %in% cols) m[, nm] else if (!is.null(alt) && alt %in% cols) m[, alt] else NA
  }
  data.table(
    atom     = get("label_atom_id"),
    res_name = get("label_comp_id"),
    chain    = get("auth_asym_id", "label_asym_id"),
    res_num  = as.integer(get("label_seq_id", "auth_seq_id")),
    x        = as.numeric(get("Cartn_x")),
    y        = as.numeric(get("Cartn_y")),
    z        = as.numeric(get("Cartn_z")),
    bfactor  = as.numeric(get("B_iso_or_equiv"))
  )
}

#' Read an AlphaFold coordinate file into a residue table
#'
#' @param model an `af_model` from [af_model()], or a path to a .cif/.pdb file
#' @param contact_cutoff distance (Angstrom) between representative atoms that
#'   counts as a contact, used for burial and 3D neighbourhoods
#' @return list of class `protein_structure`:
#'   `residues` (data.table: `pos`, `aa`, `res_name`, `x`, `y`, `z`, `plddt`,
#'   `plddt_band`, `n_neighbours`, `burial`), `atoms`, `dist` (residue-residue
#'   distance matrix), `contact_cutoff`, `model`
#' @export
read_structure <- function(model, contact_cutoff = 10) {
  path <- if (inherits(model, "af_model")) model$file else model
  fmt  <- if (grepl("\\.cif(\\.gz)?$", path, ignore.case = TRUE)) "cif" else "pdb"
  atoms <- if (fmt == "cif") .read_cif_atoms(path) else .read_pdb_atoms(path)
  atoms <- atoms[!is.na(res_num)]
  if (!nrow(atoms)) stop("No ATOM records parsed from ", path)
  chains <- unique(atoms$chain)
  if (length(chains) > 1L) atoms <- atoms[chain == chains[1]]

  # representative atom per residue: CB, or CA for glycine / missing CB
  rep_at <- atoms[atom %in% c("CA", "CB")]
  rep_at <- rep_at[order(res_num, atom != "CB")]      # CB first
  rep_at <- rep_at[!duplicated(res_num)]

  ca <- atoms[atom == "CA"][!duplicated(res_num)]
  res <- data.table(
    pos      = ca$res_num,
    res_name = ca$res_name,
    aa       = aa3to1(ca$res_name),
    plddt    = ca$bfactor
  )
  res <- merge(res, rep_at[, .(pos = res_num, x, y, z)], by = "pos", all.x = TRUE)
  setorder(res, pos)

  # offset residue numbering for fragments beyond the first
  if (inherits(model, "af_model") && !is.na(model$uniprot_start) && model$uniprot_start > 1L) {
    res[, pos := pos + (model$uniprot_start - 1L)]
  }

  xyz <- as.matrix(res[, .(x, y, z)])
  d <- as.matrix(stats::dist(xyz))
  dimnames(d) <- list(res$pos, res$pos)

  res[, n_neighbours := colSums(d <= contact_cutoff, na.rm = TRUE) - 1L]
  res[, burial := data.table::fcase(
    n_neighbours >= stats::quantile(n_neighbours, 0.66, na.rm = TRUE), "buried",
    n_neighbours >= stats::quantile(n_neighbours, 0.33, na.rm = TRUE), "intermediate",
    default = "exposed")]
  res[, plddt_band := plddt_band(plddt)]

  structure(list(
    residues = res[], atoms = atoms, dist = d,
    contact_cutoff = contact_cutoff, file = path,
    model = if (inherits(model, "af_model")) model else NULL
  ), class = "protein_structure")
}

#' @export
print.protein_structure <- function(x, ...) {
  r <- x$residues
  cat(sprintf("<protein_structure> %d residues (%d-%d), mean pLDDT %.1f\n",
              nrow(r), min(r$pos), max(r$pos), mean(r$plddt, na.rm = TRUE)))
  tb <- table(r$plddt_band)
  cat("  ", paste(sprintf("%s: %d", names(tb), as.integer(tb)), collapse = " | "), "\n", sep = "")
  invisible(x)
}

#' Residues within a distance of a given residue
#'
#' @param structure a `protein_structure`
#' @param pos residue position(s)
#' @param radius Angstrom
#' @return integer vector of neighbouring residue positions (excluding `pos`)
#' @export
residues_near <- function(structure, pos, radius = 10) {
  d <- structure$dist
  idx <- match(as.character(pos), colnames(d))
  idx <- idx[!is.na(idx)]
  if (!length(idx)) return(integer())
  hit <- which(apply(d[idx, , drop = FALSE] <= radius, 2, any))
  setdiff(as.integer(colnames(d)[hit]), pos)
}
