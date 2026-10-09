#' @importFrom r3dmol r3dmol m_add_model m_set_style m_style_cartoon m_style_sphere
#'   m_style_stick m_add_style m_sel m_zoom_to m_add_label m_style_label m_button_spin
NULL

# Interactive 3D view. Three colourings, each answering a different question:
#
#   "plddt"  - how much of this model should I believe?         (AlphaFold palette)
#   "domain" - where are the functional units?                  (categorical)
#   "burden" - where are my mutations?                          (grey -> red ramp)
#
# Mutated residues are always drawn as spheres scaled by recurrence, key
# functional sites as sticks, and low-confidence stretches as a thin translucent
# tube so they cannot be mistaken for structure.

.domain_palette <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
                     "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD")

.burden_colour <- function(n, max_n) {
  f <- if (max_n <= 1) rep(1, length(n)) else (n - 1) / (max_n - 1)
  grDevices::rgb(grDevices::colorRamp(c("#FDD0A2", "#E6550D", "#7F0000"))(pmax(0, pmin(1, f))),
                 maxColorValue = 255)
}

#' Interactive 3D view of mutations on a structure
#'
#' @param structure a `protein_structure` from [read_structure()]
#' @param residue_counts data.table with `struct_pos`, `n_mutations` (and
#'   optionally `aa_changes`), e.g. from [collapse_by_residue()]
#' @param domains table from [domain_table()]; used by `colour_by = "domain"`
#'   and to draw key sites
#' @param colour_by `"plddt"`, `"domain"` or `"burden"`
#' @param show_sites draw annotated key functional sites as sticks
#' @param label_top label the N most recurrently mutated residues
#' @param low_confidence_plddt residues below this pLDDT are drawn as a thin
#'   translucent tube
#' @param width,height widget size
#' @return an `r3dmol` htmlwidget
#' @export
view_structure_3d <- function(structure, residue_counts = NULL, domains = NULL,
                              colour_by = c("plddt", "domain", "burden"),
                              show_sites = TRUE, label_top = 5L,
                              low_confidence_plddt = 70,
                              width = NULL, height = NULL) {
  colour_by <- match.arg(colour_by)
  res <- structure$residues
  path <- structure$file %||% structure$model$file
  if (is.null(path) || !file.exists(path))
    stop("The structure object has no coordinate file to render; rebuild it with read_structure().")
  fmt  <- if (grepl("\\.cif(\\.gz)?$", path, ignore.case = TRUE)) "cif" else "pdb"

  v <- r3dmol::r3dmol(width = width, height = height)
  v <- r3dmol::m_add_model(v, data = path, format = fmt)
  v <- r3dmol::m_set_style(v, style = r3dmol::m_style_cartoon(color = "#D9D9D9", arrows = TRUE))

  # --- base colouring -------------------------------------------------------
  if (colour_by == "plddt") {
    for (b in names(plddt_palette)) {
      p <- res$pos[!is.na(res$plddt_band) & as.character(res$plddt_band) == b]
      if (!length(p)) next
      v <- r3dmol::m_add_style(v, sel = r3dmol::m_sel(resi = p),
                               style = r3dmol::m_style_cartoon(color = unname(plddt_palette[b]), arrows = TRUE))
    }
  } else if (colour_by == "domain" && !is.null(domains) && nrow(domains)) {
    d <- as.data.table(domains)[point == FALSE][order(start)]
    d <- d[!duplicated(label)]
    for (i in seq_len(nrow(d))) {
      col <- .domain_palette[((i - 1L) %% length(.domain_palette)) + 1L]
      v <- r3dmol::m_add_style(v, sel = r3dmol::m_sel(resi = seq.int(d$start[i], d$end[i])),
                               style = r3dmol::m_style_cartoon(color = col, arrows = TRUE))
    }
  }

  # --- low-confidence regions drawn thin and translucent --------------------
  lowp <- res$pos[!is.na(res$plddt) & res$plddt < low_confidence_plddt]
  if (length(lowp)) {
    v <- r3dmol::m_add_style(v, sel = r3dmol::m_sel(resi = lowp),
                             style = r3dmol::m_style_cartoon(
                               color = if (colour_by == "plddt") "#FF7D45" else "#BDBDBD",
                               style = "trace", thickness = 0.15, opacity = 0.55))
  }

  # --- key functional sites -------------------------------------------------
  if (isTRUE(show_sites) && !is.null(domains) && nrow(domains)) {
    ks <- key_sites(domains)
    if (nrow(ks)) {
      v <- r3dmol::m_add_style(v, sel = r3dmol::m_sel(resi = unique(ks$pos)),
                               style = r3dmol::m_style_stick(color = "#2B8C3E", radius = 0.18))
    }
  }

  # --- mutations ------------------------------------------------------------
  if (!is.null(residue_counts) && nrow(residue_counts)) {
    rc <- as.data.table(residue_counts)[!is.na(struct_pos)][order(-n_mutations)]
    mx <- max(rc$n_mutations)
    cols <- if (colour_by == "burden") .burden_colour(rc$n_mutations, mx) else rep("#C73E1D", nrow(rc))
    for (i in seq_len(nrow(rc))) {
      rad <- 0.6 + 0.9 * sqrt(rc$n_mutations[i] / mx)
      v <- r3dmol::m_add_style(v, sel = r3dmol::m_sel(resi = rc$struct_pos[i]),
                               style = r3dmol::m_style_sphere(color = cols[i], radius = rad))
    }
    if (label_top > 0) {
      top <- utils::head(rc, label_top)
      for (i in seq_len(nrow(top))) {
        lab <- if ("aa_changes" %in% names(top) && nzchar(top$aa_changes[i]))
          sprintf("%s (n=%d)", sub(",.*$", "", top$aa_changes[i]), top$n_mutations[i])
        else sprintf("%d (n=%d)", top$struct_pos[i], top$n_mutations[i])
        v <- r3dmol::m_add_label(v, text = lab,
                                 sel = r3dmol::m_sel(resi = top$struct_pos[i]),
                                 style = r3dmol::m_style_label(
                                   fontSize = 11, backgroundColor = "white",
                                   backgroundOpacity = 0.75, fontColor = "black"))
      }
    }
  }

  v <- r3dmol::m_zoom_to(v)
  v
}

#' Save an interactive view to a self-contained HTML file
#' @param view an `r3dmol` widget
#' @param file output path (.html)
#' @export
save_view_html <- function(view, file) {
  htmlwidgets::saveWidget(view, file, selfcontained = TRUE)
  invisible(normalizePath(file))
}

#' Save a static PNG snapshot of a 3D view (requires webshot2)
#' @param view an `r3dmol` widget
#' @param file output path (.png)
#' @param width,height pixels
#' @export
save_view_png <- function(view, file, width = 1200, height = 1000) {
  if (!requireNamespace("webshot2", quietly = TRUE))
    stop("Install webshot2 to export static 3D snapshots: install.packages('webshot2')")
  tmp <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(view, tmp, selfcontained = TRUE)
  webshot2::webshot(tmp, file, vwidth = width, vheight = height, delay = 2)
  invisible(normalizePath(file))
}
