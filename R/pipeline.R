# The one call most users will make.
#
# Deliberately mirrors xsmut::xsmut_pipeline() argument-for-argument, so a lab
# that already runs xsmut can swap the function name and get the structural view
# of the same input.

#' Project somatic mutations onto a predicted protein structure
#'
#' @param gene gene symbol, human nomenclature (the ortholog is resolved for you).
#'   May be a character vector, in which case a named list of results is returned.
#' @param mutations a `dndscv()$annotmuts` data frame (auto-detected) or any
#'   table with `gene` plus `aa_change` or `aa_pos`
#' @param species species name in any form accepted by
#'   [xsmut::ens_resolve_species()] (`"macaque"`, `"Macaca mulatta"`,
#'   `"macaca_mulatta"`); `NULL` or `"human"` for human calls
#' @param assembly expected species assembly (e.g. `"Mmul_10"`); checked against
#'   what Ensembl serves and recorded in the output
#' @param structure_species which proteome's AlphaFold model to draw on:
#'   `"auto"` (the query species' own model when one exists, else the human
#'   ortholog's), `"query"` or `"human"`
#' @param cosmic optional COSMIC mutation table or file path (see
#'   [xsmut::cosmic_load()]); adds the human COSMIC track and hotspot transfer
#' @param census_path optional COSMIC Cancer Gene Census file
#' @param human_assembly `"GRCh38"` or `"GRCh37"`; must match your COSMIC download
#' @param transfer_annotation transfer human UniProt functional annotation onto
#'   the structure across the ortholog alignment
#' @param plddt_threshold pLDDT below which residues are flagged as not
#'   interpretable (and excluded from the 3D clustering tests)
#' @param hotspot_radius sphere radius (Angstrom) for the 3D hotspot scan
#' @param n_perm permutations for the clustering tests; set to 0 to skip them
#' @param alphamissense annotate missense calls with AlphaMissense pathogenicity
#'   when AlphaFold serves it for this entry
#' @param colour_by initial colouring of the 3D view
#' @param outdir optional directory; when given, the figure (PDF), the
#'   interactive view (HTML) and the report (TSV) are written there
#' @param quiet suppress progress messages
#'
#' @return list of class `xsmut3d_result`:
#'   `gene`, `species`, `assembly`, `target`, `model`, `structure`, `domains`,
#'   `mutations` (per-call table), `residues` (per-residue table),
#'   `cosmic`, `report` (flagged, ranked mutations), `cluster`, `hotspots`,
#'   `confidence`, `plot` (patchwork), `view` (r3dmol widget), `files`
#' @export
#'
#' @examples
#' \dontrun{
#' res <- xsmut3d_pipeline(
#'   gene      = "TP53",
#'   mutations = muts,                 # dndscv(...)$annotmuts
#'   species   = "macaque",
#'   assembly  = "Mmul_10",
#'   cosmic    = "data/cosmic/Cosmic_MutantCensus_v104_GRCh38.tsv.gz",
#'   outdir    = "results/TP53"
#' )
#' res$view        # rotate in RStudio / a browser
#' res$plot        # ggsave("TP53.pdf", res$plot, width = 11, height = 9)
#' head(res$report)
#' }
xsmut3d_pipeline <- function(gene, mutations, species = NULL, assembly = NULL,
                             structure_species = c("auto", "query", "human"),
                             cosmic = NULL, census_path = NULL,
                             human_assembly = c("GRCh38", "GRCh37"),
                             transfer_annotation = TRUE,
                             plddt_threshold = 70,
                             hotspot_radius = 10,
                             n_perm = 1000L,
                             alphamissense = TRUE,
                             colour_by = c("plddt", "domain", "burden"),
                             outdir = NULL, quiet = FALSE) {
  structure_species <- match.arg(structure_species)
  human_assembly <- match.arg(human_assembly)
  colour_by <- match.arg(colour_by)
  old <- options(xsmut3d.quiet = quiet); on.exit(options(old), add = TRUE)

  if (length(gene) > 1L) {
    out <- lapply(gene, function(g) tryCatch(
      xsmut3d_pipeline(g, mutations, species, assembly, structure_species, cosmic,
                       census_path, human_assembly, transfer_annotation,
                       plddt_threshold, hotspot_radius, n_perm, alphamissense,
                       colour_by,
                       outdir = if (is.null(outdir)) NULL else file.path(outdir, g),
                       quiet = quiet),
      error = function(e) { warning(g, ": ", conditionMessage(e)); NULL }))
    names(out) <- gene
    return(out)
  }

  # ---- 1. species and assembly -------------------------------------------
  sp <- NULL; asm <- NULL
  if (!is.null(species) && !tolower(species) %in% c("human", "homo_sapiens", "homo sapiens")) {
    sp <- xsmut::ens_resolve_species(species)
    asm <- xsmut::ens_assembly(sp, assembly)
    .msg(sprintf("Species: %s (%s)", sp, asm$assembly_name))
  }

  # ---- 2. which structure? ------------------------------------------------
  target <- resolve_structure_target(gene, sp, structure_species, human_assembly)
  .msg(sprintf("Structure target: UniProt %s (%s anchor)", target$accession, target$anchor))
  model <- af_model(target$accession, alphamissense = alphamissense)
  .msg(sprintf("AlphaFold %s, %d residues, mean pLDDT %.1f",
               model$entry_id, nchar(model$sequence), model$global_plddt))
  struct <- read_structure(model, contact_cutoff = hotspot_radius)

  # ---- 3. functional annotation ------------------------------------------
  hmap <- NULL; pid <- NA_real_
  if (transfer_annotation && !is.na(target$human_accession) &&
      !identical(target$human_accession, target$accession)) {
    hseq <- uniprot_entry(target$human_accession)$sequence
    hmap <- residue_map(hseq, model$sequence)
    pid <- attr(hmap, "pid")
    .msg(sprintf("Transferring human annotation across a %.1f%% identical alignment.", pid))
  }
  domains <- domain_table(target$accession,
                          human_accession = if (transfer_annotation) target$human_accession else NULL,
                          map = hmap, transfer = transfer_annotation,
                          transfer_identity = pid)

  # ---- 4. place the mutations --------------------------------------------
  mapped <- map_mutations_to_structure(mutations, gene, model, species = sp,
                                       fallback_protein_id = target$species_protein_id)
  mapped <- annotate_confidence(mapped, struct, pae = model$pae, threshold = plddt_threshold)
  if (nrow(mapped)) {
    fl <- flag_functional(mapped$struct_pos, struct, domains, radius = 8)
    mapped <- cbind(mapped, fl[, .(domain, site_hit, site_near, site_near_dist)])
  }

  # AlphaMissense, where available and where the substitution is a real missense
  if (!is.null(model$alphamissense) && nrow(mapped)) {
    am <- model$alphamissense
    key_m <- paste0(mapped$aa_ref, mapped$struct_pos - (model$uniprot_start - 1L), mapped$aa_alt)
    i <- match(key_m, am$variant)
    mapped[, `:=`(am_pathogenicity = am$am_pathogenicity[i], am_class = am$am_class[i])]
  } else if (nrow(mapped)) {
    mapped[, `:=`(am_pathogenicity = NA_real_, am_class = NA_character_)]
  }

  residues <- collapse_by_residue(mapped)
  n_placed <- sum(!is.na(mapped$struct_pos))
  .msg(sprintf("%d of %d calls placed on the structure (%d distinct residues).",
               n_placed, nrow(mapped), nrow(residues)))

  # ---- 5. COSMIC ----------------------------------------------------------
  cos_res <- NULL; gene_check <- NULL
  if (!is.null(cosmic)) {
    # xsmut streams the file through gzip, which does not expand "~"
    cos_tab <- if (is.character(cosmic)) xsmut::cosmic_load(path.expand(cosmic), genes = gene) else cosmic
    if (!is.null(census_path)) census_path <- path.expand(census_path)
    cos_mut <- xsmut::cosmic_gene_mutations(gene, cos_tab)
    gene_check <- xsmut::cosmic_gene_listed(gene, census_path = census_path, mutations = cos_tab)
    cos_res <- cosmic_on_structure(cos_mut, model, human_accession = target$human_accession)
    if (nrow(cos_res)) {
      hot_cos <- cos_res[, .(cosmic_samples = sum(n_samples)), by = struct_pos]
      mapped[, cosmic_samples_at_residue :=
               hot_cos$cosmic_samples[match(struct_pos, hot_cos$struct_pos)]]
      mapped[is.na(cosmic_samples_at_residue), cosmic_samples_at_residue := 0L]
    }
  }
  if (!"cosmic_samples_at_residue" %in% names(mapped)) mapped[, cosmic_samples_at_residue := NA_integer_]

  # ---- 6. 3D clustering ---------------------------------------------------
  clust <- NULL; hot <- NULL
  if (n_perm > 0 && nrow(residues) >= 3L) {
    clust <- cluster_test_3d(residues, struct, n_perm = n_perm, min_plddt = plddt_threshold)
    hot <- hotspots_3d(residues, struct, radius = hotspot_radius,
                       n_perm = n_perm, min_plddt = plddt_threshold)
    if (nrow(hot)) {
      mapped[, hotspot_fdr := hot$fdr[match(struct_pos, hot$struct_pos)]]
      residues[, hotspot_fdr := hot$fdr[match(struct_pos, hot$struct_pos)]]
    }
    if (!is.na(clust$p_value)) {
      .msg(sprintf("3D clustering: mean weighted pairwise distance %.1f A vs null %.1f A (z = %.2f, p = %.4g)",
                   clust$observed, clust$null_mean, clust$z, clust$p_value))
    }
  }
  if (!"hotspot_fdr" %in% names(mapped)) mapped[, hotspot_fdr := NA_real_]

  # ---- 7. summaries, figure, view -----------------------------------------
  conf <- confidence_summary(struct, model, mapped)
  report <- .build_report(mapped)

  sp_label <- if (!is.null(sp)) sprintf("%s (%s)", tools::toTitleCase(gsub("_", " ", sp)),
                                        assembly %||% asm$assembly_name) else "Human"
  ttl <- sprintf("%s - %s on %s", gene, sp_label, model$entry_id)
  sub <- paste(
    conf$text,
    if (!is.null(clust) && !is.na(clust$p_value))
      sprintf("3D clustering of mutated residues: z = %.2f, p = %.4g (%d permutations, pLDDT >= %g).",
              clust$z, clust$p_value, n_perm, plddt_threshold),
    if (!is.na(pid)) sprintf("Functional annotation transferred from human %s at %.1f%% identity.",
                             target$human_accession, pid),
    sep = "\n")

  p <- plot_structure_tracks(struct, residues, domains, cosmic = cos_res,
                             hotspots = hot, title = ttl, subtitle = sub)
  v <- view_structure_3d(struct, residues, domains, colour_by = colour_by)

  files <- character()
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    f_pdf <- file.path(outdir, paste0(gene, "_structure_tracks.pdf"))
    ggplot2::ggsave(f_pdf, p, width = 11, height = 9.5)
    f_html <- file.path(outdir, paste0(gene, "_3D.html"))
    save_view_html(v, f_html)
    f_tsv <- file.path(outdir, paste0(gene, "_report.tsv"))
    data.table::fwrite(report, f_tsv, sep = "\t")
    f_res <- file.path(outdir, paste0(gene, "_residues.tsv"))
    data.table::fwrite(residues, f_res, sep = "\t")
    files <- c(figure = f_pdf, view = f_html, report = f_tsv, residues = f_res)
    .msg("Wrote: ", paste(basename(files), collapse = ", "), " to ", outdir)
  }

  structure(list(
    gene = gene, species = sp, assembly = asm, target = target,
    model = model, structure = struct, domains = domains,
    mutations = mapped, residues = residues,
    cosmic = cos_res, gene_check = gene_check,
    report = report, cluster = clust, hotspots = hot,
    confidence = conf, plot = p, view = v, files = files
  ), class = "xsmut3d_result")
}

# Rank the calls a reader should look at first. The ordering is deliberate:
# a confidently-modelled residue that hits a key site and recurs beats a
# high-AlphaMissense call sitting in a disordered tail.
.build_report <- function(mapped) {
  m <- as.data.table(mapped)[!is.na(struct_pos)]
  if (!nrow(m)) return(m)
  n_at <- m[, .N, by = struct_pos]
  m[, n_at_residue := n_at$N[match(struct_pos, n_at$struct_pos)]]
  m[, relevance := 0]
  m[!is.na(site_hit),      relevance := relevance + 4]
  m[!is.na(site_near),     relevance := relevance + 2]
  m[!is.na(domain),        relevance := relevance + 1]
  m[n_at_residue > 1,      relevance := relevance + pmin(3, n_at_residue - 1)]
  m[!is.na(hotspot_fdr) & hotspot_fdr < 0.1, relevance := relevance + 2]
  m[!is.na(cosmic_samples_at_residue) & cosmic_samples_at_residue > 0,
    relevance := relevance + pmin(3, log10(cosmic_samples_at_residue + 1) * 2)]
  m[!is.na(am_pathogenicity), relevance := relevance + 2 * am_pathogenicity]
  m[interpretable == FALSE, relevance := relevance * 0.4]   # discount, don't drop
  cols <- intersect(c("sample_id", "gene", "aa_change", "consequence", "struct_pos",
                      "struct_aa", "aa_ref", "aa_alt", "ref_match", "align_identity",
                      "plddt", "plddt_band", "interpretable", "confidence_note",
                      "burial", "domain", "site_hit", "site_near", "site_near_dist",
                      "n_at_residue", "hotspot_fdr", "cosmic_samples_at_residue",
                      "am_pathogenicity", "am_class", "relevance", "map_note"),
                    names(m))
  m[order(-relevance, struct_pos), ..cols]
}

#' @export
print.xsmut3d_result <- function(x, ...) {
  cat(sprintf("<xsmut3d_result> %s on %s (%s)\n", x$gene, x$model$entry_id, x$model$organism))
  cat(sprintf("  %d calls, %d placed on %d residues; mean pLDDT %.1f\n",
              nrow(x$mutations), sum(!is.na(x$mutations$struct_pos)),
              nrow(x$residues), x$model$global_plddt))
  if (!is.null(x$cluster) && !is.na(x$cluster$p_value))
    cat(sprintf("  3D clustering p = %.4g (z = %.2f)\n", x$cluster$p_value, x$cluster$z))
  n_key <- sum(!is.na(x$report$site_hit)) + sum(!is.na(x$report$site_near))
  cat(sprintf("  %d call(s) hit or neighbour an annotated key functional site\n", n_key))
  cat("  $plot  $view  $report  $residues  $confidence$text\n")
  invisible(x)
}

#' Top-ranked mutations as a printable summary
#'
#' @param x an `xsmut3d_result`
#' @param n number of rows
#' @return data.table
#' @export
top_hits <- function(x, n = 15L) {
  cols <- intersect(c("aa_change", "consequence", "struct_pos", "plddt_band",
                      "domain", "site_hit", "site_near", "n_at_residue",
                      "hotspot_fdr", "cosmic_samples_at_residue", "am_class"),
                    names(x$report))
  utils::head(unique(x$report[, ..cols]), n)
}
