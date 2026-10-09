#' @import ggplot2
NULL

# The static companion to the 3D view: everything shares one protein-position
# axis, so a reader can line up "where the mutations are", "what is there", and
# "how much of this model to believe" in a single glance.
#
#   panel 1  species mutation lollipops (coloured by consequence)
#   panel 2  human COSMIC lollipops, flipped (optional)
#   panel 3  domain architecture + key sites; transferred features hatched
#   panel 4  pLDDT profile with the AlphaFold confidence bands
#   panel 5  3D hotspot score (-log10 FDR), which has no 1D equivalent

.consequence_pal <- c(
  Missense = "#2E86AB", Nonsense = "#C73E1D", Frameshift = "#F18F01",
  "In-frame indel" = "#7B2D8E", Splice = "#3B8B5A", Synonymous = "#9E9E9E",
  Other = "#5E5E5E", Unknown = "#BDBDBD"
)

.consequence_class <- function(x) {
  x <- tolower(as.character(x))
  data.table::fcase(
    grepl("missense", x), "Missense",
    grepl("nonsense|stop_gained|stop gained", x), "Nonsense",
    grepl("frameshift", x), "Frameshift",
    grepl("inframe|in frame|deletion|insertion|indel", x), "In-frame indel",
    grepl("splice", x), "Splice",
    grepl("synonymous|silent", x), "Synonymous",
    is.na(x) | x == "", "Unknown",
    default = "Other"
  )
}

.xlim_for <- function(structure) {
  r <- structure$residues
  c(min(r$pos) - 2, max(r$pos) + 2)
}

.blank_x <- function(p) {
  p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.title.x = element_blank())
}

.lolli_panel <- function(dt, title, flip = FALSE, label_min = 2L, xlim) {
  dt <- as.data.table(dt)
  if (!nrow(dt)) {
    return(ggplot() + annotate("text", x = mean(xlim), y = 0,
                               label = "no mapped mutations", colour = "grey50", size = 3) +
             coord_cartesian(xlim = xlim) + theme_void() +
             labs(title = title))
  }
  dt[, cls := .consequence_class(top_consequence)]
  dt[, y := if (flip) -n_mutations else n_mutations]
  lab <- dt[n_mutations >= label_min]
  p <- ggplot(dt) +
    geom_segment(aes(x = struct_pos, xend = struct_pos, y = 0, yend = y),
                 colour = "grey65", linewidth = 0.3) +
    geom_point(aes(x = struct_pos, y = y, colour = cls, size = n_mutations), alpha = 0.9) +
    scale_colour_manual(values = .consequence_pal, name = NULL, drop = TRUE) +
    scale_size_area(max_size = 6, guide = "none") +
    scale_y_continuous(labels = abs) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(title = title, y = "mutations", x = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 10))
  if (nrow(lab)) {
    p <- p + geom_text(data = lab, aes(x = struct_pos, y = y, label = sub(",.*$", "", aa_changes)),
                       vjust = if (flip) 1.5 else -0.8, size = 2.5, check_overlap = TRUE)
  }
  .blank_x(p)
}

.domain_panel <- function(domains, structure, xlim, max_rows = 3L) {
  d <- as.data.table(domains)
  rng <- d[point == FALSE]
  ks  <- key_sites(domains)
  p <- ggplot() +
    annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = -0.12, ymax = 0.12,
             fill = "grey88", colour = NA)
  if (nrow(rng)) {
    rng <- rng[order(start, -(end - start))]
    # greedy row packing so overlapping domains stay readable
    rows <- integer(nrow(rng)); last_end <- rep(-Inf, max_rows)
    for (i in seq_len(nrow(rng))) {
      k <- which(last_end < rng$start[i])[1]
      if (is.na(k)) k <- which.min(last_end)
      rows[i] <- k; last_end[k] <- rng$end[i]
    }
    rng[, row := rows]
    rng[, `:=`(ymin = -0.45 - 0.55 * (row - 1), ymax = 0.45 - 0.55 * (row - 1))]
    rng[, fillc := .domain_palette[((seq_len(.N) - 1L) %% length(.domain_palette)) + 1L]]
    p <- p +
      geom_rect(data = rng, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
                fill = rng$fillc, colour = "white", alpha = 0.9) +
      geom_rect(data = rng[transferred == TRUE],
                aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
                fill = NA, colour = "black", linetype = "22", linewidth = 0.4) +
      geom_text(data = rng, aes(x = (start + end) / 2, y = (ymin + ymax) / 2, label = label),
                size = 2.4, colour = "white", fontface = "bold", check_overlap = TRUE)
  }
  if (nrow(ks)) {
    p <- p + geom_point(data = ks, aes(x = pos, y = 0.32), shape = 25, size = 1.6,
                        fill = "#2B8C3E", colour = "#2B8C3E")
  }
  p + coord_cartesian(xlim = xlim, clip = "off") +
    labs(y = NULL, x = NULL,
         caption = if (any(domains$transferred)) "dashed outline = annotation transferred from the human ortholog" else NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks = element_blank(),
          plot.caption = element_text(size = 7, colour = "grey40", hjust = 0))
}

.plddt_panel <- function(structure, xlim) {
  r <- structure$residues
  ggplot(r, aes(x = pos, y = plddt)) +
    annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = 90, ymax = 100,
             fill = plddt_palette[["Very high (>90)"]], alpha = 0.12) +
    annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = 70, ymax = 90,
             fill = plddt_palette[["Confident (70-90)"]], alpha = 0.15) +
    annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = 50, ymax = 70,
             fill = plddt_palette[["Low (50-70)"]], alpha = 0.22) +
    annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = 0, ymax = 50,
             fill = plddt_palette[["Very low (<50)"]], alpha = 0.22) +
    geom_hline(yintercept = 70, linetype = "32", colour = "grey35", linewidth = 0.3) +
    geom_line(linewidth = 0.4, colour = "grey15") +
    coord_cartesian(xlim = xlim, ylim = c(0, 100), clip = "off") +
    labs(y = "pLDDT", x = "residue (UniProt numbering)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
}

.hotspot_panel <- function(hot, xlim) {
  h <- as.data.table(hot)
  if (!nrow(h)) return(NULL)
  h[, score := -log10(fdr)]
  ggplot(h, aes(x = struct_pos, y = score)) +
    geom_col(width = 1, fill = "#7F0000") +
    geom_hline(yintercept = -log10(0.1), linetype = "32", colour = "grey35", linewidth = 0.3) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(y = expression(-log[10]~FDR), x = NULL, title = "3D neighbourhood hotspot scan") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold", size = 9))
}

#' Static multi-track figure: mutations, domains, model confidence, 3D hotspots
#'
#' @param structure a `protein_structure`
#' @param residue_counts from [collapse_by_residue()]
#' @param domains from [domain_table()]
#' @param cosmic optional COSMIC residue counts from [cosmic_on_structure()]
#' @param hotspots optional output of [hotspots_3d()]
#' @param title plot title
#' @param subtitle plot subtitle (the confidence summary is a good default)
#' @return a patchwork object; print it or pass to `ggsave()`
#' @export
plot_structure_tracks <- function(structure, residue_counts, domains,
                                  cosmic = NULL, hotspots = NULL,
                                  title = NULL, subtitle = NULL) {
  xlim <- .xlim_for(structure)
  panels <- list()
  panels[[length(panels) + 1L]] <-
    .lolli_panel(residue_counts, "Query-species mutations", flip = FALSE, xlim = xlim)
  if (!is.null(cosmic) && nrow(cosmic)) {
    cs <- as.data.table(cosmic)[, .(n_mutations = sum(n_samples),
                                    aa_changes = paste(unique(aa_change), collapse = ","),
                                    top_consequence = consequence[1]), by = struct_pos]
    panels[[length(panels) + 1L]] <-
      .lolli_panel(cs, "Human COSMIC", flip = TRUE, label_min = max(2L, ceiling(0.02 * sum(cs$n_mutations))), xlim = xlim)
  }
  panels[[length(panels) + 1L]] <- .domain_panel(domains, structure, xlim)
  hp <- if (!is.null(hotspots)) .hotspot_panel(hotspots, xlim) else NULL
  if (!is.null(hp)) panels[[length(panels) + 1L]] <- hp
  panels[[length(panels) + 1L]] <- .plddt_panel(structure, xlim)

  heights <- c(3, if (!is.null(cosmic) && nrow(cosmic)) 2.4, 1.6,
               if (!is.null(hp)) 1.2, 1.6)
  patchwork::wrap_plots(panels, ncol = 1, heights = heights) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                    plot.subtitle = element_text(size = 7.5, colour = "grey30", lineheight = 1.2))
    )
}
