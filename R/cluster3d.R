# Is the mutation pattern clustered in three dimensions?
#
# A lollipop plot can only see clustering in sequence. Selection acts on the
# folded protein, so residues that are far apart in sequence but adjacent in
# space (an active-site cleft, a dimer interface, a DNA-binding surface) are
# exactly the signal a 1D plot misses.
#
# Two complementary statistics, both permutation-based so no distributional
# assumptions are needed:
#
#   * a global test  - mean pairwise 3D distance between mutated residues,
#                      weighted by recurrence, versus the null of the same
#                      mutation load spread over the eligible residues;
#   * a per-residue  - the mutation burden within a sphere of `radius` around
#     hotspot scan     each residue, with an empirical p-value and BH FDR.
#
# Residues below the pLDDT threshold are excluded from both by default: their
# coordinates are not trustworthy, and including them inflates both tests.

#' Test whether mutated residues cluster in 3D
#'
#' @param residue_counts data.table with `struct_pos` and `n_mutations`
#'   (e.g. from [collapse_by_residue()])
#' @param structure a `protein_structure`
#' @param n_perm number of permutations
#' @param min_plddt residues below this pLDDT are excluded from both the
#'   observed set and the null (set to 0 to keep everything)
#' @param seed RNG seed for reproducibility
#' @return list: `observed`, `null_mean`, `z`, `p_value`, `n_residues`,
#'   `n_mutations`, `excluded_low_confidence`
#' @export
cluster_test_3d <- function(residue_counts, structure, n_perm = 2000L,
                            min_plddt = 70, seed = 1L) {
  set.seed(seed)
  r <- structure$residues
  eligible <- r$pos[!is.na(r$plddt) & r$plddt >= min_plddt]
  rc <- as.data.table(residue_counts)[struct_pos %in% eligible]
  excluded <- nrow(as.data.table(residue_counts)) - nrow(rc)
  if (nrow(rc) < 3L) {
    return(list(observed = NA_real_, null_mean = NA_real_, z = NA_real_,
                p_value = NA_real_, n_residues = nrow(rc),
                n_mutations = sum(rc$n_mutations), excluded_low_confidence = excluded,
                note = "fewer than 3 mutated residues with confident coordinates"))
  }
  d <- structure$dist
  idx_all <- match(as.character(eligible), colnames(d))
  stat <- function(idx, w) {
    sub <- d[idx, idx, drop = FALSE]
    ww <- outer(w, w)
    diag(ww) <- 0
    sum(sub * ww) / sum(ww)
  }
  obs_idx <- match(as.character(rc$struct_pos), colnames(d))
  obs <- stat(obs_idx, rc$n_mutations)
  k <- nrow(rc); w <- rc$n_mutations
  null <- vapply(seq_len(n_perm), function(i) {
    stat(sample(idx_all, k), w)
  }, 0)
  p <- (1 + sum(null <= obs)) / (n_perm + 1)   # one-sided: tighter than chance
  list(observed = obs, null_mean = mean(null), null_sd = stats::sd(null),
       z = (obs - mean(null)) / stats::sd(null), p_value = p,
       n_residues = k, n_mutations = sum(w), excluded_low_confidence = excluded,
       note = NA_character_)
}

#' Scan the structure for 3D mutation hotspots
#'
#' @param residue_counts data.table with `struct_pos` and `n_mutations`
#' @param structure a `protein_structure`
#' @param radius sphere radius in Angstrom
#' @param n_perm permutations for the empirical null
#' @param min_plddt pLDDT floor for eligible residues
#' @param seed RNG seed
#' @return data.table: `struct_pos`, `n_in_sphere`, `n_residues_in_sphere`,
#'   `expected`, `p_value`, `fdr`, sorted by `fdr`
#' @export
hotspots_3d <- function(residue_counts, structure, radius = 10, n_perm = 1000L,
                        min_plddt = 70, seed = 1L) {
  set.seed(seed)
  r <- structure$residues
  eligible <- r$pos[!is.na(r$plddt) & r$plddt >= min_plddt]
  rc <- as.data.table(residue_counts)[struct_pos %in% eligible]
  d <- structure$dist
  pos_all <- as.integer(colnames(d))
  keep <- pos_all %in% eligible
  dsub <- d[keep, keep, drop = FALSE]
  pos_e <- pos_all[keep]
  A <- dsub <= radius                                   # neighbourhood indicator

  cnt <- rep(0, length(pos_e))
  cnt[match(rc$struct_pos, pos_e)] <- rc$n_mutations
  obs <- as.vector(A %*% cnt)

  k <- nrow(rc); w <- rc$n_mutations
  if (!k) return(data.table(struct_pos = integer(), n_in_sphere = numeric(),
                            n_residues_in_sphere = integer(), expected = numeric(),
                            p_value = numeric(), fdr = numeric()))
  ge <- rep(0L, length(pos_e))
  for (i in seq_len(n_perm)) {
    cc <- rep(0, length(pos_e))
    cc[sample.int(length(pos_e), k)] <- w
    ge <- ge + (as.vector(A %*% cc) >= obs)
  }
  p <- (1 + ge) / (n_perm + 1)
  out <- data.table(struct_pos = pos_e,
                    n_in_sphere = obs,
                    n_residues_in_sphere = as.integer(rowSums(A)),
                    expected = sum(w) * rowSums(A) / length(pos_e),
                    p_value = p)
  out[, fdr := stats::p.adjust(p_value, "BH")]
  out[n_in_sphere > 0][order(fdr, -n_in_sphere)][]
}
