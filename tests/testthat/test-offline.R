# Tests that need no network. Anything touching Ensembl / UniProt / AlphaFold
# is skipped unless XSMUT3D_NETWORK_TESTS=true.

test_that("HGVS protein substitutions parse", {
  s <- parse_aa_sub(c("p.R175H", "p.Q331*", "p.G12D", "p.E285_K286del", NA))
  expect_equal(s$ref[1:3], c("R", "Q", "G"))
  expect_equal(s$pos[1:4], c(175L, 331L, 12L, 285L))
  expect_equal(s$alt[1:3], c("H", "*", "D"))
  expect_true(is.na(s$ref[4]))      # not a simple substitution
  expect_true(is.na(s$pos[5]))
})

test_that("pLDDT banding follows the AlphaFold DB thresholds", {
  b <- plddt_band(c(95, 80, 60, 30, NA))
  expect_equal(as.character(b),
               c("Very high (>90)", "Confident (70-90)", "Low (50-70)",
                 "Very low (<50)", NA))
  expect_true(is.ordered(b))
  expect_true(all(names(plddt_palette) %in% levels(b)))
})

test_that("identical sequences give an identity residue map", {
  s <- "MEEPQSDPSVEPPLSQETFSDLWKLLPEN"
  m <- residue_map(s, s)
  expect_equal(m, seq_len(nchar(s)))
  expect_equal(attr(m, "pid"), 100)
})

test_that("residue maps survive an internal gap", {
  skip_if_not_installed("pwalign")
  a <- "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
  b <- sub("QRQ", "", a)                       # three residues deleted in b
  m <- residue_map(b, a)
  expect_length(m, nchar(b))
  expect_equal(m[1:6], 1:6)
  expect_true(all(diff(stats::na.omit(m)) > 0))
})

test_that("a PDB block parses into residues with pLDDT in the B-factor", {
  pdb <- c(
    "ATOM      1  N   MET A   1      -8.901   4.127  -0.555  1.00 95.70           N",
    "ATOM      2  CA  MET A   1      -8.608   3.135  -1.618  1.00 95.70           C",
    "ATOM      3  CB  MET A   1      -9.437   1.886  -1.404  1.00 95.70           C",
    "ATOM      4  CA  GLY A   2      -6.503   2.776  -2.161  1.00 42.30           C",
    "ATOM      5  CA  LYS A   3      -4.101   2.100  -3.000  1.00 78.10           C",
    "ATOM      6  CB  LYS A   3      -3.900   0.700  -3.400  1.00 78.10           C"
  )
  f <- tempfile(fileext = ".pdb"); writeLines(pdb, f)
  s <- read_structure(f, contact_cutoff = 10)
  expect_s3_class(s, "protein_structure")
  expect_equal(nrow(s$residues), 3L)
  expect_equal(s$residues$aa, c("M", "G", "K"))
  expect_equal(s$residues$plddt, c(95.70, 42.30, 78.10))
  expect_equal(as.character(s$residues$plddt_band[2]), "Very low (<50)")
  expect_equal(dim(s$dist), c(3L, 3L))
})

test_that("confidence segments split high and low stretches", {
  res <- data.table::data.table(pos = 1:30, plddt = c(rep(95, 10), rep(40, 10), rep(85, 10)))
  res[, plddt_band := plddt_band(plddt)]
  s <- structure(list(residues = res), class = "protein_structure")
  seg <- confidence_segments(s, threshold = 70, min_len = 1L)
  expect_equal(nrow(seg), 3L)
  expect_equal(seg$interpretable, c(TRUE, FALSE, TRUE))
  expect_equal(seg$start, c(1L, 11L, 21L))
})

test_that("key sites expand and functional flags fire on 3D proximity", {
  dom <- data.table::data.table(
    start = c(100L, 175L), end = c(300L, 175L),
    type = c("DOMAIN", "ACT_SITE"), class = c("Domain", "Active site"),
    description = c("DNA-binding domain", "catalytic residue"),
    point = c(FALSE, TRUE), evidence = NA_character_, source = "UniProt",
    transferred = c(TRUE, TRUE), transfer_identity = 97.5, label = c("DNA-binding", "cat")
  )
  ks <- key_sites(dom)
  expect_equal(ks$pos, 175L)
  expect_equal(annotate_positions(c(150L, 500L), dom), c("DNA-binding", NA))

  # fake structure: residue 248 is 5 A from residue 175 but 73 apart in sequence
  d <- matrix(99, 3, 3, dimnames = list(c("175", "248", "400"), c("175", "248", "400")))
  diag(d) <- 0; d["175", "248"] <- d["248", "175"] <- 5
  s <- structure(list(dist = d), class = "protein_structure")
  f <- flag_functional(c(175L, 248L, 400L), s, dom, radius = 8)
  expect_false(is.na(f$site_hit[1]))
  expect_true(is.na(f$site_hit[2]))
  expect_false(is.na(f$site_near[2]))      # caught only because of the 3D test
  expect_equal(f$site_near_dist[2], 5)
  expect_true(is.na(f$site_near[3]))
})

test_that("3D clustering detects a planted cluster and ignores a spread one", {
  set.seed(42)
  n <- 120
  xyz <- matrix(stats::runif(3 * n, 0, 60), n, 3)
  xyz[1:10, ] <- matrix(stats::runif(30, 0, 4), 10, 3)   # tight cluster
  res <- data.table::data.table(pos = seq_len(n), plddt = 95,
                                x = xyz[, 1], y = xyz[, 2], z = xyz[, 3])
  res[, plddt_band := plddt_band(plddt)]
  d <- as.matrix(stats::dist(xyz)); dimnames(d) <- list(res$pos, res$pos)
  s <- structure(list(residues = res, dist = d, contact_cutoff = 10),
                 class = "protein_structure")

  clustered <- data.table::data.table(struct_pos = 1:10, n_mutations = 2L)
  spread    <- data.table::data.table(struct_pos = seq(11, 110, by = 10), n_mutations = 2L)

  a <- cluster_test_3d(clustered, s, n_perm = 300, min_plddt = 70)
  b <- cluster_test_3d(spread,    s, n_perm = 300, min_plddt = 70)
  expect_lt(a$p_value, 0.05)
  expect_gt(b$p_value, 0.05)
  expect_lt(a$observed, b$observed)
})

test_that("low-confidence residues are excluded from the clustering tests", {
  res <- data.table::data.table(pos = 1:20, plddt = c(rep(95, 10), rep(40, 10)),
                                x = 1:20, y = 0, z = 0)
  res[, plddt_band := plddt_band(plddt)]
  d <- as.matrix(stats::dist(cbind(res$x, res$y, res$z)))
  dimnames(d) <- list(res$pos, res$pos)
  s <- structure(list(residues = res, dist = d), class = "protein_structure")
  rc <- data.table::data.table(struct_pos = c(1L, 2L, 3L, 15L, 16L), n_mutations = 1L)
  out <- cluster_test_3d(rc, s, n_perm = 100, min_plddt = 70)
  expect_equal(out$n_residues, 3L)
  expect_equal(out$excluded_low_confidence, 2L)
})

test_that("collapse_by_residue aggregates calls per residue", {
  m <- data.table::data.table(
    struct_pos = c(175L, 175L, 248L, NA),
    sample_id  = c("s1", "s2", "s1", "s3"),
    aa_change  = c("p.R175H", "p.R175C", "p.R248Q", NA),
    consequence = c("missense_variant", "missense_variant", "missense_variant", "splice_site_variant")
  )
  r <- collapse_by_residue(m)
  expect_equal(r$struct_pos, c(175L, 248L))
  expect_equal(r$n_mutations, c(2L, 1L))
  expect_equal(r$n_samples, c(2L, 1L))
  expect_match(r$aa_changes[1], "p.R175C, p.R175H")
})

test_that("the dNdScv adapter from xsmut is accepted", {
  skip_if_not_installed("xsmut")
  d <- data.frame(
    sampleID = c("MQD0001d", "MQD0002d"),
    chr = c("NC_041754.1", "NC_041754.1"), pos = c(95090953, 95090999),
    ref = c("A", "C"), mut = c("G", "T"), gene = c("TP53", "TP53"), strand = c(1L, 1L),
    ref_cod = "A", mut_cod = "G", ref3_cod = "TAT", mut3_cod = "TGT",
    aachange = c("R175H", "Y220C"), ntchange = c("A524G", "A659G"),
    codonsub = "TAT>TGT", impact = c("Missense", "Missense"),
    pid = "ENSMMUP00000024831", stringsAsFactors = FALSE
  )
  expect_true(xsmut::is_dndscv(d))
  x <- xsmut::as_xsmut(d)
  expect_equal(x$aa_pos, c(175L, 220L))
  expect_equal(x$species_protein_id[1], "ENSMMUP00000024831")
})

# --- network-dependent -----------------------------------------------------

skip_offline <- function() {
  if (!identical(tolower(Sys.getenv("XSMUT3D_NETWORK_TESTS")), "true"))
    testthat::skip("set XSMUT3D_NETWORK_TESTS=true to run network tests")
}

test_that("AlphaFold serves a model for human TP53", {
  skip_offline()
  m <- af_model("P04637", pae = FALSE, alphamissense = FALSE)
  expect_s3_class(m, "af_model")
  expect_equal(m$accession, "P04637")
  expect_equal(nchar(m$sequence), 393L)
  s <- read_structure(m)
  expect_equal(nrow(s$residues), 393L)
  expect_equal(s$residues$aa[175], "R")       # the hotspot residue
})

test_that("human TP53 features include the DNA-binding region", {
  skip_offline()
  ft <- uniprot_features("P04637")
  expect_true(nrow(ft) > 0)
  expect_true(any(grepl("DNA", ft$description, ignore.case = TRUE) | ft$class == "DNA binding"))
})
