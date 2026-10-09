# Silence R CMD check notes about data.table / ggplot2 non-standard evaluation.
utils::globalVariables(c(
  ".", ".N", ".SD", "..cols", "..keep",
  "aa", "aa_alt", "aa_change", "aa_pos", "aa_ref", "aa_changes", "align_identity",
  "alt", "am_class", "am_pathogenicity", "atom", "band", "burial",
  "chain", "class", "cls", "confidence_note", "consequence", "consequences",
  "cosmic_samples_at_residue", "description", "domain", "end", "evidence",
  "fdr", "fillc", "gene", "hotspot_fdr", "impact", "interpretable",
  "label", "low", "map_note", "mean_plddt", "n_at_residue", "n_in_sphere",
  "n_mutations", "n_neighbours", "n_samples", "p_value", "pae_global", "plddt",
  "plddt_band", "point", "pos", "ref", "ref_match", "relevance", "res_num",
  "row", "sample_id", "score", "site_hit", "site_near", "site_near_dist",
  "somatic_status", "species_protein_id", "start", "struct_aa", "struct_pos",
  "transfer_identity", "transferred", "type", "variant", "x", "y", "ymax",
  "ymin", "z", "uniprot_pos"
))
