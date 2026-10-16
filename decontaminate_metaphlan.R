# Positive-control noise subtraction for the Monumentome MetaPhlAn table.
#
# The sequenced positive controls are the ZymoBIOMICS Microbial Community
# Standard (D6300). Only the two standards that matched that mix are used to
# estimate noise: BC-0247288828 (Soldier's National Monument) and
# BC-0302366907 (Orpheus). BC-0247292781 and BC-0302370816 did not match that
# mix. Those two rows were removed from the metadata and are left out of the
# transformed table.
# Any species MetaPhlAn reports in one of
# those two standards that is not one of the ten organisms in the mix is
# treated as noise. For each such species the ceiling is its highest abundance
# in either of those two standards.
# That ceiling is subtracted from the species in every sample where the species
# is present (floored at zero), the removal is carried up the lineage and down
# to the strain (t__) bins, and each sample is renormalised so its top-level
# rows (kingdoms plus UNCLASSIFIED) sum to 100 again.
#
# Outputs (written under <root>/input):
#   positive_control_noise_species.tsv   the noise list with ceilings
#   metaphlan_merged_all_taxa_decontaminated.tsv   the transformed table
#
# Usage:
#   Rscript decontaminate_metaphlan.R            # from the project directory
#   source("decontaminate_metaphlan.R"); decontaminate_metaphlan(ROOT)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# The two positive standards whose profiles matched the ZymoBIOMICS mix.
# BC-0302370816 (Seoul Tower) and BC-0247292781 (Soldier's National Monument)
# did not. Those two rows were removed from the metadata, and their columns are
# left out of the transformed table. They do not set a ceiling.
QC_POSITIVE_CONTROLS <- c(
  "BC-0247288828",
  "BC-0302366907"
)
MISLABELED_SAMPLES <- c(
  "BC-0247292781",
  "BC-0302370816"
)

ZYMO_D6300_SPECIES <- c(
  "Pseudomonas_aeruginosa",
  "Escherichia_coli",
  "Salmonella_enterica",
  "Limosilactobacillus_fermentum",  # Zymo sheet: Lactobacillus fermentum
  "Enterococcus_faecalis",
  "Staphylococcus_aureus",
  "Listeria_monocytogenes",
  "Bacillus_spizizenii",            # Zymo sheet: Bacillus subtilis
  "Saccharomyces_cerevisiae",
  "Cryptococcus_neoformans"
)

decontaminate_metaphlan <- function(root = ".",
                                    metaphlan_file = file.path(root, "input", "metaphlan_merged_all_taxa.tsv"),
                                    metadata_file = file.path(root, "input", "Monumentome_metadata.csv"),
                                    noise_out = file.path(root, "input", "positive_control_noise_species.tsv"),
                                    table_out = file.path(root, "input", "metaphlan_merged_all_taxa_decontaminated.tsv"),
                                    positive_control_type = "positive_control",
                                    positive_control_ids = QC_POSITIVE_CONTROLS,
                                    zymo_species = ZYMO_D6300_SPECIES,
                                    quiet = FALSE) {
  strip_s_suffix <- function(x) str_remove(as.character(x), "_S[0-9]+.*$")

  header_lines <- readLines(metaphlan_file, n = 5)
  header_lines <- header_lines[str_starts(header_lines, "#")]

  mp <- read_tsv(metaphlan_file, comment = "#", show_col_types = FALSE, name_repair = "minimal")
  if ("#clade_name" %in% names(mp)) names(mp)[names(mp) == "#clade_name"] <- "clade_name"
  sample_cols <- setdiff(names(mp), c("clade_name", "NCBI_tax_id"))
  sample_ids <- strip_s_suffix(sample_cols)
  # Mislabeled positive controls are not part of the analysis table.
  drop_cols <- sample_cols[sample_ids %in% MISLABELED_SAMPLES]
  if (length(drop_cols) != length(MISLABELED_SAMPLES)) {
    stop("Mislabeled samples missing from the MetaPhlAn table: ",
         paste(setdiff(MISLABELED_SAMPLES, sample_ids), collapse = ", "))
  }
  sample_cols <- setdiff(sample_cols, drop_cols)
  sample_ids <- strip_s_suffix(sample_cols)

  metadata <- read_csv(metadata_file, show_col_types = FALSE, name_repair = "unique")
  metadata_ids <- unique(strip_s_suffix(metadata$sample_name))
  # Sample columns that are not in the metadata are not paper samples.
  sample_cols <- sample_cols[sample_ids %in% metadata_ids]
  sample_ids <- strip_s_suffix(sample_cols)
  typed_positive <- metadata %>%
    filter(sample_type == positive_control_type) %>%
    pull(sample_name) %>%
    strip_s_suffix()
  positive_ids <- strip_s_suffix(positive_control_ids)
  missing_ids <- setdiff(positive_ids, typed_positive)
  if (length(missing_ids) > 0) {
    stop("These QC positive controls are not typed ", positive_control_type, " in the metadata: ",
         paste(missing_ids, collapse = ", "))
  }
  positive_cols <- sample_cols[sample_ids %in% positive_ids]
  if (length(positive_cols) != length(positive_ids)) {
    stop("These QC positive controls have no MetaPhlAn column: ",
         paste(setdiff(positive_ids, sample_ids), collapse = ", "))
  }

  mat <- as.matrix(mp[, sample_cols])
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  rownames(mat) <- mp$clade_name

  is_species_row <- str_detect(mp$clade_name, "\\|s__[^|]+$")
  species_name <- if_else(is_species_row, str_remove(str_extract(mp$clade_name, "[^|]+$"), "^s__"), NA_character_)

  # Noise list: species seen in a positive control that are not in the standard.
  pos_mat <- mat[is_species_row, positive_cols, drop = FALSE]
  rownames(pos_mat) <- species_name[is_species_row]
  seen <- rowSums(pos_mat > 0) > 0
  noise <- tibble(
    species = rownames(pos_mat)[seen],
    ceiling_abundance = apply(pos_mat[seen, , drop = FALSE], 1, max),
    n_positive_controls = rowSums(pos_mat[seen, , drop = FALSE] > 0),
    positive_controls = apply(pos_mat[seen, , drop = FALSE], 1, function(v) paste(strip_s_suffix(positive_cols)[v > 0], collapse = ";"))
  ) %>%
    filter(!species %in% zymo_species) %>%
    arrange(desc(ceiling_abundance))

  # Subtract each noise species, propagate, then renormalise.
  removed_total <- setNames(numeric(length(sample_cols)), sample_cols)
  samples_touched <- setNames(integer(nrow(noise)), noise$species)
  for (i in seq_len(nrow(noise))) {
    sp <- noise$species[i]
    row_i <- which(is_species_row & species_name == sp)
    if (length(row_i) != 1) stop("Species row not unique: ", sp)
    lineage <- mp$clade_name[row_i]
    ancestors <- which(vapply(mp$clade_name, function(cn) cn != lineage & str_starts(lineage, fixed(paste0(cn, "|"))), logical(1)))
    children <- which(str_starts(mp$clade_name, fixed(paste0(lineage, "|"))))

    ab <- mat[row_i, ]
    removed <- pmin(ab, noise$ceiling_abundance[i])
    samples_touched[sp] <- sum(removed > 0)
    removed_total <- removed_total + removed

    mat[row_i, ] <- ab - removed
    for (a in ancestors) mat[a, ] <- mat[a, ] - removed
    if (length(children)) {
      scale <- ifelse(ab > 0, (ab - removed) / ab, 1)
      mat[children, ] <- sweep(mat[children, , drop = FALSE], 2, scale, "*")
    }
  }
  mat[mat < 1e-9] <- 0

  top_rows <- !str_detect(mp$clade_name, fixed("|"))
  totals <- colSums(mat[top_rows, , drop = FALSE])
  factor <- ifelse(totals > 0, 100 / totals, 1)
  mat <- sweep(mat, 2, factor, "*")
  mat <- round(mat, 6)

  out <- bind_cols(tibble(clade_name = mp$clade_name), as_tibble(mat))
  if ("NCBI_tax_id" %in% names(mp)) out <- bind_cols(out[, 1], mp["NCBI_tax_id"], out[, -1])

  dir.create(dirname(table_out), showWarnings = FALSE, recursive = TRUE)
  writeLines(
    c(header_lines,
      paste0("#decontaminated ", format(Sys.Date()), ": ", nrow(noise),
             " positive-control noise species subtracted at their ceiling, then renormalised to 100 per sample (decontaminate_metaphlan.R)")),
    table_out
  )
  write_tsv(out, table_out, append = TRUE, col_names = TRUE)

  noise <- noise %>% mutate(n_samples_affected = unname(samples_touched[species]))
  write_tsv(noise, noise_out)

  summary <- tibble(
    sample = sample_ids,
    removed_before_renormalisation = unname(removed_total),
    renormalisation_factor = unname(factor)
  )

  if (!quiet) {
    message(sprintf("Positive controls: %s", paste(strip_s_suffix(positive_cols), collapse = ", ")))
    message(sprintf("Noise species: %d. Largest ceiling: %s (%.5f%%).", nrow(noise), noise$species[1], noise$ceiling_abundance[1]))
    message(sprintf("Samples with any subtraction: %d of %d. Largest total removed from one sample: %.4f%%.",
                    sum(removed_total > 0), length(sample_cols), max(removed_total)))
    message("Wrote ", table_out)
    message("Wrote ", noise_out)
  }

  invisible(list(noise = noise, sample_summary = summary, table_file = table_out, noise_file = noise_out))
}

if (sys.nframe() == 0) {
  decontaminate_metaphlan(".")
}
