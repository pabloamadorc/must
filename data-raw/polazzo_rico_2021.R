# Builds data/polazzo_zooplankton.rda, data/polazzo_phytoplankton.rda and
# data/polazzo_macroinvertebrates.rda from the data of Polazzo & Rico (2021).
#
# Source: Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple
#   stressors on the dimensionality of ecological stability. Dryad.
#   https://doi.org/10.5061/dryad.qfttdz0gt
#   File: Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx (version of
#   2021-06-04), sheets "Zooplankton density", "Phyto_D-5", "Phyto_D7",
#   "Phyto_D15", "Phyto_D50" and "Macroinvertebrates".
#   MD5 of the file used: ce677f670150df298b91663fb9f0feb2
# Licence: CC0 1.0 (Dryad publishes accepted datasets under CC0:
#   https://datadryad.org/publication_policy).
# Article: Polazzo, F., & Rico, A. (2021). Effects of multiple stressors on
#   the dimensionality of ecological stability. Ecology Letters, 24(8),
#   1594-1606. https://doi.org/10.1111/ele.13770
#
# Design (article, Methods): full factorial with three replicates of
# chlorpyrifos (1 ug/L or absent), diuron (18 ug/L or absent) and nutrients
# (added or not; additions started 4 weeks before the pesticides), 24
# outdoor mesocosms. Sampling days relative to the pesticide application:
# zooplankton -5, 7, 14, 30 and 50; phytoplankton -5, 7, 15 and 50;
# macroinvertebrates -32 (before the nutrient additions too), 15, 30 and 50.
#
# Decisions:
# * One data frame per community, one row per sample and one column per
#   taxon, named from the source (lowercase ASCII, punctuation as "_", a
#   final "sp." dropped); attribute "taxa" gives the source name of each
#   column. The zooplankton columns `daphnia` and `cyclopoida` keep the names
#   of the first version of the data set.
# * Treatment codes, which differ on day 50, are harmonised.
# * Phytoplankton: the sheet "Phyto_D50" says day 15 in its Day column; its
#   sheet name (day 50) is used. Each sheet lists only the taxa found that
#   day in some mesocosm; a taxon missing from a sheet was not found that
#   day, so its density is zero.
# * Macroinvertebrates: an unlabelled column ("...52", values for one
#   mesocosm per treatment only) is left out. The one missing value
#   (Hydropsychidae, mesocosm 23, day 50) stays missing, and so does the
#   total of that sample.
# * Densities and counts as in the source file, which does not state their
#   units.
#
# Packages used by this script only: readxl, usethis.
# Download the file above from Dryad into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx")
stopifnot(unname(tools::md5sum(path)) == "ce677f670150df298b91663fb9f0feb2")

read_sheet <- function(sheet) {
  x <- as.data.frame(readxl::read_excel(path, sheet = sheet))
  names(x) <- trimws(names(x))
  x
}

# ASCII version of a text (the accents of the source, explicitly, so that the
# result does not depend on the platform's iconv).
ascii <- function(x) {
  x <- chartr("\u00e1\u00e9\u00ed\u00f3\u00fa\u00f1\u00c1\u00c9\u00cd\u00d3\u00da\u00d1",
              "aeiounAEIOUN", enc2utf8(x))
  iconv(x, from = "UTF-8", to = "ASCII", sub = "")
}

# Column names from taxon names: lowercase ASCII, punctuation as "_", a final
# "sp" or "sp." dropped.
column_name <- function(x) {
  x <- tolower(ascii(trimws(x)))
  x <- sub("\\s+sp\\.?$", "", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# Treatment codes (they differ on day 50) to the three stressors.
recode <- c(
  "control" = "control", "I" = "I", "H" = "H", "N" = "N",
  "IxH" = "IxH", "IxN" = "IxN", "HxN" = "HxN", "IxHxN" = "IxHxN",
  "CPF" = "I", "DIU" = "H", "control+N" = "N", "mix" = "IxH",
  "CPF+N" = "IxN", "DIU+N" = "HxN", "mix+N" = "IxHxN"
)

design_columns <- function(raw, day) {
  stopifnot(all(raw$Treatment %in% names(recode)))
  code <- unname(recode[raw$Treatment])
  has <- function(letter) grepl(letter, code, fixed = TRUE)
  data.frame(
    mesocosm = sprintf("M%02d", as.integer(raw$Cosm)),
    day = as.integer(day),
    insecticide = ifelse(has("I"), "chlorpyrifos", "absent"),
    herbicide = ifelse(has("H"), "diuron", "absent"),
    nutrients = ifelse(has("N"), "added", "ambient"),
    stringsAsFactors = FALSE
  )
}

# One data frame: design, total, taxa (with attribute "taxa").
community <- function(design, taxa_values, source_names, total_name, renames = character()) {
  cols <- column_name(source_names)
  cols[cols %in% names(renames)] <- renames[cols[cols %in% names(renames)]]
  stopifnot(!anyDuplicated(cols), !any(cols %in% c(names(design), total_name)))
  names(taxa_values) <- cols
  out <- cbind(design, stats::setNames(data.frame(rowSums(taxa_values)), total_name), taxa_values)
  out <- out[order(out$day, out$mesocosm, method = "radix"), ]
  rownames(out) <- NULL
  # Every mesocosm keeps one treatment over the whole experiment.
  trt <- paste(out$insecticide, out$herbicide, out$nutrients)
  stopifnot(all(tapply(trt, out$mesocosm, function(x) length(unique(x))) == 1L))
  attr(out, "taxa") <- stats::setNames(ascii(source_names), cols)
  out
}

# Zooplankton ----------------------------------------------------------------
raw <- read_sheet("Zooplankton density")
taxa <- setdiff(names(raw), c("Day", "Cosm", "Treatment"))
polazzo_zooplankton <- community(
  design_columns(raw, sub("^D", "", raw$Day)), raw[taxa], taxa, "total_density",
  renames = c(ciclopoida = "cyclopoida")
)
stopifnot(
  nrow(polazzo_zooplankton) == 120L, ncol(polazzo_zooplankton) == 6L + 22L,
  all(table(polazzo_zooplankton$mesocosm) == 5L), !anyNA(polazzo_zooplankton),
  all(c("daphnia", "cyclopoida") %in% names(polazzo_zooplankton))
)

# Phytoplankton --------------------------------------------------------------
sheets <- c("Phyto_D-5" = -5L, "Phyto_D7" = 7L, "Phyto_D15" = 15L, "Phyto_D50" = 50L)
parts <- lapply(names(sheets), function(s) {
  x <- read_sheet(s)
  list(design = design_columns(x, rep(sheets[[s]], nrow(x))),
       values = x[setdiff(names(x), c("Day", "Cosm", "Treatment"))])
})
all_taxa <- unique(unlist(lapply(parts, function(p) names(p$values))))
values <- do.call(rbind, lapply(parts, function(p) {
  v <- p$values
  # A taxon missing from a sheet was not found that day.
  for (taxon in setdiff(all_taxa, names(v))) {
    v[[taxon]] <- 0
  }
  v[all_taxa]
}))
polazzo_phytoplankton <- community(
  do.call(rbind, lapply(parts, `[[`, "design")), values, all_taxa, "total_density"
)
stopifnot(
  nrow(polazzo_phytoplankton) == 96L, all(table(polazzo_phytoplankton$mesocosm) == 4L),
  !anyNA(polazzo_phytoplankton)
)

# Macroinvertebrates ---------------------------------------------------------
raw <- read_sheet("Macroinvertebrates")
taxa <- setdiff(names(raw), c("Day", "Cosm", "Treatment", "...52"))
polazzo_macroinvertebrates <- community(
  design_columns(raw, sub("^D", "", raw$Day)), raw[taxa], taxa, "total_abundance"
)
stopifnot(
  nrow(polazzo_macroinvertebrates) == 96L, all(table(polazzo_macroinvertebrates$mesocosm) == 4L),
  sum(is.na(polazzo_macroinvertebrates)) == 2L
)

usethis::use_data(polazzo_zooplankton, polazzo_phytoplankton, polazzo_macroinvertebrates,
                  overwrite = TRUE, compress = "xz")
