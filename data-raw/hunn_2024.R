# Builds data/hunn_invertebrates.rda from the data of Hunn et al. (2024).
#
# Source: Hunn, J. G. (2024). Data for Individual and combined impacts of
#   carbon dioxide enrichment, heatwaves, flow velocity variability and fine
#   sediment deposition on stream invertebrate communities (version v1).
#   Zenodo. https://doi.org/10.5281/zenodo.11064543
#   File: "Benthic data summary.txt".
#   MD5 of the file used (the same as in the record):
#   ec40c2fcf9219ed5ddd2bb9ff9a31e62
# Licence: CC BY 4.0 (record above).
# Article: Hunn, J. G., Orr, J. A., Kelly, A.-M., Piggott, J. J., & Matthaei,
#   C. D. (2024). Individual and combined impacts of carbon dioxide
#   enrichment, heatwaves, flow velocity variability, and fine sediment
#   deposition on stream invertebrate communities. Global Change Biology,
#   30(5), e17336. https://doi.org/10.1111/gcb.17336
#
# Design (article, Methods): 128 outdoor stream mesocosms (ExStream system,
# New Zealand), a full factorial of dissolved CO2 (ambient or enriched),
# temperature (ambient or two 7-day heatwaves), flow velocity (constant or
# variable) and fine sediment (none or 300 g dry sediment), eight mesocosms
# per combination, in two spatial blocks. The mesocosms were fed by eight
# header tanks of 16 mesocosms each, four per block; temperature and CO2 were
# applied to the header tanks, and flow and sediment to the mesocosms. The
# benthic community of each mesocosm was sampled on the last day (day 35);
# one quarter of each sample was counted.
#
# Decisions:
# * `header_tank` is not in the source file. It is deduced from the channel
#   numbers: channels 1-16, 17-32, ... 113-128 share their temperature and
#   CO2 and their block (checked below), as 16 mesocosms of a header tank do.
# * Counts of the counted quarter (columns "25pc_..." of the source): the
#   total of all invertebrates and the 14 taxa of the file. Annelids, which
#   fragment, are counted in halves, and so is the total. The source columns
#   scaled to whole samples (four times the quarter), richness, diversity and
#   body sizes are left out.
# * Taxon columns named from the source (lowercase); attribute "taxa" gives
#   the source name of each.
#
# Packages used by this script only: usethis.
# Download the file above from Zenodo into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "Benthic data summary.txt")
stopifnot(unname(tools::md5sum(path)) == "ec40c2fcf9219ed5ddd2bb9ff9a31e62")

raw <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(raw) == 128L, identical(sort(raw$Channel), 1:128),
          setequal(raw$Temp, c("Ambient", "Heated")),
          setequal(raw$Carbon_Dioxide, c("Ambient", "Enriched")),
          setequal(raw$Sediment, c("No", "Yes")), setequal(raw$Flow, c("Constant", "Variable")))
raw <- raw[order(raw$Channel), ]

quarter <- grep("^25pc_", names(raw), value = TRUE)
totals <- c("25pc_total_abundance", "25pc_EPT_total_abundance")
taxa_source <- sub("_abundance$", "", sub("^25pc_", "", setdiff(quarter, totals)))
stopifnot(length(taxa_source) == 14L)
taxa <- stats::setNames(taxa_source, tolower(taxa_source))

hunn_invertebrates <- data.frame(
  channel = as.integer(raw$Channel),
  block = as.integer(raw$Block),
  header_tank = as.integer(ceiling(raw$Channel / 16)),
  temperature = c(Ambient = "ambient", Heated = "heatwaves")[raw$Temp],
  co2 = c(Ambient = "ambient", Enriched = "enriched")[raw$Carbon_Dioxide],
  sediment = c(No = "none", Yes = "added")[raw$Sediment],
  flow = c(Constant = "constant", Variable = "variable")[raw$Flow],
  total_abundance = raw[["25pc_total_abundance"]],
  stringsAsFactors = FALSE
)
for (i in seq_along(taxa)) {
  hunn_invertebrates[[names(taxa)[i]]] <- raw[[paste0("25pc_", taxa[[i]], "_abundance")]]
}
rownames(hunn_invertebrates) <- NULL
attr(hunn_invertebrates, "taxa") <- taxa

# Each header tank: 16 consecutive channels in one block with one temperature
# and one CO2 level; two tanks per combination of the two, one in each block.
# Within each tank, the four combinations of flow and sediment, four channels
# each. Only annelids and the total have halves; the source scaled columns
# are four times the quarter.
tank_levels <- unique(hunn_invertebrates[c("header_tank", "block", "temperature", "co2")])
counts <- as.matrix(hunn_invertebrates[names(taxa)])
stopifnot(
  nrow(tank_levels) == 8L,
  all(table(tank_levels$temperature, tank_levels$co2) == 2L),
  all(table(tank_levels$block) == 4L),
  all(table(hunn_invertebrates$header_tank, paste(hunn_invertebrates$sediment,
                                                  hunn_invertebrates$flow)) == 4L),
  all(counts >= 0), all(counts[, colnames(counts) != "annelida"] ==
                          round(counts[, colnames(counts) != "annelida"])),
  all(2 * counts == round(2 * counts)),
  all(hunn_invertebrates$total_abundance >= rowSums(counts)),
  isTRUE(all.equal(raw$Total_abundance, 4 * raw[["25pc_total_abundance"]]))
)

usethis::use_data(hunn_invertebrates, overwrite = TRUE, compress = "xz")
