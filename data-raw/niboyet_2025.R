# Builds data/niboyet_denitrification.rda from the data of Niboyet et al.
# (2025).
#
# Source: Niboyet, A. (2025). Soil moisture and the number of simultaneous
#   stressors drive interactions among global changes on denitrification
#   (version v1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.16789113
#   File: Niboyetetal_CommunicationsEarthEnvironment_10august2025.xlsx.
#   MD5 of the file used (the same as in the record):
#   1fc0a7b3b4855b1acde5243208ae100b
# Licence: CC BY 4.0 (record above).
# Article: Niboyet, A., Le Roux, X., Chiariello, N. R., Dijkstra, P., Field,
#   C. B., & Hungate, B. A. (2025). Soil moisture and the number of
#   simultaneous stressors drive interactions among global changes on
#   denitrification. Communications Earth & Environment, 6, 704.
#   https://doi.org/10.1038/s43247-025-02703-5
#
# Design (article, Methods): Jasper Ridge Global Change Experiment
# (California, USA), an annual grassland under a full factorial of five
# treatments since 1998: atmospheric CO2 (ambient, or +275 umol/mol), warming
# (+2.0 C), precipitation (+50% and a 3-week extension of the rainy season)
# and nitrogen (+7 g N/m2/yr), plus fire (a controlled burn of half of the
# plots in July 2011). Randomized block split-plot design: CO2, temperature
# and fire on circular plots of 2 m diameter, precipitation and nitrogen on
# the four quadrants of each plot; 32 combinations with four quadrants each.
# Soil cores of every quadrant were taken 2 days after the fire and 5, 9, 21
# and 33 months later (the sampling 4 days before the fire is not in the
# file).
#
# Decisions:
# * One row per quadrant and sampling; quadrants identified by the source
#   `Id` (1 to 128) and the quadrant number within the plot.
# * Treatment codes (1 ambient or unburned, 2 elevated or burned) as labels;
#   `Time` (1 to 5) as months after the fire (0, 5, 9, 21 and 33, from the
#   header of the source column).
# * Potential denitrification (denitrifying enzyme activity) and soil water
#   content, as in the source file; the other variables are left out. The
#   source marks missing values as "NA": one denitrification and six soil
#   water values.
#
# Packages used by this script only: readxl, usethis.
# Download the file above from Zenodo into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "Niboyetetal_CommunicationsEarthEnvironment_10august2025.xlsx")
stopifnot(unname(tools::md5sum(path)) == "1fc0a7b3b4855b1acde5243208ae100b")

header <- as.data.frame(readxl::read_excel(path, col_names = FALSE, n_max = 2,
                                           .name_repair = "minimal"))
stopifnot(identical(unname(unlist(header[1, 11])), "(0, 5, 9, 21, 33 months post-fire)"),
          identical(unname(unlist(header[1, 16])), "(ng N-N2O g-1 dry soil h-1)"),
          identical(unname(unlist(header[1, 12])), "(%)"))
raw <- as.data.frame(readxl::read_excel(path, skip = 1, na = "NA"))
names(raw) <- trimws(names(raw))
stopifnot(nrow(raw) == 640L, identical(sort(unique(raw$Id)), as.numeric(1:128)),
          all(raw$Time %in% 1:5), all(table(raw$Id) == 5L),
          is.numeric(raw[["Denitrifying Enzyme Activity"]]),
          is.numeric(raw[["Soil water content"]]))

level <- function(x, labels) labels[as.integer(x)]
ambient <- c("ambient", "elevated")
niboyet_denitrification <- data.frame(
  block = as.integer(raw$Block),
  plot = as.integer(raw$Plot),
  quadrant = as.integer(raw$Quad),
  unit = as.integer(raw$Id),
  co2 = level(raw$CO2, ambient),
  temperature = level(raw$Temperature, ambient),
  water = level(raw$Water, ambient),
  nitrogen = level(raw$Nitrogen, ambient),
  fire = level(raw$Fire, c("unburned", "burned")),
  months = c(0L, 5L, 9L, 21L, 33L)[raw$Time],
  denitrification = raw[["Denitrifying Enzyme Activity"]],
  soil_water = raw[["Soil water content"]],
  stringsAsFactors = FALSE
)
niboyet_denitrification <- niboyet_denitrification[order(niboyet_denitrification$months,
                                                         niboyet_denitrification$unit,
                                                         method = "radix"), ]
rownames(niboyet_denitrification) <- NULL

# 32 plots in 4 blocks, each plot with one level of CO2, temperature and fire
# and four quadrants with the four combinations of water and nitrogen; every
# block has the 8 plot treatments and so the 32 combinations; one quadrant
# keeps its plot and treatment over the samplings.
d <- niboyet_denitrification
plots <- unique(d[c("plot", "block", "co2", "temperature", "fire")])
quads <- unique(d[c("unit", "plot", "quadrant", "water", "nitrogen")])
stopifnot(
  nrow(plots) == 32L, !anyDuplicated(plots$plot), all(table(plots$block) == 8L),
  all(table(plots$co2, plots$temperature, plots$fire) == 4L),
  nrow(quads) == 128L, !anyDuplicated(quads$unit),
  all(table(quads$plot, paste(quads$water, quads$nitrogen)) == 1L),
  all(table(do.call(paste, d[c("co2", "temperature", "water", "nitrogen", "fire")]),
            d$block) == 5L),
  sum(is.na(d$denitrification)) == 1L, sum(is.na(d$soil_water)) == 6L,
  all(d$denitrification > 0, na.rm = TRUE)
)

usethis::use_data(niboyet_denitrification, overwrite = TRUE, compress = "xz")
