# Builds data/giesler_diatoms.rda from the data of Giesler et al. (2023).
#
# Source: Giesler, J., Wohlrab, S., & Harder, T. (2023). Growth of
#   Thalassiosira gravida and Thalassiosira rotula across different levels of
#   temperature, photoperiod, and microbiome presence. Dryad.
#   https://doi.org/10.5061/dryad.x95x69pqf
#   File: plateassay_growthrates.xlsx (version of 2023-11-17).
#   MD5 of the file used: 0c6dde075f10a6d0e172da44378ba2a2
# Licence: CC0 1.0 (Dryad publishes accepted datasets under CC0:
#   https://datadryad.org/publication_policy).
# Article: Giesler, J. K., Harder, T., & Wohlrab, S. (2023). Microbiome and
#   photoperiod interactively determine thermal sensitivity of polar and
#   temperate diatoms. Biology Letters, 19(11), 20230151.
#   https://doi.org/10.1098/rsbl.2023.0151
#
# Design (article, Methods): a polar (Thalassiosira gravida) and a temperate
# (T. rotula) diatom, each as a xenic culture (with its microbiome) and an
# axenic one (microbiome removed), grown at every combination of three
# temperatures (4, 9 and 13.5 C) and three photoperiods (4, 16 and 24 h of
# light at 50 umol photons/m2/s), after 7 days of acclimation. Each
# combination had 48 replicate wells of 300 ul in 96-well plates, incubated in
# climate cabinets at the temperature of the treatment on LED tables with its
# photoperiod; fluorescence was read daily for 7 days, and the maximum growth
# rate of each well was fitted to those readings.
#
# Decisions:
# * Columns renamed: `bacteria` to `microbiome` ("xenic" as "present",
#   "axenic" as "removed"); `species` as "T. gravida" and "T. rotula";
#   `mumax` to `growth_rate`.
# * The source identifier of each row (`index`) is left out; `replicate`
#   numbers the 48 wells of each combination.
# * Growth rates as in the source file (per day), negative ones included
#   (cultures that declined).
#
# Packages used by this script only: readxl, usethis.
# Download the file above from Dryad into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "plateassay_growthrates.xlsx")
stopifnot(unname(tools::md5sum(path)) == "0c6dde075f10a6d0e172da44378ba2a2")

raw <- as.data.frame(readxl::read_excel(path))
stopifnot(nrow(raw) == 1728L, !anyNA(raw),
          setequal(raw$species, c("T_gravida", "T_rotula")),
          setequal(raw$bacteria, c("xenic", "axenic")),
          setequal(raw$temperature, c(4, 9, 13.5)), setequal(raw$photoperiod, c(4, 16, 24)),
          is.numeric(raw$mumax))

giesler_diatoms <- data.frame(
  species = c(T_gravida = "T. gravida", T_rotula = "T. rotula")[raw$species],
  microbiome = c(xenic = "present", axenic = "removed")[raw$bacteria],
  temperature = as.numeric(raw$temperature),
  photoperiod = as.integer(raw$photoperiod),
  replicate = as.integer(raw$replicate),
  growth_rate = raw$mumax,
  stringsAsFactors = FALSE
)
giesler_diatoms <- giesler_diatoms[order(giesler_diatoms$species, giesler_diatoms$microbiome,
                                         giesler_diatoms$temperature, giesler_diatoms$photoperiod,
                                         giesler_diatoms$replicate, method = "radix"), ]
rownames(giesler_diatoms) <- NULL

# 48 wells, numbered 1 to 48, in each of the 2 x 2 x 3 x 3 combinations.
cell <- do.call(paste, giesler_diatoms[c("species", "microbiome", "temperature", "photoperiod")])
stopifnot(
  length(unique(cell)) == 36L, all(table(cell) == 48L),
  all(tapply(giesler_diatoms$replicate, cell, function(r) identical(sort(r), 1:48)))
)

usethis::use_data(giesler_diatoms, overwrite = TRUE, compress = "xz")
