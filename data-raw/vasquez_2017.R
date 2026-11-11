# Builds data/vasquez_embryos.rda from the data of Vasquez et al. (2017).
#
# Source: Vasquez, M. C., Brockmann, H. J., & Julian, D. (2018). Data from:
#   Between-population differences in multi-stressor tolerance during embryo
#   development in the American horseshoe crab, Limulus polyphemus. Dryad.
#   https://doi.org/10.5061/dryad.vs1r0
#   File: Dryad+-+Embryo+Cross-Sectional+Area.xlsx (version of 2018-01-05),
#   sheet "Embryo Development Size".
#   MD5 of the file used: 9605af74035849f3cb0bba2cb2dd8cbe
# Licence: CC0 1.0 (Dryad publishes accepted datasets under CC0:
#   https://datadryad.org/publication_policy).
# Article: Vasquez, M. C., Brockmann, H. J., & Julian, D. (2017).
#   Between-population differences in multi-stressor tolerance during embryo
#   development in the American horseshoe crab, Limulus polyphemus.
#   Estuaries and Coasts, 40(5), 1487-1501.
#   https://doi.org/10.1007/s12237-017-0218-1
#
# Design (article, Materials and Methods): fertilized eggs of six mating
# pairs from Slaughter Beach (Delaware Bay, USA), May 2012, one egg per well
# of 96-well plates, exposed for 14 days to the 36 combinations of
# temperature (20, 25, 30 and 35 C), salinity (5, 15 and 34) and oxygen (5,
# 13 and 21% O2 in the air), one plate per combination, and then kept for 14
# days in the reference condition (30 C, 34, 21% O2). Every 2 days for 28
# days each egg was imaged; the cross-sectional area of eight eggs per pair
# and combination was measured, and an embryo developed when it reached the
# threshold of 5 mm2 (a late-stage embryo).
#
# Decisions:
# * One row per mating pair and combination: the number of embryos measured
#   and the number that developed (reached the threshold by day 28).
# * An embryo developed if the source column `Censor`, the event indicator of
#   the survival analysis of the article, is 1 on some day. Four embryos have
#   1 and then 0 on a later day; they count as developed. Two sizes of
#   97.93 mm2, far above all others, come after their embryos had passed the
#   threshold (5.08 and 7.21 mm2 on days 16 and 14), so they do not change
#   the result.
# * Three embryos are left out, so that their cells (pair and combination)
#   have 7 embryos: two never measured, and one (temperature 25, salinity 15,
#   oxygen 21, pair 4) whose event indicator is missing on every day; its
#   size passed 5 mm2 on days 14 to 18 and was below 2.5 mm2 afterwards.
#
# Packages used by this script only: readxl, usethis.
# Download the file above from Dryad into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "Dryad+-+Embryo+Cross-Sectional+Area.xlsx")
stopifnot(unname(tools::md5sum(path)) == "9605af74035849f3cb0bba2cb2dd8cbe")

raw <- as.data.frame(readxl::read_excel(path, sheet = "Embryo Development Size", skip = 2,
                                        na = "."))
stopifnot(identical(names(raw), c("Temperature", "Salinity", "Oxygen", "Nest", "Individual", "Row",
                                  "Column", "Day", "Size (mm)", "Censor")),
          nrow(raw) == 25920L, all(raw$Censor %in% c(0, 1, NA)))
raw$embryo <- paste(raw$Temperature, raw$Salinity, raw$Oxygen, raw$Nest, raw$Individual)
stopifnot(all(table(raw$embryo) == 15L), setequal(raw$Day, seq(0, 28, by = 2)))

embryos <- unique(raw[c("Temperature", "Salinity", "Oxygen", "Nest", "embryo")])
stopifnot(nrow(embryos) == 1728L)
by_embryo <- function(x, f) as.vector(tapply(x, raw$embryo, f)[embryos$embryo])
embryos$measured <- by_embryo(!is.na(raw$`Size (mm)`), any)
embryos$indicator <- by_embryo(!is.na(raw$Censor), all)
embryos$developed <- by_embryo(raw$Censor %in% 1, any)
# The indicator is complete or missing on every day; three embryos are left
# out (see above).
raw <- raw[order(raw$embryo, raw$Day), ]
reverts <- by_embryo(raw$Censor, function(x) any(diff(x) < 0, na.rm = TRUE))
stopifnot(identical(by_embryo(is.na(raw$Censor), any), !embryos$indicator),
          sum(!embryos$measured) == 2L, sum(embryos$measured & !embryos$indicator) == 1L,
          !any(embryos$developed[!embryos$measured]), sum(reverts) == 4L)
embryos <- embryos[embryos$measured & embryos$indicator, ]

cells <- unique(embryos[c("Temperature", "Salinity", "Oxygen", "Nest")])
key <- function(x) paste(x$Temperature, x$Salinity, x$Oxygen, x$Nest)
vasquez_embryos <- data.frame(
  temperature = as.integer(cells$Temperature),
  salinity = as.integer(cells$Salinity),
  oxygen = as.integer(cells$Oxygen),
  pair = as.integer(cells$Nest),
  embryos = as.integer(table(key(embryos))[key(cells)]),
  developed = as.integer(tapply(embryos$developed, key(embryos), sum)[key(cells)])
)
vasquez_embryos <- vasquez_embryos[order(vasquez_embryos$temperature, vasquez_embryos$salinity,
                                         vasquez_embryos$oxygen, vasquez_embryos$pair,
                                         method = "radix"), ]
rownames(vasquez_embryos) <- NULL

# Six pairs in each of the 4 x 3 x 3 combinations, eight embryos per pair
# except the three cells above.
stopifnot(
  nrow(vasquez_embryos) == 216L,
  all(table(vasquez_embryos$temperature, vasquez_embryos$salinity, vasquez_embryos$oxygen) == 6L),
  sum(vasquez_embryos$embryos == 7L) == 3L, all(vasquez_embryos$embryos %in% 7:8),
  all(vasquez_embryos$developed <= vasquez_embryos$embryos), sum(vasquez_embryos$embryos) == 1725L
)

usethis::use_data(vasquez_embryos, overwrite = TRUE, compress = "xz")
