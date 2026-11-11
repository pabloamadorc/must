# Builds data/alba_arthropods.rda from the data of Alba et al. (2026).
#
# Source: Alba, C., Lucky, A., Fahey, C., Somarriba, G., Dunbar, M., Inman,
#   W., Flory, S. L., & Angelini, C. (2025). Data from: Plant invasion
#   overrides drought in suppressing arthropod abundance and richness.
#   Dryad. https://doi.org/10.5061/dryad.02v6wwqhz
#   File: Data_For_Repository_Arthropods.xlsx (version of 2025-12-29), sheet
#   "Abundance_Raw_Data".
#   MD5 of the file used: 952592954622c19e72f19e5da65666b4
# Licence: CC0 1.0 (Dryad publishes accepted datasets under CC0:
#   https://datadryad.org/publication_policy).
# Article: Alba, C., Lucky, A., Fahey, C., Somarriba, G., Dunbar, M., Inman,
#   W., Flory, S. L., & Angelini, C. (2026). Plant invasion overrides drought
#   in suppressing arthropod abundance and richness. Oikos, 2026(5), e11920.
#   https://doi.org/10.1002/oik.11920
#
# Design (article, Material and methods): 40 plots of 4 x 4 m in an old field
# at Gainesville (Florida, USA), planted in 2012 with longleaf pine and native
# understory species. Full factorial of plant invasion (nine cogongrass
# ramets planted per plot, or none) and drought (rainout shelters, or ambient
# precipitation under control structures), applied in spring 2013 in a
# randomized block design: 10 blocks with one plot per treatment. Arthropods
# were sampled with 20 sweeps per plot in July (summer) and October (fall) of
# 2014 and 2015.
#
# Decisions:
# * One row per plot and sampling, one column per arthropod group (counts of
#   the source column `abundance`), named in English after the article, and
#   their total; attribute "taxa" gives the source code of each group.
# * Treatments from the treatment code `trt` (first letter A, ambient
#   precipitation, or D, drought; second letter N, no invasion, or C,
#   cogongrass). In the fourth sampling the ten rows of plot 38 say
#   "uninvaded" in the column `invasion`, although its code is "DC" in every
#   sampling and its invasion is "invaded" in the three earlier ones: the
#   code is used.
# * The square-root transformed counts (`sqrtabun`) and the functional groups
#   are left out; the functional groups are documented in R/data.R.
#
# Packages used by this script only: readxl, usethis.
# Download the file above from Dryad into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "Data_For_Repository_Arthropods.xlsx")
stopifnot(unname(tools::md5sum(path)) == "952592954622c19e72f19e5da65666b4")

raw <- as.data.frame(readxl::read_excel(path, sheet = "Abundance_Raw_Data"))
stopifnot(nrow(raw) == 1600L,
          !anyNA(raw[c("plot", "block", "date", "trt", "arthroID", "abundance")]))

# Groups in the order of the article (herbivores, pollinators, predators).
groups <- c(
  hemipterans = "hemtotal", caterpillars = "lepnumlarva", beetles = "coletotal",
  grasshoppers = "orthtotal", flies = "dipnumadult", moths = "lepnumadult",
  bees = "beetotal", wasps = "wasptotal", spiders = "spidnumindiv", ants = "anttotal"
)
stopifnot(setequal(raw$arthroID, groups), all(raw$abundance >= 0),
          all(raw$abundance == round(raw$abundance)))

# Only the plot 38 rows of the fourth sampling disagree with the code.
code_invaded <- substr(raw$trt, 2L, 2L) == "C"
disagree <- code_invaded != (raw$invasion == "invaded")
stopifnot(sum(disagree) == 10L, all(raw$plot[disagree] == 38), all(raw$date[disagree] == 4),
          all(raw$trt[raw$plot == 38] == "DC"),
          all(substr(raw$trt, 1L, 1L) == ifelse(raw$drought == "drought", "D", "A")))

keys <- unique(raw[c("plot", "block", "date", "year", "season", "trt")])
stopifnot(nrow(keys) == 160L, !anyDuplicated(keys[c("plot", "date")]))
counts <- vapply(groups, function(g) {
  x <- raw[raw$arthroID == g, ]
  x$abundance[match(paste(keys$plot, keys$date), paste(x$plot, x$date))]
}, numeric(nrow(keys)))
stopifnot(!anyNA(counts))

alba_arthropods <- data.frame(
  plot = as.integer(keys$plot),
  block = as.integer(keys$block),
  sampling = as.integer(keys$date),
  year = as.integer(keys$year),
  season = keys$season,
  drought = ifelse(substr(keys$trt, 1L, 1L) == "D", "drought", "ambient"),
  invasion = ifelse(substr(keys$trt, 2L, 2L) == "C", "invaded", "uninvaded"),
  stringsAsFactors = FALSE
)
storage.mode(counts) <- "integer"
alba_arthropods$total <- as.integer(rowSums(counts))
alba_arthropods <- cbind(alba_arthropods, as.data.frame(counts))
alba_arthropods <- alba_arthropods[order(alba_arthropods$sampling, alba_arthropods$plot,
                                         method = "radix"), ]
rownames(alba_arthropods) <- NULL
attr(alba_arthropods, "taxa") <- groups

# Samplings 1 to 4 are summer and fall of 2014 and 2015; each plot keeps its
# treatment and block, and every block has one plot per treatment.
stopifnot(
  identical(unique(paste(alba_arthropods$sampling, alba_arthropods$season, alba_arthropods$year)),
            c("1 summer 2014", "2 fall 2014", "3 summer 2015", "4 fall 2015")),
  all(table(alba_arthropods$plot) == 4L),
  all(tapply(paste(alba_arthropods$drought, alba_arthropods$invasion, alba_arthropods$block),
             alba_arthropods$plot, function(x) length(unique(x))) == 1L),
  all(table(alba_arthropods$block, paste(alba_arthropods$drought, alba_arthropods$invasion)) == 4L)
)

usethis::use_data(alba_arthropods, overwrite = TRUE, compress = "xz")
