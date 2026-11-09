# Builds data/polazzo_zooplankton.rda from the data of Polazzo & Rico (2021).
#
# Source: Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple
#   stressors on the dimensionality of ecological stability. Dryad.
#   https://doi.org/10.5061/dryad.qfttdz0gt
#   File: Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx,
#   sheet "Zooplankton density" (version of 2021-06-04).
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
# outdoor mesocosms. Zooplankton sampled 5 days before the pesticide
# application and on days 7, 14, 30 and 50 after it.
#
# Packages used by this script only: readxl, usethis.
# Download the file above from Dryad into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx")
stopifnot(unname(tools::md5sum(path)) == "ce677f670150df298b91663fb9f0feb2")

raw <- as.data.frame(readxl::read_excel(path, sheet = "Zooplankton density"))
names(raw) <- trimws(names(raw))
taxa <- setdiff(names(raw), c("Day", "Cosm", "Treatment"))

# Treatments are coded differently on day 50: harmonise both codings.
recode <- c(
  "control" = "control", "I" = "I", "H" = "H", "N" = "N",
  "IxH" = "IxH", "IxN" = "IxN", "HxN" = "HxN", "IxHxN" = "IxHxN",
  "CPF" = "I", "DIU" = "H", "control+N" = "N", "mix" = "IxH",
  "CPF+N" = "IxN", "DIU+N" = "HxN", "mix+N" = "IxHxN"
)
stopifnot(all(raw$Treatment %in% names(recode)))
code <- unname(recode[raw$Treatment])

# Every mesocosm keeps one treatment over the whole experiment.
stopifnot(all(tapply(code, raw$Cosm, function(x) length(unique(x))) == 1L))

has <- function(letter) grepl(letter, code, fixed = TRUE)

polazzo_zooplankton <- data.frame(
  mesocosm = sprintf("M%02d", as.integer(raw$Cosm)),
  day = as.integer(sub("^D", "", raw$Day)),
  insecticide = ifelse(has("I"), "chlorpyrifos", "absent"),
  herbicide = ifelse(has("H"), "diuron", "absent"),
  nutrients = ifelse(has("N"), "added", "ambient"),
  total_density = rowSums(raw[taxa]),
  daphnia = raw[["Daphnia sp"]],
  cyclopoida = raw[["Ciclopoida"]],
  stringsAsFactors = FALSE
)
polazzo_zooplankton <- polazzo_zooplankton[
  order(polazzo_zooplankton$day, polazzo_zooplankton$mesocosm, method = "radix"),
]
rownames(polazzo_zooplankton) <- NULL

stopifnot(
  nrow(polazzo_zooplankton) == 120L,
  all(table(polazzo_zooplankton$mesocosm) == 5L),
  !anyNA(polazzo_zooplankton)
)

usethis::use_data(polazzo_zooplankton, overwrite = TRUE)
