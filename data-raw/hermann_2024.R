# Builds the hermann_* data sets in data/ from the data of Hermann et al.
# (2024).
#
# Source: Hermann, M., Polazzo, F., Cherta Cucala, L., Crettaz Minaglia,
#   M. C., Garcia-Astillero, A., Peeters, E. T. H. M., Rico, A., & Van den
#   Brink, P. J. (2024). Data supporting: Combined stress of an insecticide
#   and heatwaves or elevated temperature induce community and food web
#   effects in a Mediterranean freshwater ecosystem (version V.1). Zenodo.
#   https://doi.org/10.5281/zenodo.11520459
#   File: CombinedStressInsecticideTemperature_Hermann_et_al.(2024).xlsx
#   MD5 of the file used (the same as in the record):
#   3a13c7c883d6f3a1a5c5d9f93d66d944
# Licence: CC BY 4.0 (record above).
# Article: Hermann, M., Polazzo, F., Cherta, L., Crettaz-Minaglia, M.,
#   Garcia-Astillero, A., Peeters, E. T. H. M., Rico, A., & Van den Brink,
#   P. J. (2024). Combined stress of an insecticide and heatwaves or elevated
#   temperature induce community and food web effects in a Mediterranean
#   freshwater ecosystem. Water Research, 260, 121903.
#   https://doi.org/10.1016/j.watres.2024.121903
#
# Design (article, section 2): 2 x 3 factorial, four replicates, 24 outdoor
# mesocosms (about 850 L) at Alcala de Henares (central Spain), 2021:
# imidacloprid (0 or 1 ug/L, applied on days 0, 14 and 28) and temperature
# (ambient; elevated, ambient + 4 C; three 7-day heatwaves, 8 C above the
# mean ambient temperature, separated by 7 days). Temperature treatments
# from day 0 to day 42 (treatment phase), then recovery to day 84.
#
# Decisions:
# * Communities: one data frame per community, one row per sample and one
#   column per taxon, named from the source (lowercase ASCII, punctuation as
#   "_", a final "sp" or "sp." dropped; "Ciclopoida" as `cyclopoida`, as in
#   polazzo_zooplankton); attribute "taxa" gives the source name of each
#   column. Sampling days from the sample labels: zooplankton and
#   phytoplankton -1, 10, 24, 38, 52, 66 and 80 (the article gives day -4
#   for the first zooplankton sample), macroinvertebrates -7, 7, 21, 35, 49,
#   63 and 77.
# * The community sheets identify mesocosms by treatment and replicate number
#   (e.g. "C0A1"), the other sheets by mesocosm number ("M14"): the source
#   gives no correspondence, so the two kinds of table cannot be joined by
#   mesocosm. A replicate number identifies one mesocosm on every sampling
#   day: within treatments, the log total abundance of a code is correlated
#   between consecutive samplings (macroinvertebrates 0.16 to 0.75,
#   zooplankton 0.14 to 0.56), and randomly paired codes are not (0.00).
# * Water quality measured in the morning and in the evening: one table with
#   `time_of_day`. Ecosystem productivity (evening rows only) is the evening
#   minus the morning dissolved oxygen of the same day in every row.
# * Leaf litter: microbial decomposition (fine mesh) and macroinvertebrate
#   feeding (coarse mesh) in one table with `process`. Two samples of
#   mesocosm M18 have no percentage loss in the source (no loss, or a gain);
#   they stay missing.
# * Other missing values stay missing: nitrate on day 52 in M07 and M09, and
#   imidacloprid in M20 on days 28, 31 and 38.
# * Values and units as in the source file; units given in the article are
#   documented, the others are not stated.
#
# Packages used by this script only: readxl, usethis.
# Download the file above from Zenodo into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "CombinedStressInsecticideTemperature_Hermann_et_al.(2024).xlsx")
stopifnot(unname(tools::md5sum(path)) == "3a13c7c883d6f3a1a5c5d9f93d66d944")

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

temperature_levels <- c(A = "ambient", E = "elevated", HW = "heatwaves")
imidacloprid_levels <- c(C0 = "absent", C1 = "applied")

# Communities ---------------------------------------------------------------
community <- function(design_sheet, abundance_sheet, total_name, renames = character()) {
  design <- read_sheet(design_sheet)
  values <- read_sheet(abundance_sheet)
  label <- design[[1L]]
  stopifnot(identical(label, values[[1L]]))
  day <- as.integer(sub("^D", "", design[["time (d)"]]))
  stopifnot(identical(as.character(design[["time (d)"]]), sub("C[01].*$", "", label)))
  chem <- design$chemtreatment
  temp <- design$temptreatment
  stopifnot(all(chem %in% names(imidacloprid_levels)), all(temp %in% names(temperature_levels)))
  mesocosm <- sub("^D-?[0-9]+", "", label)
  stopifnot(identical(mesocosm, paste0(chem, temp, sub("^.*[A-Z]", "", label))))
  taxa <- names(values)[-1L]
  cols <- column_name(taxa)
  cols[cols %in% names(renames)] <- renames[cols[cols %in% names(renames)]]
  stopifnot(!anyDuplicated(cols), !any(cols %in% c("mesocosm", "day", "imidacloprid",
                                                   "temperature", total_name)))
  v <- values[-1L]
  names(v) <- cols
  out <- data.frame(
    mesocosm = mesocosm,
    day = day,
    imidacloprid = unname(imidacloprid_levels[chem]),
    temperature = unname(temperature_levels[temp]),
    stringsAsFactors = FALSE
  )
  out[[total_name]] <- rowSums(v)
  out <- cbind(out, v)
  out <- out[order(out$day, out$mesocosm, method = "radix"), ]
  rownames(out) <- NULL
  stopifnot(all(table(out$mesocosm) == 7L), nrow(out) == 168L, !anyNA(out))
  attr(out, "taxa") <- stats::setNames(ascii(taxa), cols)
  out
}

hermann_zooplankton <- community("ZooplanktonDesign", "Zooplankton_Abundances", "total_density",
                                 renames = c(ciclopoida = "cyclopoida"))
hermann_phytoplankton <- community("PhytoplanktonDesign", "Phytoplankton_Abundances", "total_density")
hermann_macroinvertebrates <- community("MacroinvertebratesDesign", "Macroinvertebrates_Abundances",
                                        "total_abundance")

# Tables by mesocosm number --------------------------------------------------
numbered <- function(x, day, extra) {
  stopifnot(all(x$Temperature %in% names(temperature_levels)))
  out <- data.frame(
    mesocosm = sprintf("M%02d", as.integer(x$Cosm)),
    day = as.integer(day),
    imidacloprid = unname(imidacloprid_levels[x$Chemical]),
    temperature = unname(temperature_levels[x$Temperature]),
    stringsAsFactors = FALSE
  )
  stopifnot(!anyNA(out$imidacloprid))
  cbind(out, extra)
}

# Water quality, morning and evening.
morning <- read_sheet("WaterQualityMorning")
evening <- read_sheet("WaterQualityEvening")
names(evening)[names(evening) == "chla"] <- "Chla"
quality <- function(x, when) {
  numbered(x, x$Day, data.frame(
    time_of_day = when,
    water_temperature = x$Temp,
    oxygen = x$DO,
    oxygen_saturation = x$Dosat,
    conductivity = x$EC,
    ph = x$pH,
    tds = x$TDS,
    chlorophyll_a = x$Chla,
    phycoerythrin = x$PE,
    ecosystem_productivity = if ("ECOSYSP" %in% names(x)) x$ECOSYSP else NA_real_,
    stringsAsFactors = FALSE
  ))
}
hermann_water_quality <- rbind(quality(morning, "morning"), quality(evening, "evening"))
hermann_water_quality <- hermann_water_quality[
  order(hermann_water_quality$day, match(hermann_water_quality$time_of_day, c("morning", "evening")),
        hermann_water_quality$mesocosm, method = "radix"),
]
rownames(hermann_water_quality) <- NULL
stopifnot(nrow(hermann_water_quality) == 480L,
          all(table(hermann_water_quality$mesocosm) == 20L),
          !anyNA(hermann_water_quality[names(hermann_water_quality) != "ecosystem_productivity"]))
# Productivity: evening minus morning dissolved oxygen of the same day.
local({
  m <- hermann_water_quality[hermann_water_quality$time_of_day == "morning", ]
  e <- hermann_water_quality[hermann_water_quality$time_of_day == "evening", ]
  stopifnot(identical(m$mesocosm, e$mesocosm), identical(m$day, e$day),
            all(abs(e$ecosystem_productivity - (e$oxygen - m$oxygen)) < 1e-8),
            all(is.na(m$ecosystem_productivity)))
})

# Nutrients.
x <- read_sheet("Nutrients")
hermann_nutrients <- numbered(x, sub("^Day", "", x$Day), data.frame(
  nitrate = x$Nitrate,
  orthophosphate = x$Orthophosphate,
  total_phosphorus = x[["Total Phosphorus"]]
))
hermann_nutrients <- hermann_nutrients[order(hermann_nutrients$day, hermann_nutrients$mesocosm,
                                             method = "radix"), ]
rownames(hermann_nutrients) <- NULL
stopifnot(nrow(hermann_nutrients) == 240L, all(table(hermann_nutrients$mesocosm) == 10L),
          sum(is.na(hermann_nutrients)) == 2L, sum(is.na(hermann_nutrients$nitrate)) == 2L)

# Leaf litter: microbial decomposition and macroinvertebrate feeding.
litter <- function(sheet, process) {
  x <- read_sheet(sheet)
  numbered(x, x$Day, data.frame(
    process = process,
    leaf_dry_mass = x$Dryleaveweight,
    mass_loss = x$Dryleafdifference,
    mass_loss_percent = x$Degradedpercentage,
    stringsAsFactors = FALSE
  ))
}
hermann_leaf_litter <- rbind(litter("Microbial_Decomposition", "microbial"),
                             litter("Macroinvertebrate_Feeding", "macroinvertebrate feeding"))
hermann_leaf_litter <- hermann_leaf_litter[
  order(hermann_leaf_litter$process, hermann_leaf_litter$day, hermann_leaf_litter$mesocosm,
        method = "radix"),
]
rownames(hermann_leaf_litter) <- NULL
stopifnot(nrow(hermann_leaf_litter) == 264L, sum(is.na(hermann_leaf_litter)) == 2L,
          sum(is.na(hermann_leaf_litter$mass_loss_percent)) == 2L)

# Imidacloprid in water (treated mesocosms only).
x <- read_sheet("ImiConcentrationsWater")
stopifnot(all(x$Chemical == "Imi"), all(x$Temperature %in% names(temperature_levels)))
hermann_imidacloprid <- data.frame(
  mesocosm = sprintf("M%02d", as.integer(x$CosmNo)),
  day = as.integer(x$Day),
  temperature = unname(temperature_levels[x$Temperature]),
  concentration = x[["Concentration (ug/L)"]],
  stringsAsFactors = FALSE
)
hermann_imidacloprid <- hermann_imidacloprid[order(hermann_imidacloprid$day,
                                                   hermann_imidacloprid$mesocosm, method = "radix"), ]
rownames(hermann_imidacloprid) <- NULL
stopifnot(nrow(hermann_imidacloprid) == 120L, length(unique(hermann_imidacloprid$mesocosm)) == 12L,
          sum(is.na(hermann_imidacloprid)) == 3L)

# Every numbered mesocosm keeps one treatment, the same in every table.
trt <- rbind(hermann_water_quality[c("mesocosm", "imidacloprid", "temperature")],
             hermann_nutrients[c("mesocosm", "imidacloprid", "temperature")],
             hermann_leaf_litter[c("mesocosm", "imidacloprid", "temperature")])
stopifnot(all(tapply(paste(trt$imidacloprid, trt$temperature), trt$mesocosm,
                     function(v) length(unique(v))) == 1L))
treated <- unique(trt$mesocosm[trt$imidacloprid == "applied"])
stopifnot(setequal(unique(hermann_imidacloprid$mesocosm), treated))

usethis::use_data(hermann_zooplankton, hermann_phytoplankton, hermann_macroinvertebrates,
                  hermann_water_quality, hermann_nutrients, hermann_leaf_litter,
                  hermann_imidacloprid, overwrite = TRUE, compress = "xz")
