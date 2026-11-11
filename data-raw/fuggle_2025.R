# Builds data/fuggle_biofilms.rda from the data of Fuggle et al. (2025).
#
# Source: Fuggle, R., Mayer-Pinto, M., Matias, M., & Marzinelli, Z. (2025).
#   Timing of salinisation and nutrient enrichment drives freshwater
#   microbial community metabolic responses. Dryad.
#   https://doi.org/10.5061/dryad.8pk0p2nz8
#   File: univariate_data.csv (version of 2025-11-04).
#   MD5 of the file used: 5aa658f8bbb33a6e3c2caa441696222f
# Licence: CC0 1.0 (Dryad publishes accepted datasets under CC0:
#   https://datadryad.org/publication_policy).
# Article: Fuggle, R., Matias, M. G., Mayer-Pinto, M., & Marzinelli, E. M.
#   (2025). Repetition and order of exposure to salinization and nutrient
#   enrichment drive freshwater microbial community responses. Ecology,
#   106(11), e70260. https://doi.org/10.1002/ecy.70260
#
# Design (article, Methods): 1000-L outdoor ponds of the Iberian Pond Network
# at Toledo (Spain), 2023. Biofilms grown for 30 days on ceramic tiles in four
# ambient ponds were moved to treatment ponds, three or four per treatment:
# ambient (A); elevated salinity (S; 15 g NaCl/L) or nutrient enrichment (N;
# 10 mg N/L and 1 mg P/L), applied once; the same stressor applied twice (SS,
# NN); or one stressor and, a week later, the other (salinity then nutrients,
# SN; nutrients then salinity, NS). Single-stressor ponds received their
# stressor with the first application. Tiles were sampled 1, 21 and 81 days
# after the second application.
#
# Decisions:
# * Only the treatments that form an interaction of the two stressors are
#   kept: A, N, S, NS and SN. The repeated treatments (NN and SS) repeat a
#   single stressor and are left out.
# * `sequence` gives the order of the two stressors for NS and SN and
#   "single" for the ponds that share their reference cells (A, N and S).
# * Ponds from the sample names ([treatment]_[pond replicate]_[sample point]).
# * Gross and net primary production and community respiration, as in the
#   source file (its README gives nmol O2 per cm2 and hour); the oxygen
#   incubations were not done on day 81, so these values are missing then.
#   Gross production is net production plus respiration in every sample but
#   one (pond A_P3, day 1), where it equals net production minus
#   respiration; it is kept as in the source. The other responses of the
#   file are left out.
#
# Packages used by this script only: usethis.
# Download the file above from Dryad into data-raw/ (it is not kept under
# version control) and run this script from the package root.

path <- file.path("data-raw", "univariate_data.csv")
stopifnot(unname(tools::md5sum(path)) == "5aa658f8bbb33a6e3c2caa441696222f")

raw <- utils::read.csv(path, stringsAsFactors = FALSE)
stopifnot(nrow(raw) == 72L, setequal(raw$Treatment, c("A", "N", "NN", "NS", "S", "SN", "SS")),
          setequal(raw$Days, c(1L, 21L, 81L)),
          all(grepl("^[A-Z]+_P[1-4]_T[2-4]$", raw$Sample_Name)),
          all(sub("_.*$", "", raw$Sample_Name) == raw$Treatment))
raw <- raw[raw$Treatment %in% c("A", "N", "S", "NS", "SN"), ]

fuggle_biofilms <- data.frame(
  pond = sub("_T[2-4]$", "", raw$Sample_Name),
  treatment = raw$Treatment,
  sequence = ifelse(raw$Treatment == "NS", "nutrients first",
                    ifelse(raw$Treatment == "SN", "salinity first", "single")),
  nutrients = ifelse(grepl("N", raw$Treatment, fixed = TRUE), "added", "none"),
  salinity = ifelse(grepl("S", raw$Treatment, fixed = TRUE), "raised", "none"),
  day = as.integer(raw$Days),
  gross_production = raw$GPP,
  net_production = raw$NPP,
  respiration = raw$CR,
  stringsAsFactors = FALSE
)
fuggle_biofilms <- fuggle_biofilms[order(fuggle_biofilms$day, fuggle_biofilms$pond,
                                         method = "radix"), ]
rownames(fuggle_biofilms) <- NULL

# Each pond is sampled once on each day; the sample codes T2, T3 and T4 are
# days 1, 21 and 81; gross production is net production plus respiration
# except in the sample noted above.
oxygen <- c("gross_production", "net_production", "respiration")
measured <- fuggle_biofilms$day != 81L
d <- fuggle_biofilms[measured, ]
sum_ok <- abs(d$gross_production - d$net_production - d$respiration) < 1e-6
stopifnot(
  nrow(fuggle_biofilms) == 54L, all(table(fuggle_biofilms$pond) == 3L),
  identical(as.vector(table(raw$Days, sub("^.*_", "", raw$Sample_Name))),
            c(18L, 0L, 0L, 0L, 18L, 0L, 0L, 0L, 18L)),
  !anyNA(d), all(is.na(fuggle_biofilms[!measured, oxygen])),
  identical(paste(d$pond, d$day)[!sum_ok], "A_P3 1"),
  abs(d$gross_production[!sum_ok] - d$net_production[!sum_ok] + d$respiration[!sum_ok]) < 1e-6
)

usethis::use_data(fuggle_biofilms, overwrite = TRUE, compress = "xz")
