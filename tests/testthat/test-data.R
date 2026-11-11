# The data sets match their documentation (R/data.R).

# Taxon columns: those named in the attribute "taxa", after the design columns
# and the total.
expect_community <- function(d, design, total, n_taxa) {
  taxa <- attr(d, "taxa")
  expect_type(taxa, "character")
  expect_length(taxa, n_taxa)
  expect_identical(names(d), c(design, total, names(taxa)))
  expect_false(anyNA(taxa))
  expect_true(all(nzchar(taxa)))
  expect_true(all(grepl("^[a-z][a-z0-9_]*$", names(taxa))))
  values <- as.matrix(d[names(taxa)])
  expect_true(all(values >= 0, na.rm = TRUE))
  expect_equal(d[[total]], unname(rowSums(values)))
  # Sorted by day and mesocosm, with plain row names
  expect_identical(order(d$day, d$mesocosm, method = "radix"), seq_len(nrow(d)))
  expect_lt(.row_names_info(d), 0L)
}

# One treatment per mesocosm and `per_cell` mesocosms per combination.
expect_factorial <- function(d, stressors, per_cell, n_cells) {
  trt <- unique(d[c("mesocosm", stressors)])
  expect_false(anyDuplicated(trt$mesocosm) > 0L)
  counts <- table(do.call(paste, trt[stressors]))
  expect_length(counts, n_cells)
  expect_true(all(counts == per_cell))
}

polazzo_design <- c("mesocosm", "day", "insecticide", "herbicide", "nutrients")
polazzo_stressors <- c("insecticide", "herbicide", "nutrients")

test_that("polazzo_zooplankton matches its documentation", {
  d <- polazzo_zooplankton
  expect_equal(dim(d), c(120L, 28L))
  expect_community(d, polazzo_design, "total_density", 22L)
  # The columns of the first version of the data set are kept
  expect_true(all(c("daphnia", "cyclopoida") %in% names(d)))
  expect_identical(attr(d, "taxa")[["cyclopoida"]], "Ciclopoida")
  expect_false(anyNA(d))
  expect_setequal(unique(d$day), c(-5L, 7L, 14L, 30L, 50L))
  expect_true(all(table(d$mesocosm) == 5L))
  expect_setequal(unique(d$mesocosm), sprintf("M%02d", 1:24))
  expect_setequal(unique(d$insecticide), c("chlorpyrifos", "absent"))
  expect_setequal(unique(d$herbicide), c("diuron", "absent"))
  expect_setequal(unique(d$nutrients), c("added", "ambient"))
  expect_factorial(d, polazzo_stressors, 3L, 8L)
})

test_that("polazzo_phytoplankton matches its documentation", {
  d <- polazzo_phytoplankton
  expect_equal(dim(d), c(96L, 77L))
  expect_community(d, polazzo_design, "total_density", 71L)
  expect_false(anyNA(d))
  expect_setequal(unique(d$day), c(-5L, 7L, 15L, 50L))
  expect_true(all(table(d$mesocosm) == 4L))
  expect_factorial(d, polazzo_stressors, 3L, 8L)
})

test_that("polazzo_macroinvertebrates matches its documentation", {
  d <- polazzo_macroinvertebrates
  expect_equal(dim(d), c(96L, 54L))
  expect_community(d, polazzo_design, "total_abundance", 48L)
  expect_setequal(unique(d$day), c(-32L, 15L, 30L, 50L))
  expect_true(all(table(d$mesocosm) == 4L))
  expect_factorial(d, polazzo_stressors, 3L, 8L)
  # One missing count, and so the total of that sample
  expect_equal(sum(is.na(d)), 2L)
  missing <- d[is.na(d$total_abundance), ]
  expect_identical(c(missing$mesocosm, missing$day), c("M23", "50"))
  expect_true(is.na(missing$hydropsychidae))
  values <- as.matrix(d[names(attr(d, "taxa"))])
  expect_true(all(values == round(values), na.rm = TRUE))
})

test_that("the Polazzo tables share their mesocosms and treatments", {
  trt <- function(d) unique(d[c("mesocosm", polazzo_stressors)])
  sorted <- function(x) x[order(x$mesocosm), , drop = FALSE]
  zoo <- sorted(trt(polazzo_zooplankton))
  rownames(zoo) <- NULL
  for (d in list(polazzo_phytoplankton, polazzo_macroinvertebrates)) {
    other <- sorted(trt(d))
    rownames(other) <- NULL
    expect_identical(other, zoo)
  }
})

hermann_design <- c("mesocosm", "day", "imidacloprid", "temperature")
hermann_stressors <- c("imidacloprid", "temperature")
hermann_days <- c(-4L, 3L, 10L, 17L, 24L, 31L, 38L, 52L, 66L, 80L)

test_that("the Hermann community tables match their documentation", {
  cases <- list(
    list(d = hermann_zooplankton, total = "total_density", taxa = 21L,
         days = c(-1L, 10L, 24L, 38L, 52L, 66L, 80L)),
    list(d = hermann_phytoplankton, total = "total_density", taxa = 91L,
         days = c(-1L, 10L, 24L, 38L, 52L, 66L, 80L)),
    list(d = hermann_macroinvertebrates, total = "total_abundance", taxa = 29L,
         days = c(-7L, 7L, 21L, 35L, 49L, 63L, 77L))
  )
  for (case in cases) {
    d <- case$d
    expect_equal(dim(d), c(168L, 4L + 1L + case$taxa))
    expect_community(d, hermann_design, case$total, case$taxa)
    expect_false(anyNA(d))
    expect_setequal(unique(d$day), case$days)
    expect_true(all(table(d$mesocosm) == 7L))
    expect_setequal(unique(d$imidacloprid), c("absent", "applied"))
    expect_setequal(unique(d$temperature), c("ambient", "elevated", "heatwaves"))
    expect_factorial(d, hermann_stressors, 4L, 6L)
    # The code is the treatment and the replicate number
    code <- paste0(c(absent = "C0", applied = "C1")[d$imidacloprid],
                   c(ambient = "A", elevated = "E", heatwaves = "HW")[d$temperature])
    expect_identical(sub("[0-9]+$", "", d$mesocosm), unname(code))
    expect_setequal(sub("^[A-Z0-9]*[A-Z]", "", d$mesocosm), as.character(1:4))
  }
  expect_identical(attr(hermann_zooplankton, "taxa")[["cyclopoida"]], "Ciclopoida")
  values <- as.matrix(hermann_macroinvertebrates[names(attr(hermann_macroinvertebrates, "taxa"))])
  expect_true(all(values == round(values)))
  # The three communities code the same 24 mesocosms
  expect_setequal(unique(hermann_phytoplankton$mesocosm), unique(hermann_zooplankton$mesocosm))
  expect_setequal(unique(hermann_macroinvertebrates$mesocosm), unique(hermann_zooplankton$mesocosm))
})

test_that("hermann_water_quality matches its documentation", {
  d <- hermann_water_quality
  expect_equal(dim(d), c(480L, 14L))
  expect_identical(names(d), c(hermann_design, "time_of_day", "water_temperature", "oxygen",
                               "oxygen_saturation", "conductivity", "ph", "tds",
                               "chlorophyll_a", "phycoerythrin", "ecosystem_productivity"))
  expect_setequal(unique(d$mesocosm), sprintf("M%02d", 1:24))
  expect_setequal(unique(d$day), hermann_days)
  expect_identical(as.vector(table(d$time_of_day)), c(240L, 240L))
  expect_factorial(d, hermann_stressors, 4L, 6L)
  # Productivity: evening rows only, evening minus morning oxygen
  morning <- d[d$time_of_day == "morning", ]
  evening <- d[d$time_of_day == "evening", ]
  expect_true(all(is.na(morning$ecosystem_productivity)))
  expect_identical(evening$mesocosm, morning$mesocosm)
  expect_identical(evening$day, morning$day)
  expect_equal(evening$ecosystem_productivity, evening$oxygen - morning$oxygen)
  expect_false(anyNA(d[names(d) != "ecosystem_productivity"]))
  expect_true(all(d$ph > 6 & d$ph < 11))
})

test_that("hermann_nutrients matches its documentation", {
  d <- hermann_nutrients
  expect_equal(dim(d), c(240L, 7L))
  expect_identical(names(d), c(hermann_design, "nitrate", "orthophosphate", "total_phosphorus"))
  expect_setequal(unique(d$day), hermann_days)
  expect_true(all(table(d$mesocosm) == 10L))
  expect_factorial(d, hermann_stressors, 4L, 6L)
  missing <- d[is.na(d$nitrate), ]
  expect_identical(missing$mesocosm, c("M07", "M09"))
  expect_identical(missing$day, c(52L, 52L))
  expect_equal(sum(is.na(d)), 2L)
})

test_that("hermann_leaf_litter matches its documentation", {
  d <- hermann_leaf_litter
  expect_equal(dim(d), c(264L, 8L))
  expect_identical(names(d), c(hermann_design, "process", "leaf_dry_mass", "mass_loss",
                               "mass_loss_percent"))
  expect_setequal(unique(d$day[d$process == "microbial"]), c(11L, 25L, 39L, 53L, 67L, 81L))
  expect_setequal(unique(d$day[d$process == "macroinvertebrate feeding"]),
                  c(11L, 25L, 39L, 53L, 81L))
  expect_identical(as.vector(table(d$process)), c(120L, 144L))
  expect_factorial(d, hermann_stressors, 4L, 6L)
  # Each bag held 2 g
  expect_true(all(abs(d$leaf_dry_mass + d$mass_loss - 2) < 0.005))
  ok <- !is.na(d$mass_loss_percent)
  expect_equal(d$mass_loss_percent[ok], 100 * d$mass_loss[ok] / 2, tolerance = 1e-3)
  missing <- d[!ok, ]
  expect_identical(missing$mesocosm, c("M18", "M18"))
  expect_equal(sum(is.na(d)), 2L)
})

test_that("hermann_imidacloprid matches its documentation", {
  d <- hermann_imidacloprid
  expect_equal(dim(d), c(120L, 4L))
  expect_identical(names(d), c("mesocosm", "day", "temperature", "concentration"))
  expect_setequal(unique(d$day), c(0L, 3L, 10L, 14L, 17L, 24L, 28L, 31L, 38L, 45L))
  expect_true(all(table(d$mesocosm) == 10L))
  # Only the treated mesocosms, with the treatments of the other tables
  treated <- unique(hermann_water_quality[hermann_water_quality$imidacloprid == "applied",
                                          c("mesocosm", "temperature")])
  expect_setequal(unique(d$mesocosm), treated$mesocosm)
  expect_identical(unname(d$temperature),
                   unname(treated$temperature[match(d$mesocosm, treated$mesocosm)]))
  missing <- d[is.na(d$concentration), ]
  expect_identical(missing$mesocosm, rep("M20", 3L))
  expect_identical(missing$day, c(28L, 31L, 38L))
  expect_true(all(d$concentration > 0, na.rm = TRUE))
})

test_that("the Hermann tables numbered by mesocosm share their treatments", {
  trt <- function(d) {
    x <- unique(d[c("mesocosm", hermann_stressors)])
    x <- x[order(x$mesocosm), ]
    rownames(x) <- NULL
    x
  }
  reference <- trt(hermann_water_quality)
  expect_identical(trt(hermann_nutrients), reference)
  expect_identical(trt(hermann_leaf_litter), reference)
})

# Data of the case studies ---------------------------------------------------

# Each unit keeps the same value of `columns` in all its rows.
expect_constant_within <- function(d, unit, columns) {
  values <- do.call(paste, d[columns])
  expect_true(all(tapply(values, d[[unit]], function(x) length(unique(x))) == 1L))
}

test_that("alba_arthropods matches its documentation", {
  d <- alba_arthropods
  groups <- c("hemipterans", "caterpillars", "beetles", "grasshoppers", "flies", "moths",
              "bees", "wasps", "spiders", "ants")
  expect_equal(dim(d), c(160L, 18L))
  expect_identical(names(d), c("plot", "block", "sampling", "year", "season", "drought",
                               "invasion", "total", groups))
  expect_identical(names(attr(d, "taxa")), groups)
  expect_true(all(vapply(d[c("plot", "block", "sampling", "year", "total", groups)],
                         is.integer, logical(1))))
  expect_equal(d$total, unname(rowSums(d[groups])))
  expect_false(anyNA(d))
  expect_identical(order(d$sampling, d$plot, method = "radix"), seq_len(nrow(d)))
  expect_identical(unique(paste(d$sampling, d$season, d$year)),
                   c("1 summer 2014", "2 fall 2014", "3 summer 2015", "4 fall 2015"))
  expect_constant_within(d, "plot", c("block", "drought", "invasion"))
  # Ten complete blocks of the 2 x 2 factorial, four samplings
  first <- d[d$sampling == 1L, ]
  expect_true(all(table(first$block, paste(first$drought, first$invasion)) == 1L))
  expect_true(all(table(d$plot) == 4L))
  expect_identical(d$invasion[d$plot == 38L], rep("invaded", 4L))
})

test_that("fuggle_biofilms matches its documentation", {
  d <- fuggle_biofilms
  expect_equal(dim(d), c(54L, 9L))
  expect_identical(names(d), c("pond", "treatment", "sequence", "nutrients", "salinity", "day",
                               "gross_production", "net_production", "respiration"))
  expect_setequal(d$day, c(1L, 21L, 81L))
  expect_true(all(table(d$pond) == 3L))
  expect_constant_within(d, "pond", c("treatment", "sequence", "nutrients", "salinity"))
  ponds <- unique(d[c("pond", "treatment", "sequence", "nutrients", "salinity")])
  expect_identical(as.vector(table(ponds$treatment)[c("A", "N", "S", "NS", "SN")]),
                   c(4L, 3L, 3L, 4L, 4L))
  expect_identical(unique(ponds$sequence[ponds$treatment %in% c("A", "N", "S")]), "single")
  expect_identical(unique(ponds$sequence[ponds$treatment == "NS"]), "nutrients first")
  expect_identical(unique(ponds$sequence[ponds$treatment == "SN"]), "salinity first")
  expect_identical(ponds$nutrients == "added", grepl("N", ponds$treatment))
  expect_identical(ponds$salinity == "raised", grepl("S", ponds$treatment))
  # The oxygen rates are missing on day 81 only
  rates <- c("gross_production", "net_production", "respiration")
  expect_true(all(is.na(d[d$day == 81L, rates])))
  expect_false(anyNA(d[d$day != 81L, ]))
  ok <- d[d$day != 81L, ]
  sum_ok <- abs(ok$gross_production - ok$net_production - ok$respiration) < 1e-6
  expect_identical(paste(ok$pond, ok$day)[!sum_ok], "A_P3 1")
})

test_that("giesler_diatoms matches its documentation", {
  d <- giesler_diatoms
  expect_equal(dim(d), c(1728L, 6L))
  expect_identical(names(d), c("species", "microbiome", "temperature", "photoperiod",
                               "replicate", "growth_rate"))
  expect_setequal(d$species, c("T. gravida", "T. rotula"))
  expect_setequal(d$microbiome, c("present", "removed"))
  expect_setequal(d$temperature, c(4, 9, 13.5))
  expect_setequal(d$photoperiod, c(4L, 16L, 24L))
  cells <- table(do.call(paste, d[c("species", "microbiome", "temperature", "photoperiod")]))
  expect_length(cells, 36L)
  expect_true(all(cells == 48L))
  expect_false(anyNA(d))
  expect_true(any(d$growth_rate < 0))
})

test_that("vasquez_embryos matches its documentation", {
  d <- vasquez_embryos
  expect_equal(dim(d), c(216L, 6L))
  expect_identical(names(d), c("temperature", "salinity", "oxygen", "pair", "embryos",
                               "developed"))
  expect_true(all(vapply(d, is.integer, logical(1))))
  expect_true(all(table(d$temperature, d$salinity, d$oxygen) == 6L))
  expect_setequal(d$pair, 1:6)
  expect_identical(sort(unique(d$embryos)), c(7L, 8L))
  expect_identical(sum(d$embryos == 7L), 3L)
  expect_true(all(d$developed >= 0L & d$developed <= d$embryos))
  expect_false(anyNA(d))
  # Reference condition: 30 C, 34 ppt and 21% O2
  reference <- d[d$temperature == 30L & d$salinity == 34L & d$oxygen == 21L, ]
  expect_identical(nrow(reference), 6L)
})

test_that("hunn_invertebrates matches its documentation", {
  d <- hunn_invertebrates
  expect_equal(dim(d), c(128L, 22L))
  taxa <- attr(d, "taxa")
  expect_length(taxa, 14L)
  expect_identical(names(d), c("channel", "block", "header_tank", "temperature", "co2",
                               "sediment", "flow", "total_abundance", names(taxa)))
  expect_identical(d$channel, 1:128)
  expect_false(anyNA(d))
  # Eight header tanks of 16 channels, one temperature and CO2 level each,
  # two per combination of the two, one in each block
  expect_true(all(table(d$header_tank) == 16L))
  expect_constant_within(d, "header_tank", c("block", "temperature", "co2"))
  tanks <- unique(d[c("header_tank", "block", "temperature", "co2")])
  expect_true(all(table(tanks$temperature, tanks$co2) == 2L))
  expect_true(all(table(tanks$block) == 4L))
  cells <- table(do.call(paste, d[c("temperature", "co2", "sediment", "flow")]))
  expect_length(cells, 16L)
  expect_true(all(cells == 8L))
  counts <- as.matrix(d[names(taxa)])
  expect_true(all(counts >= 0))
  expect_true(all(2 * counts == round(2 * counts)))
  expect_true(all(d$total_abundance >= rowSums(counts)))
})

test_that("niboyet_denitrification matches its documentation", {
  d <- niboyet_denitrification
  expect_equal(dim(d), c(640L, 12L))
  expect_identical(names(d), c("block", "plot", "quadrant", "unit", "co2", "temperature",
                               "water", "nitrogen", "fire", "months", "denitrification",
                               "soil_water"))
  expect_setequal(d$months, c(0L, 5L, 9L, 21L, 33L))
  expect_true(all(table(d$unit) == 5L))
  # Split plots: CO2, temperature and fire per plot, water and nitrogen per
  # quadrant; 32 combinations with four quadrants each
  expect_constant_within(d, "plot", c("block", "co2", "temperature", "fire"))
  expect_constant_within(d, "unit", c("plot", "quadrant", "water", "nitrogen"))
  plots <- unique(d[c("plot", "block", "co2", "temperature", "fire")])
  expect_identical(nrow(plots), 32L)
  expect_true(all(table(plots$block) == 8L))
  cells <- table(do.call(paste, d[d$months == 0L, c("co2", "temperature", "water", "nitrogen",
                                                    "fire")]))
  expect_length(cells, 32L)
  expect_true(all(cells == 4L))
  expect_identical(sum(is.na(d$denitrification)), 1L)
  expect_identical(sum(is.na(d$soil_water)), 6L)
  expect_true(all(d$denitrification > 0, na.rm = TRUE))
})
