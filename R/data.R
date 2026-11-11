# Rd items for the taxon columns of a community table: one per column, with the
# name of the taxon in the source file (attribute "taxa"). Called by the inline
# code of the @format sections below when the documentation is generated.
rd_taxa <- function(data) {
  taxa <- attr(data, "taxa")
  escape <- function(x) gsub("([%{}\\\\])", "\\\\\\1", x)
  paste(sprintf("\\item{%s}{%s}", escape(names(taxa)), escape(unname(taxa))), collapse = "\n")
}

# Polazzo & Rico (2021) --------------------------------------------------------

#' Zooplankton in a three-stressor mesocosm experiment
#'
#' Zooplankton densities from an outdoor mesocosm experiment with a full
#' factorial of an insecticide, a herbicide and nutrient enrichment (Polazzo &
#' Rico 2021), sampled once before the pesticide application and four times
#' after it.
#'
#' @section Experimental design:
#' Twenty-four outdoor mesocosms, three per combination of a full factorial of
#' chlorpyrifos (insecticide, 1 µg/L, or absent), diuron (herbicide, 18 µg/L,
#' or absent) and nutrients (750 µg/L N and 75 µg/L P added twice a week, or
#' ambient). The nutrient additions started four weeks before the pesticide
#' application (day 0). The community tables [polazzo_zooplankton],
#' [polazzo_phytoplankton] and [polazzo_macroinvertebrates] share the
#' mesocosm codes and the design columns, so they can be joined by `mesocosm`
#' and, where the sampling days coincide, `day`.
#'
#' @format A data frame with 120 rows (24 mesocosms x 5 sampling days) and 28
#'   columns:
#'   \describe{
#'     \item{mesocosm}{Mesocosm code (text, `"M01"` to `"M24"`), the
#'       experimental unit.}
#'     \item{day}{Day relative to the pesticide application (integer: -5, 7,
#'       14, 30, 50). Day -5 is before the pesticides but not before the
#'       nutrients.}
#'     \item{insecticide}{`"chlorpyrifos"` or `"absent"`.}
#'     \item{herbicide}{`"diuron"` or `"absent"`.}
#'     \item{nutrients}{`"added"` or `"ambient"`.}
#'     \item{total_density}{Sum of the 22 taxon columns.}
#'     `r rd_taxa(polazzo_zooplankton)`
#'   }
#'   The taxon columns hold the density of each taxon of the source file,
#'   whose name for the taxon is given above and in the attribute `"taxa"`.
#'   Densities are as reported in the source file, which does not state
#'   their unit. Zeros are taxa not found in the sample; there are no
#'   missing values.
#'
#' @source Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple
#'   stressors on the dimensionality of ecological stability. Dryad,
#'   \doi{10.5061/dryad.qfttdz0gt}. Licence CC0 1.0. Sheet "Zooplankton
#'   density" of `Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx`; the
#'   treatment codes, which differ on day 50, were harmonised (script
#'   `data-raw/polazzo_rico_2021.R` of the source repository).
#'
#' @references
#' Polazzo, F., & Rico, A. (2021). Effects of multiple stressors on the
#' dimensionality of ecological stability. *Ecology Letters*, 24(8),
#' 1594–1606. \doi{10.1111/ele.13770}
#'
#' @family data of Polazzo and Rico (2021)
#' @examples
#' str(polazzo_zooplankton[1:8])
#' # Three mesocosms per treatment combination
#' with(polazzo_zooplankton[polazzo_zooplankton$day == 7, ],
#'      table(insecticide, herbicide, nutrients))
#' head(attr(polazzo_zooplankton, "taxa"))
"polazzo_zooplankton"

#' Phytoplankton in a three-stressor mesocosm experiment
#'
#' Phytoplankton densities from the mesocosm experiment of Polazzo & Rico
#' (2021) with an insecticide, a herbicide and nutrient enrichment, sampled
#' once before the pesticide application and three times after it.
#'
#' @inheritSection polazzo_zooplankton Experimental design
#'
#' @format A data frame with 96 rows (24 mesocosms x 4 sampling days) and 77
#'   columns:
#'   \describe{
#'     \item{mesocosm, insecticide, herbicide, nutrients}{As in
#'       [polazzo_zooplankton].}
#'     \item{day}{Day relative to the pesticide application (integer: -5, 7,
#'       15, 50).}
#'     \item{total_density}{Sum of the 71 taxon columns.}
#'     `r rd_taxa(polazzo_phytoplankton)`
#'   }
#'   The taxon columns hold the density of each taxon found in the source
#'   file on any sampling day, whose name for the taxon is given above and in
#'   the attribute `"taxa"`. Densities are as reported in the source file,
#'   which does not state their unit. The source has one sheet per sampling
#'   day, each with the taxa found that day in some mesocosm; a taxon missing
#'   from a sheet was not found that day, so its density is zero. Zeros are
#'   taxa not found in the sample; there are no missing values.
#'
#' @source Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple
#'   stressors on the dimensionality of ecological stability. Dryad,
#'   \doi{10.5061/dryad.qfttdz0gt}. Licence CC0 1.0. Sheets "Phyto_D-5",
#'   "Phyto_D7", "Phyto_D15" and "Phyto_D50" of
#'   `Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx` (the day of the
#'   last sheet is taken from its name: its day column says 15); script
#'   `data-raw/polazzo_rico_2021.R` of the source repository.
#'
#' @inherit polazzo_zooplankton references
#' @family data of Polazzo and Rico (2021)
#' @examples
#' range(polazzo_phytoplankton$total_density)
#' # Taxa found in each sample
#' taxa <- names(attr(polazzo_phytoplankton, "taxa"))
#' summary(rowSums(polazzo_phytoplankton[taxa] > 0))
"polazzo_phytoplankton"

#' Macroinvertebrates in a three-stressor mesocosm experiment
#'
#' Macroinvertebrate counts from the mesocosm experiment of Polazzo & Rico
#' (2021) with an insecticide, a herbicide and nutrient enrichment, sampled
#' once before the pesticide application and three times after it.
#'
#' @inheritSection polazzo_zooplankton Experimental design
#'
#' @format A data frame with 96 rows (24 mesocosms x 4 sampling days) and 54
#'   columns:
#'   \describe{
#'     \item{mesocosm, insecticide, herbicide, nutrients}{As in
#'       [polazzo_zooplankton].}
#'     \item{day}{Day relative to the pesticide application (integer: -32,
#'       15, 30, 50). Day -32 is before the nutrient additions too.}
#'     \item{total_abundance}{Sum of the 48 taxon columns.}
#'     `r rd_taxa(polazzo_macroinvertebrates)`
#'   }
#'   The taxon columns hold the number of individuals of each taxon of the
#'   source file, whose name for the taxon is given above and in the
#'   attribute `"taxa"`. Counts are as reported in the source file. Zeros are
#'   taxa not found in the sample. One count is missing (`hydropsychidae`,
#'   mesocosm `"M23"`, day 50), and so is the total of that sample.
#'
#' @source Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple
#'   stressors on the dimensionality of ecological stability. Dryad,
#'   \doi{10.5061/dryad.qfttdz0gt}. Licence CC0 1.0. Sheet
#'   "Macroinvertebrates" of
#'   `Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx`, without an
#'   unlabelled column that has values for one mesocosm per treatment only;
#'   script `data-raw/polazzo_rico_2021.R` of the source repository.
#'
#' @inherit polazzo_zooplankton references
#' @family data of Polazzo and Rico (2021)
#' @examples
#' range(polazzo_macroinvertebrates$total_abundance, na.rm = TRUE)
#' polazzo_macroinvertebrates[is.na(polazzo_macroinvertebrates$total_abundance), 1:6]
"polazzo_macroinvertebrates"

# Hermann et al. (2024) --------------------------------------------------------

#' Zooplankton in an insecticide and temperature mesocosm experiment
#'
#' Zooplankton densities from an outdoor mesocosm experiment with an
#' insecticide (imidacloprid) and three temperature regimes: ambient,
#' elevated and repeated heatwaves (Hermann et al. 2024), sampled once before
#' the first insecticide application and six times after it.
#'
#' @section Experimental design:
#' Twenty-four outdoor mesocosms (about 850 L of water over a sediment layer)
#' at Alcalá de Henares (central Spain), from 19 April (day -7) to 19 July
#' 2021, four per combination of a 2 x 3 factorial: imidacloprid (nominal
#' 1 µg/L, applied on days 0, 14 and 28, or absent) and temperature (ambient;
#' elevated, 4 °C above ambient; or three heatwaves of 7 days, separated by
#' 7 days at ambient temperature, with targets of 24, 28 and 32 °C). The
#' temperature treatments lasted from day 0 to day 42; days 42 to 84 were a
#' recovery phase.
#'
#' The community tables [hermann_zooplankton], [hermann_phytoplankton] and
#' [hermann_macroinvertebrates] code the mesocosms by treatment and
#' replicate number (`"C0A1"`); [hermann_water_quality],
#' [hermann_nutrients], [hermann_leaf_litter] and [hermann_imidacloprid] by
#' mesocosm number (`"M01"`). The source gives no correspondence between the
#' two codings, so the two kinds of table cannot be joined by mesocosm.
#'
#' @format A data frame with 168 rows (24 mesocosms x 7 sampling days) and 26
#'   columns:
#'   \describe{
#'     \item{mesocosm}{Mesocosm code (text) from the sample labels of the
#'       source: imidacloprid (`C0` absent, `C1` applied), temperature (`A`
#'       ambient, `E` elevated, `HW` heatwaves) and replicate number, for
#'       example `"C1HW3"`. The codes appear to identify the same mesocosm
#'       on every sampling day: within a treatment, the abundances of a code
#'       are correlated between consecutive samplings, and those of
#'       different codes are not.}
#'     \item{day}{Day relative to the first imidacloprid application
#'       (integer: -1, 10, 24, 38, 52, 66, 80), from the sample labels; the
#'       article gives day -4 for the first sampling.}
#'     \item{imidacloprid}{`"applied"` or `"absent"`.}
#'     \item{temperature}{`"ambient"`, `"elevated"` or `"heatwaves"`.}
#'     \item{total_density}{Sum of the 21 taxon columns.}
#'     `r rd_taxa(hermann_zooplankton)`
#'   }
#'   The taxon columns hold the density of each taxon of the source file
#'   (individuals per litre), whose name for the taxon is given above and in
#'   the attribute `"taxa"`. Zeros are taxa not found in the sample; there are
#'   no missing values.
#'
#' @source Hermann, M., Polazzo, F., Cherta Cucala, L., Crettaz Minaglia,
#'   M. C., García-Astillero, A., Peeters, E. T. H. M., Rico, A., & Van den
#'   Brink, P. J. (2024). Data supporting: Combined stress of an insecticide
#'   and heatwaves or elevated temperature induce community and food web
#'   effects in a Mediterranean freshwater ecosystem (version V.1). Zenodo,
#'   \doi{10.5281/zenodo.11520459}. Licence CC BY 4.0. Sheets
#'   "ZooplanktonDesign" and "Zooplankton_Abundances" of
#'   `CombinedStressInsecticideTemperature_Hermann_et_al.(2024).xlsx`; script
#'   `data-raw/hermann_2024.R` of the source repository.
#'
#' @references
#' Hermann, M., Polazzo, F., Cherta, L., Crettaz-Minaglia, M.,
#' García-Astillero, A., Peeters, E. T. H. M., Rico, A., & Van den Brink,
#' P. J. (2024). Combined stress of an insecticide and heatwaves or elevated
#' temperature induce community and food web effects in a Mediterranean
#' freshwater ecosystem. *Water Research*, 260, 121903.
#' \doi{10.1016/j.watres.2024.121903}
#'
#' @family data of Hermann et al. (2024)
#' @examples
#' str(hermann_zooplankton[1:9])
#' # Four mesocosms per treatment combination on each sampling day
#' with(hermann_zooplankton, table(imidacloprid, temperature, day))[, , "24"]
#' # The design as must sees it
#' describe_design(hermann_zooplankton[hermann_zooplankton$day > 0, ],
#'                 stressors = c("imidacloprid", "temperature"),
#'                 control = c(imidacloprid = "absent", temperature = "ambient"),
#'                 unit = "mesocosm", time = "day")
"hermann_zooplankton"

#' Phytoplankton in an insecticide and temperature mesocosm experiment
#'
#' Phytoplankton densities from the mesocosm experiment of Hermann et al.
#' (2024) with an insecticide and three temperature regimes, sampled once
#' before the first insecticide application and six times after it.
#'
#' @inheritSection hermann_zooplankton Experimental design
#'
#' @format A data frame with 168 rows (24 mesocosms x 7 sampling days) and 96
#'   columns:
#'   \describe{
#'     \item{mesocosm, imidacloprid, temperature}{As in
#'       [hermann_zooplankton].}
#'     \item{day}{Day relative to the first imidacloprid application
#'       (integer: -1, 10, 24, 38, 52, 66, 80).}
#'     \item{total_density}{Sum of the 91 taxon columns.}
#'     `r rd_taxa(hermann_phytoplankton)`
#'   }
#'   The taxon columns hold the density of each taxon of the source file
#'   (ciliates included), whose name for the taxon is given above and in the
#'   attribute `"taxa"`. Densities are as reported in the source file, which
#'   does not state their unit. Zeros are taxa not found in the sample; there
#'   are no missing values.
#'
#' @source Hermann et al. (2024), Zenodo, \doi{10.5281/zenodo.11520459}
#'   (full reference in [hermann_zooplankton]). Licence CC BY 4.0. Sheets
#'   "PhytoplanktonDesign" and "Phytoplankton_Abundances"; script
#'   `data-raw/hermann_2024.R` of the source repository.
#'
#' @inherit hermann_zooplankton references
#' @family data of Hermann et al. (2024)
#' @examples
#' # Total density, median by treatment combination and day
#' with(hermann_phytoplankton,
#'      tapply(total_density, list(paste(imidacloprid, temperature), day), median))
"hermann_phytoplankton"

#' Macroinvertebrates in an insecticide and temperature mesocosm experiment
#'
#' Macroinvertebrate counts from the mesocosm experiment of Hermann et al.
#' (2024) with an insecticide and three temperature regimes, sampled once
#' before the first insecticide application and six times after it.
#'
#' @inheritSection hermann_zooplankton Experimental design
#'
#' @format A data frame with 168 rows (24 mesocosms x 7 sampling days) and 34
#'   columns:
#'   \describe{
#'     \item{mesocosm, imidacloprid, temperature}{As in
#'       [hermann_zooplankton].}
#'     \item{day}{Day relative to the first imidacloprid application
#'       (integer: -7, 7, 21, 35, 49, 63, 77).}
#'     \item{total_abundance}{Sum of the 29 taxon columns.}
#'     `r rd_taxa(hermann_macroinvertebrates)`
#'   }
#'   The taxon columns hold the number of individuals of each taxon of the
#'   source file, whose name for the taxon is given above and in the
#'   attribute `"taxa"`; the article describes four bottom traps per
#'   mesocosm. Counts are as reported in the source file. Zeros are taxa not
#'   found in the sample; there are no missing values.
#'
#' @source Hermann et al. (2024), Zenodo, \doi{10.5281/zenodo.11520459}
#'   (full reference in [hermann_zooplankton]). Licence CC BY 4.0. Sheets
#'   "MacroinvertebratesDesign" and "Macroinvertebrates_Abundances"; script
#'   `data-raw/hermann_2024.R` of the source repository.
#'
#' @inherit hermann_zooplankton references
#' @family data of Hermann et al. (2024)
#' @examples
#' range(hermann_macroinvertebrates$total_abundance)
#' # Samples in which each taxon was found
#' taxa <- names(attr(hermann_macroinvertebrates, "taxa"))
#' sort(colSums(hermann_macroinvertebrates[taxa] > 0), decreasing = TRUE)
"hermann_macroinvertebrates"

#' Water quality in an insecticide and temperature mesocosm experiment
#'
#' Physico-chemical variables measured in the morning and in the evening in
#' the mesocosm experiment of Hermann et al. (2024) with an insecticide and
#' three temperature regimes.
#'
#' @inheritSection hermann_zooplankton Experimental design
#'
#' @format A data frame with 480 rows (24 mesocosms x 10 days x 2 times of
#'   day) and 14 columns:
#'   \describe{
#'     \item{mesocosm}{Mesocosm number (text, `"M01"` to `"M24"`), the
#'       experimental unit; the same coding in [hermann_nutrients],
#'       [hermann_leaf_litter] and [hermann_imidacloprid].}
#'     \item{day}{Day relative to the first imidacloprid application
#'       (integer: -4, 3, 10, 17, 24, 31, 38, 52, 66, 80).}
#'     \item{imidacloprid}{`"applied"` or `"absent"`.}
#'     \item{temperature}{`"ambient"`, `"elevated"` or `"heatwaves"`.}
#'     \item{time_of_day}{`"morning"` or `"evening"`.}
#'     \item{water_temperature}{Water temperature (°C).}
#'     \item{oxygen}{Dissolved oxygen (mg/L).}
#'     \item{oxygen_saturation}{Oxygen saturation (%).}
#'     \item{conductivity}{Electrical conductivity (µS/cm).}
#'     \item{ph}{pH.}
#'     \item{tds}{Total dissolved solids (ppm).}
#'     \item{chlorophyll_a, phycoerythrin}{Chlorophyll-a and phycoerythrin,
#'       as reported in the source file, which does not state their unit.}
#'     \item{ecosystem_productivity}{Evening minus morning dissolved oxygen
#'       of the same day (mg/L), as given in the source; missing in the
#'       morning rows.}
#'   }
#'   There are no other missing values.
#'
#' @source Hermann et al. (2024), Zenodo, \doi{10.5281/zenodo.11520459}
#'   (full reference in [hermann_zooplankton]). Licence CC BY 4.0. Sheets
#'   "WaterQualityMorning" and "WaterQualityEvening"; script
#'   `data-raw/hermann_2024.R` of the source repository.
#'
#' @inherit hermann_zooplankton references
#' @family data of Hermann et al. (2024)
#' @examples
#' # Mean evening water temperature by temperature regime and day
#' with(hermann_water_quality[hermann_water_quality$time_of_day == "evening", ],
#'      round(tapply(water_temperature, list(temperature, day), mean), 1))
"hermann_water_quality"

#' Nutrients in an insecticide and temperature mesocosm experiment
#'
#' Nitrate, orthophosphate and total phosphorus in the mesocosm experiment of
#' Hermann et al. (2024) with an insecticide and three temperature regimes.
#'
#' @inheritSection hermann_zooplankton Experimental design
#'
#' @format A data frame with 240 rows (24 mesocosms x 10 days) and 7 columns:
#'   \describe{
#'     \item{mesocosm, imidacloprid, temperature}{As in
#'       [hermann_water_quality].}
#'     \item{day}{Day relative to the first imidacloprid application
#'       (integer: -4, 3, 10, 17, 24, 31, 38, 52, 66, 80).}
#'     \item{nitrate, orthophosphate, total_phosphorus}{Concentrations as
#'       reported in the source file, which does not state their units.}
#'   }
#'   Nitrate is missing on day 52 in mesocosms `"M07"` and `"M09"`; there are
#'   no other missing values.
#'
#' @source Hermann et al. (2024), Zenodo, \doi{10.5281/zenodo.11520459}
#'   (full reference in [hermann_zooplankton]). Licence CC BY 4.0. Sheet
#'   "Nutrients"; script `data-raw/hermann_2024.R` of the source repository.
#'
#' @inherit hermann_zooplankton references
#' @family data of Hermann et al. (2024)
#' @examples
#' with(hermann_nutrients, round(tapply(total_phosphorus, list(temperature, day), median)))
"hermann_nutrients"

#' Leaf litter decomposition in an insecticide and temperature mesocosm
#' experiment
#'
#' Dry mass lost by poplar leaves in litter bags, by microbial decomposition
#' and by macroinvertebrate feeding, in the mesocosm experiment of Hermann et
#' al. (2024) with an insecticide and three temperature regimes.
#'
#' @inheritSection hermann_zooplankton Experimental design
#'
#' @format A data frame with 264 rows (24 mesocosms x 6 days of microbial
#'   decomposition and 5 of macroinvertebrate feeding) and 8 columns:
#'   \describe{
#'     \item{mesocosm, imidacloprid, temperature}{As in
#'       [hermann_water_quality].}
#'     \item{day}{Day of collection relative to the first imidacloprid
#'       application (integer: 11, 25, 39, 53, 67 and 81 for microbial
#'       decomposition, and the same days except 67 for macroinvertebrate
#'       feeding). The article reports collections after 14 to 84 days of
#'       exposure.}
#'     \item{process}{`"microbial"` (fine-mesh bags, 500 µm) or
#'       `"macroinvertebrate feeding"` (coarse-mesh bags, up to 0.5 cm).}
#'     \item{leaf_dry_mass}{Dry mass of the leaves left in the bag (g); each
#'       bag held 2 g.}
#'     \item{mass_loss}{Dry mass lost (g).}
#'     \item{mass_loss_percent}{Dry mass lost as a percentage of the initial
#'       2 g.}
#'   }
#'   Values are as reported in the source file. Two samples of mesocosm
#'   `"M18"` have no percentage in the source and are missing: microbial
#'   decomposition on day 11 (no loss) and macroinvertebrate feeding on day
#'   25 (a gain).
#'
#' @source Hermann et al. (2024), Zenodo, \doi{10.5281/zenodo.11520459}
#'   (full reference in [hermann_zooplankton]). Licence CC BY 4.0. Sheets
#'   "Microbial_Decomposition" and "Macroinvertebrate_Feeding"; script
#'   `data-raw/hermann_2024.R` of the source repository.
#'
#' @inherit hermann_zooplankton references
#' @family data of Hermann et al. (2024)
#' @examples
#' with(hermann_leaf_litter,
#'      round(tapply(mass_loss_percent, list(process, day), mean, na.rm = TRUE), 1))
"hermann_leaf_litter"

#' Imidacloprid in the water of an insecticide and temperature mesocosm
#' experiment
#'
#' Measured imidacloprid concentrations in the treated mesocosms of the
#' experiment of Hermann et al. (2024) with an insecticide and three
#' temperature regimes.
#'
#' @inheritSection hermann_zooplankton Experimental design
#'
#' @format A data frame with 120 rows (12 treated mesocosms x 10 days) and 4
#'   columns:
#'   \describe{
#'     \item{mesocosm}{Mesocosm number (text), as in
#'       [hermann_water_quality]; only those with imidacloprid.}
#'     \item{day}{Day relative to the first imidacloprid application
#'       (integer: 0, 3, 10, 14, 17, 24, 28, 31, 38, 45). The article
#'       describes the samples of day 0 as taken 4 hours after the
#'       application.}
#'     \item{temperature}{`"ambient"`, `"elevated"` or `"heatwaves"`.}
#'     \item{concentration}{Imidacloprid in water (µg/L).}
#'   }
#'   The concentration is missing in mesocosm `"M20"` on days 28, 31 and 38;
#'   there are no other missing values.
#'
#' @source Hermann et al. (2024), Zenodo, \doi{10.5281/zenodo.11520459}
#'   (full reference in [hermann_zooplankton]). Licence CC BY 4.0. Sheet
#'   "ImiConcentrationsWater"; script `data-raw/hermann_2024.R` of the source
#'   repository.
#'
#' @inherit hermann_zooplankton references
#' @family data of Hermann et al. (2024)
#' @examples
#' with(hermann_imidacloprid, round(tapply(concentration, list(temperature, day),
#'                                         mean, na.rm = TRUE), 2))
"hermann_imidacloprid"

# Case studies (article of the package website) -------------------------------

#' Arthropods in a drought and plant invasion field experiment
#'
#' Arthropod counts from a field experiment with a full factorial of
#' simulated drought and invasion by cogongrass (*Imperata cylindrica*) in
#' regenerating longleaf pine understory (Alba et al. 2026), sampled twice a
#' year for two years.
#'
#' @section Experimental design:
#' Forty plots of 4 x 4 m in an old field at Gainesville (Florida, USA),
#' planted in 2012 with longleaf pine and native understory species, in 10
#' blocks with one plot of each combination of drought (rainout shelters, or
#' ambient precipitation under control structures) and invasion (nine
#' cogongrass ramets planted, or none). Both treatments started in spring
#' 2013. Arthropods were sampled with 20 sweeps of a net per plot in July
#' (summer) and October (fall) of 2014 and 2015.
#'
#' @format A data frame with 160 rows (40 plots x 4 samplings) and 18
#'   columns:
#'   \describe{
#'     \item{plot}{Plot number (integer, 1 to 40), the experimental unit.}
#'     \item{block}{Block number (integer, 1 to 10).}
#'     \item{sampling}{Sampling number (integer, 1 to 4) in chronological
#'       order: summer 2014, fall 2014, summer 2015 and fall 2015.}
#'     \item{year}{Year of the sampling (integer, 2014 or 2015).}
#'     \item{season}{`"summer"` (July) or `"fall"` (October).}
#'     \item{drought}{`"drought"` or `"ambient"`.}
#'     \item{invasion}{`"invaded"` or `"uninvaded"`.}
#'     \item{total}{Sum of the 10 group columns.}
#'     \item{hemipterans}{Hemiptera (source code `hemtotal`).}
#'     \item{caterpillars}{Lepidoptera larvae (`lepnumlarva`).}
#'     \item{beetles}{Coleoptera (`coletotal`).}
#'     \item{grasshoppers}{Orthoptera (`orthtotal`).}
#'     \item{flies}{Adult Diptera (`dipnumadult`).}
#'     \item{moths}{Adult Lepidoptera (`lepnumadult`).}
#'     \item{bees}{Bees (`beetotal`).}
#'     \item{wasps}{Wasps (`wasptotal`).}
#'     \item{spiders}{Spiders (`spidnumindiv`).}
#'     \item{ants}{Ants (`anttotal`).}
#'   }
#'   The group columns, named after the article, count the individuals of
#'   each group in the sample of a plot; the attribute `"taxa"` gives the
#'   code of each group in the source file. The source assigns hemipterans,
#'   caterpillars, beetles and grasshoppers to herbivores, flies, moths and
#'   bees to pollinators, and wasps, spiders and ants to predators. Zeros are
#'   groups not found in the sample; there are no missing values.
#'
#' @source Alba, C., Lucky, A., Fahey, C., Somarriba, G., Dunbar, M., Inman,
#'   W., Flory, S. L., & Angelini, C. (2025). Data from: Plant invasion
#'   overrides drought in suppressing arthropod abundance and richness.
#'   Dryad, \doi{10.5061/dryad.02v6wwqhz}. Licence CC0 1.0. Sheet
#'   "Abundance_Raw_Data" of `Data_For_Repository_Arthropods.xlsx`, one
#'   column per group; treatments from the treatment code, which in the
#'   fourth sampling disagrees with the invasion column for plot 38 (script
#'   `data-raw/alba_2026.R` of the source repository).
#'
#' @references
#' Alba, C., Lucky, A., Fahey, C., Somarriba, G., Dunbar, M., Inman, W.,
#' Flory, S. L., & Angelini, C. (2026). Plant invasion overrides drought in
#' suppressing arthropod abundance and richness. *Oikos*, 2026(5), e11920.
#' \doi{10.1002/oik.11920}
#'
#' @family case-study data
#' @examples
#' str(alba_arthropods[1:9])
#' # One plot of each treatment in every block
#' with(alba_arthropods[alba_arthropods$sampling == 1, ], table(block, paste(drought, invasion)))
#' describe_design(alba_arthropods, stressors = c("drought", "invasion"),
#'                 control = c(drought = "ambient", invasion = "uninvaded"),
#'                 unit = "plot", time = "sampling", groups = "block")
"alba_arthropods"

#' Biofilm metabolism under salinisation and nutrient enrichment in two orders
#'
#' Primary production and respiration of benthic biofilms in outdoor ponds
#' exposed to elevated salinity, nutrient enrichment, or both in either
#' order (Fuggle et al. 2025).
#'
#' @section Experimental design:
#' Outdoor ponds of 1000 L of the Iberian Pond Network at Toledo (Spain),
#' 2023. Biofilms grown for 30 days on ceramic tiles in four ambient ponds
#' were moved to treatment ponds, three or four per treatment: ambient
#' (`"A"`); elevated salinity (`"S"`, 15 g NaCl/L) or nutrient enrichment
#' (`"N"`, 10 mg N/L and 1 mg P/L), applied once; or one of the two stressors
#' and, a week later, the other (`"SN"`, salinity then nutrients; `"NS"`,
#' nutrients then salinity). The ponds with a single stressor received it
#' with the first application. Tiles were sampled 1, 21 and 81 days after the
#' second application. The article also has ponds with the same stressor
#' applied twice, which do not form an interaction and are not included.
#'
#' @format A data frame with 54 rows (18 ponds x 3 sampling days) and 9
#'   columns:
#'   \describe{
#'     \item{pond}{Pond code (text), from the treatment and the replicate
#'       number, such as `"NS_P2"`: the experimental unit.}
#'     \item{treatment}{Treatment code of the source: `"A"`, `"N"`, `"S"`,
#'       `"NS"` or `"SN"`.}
#'     \item{sequence}{Order of the two stressors: `"nutrients first"`
#'       (`"NS"`), `"salinity first"` (`"SN"`) or `"single"` for the ponds
#'       that both orders share (`"A"`, `"N"` and `"S"`).}
#'     \item{nutrients}{`"added"` or `"none"`.}
#'     \item{salinity}{`"raised"` or `"none"`.}
#'     \item{day}{Days after the second application of stressors (integer:
#'       1, 21, 81).}
#'     \item{gross_production}{Gross primary production of the biofilm of a
#'       tile.}
#'     \item{net_production}{Net primary production.}
#'     \item{respiration}{Community respiration.}
#'   }
#'   The three rates are as in the source file, whose README gives them in
#'   nmol O2 per cm² and hour. They were not measured on day 81 (the pond
#'   water had become anoxic), where they are missing; there are no other
#'   missing values. Gross production is net production plus respiration in
#'   every sample but one (pond `"A_P3"`, day 1), where the source has net
#'   production minus respiration.
#'
#' @source Fuggle, R., Mayer-Pinto, M., Matias, M., & Marzinelli, Z. (2025).
#'   Timing of salinisation and nutrient enrichment drives freshwater
#'   microbial community metabolic responses. Dryad,
#'   \doi{10.5061/dryad.8pk0p2nz8}. Licence CC0 1.0. File
#'   `univariate_data.csv`, the treatments that form an interaction of the
#'   two stressors and the three oxygen rates (script
#'   `data-raw/fuggle_2025.R` of the source repository).
#'
#' @references
#' Fuggle, R., Matias, M. G., Mayer-Pinto, M., & Marzinelli, E. M. (2025).
#' Repetition and order of exposure to salinization and nutrient enrichment
#' drive freshwater microbial community responses. *Ecology*, 106(11),
#' e70260. \doi{10.1002/ecy.70260}
#'
#' @family case-study data
#' @examples
#' str(fuggle_biofilms)
#' with(fuggle_biofilms[fuggle_biofilms$day == 1, ], table(treatment, sequence))
#' # Each order is a `by` group; the ponds of both orders are shared
#' describe_design(fuggle_biofilms, stressors = c("nutrients", "salinity"),
#'                 control = c(nutrients = "none", salinity = "none"),
#'                 unit = "pond", time = "day", by = "sequence",
#'                 shared = c(sequence = "single"))
"fuggle_biofilms"

#' Growth of two diatoms under temperature, photoperiod and microbiome removal
#'
#' Maximum growth rates of a polar and a temperate diatom at three
#' temperatures and three photoperiods, with and without their microbiome
#' (Giesler et al. 2023).
#'
#' @section Experimental design:
#' Laboratory cultures of the polar diatom *Thalassiosira gravida* and the
#' temperate *T. rotula*, each as a xenic culture (with its microbiome) and an
#' axenic one (microbiome removed), grown at each combination of three
#' temperatures (4, 9 and 13.5 °C) and three photoperiods (4, 16 and 24 h of
#' light at 50 µmol of photons per m² and second) after 7 days of
#' acclimation. Each combination had 48 replicate wells of 300 µL in 96-well
#' plates, incubated in climate cabinets at the temperature of the treatment,
#' on LED tables with its photoperiod. The fluorescence of each well was read
#' daily for 7 days, and its maximum growth rate was fitted to those
#' readings. The data do not identify plates, cabinets or tables: the wells
#' of a combination shared its incubation conditions.
#'
#' @format A data frame with 1728 rows (2 species x 2 microbiome treatments x
#'   3 temperatures x 3 photoperiods x 48 wells) and 6 columns:
#'   \describe{
#'     \item{species}{`"T. gravida"` (polar) or `"T. rotula"` (temperate).}
#'     \item{microbiome}{`"present"` (xenic culture) or `"removed"`
#'       (axenic).}
#'     \item{temperature}{Temperature (°C: 4, 9 or 13.5).}
#'     \item{photoperiod}{Hours of light per day (integer: 4, 16 or 24).}
#'     \item{replicate}{Well number within the combination (integer, 1 to
#'       48).}
#'     \item{growth_rate}{Maximum growth rate (per day). Negative values are
#'       wells whose fluorescence declined.}
#'   }
#'   There are no missing values.
#'
#' @source Giesler, J., Wohlrab, S., & Harder, T. (2023). Growth of
#'   *Thalassiosira gravida* and *Thalassiosira rotula* across different
#'   levels of temperature, photoperiod, and microbiome presence. Dryad,
#'   \doi{10.5061/dryad.x95x69pqf}. Licence CC0 1.0. File
#'   `plateassay_growthrates.xlsx`, with the columns renamed (script
#'   `data-raw/giesler_2023.R` of the source repository).
#'
#' @references
#' Giesler, J. K., Harder, T., & Wohlrab, S. (2023). Microbiome and
#' photoperiod interactively determine thermal sensitivity of polar and
#' temperate diatoms. *Biology Letters*, 19(11), 20230151.
#' \doi{10.1098/rsbl.2023.0151}
#'
#' @family case-study data
#' @examples
#' str(giesler_diatoms)
#' with(giesler_diatoms[giesler_diatoms$species == "T. gravida", ],
#'      round(tapply(growth_rate, list(temperature, photoperiod, microbiome), mean), 2))
"giesler_diatoms"

#' Development of horseshoe crab embryos under temperature, salinity and
#' oxygen
#'
#' Embryos of the American horseshoe crab (*Limulus polyphemus*) that
#' developed in each of 36 combinations of temperature, salinity and oxygen,
#' by mating pair (Vasquez et al. 2017).
#'
#' @section Experimental design:
#' Fertilized eggs of six mating pairs from Slaughter Beach (Delaware Bay,
#' USA), May 2012, one egg per well of 96-well plates, were exposed for 14
#' days to each of the 36 combinations of temperature (20, 25, 30 and 35 °C),
#' salinity (5, 15 and 34 ppt) and oxygen (5, 13 and 21% O2 in the air), one
#' plate per combination, and then kept for 14 days in the reference
#' condition (30 °C, 34 ppt and 21% O2). Every 2 days the eggs were imaged,
#' and the cross-sectional area of eight eggs per pair and combination was
#' measured. An embryo developed when it reached 5 mm², the size of a
#' late-stage embryo. The plate of each combination is confounded with it.
#'
#' @format A data frame with 216 rows (36 combinations x 6 mating pairs) and
#'   6 columns:
#'   \describe{
#'     \item{temperature}{Temperature (°C, integer: 20, 25, 30, 35).}
#'     \item{salinity}{Salinity (ppt, integer: 5, 15, 34).}
#'     \item{oxygen}{Oxygen in the air (% O2, integer: 5, 13, 21).}
#'     \item{pair}{Mating pair (integer, 1 to 6), whose eggs were spread
#'       over all combinations.}
#'     \item{embryos}{Embryos measured (integer, 8, or 7 in three cells).}
#'     \item{developed}{Embryos that reached 5 mm² by day 28 (integer).}
#'   }
#'   An embryo counts as developed when the event indicator of the source,
#'   used by the survival analysis of the article, became 1. Three embryos
#'   are left out: two never measured and one without that indicator. There
#'   are no missing values.
#'
#' @source Vasquez, M. C., Brockmann, H. J., & Julian, D. (2018). Data from:
#'   Between-population differences in multi-stressor tolerance during
#'   embryo development in the American horseshoe crab, *Limulus
#'   polyphemus*. Dryad, \doi{10.5061/dryad.vs1r0}. Licence CC0 1.0. Sheet
#'   "Embryo Development Size" of
#'   `Dryad+-+Embryo+Cross-Sectional+Area.xlsx`, one embryo per row and day,
#'   summed by pair and combination (script `data-raw/vasquez_2017.R` of the
#'   source repository).
#'
#' @references
#' Vasquez, M. C., Brockmann, H. J., & Julian, D. (2017). Between-population
#' differences in multi-stressor tolerance during embryo development in the
#' American horseshoe crab, *Limulus polyphemus*. *Estuaries and Coasts*,
#' 40(5), 1487–1501. \doi{10.1007/s12237-017-0218-1}
#'
#' @family case-study data
#' @examples
#' str(vasquez_embryos)
#' # Share of embryos that developed in each combination (all pairs)
#' n <- with(vasquez_embryos, tapply(developed, list(temperature, salinity, oxygen), sum))
#' d <- with(vasquez_embryos, tapply(embryos, list(temperature, salinity, oxygen), sum))
#' round(n / d, 2)
"vasquez_embryos"

#' Stream invertebrates under CO2 enrichment, heatwaves, flow variability
#' and fine sediment
#'
#' Benthic invertebrate counts at the end of a stream mesocosm experiment
#' with a full factorial of four stressors (Hunn et al. 2024).
#'
#' @section Experimental design:
#' One hundred and twenty-eight outdoor stream mesocosms (ExStream system,
#' New Zealand), eight per combination of dissolved CO2 (ambient or
#' enriched), temperature (ambient or two heatwaves of 7 days), flow velocity
#' (constant or variable) and fine sediment (none or 300 g of dry sediment),
#' in two spatial blocks. Eight header tanks, four per block, fed 16
#' mesocosms each: temperature and CO2 were applied to the header tanks, and
#' flow and sediment to the mesocosms. After 7 weeks, on day 35, the whole
#' benthic community of each mesocosm was sampled, and one quarter of each
#' sample was counted.
#'
#' @format A data frame with 128 rows (one per mesocosm) and 22 columns:
#'   \describe{
#'     \item{channel}{Mesocosm number (integer, 1 to 128).}
#'     \item{block}{Spatial block (integer, 1 or 2).}
#'     \item{header_tank}{Header tank (integer, 1 to 8), which is not in the
#'       source file: deduced from the channel numbers, since channels 1 to
#'       16, 17 to 32, and so on, share their block, temperature and CO2,
#'       as the 16 mesocosms of a header tank do.}
#'     \item{temperature}{`"ambient"` or `"heatwaves"`.}
#'     \item{co2}{`"ambient"` or `"enriched"`.}
#'     \item{sediment}{`"none"` or `"added"`.}
#'     \item{flow}{`"constant"` or `"variable"`.}
#'     \item{total_abundance}{Invertebrates counted in the quarter of the
#'       sample, all taxa.}
#'     `r rd_taxa(hunn_invertebrates)`
#'   }
#'   The taxon columns count the 14 taxa of the source file in the counted
#'   quarter; the attribute `"taxa"` gives their source names. Annelids,
#'   which fragment, are counted in halves, and so is the total. Zeros are
#'   taxa not found; there are no missing values.
#'
#' @source Hunn, J. G. (2024). Data for Individual and combined impacts of
#'   carbon dioxide enrichment, heatwaves, flow velocity variability and fine
#'   sediment deposition on stream invertebrate communities (version v1).
#'   Zenodo, \doi{10.5281/zenodo.11064543}. Licence CC BY 4.0. File
#'   `Benthic data summary.txt`: the counts of the quarter (columns
#'   `25pc_...`), with the header tanks added and the columns renamed
#'   (script `data-raw/hunn_2024.R` of the source repository).
#'
#' @references
#' Hunn, J. G., Orr, J. A., Kelly, A.-M., Piggott, J. J., & Matthaei, C. D.
#' (2024). Individual and combined impacts of carbon dioxide enrichment,
#' heatwaves, flow velocity variability, and fine sediment deposition on
#' stream invertebrate communities. *Global Change Biology*, 30(5), e17336.
#' \doi{10.1111/gcb.17336}
#'
#' @family case-study data
#' @examples
#' str(hunn_invertebrates[1:9])
#' # Temperature and CO2 are constant within each header tank
#' unique(hunn_invertebrates[c("header_tank", "block", "temperature", "co2")])
"hunn_invertebrates"

#' Potential denitrification under five global changes
#'
#' Potential denitrification and soil water in a grassland under a full
#' factorial of elevated CO2, warming, increased precipitation, nitrogen
#' addition and fire, sampled five times after the fire (Niboyet et al.
#' 2025).
#'
#' @section Experimental design:
#' Jasper Ridge Global Change Experiment (California, USA), an annual
#' grassland under CO2 enrichment (+275 µmol/mol), warming (+2.0 °C),
#' increased precipitation (+50% and a rainy season 3 weeks longer) and
#' nitrogen addition (+7 g N per m² and year) since 1998, and a controlled
#' fire on half of the plots in July 2011. Randomized block split-plot
#' design: CO2, temperature and fire were applied to circular plots of 2 m
#' diameter, and precipitation and nitrogen to the four quadrants of each
#' plot. With 4 blocks of 8 plots, each of the 32 combinations has four
#' quadrants. Soil cores of every quadrant were taken 2 days after the fire
#' and 5, 9, 21 and 33 months later.
#'
#' @format A data frame with 640 rows (128 quadrants x 5 samplings) and 12
#'   columns:
#'   \describe{
#'     \item{block}{Block (integer, 1 to 4).}
#'     \item{plot}{Plot (integer, 1 to 32), to which CO2, temperature and fire
#'       were applied.}
#'     \item{quadrant}{Quadrant within the plot (integer, 1 to 4).}
#'     \item{unit}{Quadrant identifier (integer, 1 to 128), the experimental
#'       unit of precipitation and nitrogen.}
#'     \item{co2}{`"ambient"` or `"elevated"`.}
#'     \item{temperature}{`"ambient"` or `"elevated"`.}
#'     \item{water}{Precipitation: `"ambient"` or `"elevated"`.}
#'     \item{nitrogen}{`"ambient"` or `"elevated"`.}
#'     \item{fire}{`"unburned"` or `"burned"`.}
#'     \item{months}{Months after the fire (integer: 0, 5, 9, 21, 33).}
#'     \item{denitrification}{Potential denitrification, measured as
#'       denitrifying enzyme activity (ng N2O-N per g of dry soil and hour).}
#'     \item{soil_water}{Soil water content (%).}
#'   }
#'   Values and units as in the source file. One denitrification value and
#'   six soil water values are missing in the source.
#'
#' @source Niboyet, A. (2025). Soil moisture and the number of simultaneous
#'   stressors drive interactions among global changes on denitrification
#'   (version v1). Zenodo, \doi{10.5281/zenodo.16789113}.
#'   Licence CC BY 4.0. File
#'   `Niboyetetal_CommunicationsEarthEnvironment_10august2025.xlsx`, with the
#'   treatment codes as labels and two of its variables (script
#'   `data-raw/niboyet_2025.R` of the source repository).
#'
#' @references
#' Niboyet, A., Le Roux, X., Chiariello, N. R., Dijkstra, P., Field, C. B.,
#' & Hungate, B. A. (2025). Soil moisture and the number of simultaneous
#' stressors drive interactions among global changes on denitrification.
#' *Communications Earth & Environment*, 6, 704.
#' \doi{10.1038/s43247-025-02703-5}
#'
#' @family case-study data
#' @examples
#' str(niboyet_denitrification)
#' # CO2, temperature and fire are constant within each plot
#' head(unique(niboyet_denitrification[c("plot", "block", "co2", "temperature", "fire")]))
#' with(niboyet_denitrification, round(tapply(soil_water, months, mean, na.rm = TRUE), 1))
"niboyet_denitrification"
