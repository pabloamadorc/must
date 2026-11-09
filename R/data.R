#' Zooplankton in a three-stressor mesocosm experiment
#'
#' Zooplankton densities from an outdoor mesocosm experiment with a full
#' factorial of an insecticide (chlorpyrifos), a herbicide (diuron) and
#' nutrient enrichment, three mesocosms per treatment combination, sampled
#' once before the pesticide application and four times after it (Polazzo &
#' Rico 2021).
#'
#' @format A data frame with 120 rows (24 mesocosms x 5 sampling days) and 8
#'   columns:
#'   \describe{
#'     \item{mesocosm}{Mesocosm identifier (text, `"M01"` to `"M24"`), the
#'       experimental unit.}
#'     \item{day}{Day relative to the pesticide application (integer: -5, 7,
#'       14, 30, 50). Nutrient additions started four weeks before the
#'       application, so day -5 is before the pesticides but not before the
#'       nutrients.}
#'     \item{insecticide}{`"chlorpyrifos"` (1 ug/L) or `"absent"`.}
#'     \item{herbicide}{`"diuron"` (18 ug/L) or `"absent"`.}
#'     \item{nutrients}{`"added"` (750 ug/L N and 75 ug/L P, twice a week) or
#'       `"ambient"`.}
#'     \item{total_density}{Sum of the densities of the 22 zooplankton taxa of
#'       the source file.}
#'     \item{daphnia}{Density of *Daphnia* sp.}
#'     \item{cyclopoida}{Density of cyclopoid copepods.}
#'   }
#'   Densities are as reported in the source file, which does not state
#'   their unit. Zeros are taxa not found in the sample; there are no
#'   missing values.
#'
#' @source Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple
#'   stressors on the dimensionality of ecological stability. Dryad,
#'   \doi{10.5061/dryad.qfttdz0gt}. Licence CC0 1.0. Sheet "Zooplankton
#'   density" of `Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx`;
#'   the treatment codes, which differ on day 50, were harmonised (script in
#'   `data-raw/` of the source repository).
#'
#' @references
#' Polazzo, F., & Rico, A. (2021). Effects of multiple stressors on the
#' dimensionality of ecological stability. *Ecology Letters*, 24(8),
#' 1594–1606. \doi{10.1111/ele.13770}
#'
#' @examples
#' str(polazzo_zooplankton)
#' table(polazzo_zooplankton$insecticide, polazzo_zooplankton$herbicide,
#'       polazzo_zooplankton$nutrients) / 5
"polazzo_zooplankton"
