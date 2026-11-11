# Zooplankton in a three-stressor mesocosm experiment

Zooplankton densities from an outdoor mesocosm experiment with a full
factorial of an insecticide, a herbicide and nutrient enrichment
(Polazzo & Rico 2021), sampled once before the pesticide application and
four times after it.

## Usage

``` r
polazzo_zooplankton
```

## Format

A data frame with 120 rows (24 mesocosms x 5 sampling days) and 28
columns:

- mesocosm:

  Mesocosm code (text, `"M01"` to `"M24"`), the experimental unit.

- day:

  Day relative to the pesticide application (integer: -5, 7, 14, 30,
  50). Day -5 is before the pesticides but not before the nutrients.

- insecticide:

  `"chlorpyrifos"` or `"absent"`.

- herbicide:

  `"diuron"` or `"absent"`.

- nutrients:

  `"added"` or `"ambient"`.

- total_density:

  Sum of the 22 taxon columns.

- daphnia:

  Daphnia sp

- simocephalus:

  Simocephalus sp

- cyclopoida:

  Ciclopoida

- alonella:

  Alonella

- alona:

  Alona

- ceriodafnia:

  Ceriodafnia

- ostracoda:

  Ostracoda

- nauplii:

  Nauplii

- keratella_quadrata:

  Keratella quadrata

- lecanidae:

  Lecanidae

- cephalodella:

  Cephalodella

- lepadella:

  Lepadella

- tricocerca:

  Tricocerca

- mytilinia:

  Mytilinia

- polyarthra:

  Polyarthra

- colurella:

  Colurella

- tricotria:

  Tricotria

- ascomorfa:

  Ascomorfa

- squatinella:

  Squatinella

- testudinella:

  Testudinella

- notholca:

  Notholca

- proales:

  Proales

The taxon columns hold the density of each taxon of the source file,
whose name for the taxon is given above and in the attribute `"taxa"`.
Densities are as reported in the source file, which does not state their
unit. Zeros are taxa not found in the sample; there are no missing
values.

## Source

Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple stressors
on the dimensionality of ecological stability. Dryad,
[doi:10.5061/dryad.qfttdz0gt](https://doi.org/10.5061/dryad.qfttdz0gt) .
Licence CC0 1.0. Sheet "Zooplankton density" of
`Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx`; the treatment
codes, which differ on day 50, were harmonised (script
`data-raw/polazzo_rico_2021.R` of the source repository).

## Experimental design

Twenty-four outdoor mesocosms, three per combination of a full factorial
of chlorpyrifos (insecticide, 1 µg/L, or absent), diuron (herbicide, 18
µg/L, or absent) and nutrients (750 µg/L N and 75 µg/L P added twice a
week, or ambient). The nutrient additions started four weeks before the
pesticide application (day 0). The community tables polazzo_zooplankton,
[polazzo_phytoplankton](https://pabloamadorc.github.io/must/reference/polazzo_phytoplankton.md)
and
[polazzo_macroinvertebrates](https://pabloamadorc.github.io/must/reference/polazzo_macroinvertebrates.md)
share the mesocosm codes and the design columns, so they can be joined
by `mesocosm` and, where the sampling days coincide, `day`.

## References

Polazzo, F., & Rico, A. (2021). Effects of multiple stressors on the
dimensionality of ecological stability. *Ecology Letters*, 24(8),
1594–1606. [doi:10.1111/ele.13770](https://doi.org/10.1111/ele.13770)

## See also

Other data of Polazzo and Rico (2021):
[`polazzo_macroinvertebrates`](https://pabloamadorc.github.io/must/reference/polazzo_macroinvertebrates.md),
[`polazzo_phytoplankton`](https://pabloamadorc.github.io/must/reference/polazzo_phytoplankton.md)

## Examples

``` r
str(polazzo_zooplankton[1:8])
#> 'data.frame':    120 obs. of  8 variables:
#>  $ mesocosm     : chr  "M01" "M02" "M03" "M04" ...
#>  $ day          : int  -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 ...
#>  $ insecticide  : chr  "chlorpyrifos" "absent" "absent" "chlorpyrifos" ...
#>  $ herbicide    : chr  "absent" "absent" "diuron" "diuron" ...
#>  $ nutrients    : chr  "added" "added" "ambient" "added" ...
#>  $ total_density: num  805 485 510 739 999 ...
#>  $ daphnia      : num  44.2 33.4 2.2 11.8 3.2 15.2 0 2 0.6 13 ...
#>  $ simocephalus : num  1.8 0.6 1.2 0.8 1.6 0.4 0 0.2 0.4 0.2 ...
# Three mesocosms per treatment combination
with(polazzo_zooplankton[polazzo_zooplankton$day == 7, ],
     table(insecticide, herbicide, nutrients))
#> , , nutrients = added
#> 
#>               herbicide
#> insecticide    absent diuron
#>   absent            3      3
#>   chlorpyrifos      3      3
#> 
#> , , nutrients = ambient
#> 
#>               herbicide
#> insecticide    absent diuron
#>   absent            3      3
#>   chlorpyrifos      3      3
#> 
head(attr(polazzo_zooplankton, "taxa"))
#>           daphnia      simocephalus        cyclopoida          alonella 
#>      "Daphnia sp" "Simocephalus sp"      "Ciclopoida"        "Alonella" 
#>             alona       ceriodafnia 
#>           "Alona"     "Ceriodafnia" 
```
