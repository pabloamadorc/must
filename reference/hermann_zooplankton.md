# Zooplankton in an insecticide and temperature mesocosm experiment

Zooplankton densities from an outdoor mesocosm experiment with an
insecticide (imidacloprid) and three temperature regimes: ambient,
elevated and repeated heatwaves (Hermann et al. 2024), sampled once
before the first insecticide application and six times after it.

## Usage

``` r
hermann_zooplankton
```

## Format

A data frame with 168 rows (24 mesocosms x 7 sampling days) and 26
columns:

- mesocosm:

  Mesocosm code (text) from the sample labels of the source:
  imidacloprid (`C0` absent, `C1` applied), temperature (`A` ambient,
  `E` elevated, `HW` heatwaves) and replicate number, for example
  `"C1HW3"`. The codes appear to identify the same mesocosm on every
  sampling day: within a treatment, the abundances of a code are
  correlated between consecutive samplings, and those of different codes
  are not.

- day:

  Day relative to the first imidacloprid application (integer: -1, 10,
  24, 38, 52, 66, 80), from the sample labels; the article gives day -4
  for the first sampling.

- imidacloprid:

  `"applied"` or `"absent"`.

- temperature:

  `"ambient"`, `"elevated"` or `"heatwaves"`.

- total_density:

  Sum of the 21 taxon columns.

- daphnia_pulex:

  Daphnia pulex

- simocephalus:

  Simocephalus sp

- chydorus:

  Chydorus sp

- cyclopoida:

  Ciclopoida

- calanoida:

  Calanoida

- ostracoda:

  Ostracoda

- daphnia_longispina:

  Daphnia longispina

- alona:

  Alona

- alonella:

  Alonella

- ceriodaphnia:

  Ceriodaphnia

- nauplii:

  Nauplii

- keratella_quadrata:

  Keratella quadrata

- lecanidae:

  Lecanidae

- lepadella:

  Lepadella

- testudinella:

  Testudinella

- polyarthra:

  Polyarthra

- brachionus:

  Brachionus

- squatinella:

  Squatinella

- ascomorpha:

  Ascomorpha

- trichocerca:

  Trichocerca

- cephalodella:

  Cephalodella

The taxon columns hold the density of each taxon of the source file
(individuals per litre), whose name for the taxon is given above and in
the attribute `"taxa"`. Zeros are taxa not found in the sample; there
are no missing values.

## Source

Hermann, M., Polazzo, F., Cherta Cucala, L., Crettaz Minaglia, M. C.,
García-Astillero, A., Peeters, E. T. H. M., Rico, A., & Van den Brink,
P. J. (2024). Data supporting: Combined stress of an insecticide and
heatwaves or elevated temperature induce community and food web effects
in a Mediterranean freshwater ecosystem (version V.1). Zenodo,
[doi:10.5281/zenodo.11520459](https://doi.org/10.5281/zenodo.11520459) .
Licence CC BY 4.0. Sheets "ZooplanktonDesign" and
"Zooplankton_Abundances" of
`CombinedStressInsecticideTemperature_Hermann_et_al.(2024).xlsx`; script
`data-raw/hermann_2024.R` of the source repository.

## Experimental design

Twenty-four outdoor mesocosms (about 850 L of water over a sediment
layer) at Alcalá de Henares (central Spain), from 19 April (day -7) to
19 July 2021, four per combination of a 2 x 3 factorial: imidacloprid
(nominal 1 µg/L, applied on days 0, 14 and 28, or absent) and
temperature (ambient; elevated, 4 °C above ambient; or three heatwaves
of 7 days, separated by 7 days at ambient temperature, with targets of
24, 28 and 32 °C). The temperature treatments lasted from day 0 to day
42; days 42 to 84 were a recovery phase.

The community tables hermann_zooplankton,
[hermann_phytoplankton](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md)
and
[hermann_macroinvertebrates](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md)
code the mesocosms by treatment and replicate number (`"C0A1"`);
[hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[hermann_nutrients](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
[hermann_leaf_litter](https://pabloamadorc.github.io/must/reference/hermann_leaf_litter.md)
and
[hermann_imidacloprid](https://pabloamadorc.github.io/must/reference/hermann_imidacloprid.md)
by mesocosm number (`"M01"`). The source gives no correspondence between
the two codings, so the two kinds of table cannot be joined by mesocosm.

## References

Hermann, M., Polazzo, F., Cherta, L., Crettaz-Minaglia, M.,
García-Astillero, A., Peeters, E. T. H. M., Rico, A., & Van den Brink,
P. J. (2024). Combined stress of an insecticide and heatwaves or
elevated temperature induce community and food web effects in a
Mediterranean freshwater ecosystem. *Water Research*, 260, 121903.
[doi:10.1016/j.watres.2024.121903](https://doi.org/10.1016/j.watres.2024.121903)

## See also

Other data of Hermann et al. (2024):
[`hermann_imidacloprid`](https://pabloamadorc.github.io/must/reference/hermann_imidacloprid.md),
[`hermann_leaf_litter`](https://pabloamadorc.github.io/must/reference/hermann_leaf_litter.md),
[`hermann_macroinvertebrates`](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md),
[`hermann_nutrients`](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
[`hermann_phytoplankton`](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md),
[`hermann_water_quality`](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md)

## Examples

``` r
str(hermann_zooplankton[1:9])
#> 'data.frame':    168 obs. of  9 variables:
#>  $ mesocosm     : chr  "C0A1" "C0A2" "C0A3" "C0A4" ...
#>  $ day          : int  -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#>  $ imidacloprid : chr  "absent" "absent" "absent" "absent" ...
#>  $ temperature  : chr  "ambient" "ambient" "ambient" "ambient" ...
#>  $ total_density: num  397 617 498 571 316 ...
#>  $ daphnia_pulex: num  2.2 3.4 21.4 10 3.2 8.6 2.4 9.8 3.6 7.2 ...
#>  $ simocephalus : num  0.2 0.4 0.8 1.8 0.6 0.4 0.4 4 0.4 0.6 ...
#>  $ chydorus     : num  1.8 5.2 1 2.8 0.8 1.6 1.6 0.8 0 1 ...
#>  $ cyclopoida   : num  13.2 80.6 9.4 12.8 72.2 41.4 32.8 16.6 52.8 93.8 ...
# Four mesocosms per treatment combination on each sampling day
with(hermann_zooplankton, table(imidacloprid, temperature, day))[, , "24"]
#>             temperature
#> imidacloprid ambient elevated heatwaves
#>      absent        4        4         4
#>      applied       4        4         4
# The design as must sees it
describe_design(hermann_zooplankton[hermann_zooplankton$day > 0, ],
                stressors = c("imidacloprid", "temperature"),
                control = c(imidacloprid = "absent", temperature = "ambient"),
                unit = "mesocosm", time = "day")
#> <must_design> 144 rows, 2 stressors, 6 observed treatment combinations
#> Stressors (control): imidacloprid (absent), temperature (ambient)
#> Unit: mesocosm | Time: day | By: none | Groups: none
#> Contrast sets: 2 candidate(s), 2 identifiable
#> Checks:
#>   info     Experimental unit: 'mesocosm' (24 units, each in a single
#>            treatment combination).
#>   info     Units per treatment combination: 4.
#>   info     Units are measured repeatedly (6 sampling times; 6 per unit).
#>            Samplings of the same unit are not independent: by default
#>            fit_responses() adds a random intercept per unit, which assumes
#>            the same correlation between any two samplings of a unit;
#>            temporal autocorrelation can be modelled through its `random`
#>            argument.
#>   info     Identifiable contrast sets (net deviation): 2 with 2 stressors.
#>   info     2 identifiable set(s) share cells (one cell is used by up to 2
#>            sets, typically the control): their estimates are correlated,
#>            not independent replicates of an interaction.
```
