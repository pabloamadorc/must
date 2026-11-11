# Macroinvertebrates in an insecticide and temperature mesocosm experiment

Macroinvertebrate counts from the mesocosm experiment of Hermann et al.
(2024) with an insecticide and three temperature regimes, sampled once
before the first insecticide application and six times after it.

## Usage

``` r
hermann_macroinvertebrates
```

## Format

A data frame with 168 rows (24 mesocosms x 7 sampling days) and 34
columns:

- mesocosm, imidacloprid, temperature:

  As in
  [hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md).

- day:

  Day relative to the first imidacloprid application (integer: -7, 7,
  21, 35, 49, 63, 77).

- total_abundance:

  Sum of the 29 taxon columns.

- erpobdella:

  Erpobdella sp.

- lumbriculidae:

  Lumbriculidae

- hydracarina:

  Hydracarina

- sphaeriidae:

  Sphaeriidae

- elmidae_adult:

  Elmidae adult

- elmidae_larvae:

  Elmidae larvae

- echinogammarus:

  Echinogammarus

- chironomini:

  Chironomini

- culicidae:

  Culicidae

- diptera_pupa:

  Diptera pupa

- orthocladiinae:

  Orthocladiinae

- caenidae:

  Caenidae

- cloeon:

  Cloeon sp.

- heptageniidae:

  Heptageniidae

- ancylidae:

  Ancylidae

- hydrobiidae:

  Hydrobiidae

- physella_acuta:

  Physella acuta

- planorbidae:

  Planorbidae

- notonectidae:

  Notonectidae

- anisoptera_i_green:

  Anisoptera I (green)

- zygoptera:

  Zygoptera

- dugesia:

  Dugesia sp.

- polycelis:

  Polycelis sp.

- leptoceridae_stone_house:

  Leptoceridae (stone house)

- leptoceridae_wood_house:

  Leptoceridae (wood house)

- lymnaeidae:

  Lymnaeidae

- tanypodinae:

  Tanypodinae

- chaoborus:

  Chaoborus sp.

- beetle_big:

  Beetle (big)

The taxon columns hold the number of individuals of each taxon of the
source file, whose name for the taxon is given above and in the
attribute `"taxa"`; the article describes four bottom traps per
mesocosm. Counts are as reported in the source file. Zeros are taxa not
found in the sample; there are no missing values.

## Source

Hermann et al. (2024), Zenodo,
[doi:10.5281/zenodo.11520459](https://doi.org/10.5281/zenodo.11520459)
(full reference in
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)).
Licence CC BY 4.0. Sheets "MacroinvertebratesDesign" and
"Macroinvertebrates_Abundances"; script `data-raw/hermann_2024.R` of the
source repository.

## Experimental design

Twenty-four outdoor mesocosms (about 850 L of water over a sediment
layer) at Alcalá de Henares (central Spain), from 19 April (day -7) to
19 July 2021, four per combination of a 2 x 3 factorial: imidacloprid
(nominal 1 µg/L, applied on days 0, 14 and 28, or absent) and
temperature (ambient; elevated, 4 °C above ambient; or three heatwaves
of 7 days, separated by 7 days at ambient temperature, with targets of
24, 28 and 32 °C). The temperature treatments lasted from day 0 to day
42; days 42 to 84 were a recovery phase.

The community tables
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md),
[hermann_phytoplankton](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md)
and hermann_macroinvertebrates code the mesocosms by treatment and
replicate number (`"C0A1"`);
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
[`hermann_nutrients`](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
[`hermann_phytoplankton`](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md),
[`hermann_water_quality`](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[`hermann_zooplankton`](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)

## Examples

``` r
range(hermann_macroinvertebrates$total_abundance)
#> [1]  48 635
# Samples in which each taxon was found
taxa <- names(attr(hermann_macroinvertebrates, "taxa"))
sort(colSums(hermann_macroinvertebrates[taxa] > 0), decreasing = TRUE)
#>                   cloeon                zygoptera       anisoptera_i_green 
#>                      167                      163                      160 
#>           physella_acuta             notonectidae            lumbriculidae 
#>                      150                      150                      144 
#>                 caenidae              chironomini              tanypodinae 
#>                      124                      120                      111 
#>               lymnaeidae               erpobdella             diptera_pupa 
#>                      107                       99                       93 
#>           orthocladiinae                culicidae           echinogammarus 
#>                       89                       60                       58 
#>              hydrobiidae              hydracarina                chaoborus 
#>                       53                       38                       30 
#>              sphaeriidae  leptoceridae_wood_house leptoceridae_stone_house 
#>                       29                       27                       24 
#>            elmidae_adult              planorbidae                  dugesia 
#>                       22                       17                       16 
#>            heptageniidae                polycelis           elmidae_larvae 
#>                       12                        7                        4 
#>                ancylidae               beetle_big 
#>                        2                        2 
```
