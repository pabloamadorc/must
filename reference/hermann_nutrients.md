# Nutrients in an insecticide and temperature mesocosm experiment

Nitrate, orthophosphate and total phosphorus in the mesocosm experiment
of Hermann et al. (2024) with an insecticide and three temperature
regimes.

## Usage

``` r
hermann_nutrients
```

## Format

A data frame with 240 rows (24 mesocosms x 10 days) and 7 columns:

- mesocosm, imidacloprid, temperature:

  As in
  [hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md).

- day:

  Day relative to the first imidacloprid application (integer: -4, 3,
  10, 17, 24, 31, 38, 52, 66, 80).

- nitrate, orthophosphate, total_phosphorus:

  Concentrations as reported in the source file, which does not state
  their units.

Nitrate is missing on day 52 in mesocosms `"M07"` and `"M09"`; there are
no other missing values.

## Source

Hermann et al. (2024), Zenodo,
[doi:10.5281/zenodo.11520459](https://doi.org/10.5281/zenodo.11520459)
(full reference in
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)).
Licence CC BY 4.0. Sheet "Nutrients"; script `data-raw/hermann_2024.R`
of the source repository.

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
and
[hermann_macroinvertebrates](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md)
code the mesocosms by treatment and replicate number (`"C0A1"`);
[hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
hermann_nutrients,
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
[`hermann_phytoplankton`](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md),
[`hermann_water_quality`](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[`hermann_zooplankton`](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)

## Examples

``` r
with(hermann_nutrients, round(tapply(total_phosphorus, list(temperature, day), median)))
#>           -4  3 10 17 24 31 38 52 66 80
#> ambient   21 30 30 26 29 37 35 44 43 27
#> elevated  31 31 34 34 27 36 37 38 56 30
#> heatwaves 27 32 32 27 37 38 40 41 38 38
```
