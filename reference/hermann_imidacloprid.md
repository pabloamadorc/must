# Imidacloprid in the water of an insecticide and temperature mesocosm experiment

Measured imidacloprid concentrations in the treated mesocosms of the
experiment of Hermann et al. (2024) with an insecticide and three
temperature regimes.

## Usage

``` r
hermann_imidacloprid
```

## Format

A data frame with 120 rows (12 treated mesocosms x 10 days) and 4
columns:

- mesocosm:

  Mesocosm number (text), as in
  [hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md);
  only those with imidacloprid.

- day:

  Day relative to the first imidacloprid application (integer: 0, 3, 10,
  14, 17, 24, 28, 31, 38, 45). The article describes the samples of day
  0 as taken 4 hours after the application.

- temperature:

  `"ambient"`, `"elevated"` or `"heatwaves"`.

- concentration:

  Imidacloprid in water (µg/L).

The concentration is missing in mesocosm `"M20"` on days 28, 31 and 38;
there are no other missing values.

## Source

Hermann et al. (2024), Zenodo,
[doi:10.5281/zenodo.11520459](https://doi.org/10.5281/zenodo.11520459)
(full reference in
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)).
Licence CC BY 4.0. Sheet "ImiConcentrationsWater"; script
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

The community tables
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md),
[hermann_phytoplankton](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md)
and
[hermann_macroinvertebrates](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md)
code the mesocosms by treatment and replicate number (`"C0A1"`);
[hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[hermann_nutrients](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
[hermann_leaf_litter](https://pabloamadorc.github.io/must/reference/hermann_leaf_litter.md)
and hermann_imidacloprid by mesocosm number (`"M01"`). The source gives
no correspondence between the two codings, so the two kinds of table
cannot be joined by mesocosm.

## References

Hermann, M., Polazzo, F., Cherta, L., Crettaz-Minaglia, M.,
García-Astillero, A., Peeters, E. T. H. M., Rico, A., & Van den Brink,
P. J. (2024). Combined stress of an insecticide and heatwaves or
elevated temperature induce community and food web effects in a
Mediterranean freshwater ecosystem. *Water Research*, 260, 121903.
[doi:10.1016/j.watres.2024.121903](https://doi.org/10.1016/j.watres.2024.121903)

## See also

Other data of Hermann et al. (2024):
[`hermann_leaf_litter`](https://pabloamadorc.github.io/must/reference/hermann_leaf_litter.md),
[`hermann_macroinvertebrates`](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md),
[`hermann_nutrients`](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
[`hermann_phytoplankton`](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md),
[`hermann_water_quality`](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[`hermann_zooplankton`](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)

## Examples

``` r
with(hermann_imidacloprid, round(tapply(concentration, list(temperature, day),
                                        mean, na.rm = TRUE), 2))
#>              0    3   10   14   17   24   28   31   38   45
#> ambient   1.14 0.79 0.23 1.35 0.88 0.20 1.19 0.62 0.15 0.11
#> elevated  1.17 0.74 0.18 1.44 0.87 0.19 1.33 0.58 0.10 0.03
#> heatwaves 1.23 0.82 0.22 1.43 0.81 0.17 1.33 0.54 0.12 0.03
```
