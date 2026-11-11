# Water quality in an insecticide and temperature mesocosm experiment

Physico-chemical variables measured in the morning and in the evening in
the mesocosm experiment of Hermann et al. (2024) with an insecticide and
three temperature regimes.

## Usage

``` r
hermann_water_quality
```

## Format

A data frame with 480 rows (24 mesocosms x 10 days x 2 times of day) and
14 columns:

- mesocosm:

  Mesocosm number (text, `"M01"` to `"M24"`), the experimental unit; the
  same coding in
  [hermann_nutrients](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
  [hermann_leaf_litter](https://pabloamadorc.github.io/must/reference/hermann_leaf_litter.md)
  and
  [hermann_imidacloprid](https://pabloamadorc.github.io/must/reference/hermann_imidacloprid.md).

- day:

  Day relative to the first imidacloprid application (integer: -4, 3,
  10, 17, 24, 31, 38, 52, 66, 80).

- imidacloprid:

  `"applied"` or `"absent"`.

- temperature:

  `"ambient"`, `"elevated"` or `"heatwaves"`.

- time_of_day:

  `"morning"` or `"evening"`.

- water_temperature:

  Water temperature (°C).

- oxygen:

  Dissolved oxygen (mg/L).

- oxygen_saturation:

  Oxygen saturation (%).

- conductivity:

  Electrical conductivity (µS/cm).

- ph:

  pH.

- tds:

  Total dissolved solids (ppm).

- chlorophyll_a, phycoerythrin:

  Chlorophyll-a and phycoerythrin, as reported in the source file, which
  does not state their unit.

- ecosystem_productivity:

  Evening minus morning dissolved oxygen of the same day (mg/L), as
  given in the source; missing in the morning rows.

There are no other missing values.

## Source

Hermann et al. (2024), Zenodo,
[doi:10.5281/zenodo.11520459](https://doi.org/10.5281/zenodo.11520459)
(full reference in
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)).
Licence CC BY 4.0. Sheets "WaterQualityMorning" and
"WaterQualityEvening"; script `data-raw/hermann_2024.R` of the source
repository.

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
hermann_water_quality,
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
[`hermann_zooplankton`](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)

## Examples

``` r
# Mean evening water temperature by temperature regime and day
with(hermann_water_quality[hermann_water_quality$time_of_day == "evening", ],
     round(tapply(water_temperature, list(temperature, day), mean), 1))
#>           -4    3   10   17   24   31   38   52   66   80
#> ambient   17 17.0 19.9 17.4 23.1 22.4 24.1 22.6 26.0 25.2
#> elevated  17 21.1 24.5 21.6 27.6 26.6 28.2 22.7 26.0 25.1
#> heatwaves 17 24.2 20.5 27.9 23.8 32.0 24.8 22.7 25.9 25.0
```
