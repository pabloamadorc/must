# Leaf litter decomposition in an insecticide and temperature mesocosm experiment

Dry mass lost by poplar leaves in litter bags, by microbial
decomposition and by macroinvertebrate feeding, in the mesocosm
experiment of Hermann et al. (2024) with an insecticide and three
temperature regimes.

## Usage

``` r
hermann_leaf_litter
```

## Format

A data frame with 264 rows (24 mesocosms x 6 days of microbial
decomposition and 5 of macroinvertebrate feeding) and 8 columns:

- mesocosm, imidacloprid, temperature:

  As in
  [hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md).

- day:

  Day of collection relative to the first imidacloprid application
  (integer: 11, 25, 39, 53, 67 and 81 for microbial decomposition, and
  the same days except 67 for macroinvertebrate feeding). The article
  reports collections after 14 to 84 days of exposure.

- process:

  `"microbial"` (fine-mesh bags, 500 µm) or
  `"macroinvertebrate feeding"` (coarse-mesh bags, up to 0.5 cm).

- leaf_dry_mass:

  Dry mass of the leaves left in the bag (g); each bag held 2 g.

- mass_loss:

  Dry mass lost (g).

- mass_loss_percent:

  Dry mass lost as a percentage of the initial 2 g.

Values are as reported in the source file. Two samples of mesocosm
`"M18"` have no percentage in the source and are missing: microbial
decomposition on day 11 (no loss) and macroinvertebrate feeding on day
25 (a gain).

## Source

Hermann et al. (2024), Zenodo,
[doi:10.5281/zenodo.11520459](https://doi.org/10.5281/zenodo.11520459)
(full reference in
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)).
Licence CC BY 4.0. Sheets "Microbial_Decomposition" and
"Macroinvertebrate_Feeding"; script `data-raw/hermann_2024.R` of the
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
and
[hermann_macroinvertebrates](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md)
code the mesocosms by treatment and replicate number (`"C0A1"`);
[hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[hermann_nutrients](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
hermann_leaf_litter and
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
[`hermann_macroinvertebrates`](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md),
[`hermann_nutrients`](https://pabloamadorc.github.io/must/reference/hermann_nutrients.md),
[`hermann_phytoplankton`](https://pabloamadorc.github.io/must/reference/hermann_phytoplankton.md),
[`hermann_water_quality`](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[`hermann_zooplankton`](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)

## Examples

``` r
with(hermann_leaf_litter,
     round(tapply(mass_loss_percent, list(process, day), mean, na.rm = TRUE), 1))
#>                             11   25   39   53   67   81
#> macroinvertebrate feeding 12.9 13.3 17.5 22.3   NA 25.5
#> microbial                 12.3 10.5 20.1 21.7 26.8 28.7
```
