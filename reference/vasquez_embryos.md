# Development of horseshoe crab embryos under temperature, salinity and oxygen

Embryos of the American horseshoe crab (*Limulus polyphemus*) that
developed in each of 36 combinations of temperature, salinity and
oxygen, by mating pair (Vasquez et al. 2017).

## Usage

``` r
vasquez_embryos
```

## Format

A data frame with 216 rows (36 combinations x 6 mating pairs) and 6
columns:

- temperature:

  Temperature (°C, integer: 20, 25, 30, 35).

- salinity:

  Salinity (ppt, integer: 5, 15, 34).

- oxygen:

  Oxygen in the air (% O2, integer: 5, 13, 21).

- pair:

  Mating pair (integer, 1 to 6), whose eggs were spread over all
  combinations.

- embryos:

  Embryos measured (integer, 8, or 7 in three cells).

- developed:

  Embryos that reached 5 mm² by day 28 (integer).

An embryo counts as developed when the event indicator of the source,
used by the survival analysis of the article, became 1. Three embryos
are left out: two never measured and one without that indicator. There
are no missing values.

## Source

Vasquez, M. C., Brockmann, H. J., & Julian, D. (2018). Data from:
Between-population differences in multi-stressor tolerance during embryo
development in the American horseshoe crab, *Limulus polyphemus*. Dryad,
[doi:10.5061/dryad.vs1r0](https://doi.org/10.5061/dryad.vs1r0) . Licence
CC0 1.0. Sheet "Embryo Development Size" of
`Dryad+-+Embryo+Cross-Sectional+Area.xlsx`, one embryo per row and day,
summed by pair and combination (script `data-raw/vasquez_2017.R` of the
source repository).

## Experimental design

Fertilized eggs of six mating pairs from Slaughter Beach (Delaware Bay,
USA), May 2012, one egg per well of 96-well plates, were exposed for 14
days to each of the 36 combinations of temperature (20, 25, 30 and 35
°C), salinity (5, 15 and 34 ppt) and oxygen (5, 13 and 21% O2 in the
air), one plate per combination, and then kept for 14 days in the
reference condition (30 °C, 34 ppt and 21% O2). Every 2 days the eggs
were imaged, and the cross-sectional area of eight eggs per pair and
combination was measured. An embryo developed when it reached 5 mm², the
size of a late-stage embryo. The plate of each combination is confounded
with it.

## References

Vasquez, M. C., Brockmann, H. J., & Julian, D. (2017).
Between-population differences in multi-stressor tolerance during embryo
development in the American horseshoe crab, *Limulus polyphemus*.
*Estuaries and Coasts*, 40(5), 1487–1501.
[doi:10.1007/s12237-017-0218-1](https://doi.org/10.1007/s12237-017-0218-1)

## See also

Other case-study data:
[`alba_arthropods`](https://pabloamadorc.github.io/must/reference/alba_arthropods.md),
[`fuggle_biofilms`](https://pabloamadorc.github.io/must/reference/fuggle_biofilms.md),
[`giesler_diatoms`](https://pabloamadorc.github.io/must/reference/giesler_diatoms.md),
[`hunn_invertebrates`](https://pabloamadorc.github.io/must/reference/hunn_invertebrates.md),
[`niboyet_denitrification`](https://pabloamadorc.github.io/must/reference/niboyet_denitrification.md)

## Examples

``` r
str(vasquez_embryos)
#> 'data.frame':    216 obs. of  6 variables:
#>  $ temperature: int  20 20 20 20 20 20 20 20 20 20 ...
#>  $ salinity   : int  5 5 5 5 5 5 5 5 5 5 ...
#>  $ oxygen     : int  5 5 5 5 5 5 13 13 13 13 ...
#>  $ pair       : int  1 2 3 4 5 6 1 2 3 4 ...
#>  $ embryos    : int  8 8 8 8 8 8 8 8 8 8 ...
#>  $ developed  : int  6 0 0 1 0 5 0 0 2 0 ...
# Share of embryos that developed in each combination (all pairs)
n <- with(vasquez_embryos, tapply(developed, list(temperature, salinity, oxygen), sum))
d <- with(vasquez_embryos, tapply(embryos, list(temperature, salinity, oxygen), sum))
round(n / d, 2)
#> , , 5
#> 
#>       5   15   34
#> 20 0.25 0.94 0.85
#> 25 0.10 0.88 0.92
#> 30 0.42 0.96 0.83
#> 35 0.00 0.40 0.00
#> 
#> , , 13
#> 
#>       5   15   34
#> 20 0.06 0.96 0.90
#> 25 0.31 0.92 0.85
#> 30 0.04 0.94 0.81
#> 35 0.00 0.81 0.65
#> 
#> , , 21
#> 
#>       5   15   34
#> 20 0.77 0.90 0.94
#> 25 0.27 0.85 0.79
#> 30 0.00 0.88 0.88
#> 35 0.00 0.04 0.77
#> 
```
