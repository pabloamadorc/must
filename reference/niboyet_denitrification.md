# Potential denitrification under five global changes

Potential denitrification and soil water in a grassland under a full
factorial of elevated CO2, warming, increased precipitation, nitrogen
addition and fire, sampled five times after the fire (Niboyet et al.
2025).

## Usage

``` r
niboyet_denitrification
```

## Format

A data frame with 640 rows (128 quadrants x 5 samplings) and 12 columns:

- block:

  Block (integer, 1 to 4).

- plot:

  Plot (integer, 1 to 32), to which CO2, temperature and fire were
  applied.

- quadrant:

  Quadrant within the plot (integer, 1 to 4).

- unit:

  Quadrant identifier (integer, 1 to 128), the experimental unit of
  precipitation and nitrogen.

- co2:

  `"ambient"` or `"elevated"`.

- temperature:

  `"ambient"` or `"elevated"`.

- water:

  Precipitation: `"ambient"` or `"elevated"`.

- nitrogen:

  `"ambient"` or `"elevated"`.

- fire:

  `"unburned"` or `"burned"`.

- months:

  Months after the fire (integer: 0, 5, 9, 21, 33).

- denitrification:

  Potential denitrification, measured as denitrifying enzyme activity
  (ng N2O-N per g of dry soil and hour).

- soil_water:

  Soil water content (%).

Values and units as in the source file. One denitrification value and
six soil water values are missing in the source.

## Source

Niboyet, A. (2025). Soil moisture and the number of simultaneous
stressors drive interactions among global changes on denitrification
(version v1). Zenodo,
[doi:10.5281/zenodo.16789113](https://doi.org/10.5281/zenodo.16789113) .
Licence CC BY 4.0. File
`Niboyetetal_CommunicationsEarthEnvironment_10august2025.xlsx`, with the
treatment codes as labels and two of its variables (script
`data-raw/niboyet_2025.R` of the source repository).

## Experimental design

Jasper Ridge Global Change Experiment (California, USA), an annual
grassland under CO2 enrichment (+275 µmol/mol), warming (+2.0 °C),
increased precipitation (+50% and a rainy season 3 weeks longer) and
nitrogen addition (+7 g N per m² and year) since 1998, and a controlled
fire on half of the plots in July 2011. Randomized block split-plot
design: CO2, temperature and fire were applied to circular plots of 2 m
diameter, and precipitation and nitrogen to the four quadrants of each
plot. With 4 blocks of 8 plots, each of the 32 combinations has four
quadrants. Soil cores of every quadrant were taken 2 days after the fire
and 5, 9, 21 and 33 months later.

## References

Niboyet, A., Le Roux, X., Chiariello, N. R., Dijkstra, P., Field, C. B.,
& Hungate, B. A. (2025). Soil moisture and the number of simultaneous
stressors drive interactions among global changes on denitrification.
*Communications Earth & Environment*, 6, 704.
[doi:10.1038/s43247-025-02703-5](https://doi.org/10.1038/s43247-025-02703-5)

## See also

Other case-study data:
[`alba_arthropods`](https://pabloamadorc.github.io/must/reference/alba_arthropods.md),
[`fuggle_biofilms`](https://pabloamadorc.github.io/must/reference/fuggle_biofilms.md),
[`giesler_diatoms`](https://pabloamadorc.github.io/must/reference/giesler_diatoms.md),
[`hunn_invertebrates`](https://pabloamadorc.github.io/must/reference/hunn_invertebrates.md),
[`vasquez_embryos`](https://pabloamadorc.github.io/must/reference/vasquez_embryos.md)

## Examples

``` r
str(niboyet_denitrification)
#> 'data.frame':    640 obs. of  12 variables:
#>  $ block          : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ plot           : int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ quadrant       : int  1 2 3 4 1 2 3 4 1 2 ...
#>  $ unit           : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ co2            : chr  "ambient" "ambient" "ambient" "ambient" ...
#>  $ temperature    : chr  "ambient" "ambient" "ambient" "ambient" ...
#>  $ water          : chr  "elevated" "ambient" "ambient" "elevated" ...
#>  $ nitrogen       : chr  "ambient" "elevated" "ambient" "elevated" ...
#>  $ fire           : chr  "unburned" "unburned" "unburned" "unburned" ...
#>  $ months         : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ denitrification: num  537 503 299 566 326 ...
#>  $ soil_water     : num  2.78 3.33 3.55 3.37 3.03 ...
# CO2, temperature and fire are constant within each plot
head(unique(niboyet_denitrification[c("plot", "block", "co2", "temperature", "fire")]))
#>    plot block      co2 temperature     fire
#> 1     1     1  ambient     ambient unburned
#> 5     2     1 elevated    elevated unburned
#> 9     3     1  ambient    elevated unburned
#> 13    4     1 elevated     ambient unburned
#> 17    5     2  ambient     ambient unburned
#> 21    6     2  ambient    elevated unburned
with(niboyet_denitrification, round(tapply(soil_water, months, mean, na.rm = TRUE), 1))
#>    0    5    9   21   33 
#>  2.8 11.8 19.1  7.9 13.9 
```
