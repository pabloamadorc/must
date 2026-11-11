# Growth of two diatoms under temperature, photoperiod and microbiome removal

Maximum growth rates of a polar and a temperate diatom at three
temperatures and three photoperiods, with and without their microbiome
(Giesler et al. 2023).

## Usage

``` r
giesler_diatoms
```

## Format

A data frame with 1728 rows (2 species x 2 microbiome treatments x 3
temperatures x 3 photoperiods x 48 wells) and 6 columns:

- species:

  `"T. gravida"` (polar) or `"T. rotula"` (temperate).

- microbiome:

  `"present"` (xenic culture) or `"removed"` (axenic).

- temperature:

  Temperature (°C: 4, 9 or 13.5).

- photoperiod:

  Hours of light per day (integer: 4, 16 or 24).

- replicate:

  Well number within the combination (integer, 1 to 48).

- growth_rate:

  Maximum growth rate (per day). Negative values are wells whose
  fluorescence declined.

There are no missing values.

## Source

Giesler, J., Wohlrab, S., & Harder, T. (2023). Growth of *Thalassiosira
gravida* and *Thalassiosira rotula* across different levels of
temperature, photoperiod, and microbiome presence. Dryad,
[doi:10.5061/dryad.x95x69pqf](https://doi.org/10.5061/dryad.x95x69pqf) .
Licence CC0 1.0. File `plateassay_growthrates.xlsx`, with the columns
renamed (script `data-raw/giesler_2023.R` of the source repository).

## Experimental design

Laboratory cultures of the polar diatom *Thalassiosira gravida* and the
temperate *T. rotula*, each as a xenic culture (with its microbiome) and
an axenic one (microbiome removed), grown at each combination of three
temperatures (4, 9 and 13.5 °C) and three photoperiods (4, 16 and 24 h
of light at 50 µmol of photons per m² and second) after 7 days of
acclimation. Each combination had 48 replicate wells of 300 µL in
96-well plates, incubated in climate cabinets at the temperature of the
treatment, on LED tables with its photoperiod. The fluorescence of each
well was read daily for 7 days, and its maximum growth rate was fitted
to those readings. The data do not identify plates, cabinets or tables:
the wells of a combination shared its incubation conditions.

## References

Giesler, J. K., Harder, T., & Wohlrab, S. (2023). Microbiome and
photoperiod interactively determine thermal sensitivity of polar and
temperate diatoms. *Biology Letters*, 19(11), 20230151.
[doi:10.1098/rsbl.2023.0151](https://doi.org/10.1098/rsbl.2023.0151)

## See also

Other case-study data:
[`alba_arthropods`](https://pabloamadorc.github.io/must/reference/alba_arthropods.md),
[`fuggle_biofilms`](https://pabloamadorc.github.io/must/reference/fuggle_biofilms.md),
[`hunn_invertebrates`](https://pabloamadorc.github.io/must/reference/hunn_invertebrates.md),
[`niboyet_denitrification`](https://pabloamadorc.github.io/must/reference/niboyet_denitrification.md),
[`vasquez_embryos`](https://pabloamadorc.github.io/must/reference/vasquez_embryos.md)

## Examples

``` r
str(giesler_diatoms)
#> 'data.frame':    1728 obs. of  6 variables:
#>  $ species    : chr  "T. gravida" "T. gravida" "T. gravida" "T. gravida" ...
#>  $ microbiome : chr  "present" "present" "present" "present" ...
#>  $ temperature: num  4 4 4 4 4 4 4 4 4 4 ...
#>  $ photoperiod: int  4 4 4 4 4 4 4 4 4 4 ...
#>  $ replicate  : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ growth_rate: num  0.184 0.183 0.179 0.168 0.158 ...
with(giesler_diatoms[giesler_diatoms$species == "T. gravida", ],
     round(tapply(growth_rate, list(temperature, photoperiod, microbiome), mean), 2))
#> , , present
#> 
#>         4   16   24
#> 4    0.17 0.42 0.43
#> 9    0.39 0.68 0.63
#> 13.5 0.39 0.61 0.34
#> 
#> , , removed
#> 
#>          4   16   24
#> 4     0.09 0.44 0.41
#> 9     0.03 0.44 0.46
#> 13.5 -0.02 0.27 0.51
#> 
```
