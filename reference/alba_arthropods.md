# Arthropods in a drought and plant invasion field experiment

Arthropod counts from a field experiment with a full factorial of
simulated drought and invasion by cogongrass (*Imperata cylindrica*) in
regenerating longleaf pine understory (Alba et al. 2026), sampled twice
a year for two years.

## Usage

``` r
alba_arthropods
```

## Format

A data frame with 160 rows (40 plots x 4 samplings) and 18 columns:

- plot:

  Plot number (integer, 1 to 40), the experimental unit.

- block:

  Block number (integer, 1 to 10).

- sampling:

  Sampling number (integer, 1 to 4) in chronological order: summer 2014,
  fall 2014, summer 2015 and fall 2015.

- year:

  Year of the sampling (integer, 2014 or 2015).

- season:

  `"summer"` (July) or `"fall"` (October).

- drought:

  `"drought"` or `"ambient"`.

- invasion:

  `"invaded"` or `"uninvaded"`.

- total:

  Sum of the 10 group columns.

- hemipterans:

  Hemiptera (source code `hemtotal`).

- caterpillars:

  Lepidoptera larvae (`lepnumlarva`).

- beetles:

  Coleoptera (`coletotal`).

- grasshoppers:

  Orthoptera (`orthtotal`).

- flies:

  Adult Diptera (`dipnumadult`).

- moths:

  Adult Lepidoptera (`lepnumadult`).

- bees:

  Bees (`beetotal`).

- wasps:

  Wasps (`wasptotal`).

- spiders:

  Spiders (`spidnumindiv`).

- ants:

  Ants (`anttotal`).

The group columns, named after the article, count the individuals of
each group in the sample of a plot; the attribute `"taxa"` gives the
code of each group in the source file. The source assigns hemipterans,
caterpillars, beetles and grasshoppers to herbivores, flies, moths and
bees to pollinators, and wasps, spiders and ants to predators. Zeros are
groups not found in the sample; there are no missing values.

## Source

Alba, C., Lucky, A., Fahey, C., Somarriba, G., Dunbar, M., Inman, W.,
Flory, S. L., & Angelini, C. (2025). Data from: Plant invasion overrides
drought in suppressing arthropod abundance and richness. Dryad,
[doi:10.5061/dryad.02v6wwqhz](https://doi.org/10.5061/dryad.02v6wwqhz) .
Licence CC0 1.0. Sheet "Abundance_Raw_Data" of
`Data_For_Repository_Arthropods.xlsx`, one column per group; treatments
from the treatment code, which in the fourth sampling disagrees with the
invasion column for plot 38 (script `data-raw/alba_2026.R` of the source
repository).

## Experimental design

Forty plots of 4 x 4 m in an old field at Gainesville (Florida, USA),
planted in 2012 with longleaf pine and native understory species, in 10
blocks with one plot of each combination of drought (rainout shelters,
or ambient precipitation under control structures) and invasion (nine
cogongrass ramets planted, or none). Both treatments started in spring
2013. Arthropods were sampled with 20 sweeps of a net per plot in July
(summer) and October (fall) of 2014 and 2015.

## References

Alba, C., Lucky, A., Fahey, C., Somarriba, G., Dunbar, M., Inman, W.,
Flory, S. L., & Angelini, C. (2026). Plant invasion overrides drought in
suppressing arthropod abundance and richness. *Oikos*, 2026(5), e11920.
[doi:10.1002/oik.11920](https://doi.org/10.1002/oik.11920)

## See also

Other case-study data:
[`fuggle_biofilms`](https://pabloamadorc.github.io/must/reference/fuggle_biofilms.md),
[`giesler_diatoms`](https://pabloamadorc.github.io/must/reference/giesler_diatoms.md),
[`hunn_invertebrates`](https://pabloamadorc.github.io/must/reference/hunn_invertebrates.md),
[`niboyet_denitrification`](https://pabloamadorc.github.io/must/reference/niboyet_denitrification.md),
[`vasquez_embryos`](https://pabloamadorc.github.io/must/reference/vasquez_embryos.md)

## Examples

``` r
str(alba_arthropods[1:9])
#> 'data.frame':    160 obs. of  9 variables:
#>  $ plot       : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ block      : int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ sampling   : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ year       : int  2014 2014 2014 2014 2014 2014 2014 2014 2014 2014 ...
#>  $ season     : chr  "summer" "summer" "summer" "summer" ...
#>  $ drought    : chr  "drought" "ambient" "drought" "ambient" ...
#>  $ invasion   : chr  "invaded" "uninvaded" "uninvaded" "invaded" ...
#>  $ total      : int  22 7 79 58 121 65 31 59 128 44 ...
#>  $ hemipterans: int  16 3 49 31 53 42 13 18 94 24 ...
# One plot of each treatment in every block
with(alba_arthropods[alba_arthropods$sampling == 1, ], table(block, paste(drought, invasion)))
#>      
#> block ambient invaded ambient uninvaded drought invaded drought uninvaded
#>    1                1                 1               1                 1
#>    2                1                 1               1                 1
#>    3                1                 1               1                 1
#>    4                1                 1               1                 1
#>    5                1                 1               1                 1
#>    6                1                 1               1                 1
#>    7                1                 1               1                 1
#>    8                1                 1               1                 1
#>    9                1                 1               1                 1
#>    10               1                 1               1                 1
describe_design(alba_arthropods, stressors = c("drought", "invasion"),
                control = c(drought = "ambient", invasion = "uninvaded"),
                unit = "plot", time = "sampling", groups = "block")
#> <must_design> 160 rows, 2 stressors, 4 observed treatment combinations
#> Stressors (control): drought (ambient), invasion (uninvaded)
#> Unit: plot | Time: sampling | By: none | Groups: block
#> Contrast sets: 1 candidate(s), 1 identifiable
#> Checks:
#>   info     Experimental unit: 'plot' (40 units, each in a single treatment
#>            combination).
#>   info     Units per treatment combination: 10.
#>   info     Units are measured repeatedly (4 sampling times; 4 per unit).
#>            Samplings of the same unit are not independent: by default
#>            fit_responses() adds a random intercept per unit, which assumes
#>            the same correlation between any two samplings of a unit;
#>            temporal autocorrelation can be modelled through its `random`
#>            argument.
#>   info     Identifiable contrast sets (net deviation): 1 with 2 stressors.
#>   info     'block': 10 groups with 4 treatment combination(s) each
#>            (complete blocks); no stressor is applied to whole groups.
#>            fit_responses() adds a random intercept for 'block'.
```
