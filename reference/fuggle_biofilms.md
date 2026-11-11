# Biofilm metabolism under salinisation and nutrient enrichment in two orders

Primary production and respiration of benthic biofilms in outdoor ponds
exposed to elevated salinity, nutrient enrichment, or both in either
order (Fuggle et al. 2025).

## Usage

``` r
fuggle_biofilms
```

## Format

A data frame with 54 rows (18 ponds x 3 sampling days) and 9 columns:

- pond:

  Pond code (text), from the treatment and the replicate number, such as
  `"NS_P2"`: the experimental unit.

- treatment:

  Treatment code of the source: `"A"`, `"N"`, `"S"`, `"NS"` or `"SN"`.

- sequence:

  Order of the two stressors: `"nutrients first"` (`"NS"`),
  `"salinity first"` (`"SN"`) or `"single"` for the ponds that both
  orders share (`"A"`, `"N"` and `"S"`).

- nutrients:

  `"added"` or `"none"`.

- salinity:

  `"raised"` or `"none"`.

- day:

  Days after the second application of stressors (integer: 1, 21, 81).

- gross_production:

  Gross primary production of the biofilm of a tile.

- net_production:

  Net primary production.

- respiration:

  Community respiration.

The three rates are as in the source file, whose README gives them in
nmol O2 per cm² and hour. They were not measured on day 81 (the pond
water had become anoxic), where they are missing; there are no other
missing values. Gross production is net production plus respiration in
every sample but one (pond `"A_P3"`, day 1), where the source has net
production minus respiration.

## Source

Fuggle, R., Mayer-Pinto, M., Matias, M., & Marzinelli, Z. (2025). Timing
of salinisation and nutrient enrichment drives freshwater microbial
community metabolic responses. Dryad,
[doi:10.5061/dryad.8pk0p2nz8](https://doi.org/10.5061/dryad.8pk0p2nz8) .
Licence CC0 1.0. File `univariate_data.csv`, the treatments that form an
interaction of the two stressors and the three oxygen rates (script
`data-raw/fuggle_2025.R` of the source repository).

## Experimental design

Outdoor ponds of 1000 L of the Iberian Pond Network at Toledo (Spain),
2023. Biofilms grown for 30 days on ceramic tiles in four ambient ponds
were moved to treatment ponds, three or four per treatment: ambient
(`"A"`); elevated salinity (`"S"`, 15 g NaCl/L) or nutrient enrichment
(`"N"`, 10 mg N/L and 1 mg P/L), applied once; or one of the two
stressors and, a week later, the other (`"SN"`, salinity then nutrients;
`"NS"`, nutrients then salinity). The ponds with a single stressor
received it with the first application. Tiles were sampled 1, 21 and 81
days after the second application. The article also has ponds with the
same stressor applied twice, which do not form an interaction and are
not included.

## References

Fuggle, R., Matias, M. G., Mayer-Pinto, M., & Marzinelli, E. M. (2025).
Repetition and order of exposure to salinization and nutrient enrichment
drive freshwater microbial community responses. *Ecology*, 106(11),
e70260. [doi:10.1002/ecy.70260](https://doi.org/10.1002/ecy.70260)

## See also

Other case-study data:
[`alba_arthropods`](https://pabloamadorc.github.io/must/reference/alba_arthropods.md),
[`giesler_diatoms`](https://pabloamadorc.github.io/must/reference/giesler_diatoms.md),
[`hunn_invertebrates`](https://pabloamadorc.github.io/must/reference/hunn_invertebrates.md),
[`niboyet_denitrification`](https://pabloamadorc.github.io/must/reference/niboyet_denitrification.md),
[`vasquez_embryos`](https://pabloamadorc.github.io/must/reference/vasquez_embryos.md)

## Examples

``` r
str(fuggle_biofilms)
#> 'data.frame':    54 obs. of  9 variables:
#>  $ pond            : chr  "A_P1" "A_P2" "A_P3" "A_P4" ...
#>  $ treatment       : chr  "A" "A" "A" "A" ...
#>  $ sequence        : chr  "single" "single" "single" "single" ...
#>  $ nutrients       : chr  "none" "none" "none" "none" ...
#>  $ salinity        : chr  "none" "none" "none" "none" ...
#>  $ day             : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ gross_production: num  28242 14844 4277 52871 19980 ...
#>  $ net_production  : num  20801 14844 5449 37695 18242 ...
#>  $ respiration     : num  7441 0 1172 15176 1738 ...
with(fuggle_biofilms[fuggle_biofilms$day == 1, ], table(treatment, sequence))
#>          sequence
#> treatment nutrients first salinity first single
#>        A                0              0      4
#>        N                0              0      3
#>        NS               4              0      0
#>        S                0              0      3
#>        SN               0              4      0
# Each order is a `by` group; the ponds of both orders are shared
describe_design(fuggle_biofilms, stressors = c("nutrients", "salinity"),
                control = c(nutrients = "none", salinity = "none"),
                unit = "pond", time = "day", by = "sequence",
                shared = c(sequence = "single"))
#> <must_design> 54 rows, 2 stressors, 5 observed treatment combinations
#> Stressors (control): nutrients (none), salinity (none)
#> Unit: pond | Time: day | By: sequence (shared: sequence = single) | Groups: none
#> Contrast sets: 2 candidate(s), 2 identifiable
#> Checks:
#>   info     Experimental unit: 'pond' (18 units, each in a single treatment
#>            combination).
#>   info     Units per treatment combination: 3-4.
#>   info     Units are measured repeatedly (3 sampling times; 3 per unit).
#>            Samplings of the same unit are not independent: by default
#>            fit_responses() adds a random intercept per unit, which assumes
#>            the same correlation between any two samplings of a unit;
#>            temporal autocorrelation can be modelled through its `random`
#>            argument.
#>   info     Identifiable contrast sets (net deviation): 2 with 2 stressors.
#>   info     2 identifiable set(s) share cells (one cell is used by up to 2
#>            sets, typically the control): their estimates are correlated,
#>            not independent replicates of an interaction.
#>   info     Rows of 'sequence=single' are shared: their 3 treatment
#>            combination(s) ('control', 'salinity=raised', 'nutrients=added')
#>            complete the contrast sets of 2 other `by` group(s), whose
#>            estimates are therefore correlated.
```
