# Macroinvertebrates in a three-stressor mesocosm experiment

Macroinvertebrate counts from the mesocosm experiment of Polazzo & Rico
(2021) with an insecticide, a herbicide and nutrient enrichment, sampled
once before the pesticide application and three times after it.

## Usage

``` r
polazzo_macroinvertebrates
```

## Format

A data frame with 96 rows (24 mesocosms x 4 sampling days) and 54
columns:

- mesocosm, insecticide, herbicide, nutrients:

  As in
  [polazzo_zooplankton](https://pabloamadorc.github.io/must/reference/polazzo_zooplankton.md).

- day:

  Day relative to the pesticide application (integer: -32, 15, 30, 50).
  Day -32 is before the nutrient additions too.

- total_abundance:

  Sum of the 48 taxon columns.

- erpobdella:

  Erpobdella sp.

- lumbriculidae:

  Lumbriculidae

- hydracarina:

  Hydracarina

- sphaeriidae:

  Sphaeriidae

- dytiscidae_larvae:

  Dytiscidae larvae

- elmidae_adult:

  Elmidae adult

- elmidae_larvae:

  Elmidae larvae

- echinogammarus:

  Echinogammarus

- chironomini:

  Chironomini

- corynoneura:

  Corynoneura

- culicidae:

  Culicidae

- diptera_pupa:

  Diptera pupa

- orthocladiinae:

  Orthocladiinae

- sciomyzidae:

  Sciomyzidae

- tanytarsini:

  Tanytarsini

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

- anisoptera_ii_stripped:

  Anisoptera II (stripped)

- zygoptera:

  Zygoptera

- dugesia:

  Dugesia sp.

- polycelis:

  Polycelis sp.

- nemouridae:

  Nemouridae

- leptoceratidae_stone_house:

  Leptoceratidae (stone house)

- leptoceridae_wood_house:

  Leptoceridae (wood house)

- lymnaeidae:

  Lymnaeidae

- hydropsychidae:

  Hydropsychidae

- haliplidae_adult:

  Haliplidae adult

- polypocentridae:

  Polypocentridae

- limnius:

  Limnius sp.

- hygrobiidae:

  Hygrobiidae

- simulidae:

  Simulidae

- enchytraeidae:

  Enchytraeidae

- chaoborus:

  Chaoborus sp.

- dytiscidae_adult:

  Dytiscidae adult

- potamanthidae:

  Potamanthidae

- plea:

  Plea sp.

- hydrophilidae:

  Hydrophilidae

- salix:

  Salix sp.

- lumbriciadae:

  Lumbriciadae

- tabanidae_larvae:

  Tabanidae larvae

- diptera_larvae:

  Diptera larvae

The taxon columns hold the number of individuals of each taxon of the
source file, whose name for the taxon is given above and in the
attribute `"taxa"`. Counts are as reported in the source file. Zeros are
taxa not found in the sample. One count is missing (`hydropsychidae`,
mesocosm `"M23"`, day 50), and so is the total of that sample.

## Source

Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple stressors
on the dimensionality of ecological stability. Dryad,
[doi:10.5061/dryad.qfttdz0gt](https://doi.org/10.5061/dryad.qfttdz0gt) .
Licence CC0 1.0. Sheet "Macroinvertebrates" of
`Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx`, without an
unlabelled column that has values for one mesocosm per treatment only;
script `data-raw/polazzo_rico_2021.R` of the source repository.

## Experimental design

Twenty-four outdoor mesocosms, three per combination of a full factorial
of chlorpyrifos (insecticide, 1 µg/L, or absent), diuron (herbicide, 18
µg/L, or absent) and nutrients (750 µg/L N and 75 µg/L P added twice a
week, or ambient). The nutrient additions started four weeks before the
pesticide application (day 0). The community tables
[polazzo_zooplankton](https://pabloamadorc.github.io/must/reference/polazzo_zooplankton.md),
[polazzo_phytoplankton](https://pabloamadorc.github.io/must/reference/polazzo_phytoplankton.md)
and polazzo_macroinvertebrates share the mesocosm codes and the design
columns, so they can be joined by `mesocosm` and, where the sampling
days coincide, `day`.

## References

Polazzo, F., & Rico, A. (2021). Effects of multiple stressors on the
dimensionality of ecological stability. *Ecology Letters*, 24(8),
1594–1606. [doi:10.1111/ele.13770](https://doi.org/10.1111/ele.13770)

## See also

Other data of Polazzo and Rico (2021):
[`polazzo_phytoplankton`](https://pabloamadorc.github.io/must/reference/polazzo_phytoplankton.md),
[`polazzo_zooplankton`](https://pabloamadorc.github.io/must/reference/polazzo_zooplankton.md)

## Examples

``` r
range(polazzo_macroinvertebrates$total_abundance, na.rm = TRUE)
#> [1]   65 1207
polazzo_macroinvertebrates[is.na(polazzo_macroinvertebrates$total_abundance), 1:6]
#>    mesocosm day  insecticide herbicide nutrients total_abundance
#> 95      M23  50 chlorpyrifos    diuron     added              NA
```
