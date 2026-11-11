# Phytoplankton in a three-stressor mesocosm experiment

Phytoplankton densities from the mesocosm experiment of Polazzo & Rico
(2021) with an insecticide, a herbicide and nutrient enrichment, sampled
once before the pesticide application and three times after it.

## Usage

``` r
polazzo_phytoplankton
```

## Format

A data frame with 96 rows (24 mesocosms x 4 sampling days) and 77
columns:

- mesocosm, insecticide, herbicide, nutrients:

  As in
  [polazzo_zooplankton](https://pabloamadorc.github.io/must/reference/polazzo_zooplankton.md).

- day:

  Day relative to the pesticide application (integer: -5, 7, 15, 50).

- total_density:

  Sum of the 71 taxon columns.

- ankyra:

  Ankyra sp.

- chamydomonas_sp1:

  Chamydomonas sp1.

- chamydomonas_sp3:

  Chamydomonas sp3.

- chlamydomonas_sp4:

  Chlamydomonas sp4.

- chlorella:

  Chlorella sp.

- coelastrum_cf_microporum:

  Coelastrum cf. microporum

- monoraphidium_contortum:

  Monoraphidium contortum

- monoraphidium_minutum:

  Monoraphidium minutum

- neochloris_treuboxia:

  Neochloris/Treuboxia

- oocystis:

  Oocystis sp.

- scenedesmus_quadricauda:

  Scenedesmus quadricauda

- tetraedron_caudatum:

  Tetraedron caudatum

- tetraedron_minutum:

  Tetraedron minutum

- tetraedron_regulare:

  Tetraedron regulare

- closterium:

  Closterium sp.

- cosmarium_sp1:

  Cosmarium sp1.

- cosmarium_sp2:

  Cosmarium sp2.

- cosmarium_sp3:

  Cosmarium sp3.

- staurastrum:

  Staurastrum sp.

- peridinium:

  Peridinium sp.

- cryptophyta_sp1:

  Cryptophyta sp1.

- cryptophyta_sp2:

  Cryptophyta sp2.

- rhodomonas:

  Rhodomonas sp.

- achanthes:

  Achanthes sp.

- amphora:

  Amphora sp.

- cyclotella:

  Cyclotella sp.

- cymbella_minutissima:

  Cymbella minutissima

- epithemia:

  Epithemia sp.

- nitzschia_sigmoidea:

  Nitzschia sigmoidea

- nitzschia:

  Nitzschia sp.

- rhoicosphenia_cf_curvata:

  Rhoicosphenia cf. curvata

- anabaena:

  Anabaena sp.

- chroococcus_cf_limneticus:

  Chroococcus cf. limneticus

- chroococcus:

  Chroococcus sp.

- jaeginema:

  Jaeginema sp.

- microcystis:

  Microcystis sp.

- oscillatoria:

  Oscillatoria sp.

- pseudanabaena:

  Pseudanabaena sp.

- spirulina:

  Spirulina sp.

- monoraphidium_nanum:

  Monoraphidium nanum

- nephrocytium_cf_agardhianum:

  Nephrocytium cf. agardhianum

- pediastrum_cf_boryanum:

  Pediastrum cf. boryanum

- scenedesmus_cf_acunae:

  Scenedesmus cf. acunae

- scenedesmus_cf_ecornis:

  Scenedesmus cf ecornis

- scenedesmus_obliquus:

  Scenedesmus obliquus

- scenedesmus:

  Scenedesmus sp.

- cymbella:

  Cymbella sp.

- eunotia:

  Eunotia sp.

- navicula:

  Navicula sp.

- rhopalodia:

  Rhopalodia sp.

- synedra_acus:

  Synedra acus

- merismopedia_puntacta:

  Merismopedia puntacta

- merismopedia_tenuissima:

  Merismopedia tenuissima

- planktolyngbya_cf_limnetica:

  Planktolyngbya cf. limnetica

- planktothrix:

  Planktothrix sp.

- botryococcus:

  Botryococcus sp.

- pediastrum_tetras:

  Pediastrum tetras

- volvox:

  Volvox sp.

- closterium_sp_2:

  Closterium sp. 2

- cocconeis:

  Cocconeis sp.

- epithemia_sp_1:

  Epithemia sp. 1

- epithemia_sp_2_pequena:

  Epithemia sp. 2 (pequena)

- gomphonema:

  Gomphonema sp

- tryblionella:

  Tryblionella sp.

- komvophoron:

  Komvophoron sp.

- euglena:

  Euglena sp.

- lepocinclis:

  Lepocinclis sp.

- phacus_curvicauda:

  Phacus curvicauda

- trachelomonas_hispida:

  Trachelomonas hispida

- unicellular_flagellate_2_3_4_micras_size:

  Unicellular flagellate 2-3 (4) micras-size

- ciliates:

  Ciliates

The taxon columns hold the density of each taxon found in the source
file on any sampling day, whose name for the taxon is given above and in
the attribute `"taxa"`. Densities are as reported in the source file,
which does not state their unit. The source has one sheet per sampling
day, each with the taxa found that day in some mesocosm; a taxon missing
from a sheet was not found that day, so its density is zero. Zeros are
taxa not found in the sample; there are no missing values.

## Source

Polazzo, F., & Rico, A. (2021). Data from: Effects of multiple stressors
on the dimensionality of ecological stability. Dryad,
[doi:10.5061/dryad.qfttdz0gt](https://doi.org/10.5061/dryad.qfttdz0gt) .
Licence CC0 1.0. Sheets "Phyto_D-5", "Phyto_D7", "Phyto_D15" and
"Phyto_D50" of `Polazzo_and_Rico_2021_EcologyLetters_Full_data.xlsx`
(the day of the last sheet is taken from its name: its day column says
15); script `data-raw/polazzo_rico_2021.R` of the source repository.

## Experimental design

Twenty-four outdoor mesocosms, three per combination of a full factorial
of chlorpyrifos (insecticide, 1 µg/L, or absent), diuron (herbicide, 18
µg/L, or absent) and nutrients (750 µg/L N and 75 µg/L P added twice a
week, or ambient). The nutrient additions started four weeks before the
pesticide application (day 0). The community tables
[polazzo_zooplankton](https://pabloamadorc.github.io/must/reference/polazzo_zooplankton.md),
polazzo_phytoplankton and
[polazzo_macroinvertebrates](https://pabloamadorc.github.io/must/reference/polazzo_macroinvertebrates.md)
share the mesocosm codes and the design columns, so they can be joined
by `mesocosm` and, where the sampling days coincide, `day`.

## References

Polazzo, F., & Rico, A. (2021). Effects of multiple stressors on the
dimensionality of ecological stability. *Ecology Letters*, 24(8),
1594–1606. [doi:10.1111/ele.13770](https://doi.org/10.1111/ele.13770)

## See also

Other data of Polazzo and Rico (2021):
[`polazzo_macroinvertebrates`](https://pabloamadorc.github.io/must/reference/polazzo_macroinvertebrates.md),
[`polazzo_zooplankton`](https://pabloamadorc.github.io/must/reference/polazzo_zooplankton.md)

## Examples

``` r
range(polazzo_phytoplankton$total_density)
#> [1]     16 181676
# Taxa found in each sample
taxa <- names(attr(polazzo_phytoplankton, "taxa"))
summary(rowSums(polazzo_phytoplankton[taxa] > 0))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    2.00    7.00   10.00   10.02   12.00   21.00 
```
