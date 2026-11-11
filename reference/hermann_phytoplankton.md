# Phytoplankton in an insecticide and temperature mesocosm experiment

Phytoplankton densities from the mesocosm experiment of Hermann et al.
(2024) with an insecticide and three temperature regimes, sampled once
before the first insecticide application and six times after it.

## Usage

``` r
hermann_phytoplankton
```

## Format

A data frame with 168 rows (24 mesocosms x 7 sampling days) and 96
columns:

- mesocosm, imidacloprid, temperature:

  As in
  [hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md).

- day:

  Day relative to the first imidacloprid application (integer: -1, 10,
  24, 38, 52, 66, 80).

- total_density:

  Sum of the 91 taxon columns.

- rhodomonas:

  Rhodomonas sp

- cryptomonas:

  Cryptomonas sp

- chlamydomonas:

  Chlamydomonas

- scenedesmus:

  Scenedesmus

- scenedesmus_obliquus:

  Scenedesmus obliquus

- desmodesmus_quadricauda:

  Desmodesmus quadricauda

- desmodesmus_2_cells:

  Desmodesmus (2 cells)

- scenedesmus_accuminatum:

  Scenedesmus accuminatum

- scenedesmus_arcuatus:

  Scenedesmus arcuatus

- ankyra:

  Ankyra sp

- ankyra_sp_2_lanceolada:

  Ankyra sp 2 (lanceolada)

- staurodesmus:

  Staurodesmus

- pediastrum_simplex:

  Pediastrum simplex

- pediastrum_tetras:

  Pediastrum tetras

- pediastrum_boryanum:

  Pediastrum boryanum

- pediastrum_duplex:

  Pediastrum duplex

- ankystrodesmus:

  Ankystrodesmus

- eudorina:

  Eudorina

- coelastrum:

  Coelastrum

- monoraphidium_sp1:

  Monoraphidium sp1

- monoraphidium_sp2:

  Monoraphidium sp2

- tetraedron:

  Tetraedron

- oocystis:

  Oocystis

- volvox:

  Volvox

- coccal_cells_5_um:

  Coccal cells \< 5 um

- nephrocytium:

  Nephrocytium

- oedogonium:

  Oedogonium

- peridinium:

  Peridinium

- closterium:

  Closterium

- closterium_lineal:

  Closterium lineal

- cosmarium_sp_1:

  Cosmarium sp 1

- cosmarium_sp_2_small:

  Cosmarium sp 2 (small)

- staurastrum:

  Staurastrum sp

- pleurotaenium:

  Pleurotaenium

- desmidium:

  Desmidium

- coleochaete:

  Coleochaete

- gyrosigma:

  Gyrosigma

- cymatopleura_solea:

  Cymatopleura solea

- pinnularia:

  Pinnularia

- navicula_sp_1:

  Navicula sp 1

- navicula_sp_2_10_um:

  Navicula sp 2 (10 um)

- nitzschia_sp_3:

  Nitzschia sp 3

- nitszchia_acicularis_sp_1:

  Nitszchia (acicularis) sp 1

- nitzschia_recta:

  Nitzschia recta

- nitzschia_sp_2_dubia:

  Nitzschia sp 2 (dubia)

- pseudo_nizschia:

  Pseudo-nizschia

- fragilaria:

  Fragilaria

- cymbella:

  Cymbella

- synedra:

  Synedra sp

- rhopalodia:

  Rhopalodia sp

- cocconeis:

  Cocconeis

- amphora:

  Amphora

- eunotia:

  Eunotia

- epithemia:

  Epithemia

- diatoma_sp_2:

  Diatoma sp 2

- achnanthes:

  Achnanthes

- achnanthes_5_um:

  Achnanthes \< 5 um

- rhoicosphenia_cf_curvata:

  Rhoicosphenia cf curvata

- diatoma_sp1_vulgaris:

  Diatoma sp1 (vulgaris)

- central_diatom:

  Central diatom

- reimeria:

  Reimeria

- surirella:

  Surirella

- frustulia:

  Frustulia

- diploneis:

  Diploneis

- pennal_diatom_si_30_um:

  Pennal diatom si \> 30 um

- gomphonema:

  Gomphonema

- synura:

  Synura

- trachelomonas_armata:

  Trachelomonas armata

- trachelomonas_hispida:

  Trachelomonas hispida

- trachelomonas_sp_1_volvocina:

  Trachelomonas sp 1 (volvocina)

- lepocinclis_acus:

  Lepocinclis acus

- euglenophyta:

  Euglenophyta

- phacus_longicauda:

  Phacus longicauda

- euglena_ehrenbergii:

  Euglena ehrenbergii

- phacus_sp_triqueter:

  Phacus sp (triqueter)

- phacus_contortum:

  Phacus contortum

- strombomonas:

  Strombomonas

- nostoc:

  Nostoc

- merismopedia:

  Merismopedia

- chroococcus:

  Chroococcus

- oscillatoria_in_filaments:

  Oscillatoria (in filaments)

- komvophoron:

  Komvophoron

- aphanocapsa:

  Aphanocapsa

- filamentous_cyanobacteria:

  Filamentous cyanobacteria

- oscillatorial:

  Oscillatorial

- spirulina:

  Spirulina

- gloeocapsa:

  Gloeocapsa

- dolichospermum:

  Dolichospermum

- ciliates_1:

  Ciliates 1

- ciliates_2:

  Ciliates 2

- vorticella:

  Vorticella

The taxon columns hold the density of each taxon of the source file
(ciliates included), whose name for the taxon is given above and in the
attribute `"taxa"`. Densities are as reported in the source file, which
does not state their unit. Zeros are taxa not found in the sample; there
are no missing values.

## Source

Hermann et al. (2024), Zenodo,
[doi:10.5281/zenodo.11520459](https://doi.org/10.5281/zenodo.11520459)
(full reference in
[hermann_zooplankton](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)).
Licence CC BY 4.0. Sheets "PhytoplanktonDesign" and
"Phytoplankton_Abundances"; script `data-raw/hermann_2024.R` of the
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
hermann_phytoplankton and
[hermann_macroinvertebrates](https://pabloamadorc.github.io/must/reference/hermann_macroinvertebrates.md)
code the mesocosms by treatment and replicate number (`"C0A1"`);
[hermann_water_quality](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
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
[`hermann_water_quality`](https://pabloamadorc.github.io/must/reference/hermann_water_quality.md),
[`hermann_zooplankton`](https://pabloamadorc.github.io/must/reference/hermann_zooplankton.md)

## Examples

``` r
# Total density, median by treatment combination and day
with(hermann_phytoplankton,
     tapply(total_density, list(paste(imidacloprid, temperature), day), median))
#>                          -1        10        24         38         52
#> absent ambient    293692000 216558750 406389250 1251727000 1110223250
#> absent elevated   175444250 179265000 127172750  380268750  465200750
#> absent heatwaves  435854500  67345500 216240000   65012250  241595500
#> applied ambient   242738750 359341750 524224750  610691000 1032444000
#> applied elevated  779182250  94826000 125927500  548165000  509532500
#> applied heatwaves 437053000 120164500 184849500  149383250  541352250
#>                           66        80
#> absent ambient     685724750 402730000
#> absent elevated    721433250 332630500
#> absent heatwaves  2515320000 329876500
#> applied ambient   1095150625 600002250
#> applied elevated  1270401500 188260125
#> applied heatwaves  823909250 307100750
```
