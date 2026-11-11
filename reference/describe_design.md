# Describe what a multiple-stressor design can identify

Inspects the treatment structure of an experiment before any model is
fitted: which combinations of stressor levels were observed and which
are missing, which interactions each group of cells identifies, which
cells are shared by several of those groups, how the experimental units
are replicated and whether they were measured repeatedly.
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
refuses designs with problems, and
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
and
[`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)
only work on the contrast sets found here.

## Usage

``` r
describe_design(
  data,
  stressors,
  control,
  unit = NULL,
  time = NULL,
  by = NULL,
  groups = NULL,
  max_order = NULL,
  baseline = NULL,
  background = "control",
  shared = NULL
)
```

## Arguments

- data:

  A data frame with one row per observation. Tibbles are accepted;
  grouped data frames are not.

- stressors:

  Character vector with the names of at least two stressor columns. A
  column may hold factors, text or numbers; every distinct value is
  treated as a level (a category, not a gradient).

- control:

  Named character vector (or list) with, for every stressor, the level
  that represents its absence (control or ambient), e.g.
  `c(temperature = "20", nitrate = "0")`. There is no default: which
  level means "no stressor" is a scientific decision.

- unit:

  Name of the column that identifies the experimental (replication)
  unit, or `NULL` (default) if every row is an independent unit.
  Required when `time` is given.

- time:

  Name of the column with the sampling time, or `NULL` (default) for a
  single sampling.

- by:

  Optional names of columns with other factors (population, species,
  experiment). Contrast sets are formed within each combination of their
  levels.

- groups:

  Optional names of columns that group experimental units: blocks, whole
  plots of a split-plot design, rooms, incubators, header tanks or
  batches. For each one the checks report which stressors were applied
  to whole groups and whether their effects can be separated from
  differences among groups, and
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
  adds a random intercept for it when it can be estimated. Group labels
  must be unique within each `by` group (number plots 1-32, not 1-8
  within each block).

- max_order:

  Largest number of stressors in a contrast set. `NULL` (default) means
  all stressors.

- baseline:

  Optional values of `time` taken before the stressors were applied
  (pre-treatment samplings), e.g. `baseline = 0`. They are kept in the
  models fitted by
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md),
  where they inform the differences among units (with
  `time_structure = "additive"` they are coded as untreated, the control
  combination), but
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
  does not estimate or test interactions at them, so they do not enter
  the multiple-testing adjustment either. If the stressors were applied
  at different times, give a named list with the times before each one,
  e.g. `baseline = list(drying = c(-36, -1))` when drying started on day
  0 and the other stressors were applied from the start: interactions
  that involve drying are then not estimated at those times, and those
  among the other stressors are (with `time_structure = "additive"`,
  drying is coded at its control there). Stressors left out of the list
  were applied before the first sampling. Each unit keeps the treatment
  it was assigned in the stressor columns. Requires `time`.

- background:

  `"control"` (default) or `"all"`. With `"control"`, interactions are
  estimated with the stressors outside each contrast set at their
  control (interactions "in the absence of the other stressors"). With
  `"all"`, sets are also formed with the other stressors held at each of
  their observed non-control levels (for example, the interaction
  between salinity and warming in the presence of nutrients); their
  reference cell is then that background combination instead of the
  control.

- shared:

  Optional named character vector that marks rows shared by every group
  of a `by` column, e.g. `c(regime = "both")`: the rows whose `regime`
  is `"both"` (usually the control, sometimes also the single-stressor
  treatments) belong to every other level of `regime`. Use it when one
  control serves several `by` groups (two exposure regimes, two orders
  of application, upstream and downstream ponds), so that it enters the
  analysis once instead of being copied into each group. Contrast sets
  are still formed within each `by` group, with the shared rows
  supplying the cells that the group lacks; sets of different groups
  that use them are correlated. Each name must be a `by` column and each
  value one of its levels.

## Value

An object of class `must_design`, a list with:

- `spec`: the arguments used, with the observed levels of each stressor
  (control first).

- `cells`: data frame with one row per observed treatment combination
  within each `by` group: `context` (label of the `by` group, `""`
  without `by`), the `by` columns and the stressor columns (as text),
  `cell` (label), `n_stressed` (number of stressors not at their
  control), `n_units`, `n_obs` and, with `time`, `n_times`. Ordered by
  `context`, `n_stressed` and stressor levels.

- `missing_cells`: combinations of the observed levels that are absent,
  within each `by` group (same columns as `cells` without the counts).

- `sets`: one row per candidate contrast set, that is, per observed cell
  with two or more stressors away from control (and, with
  `background = "all"`, per subset of two or more of those stressors):
  `set_id`, `context`, `order` (number of stressors), `stressors`,
  `levels` (label of the combined cell), `background` (levels of the
  stressors outside the set that are not at control, `""` if none),
  `net_identifiable`, `pure_identifiable` (`NA` for two stressors, where
  net and pure coincide), `missing_for_net`, `missing_for_pure`,
  `min_units` (fewest units among the cells needed for the net
  deviation), `shared` (`TRUE` if the set uses a cell that another
  identifiable set also uses) and `confounded` (`TRUE` if the set
  involves a stressor applied to whole groups of a `groups` column and
  some of its combinations were given to a single group).

- `groups`: one row per `groups` column: `group`, `n_groups`, `between`
  (stressors applied to whole groups), `structure` (`"blocks"`,
  `"whole_units"`, `"confounded"` or `"units"`, see Details),
  `groups_per_combination`, `units_per_group`, `cells_per_group` and
  `random` (whether
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
  adds a random intercept for it).

- `set_cells`: for each set, the cells it uses, with the `context` in
  which each one was observed (used by
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)).

- `cell_use`: number of identifiable sets that use each cell.

- `checks`: data frame with `check`, `severity` (`"problem"`,
  `"warning"` or `"info"`) and `message`.

## Details

A contrast set is the group of cells needed to estimate how a
combination of stressors deviates from a null model. For a combination
of the stressors in S at given levels (all other stressors at control):

- the **net** deviation needs the control, each stressor of S alone and
  their combination. It compares the combined response with the
  combination of the single-stressor effects predicted by the null
  model.

- the **pure** interaction of order \|S\| needs all 2^\|S\| combinations
  of S. It is what remains once all lower-order interactions are
  removed.

With two stressors both are the same. With three or more, the net
deviation of ABC includes the pairwise interactions AB, AC and BC, so it
is not a pure third-order interaction; designs that lack the pairwise
cells cannot separate them.

Cells shared by several sets (typically the control in designs with
several doses, or the single-stressor cells of a factorial with several
levels) make the estimates of those sets correlated.

**Rows shared by several `by` groups.** When a control (or another cell)
is common to several groups, copying it into each group would count its
units several times. With `shared`, those rows form their own group
(`context` such as `"regime=both"`), and each of the other groups uses
its cells when it lacks them: a set of `"regime=pulse"` takes its
control from `"regime=both"`. The shared rows supply cells, not sets:
sets are formed from the combinations observed within each group. A
design with the order of application as a `by` column (nutrients then
salinity, salinity then nutrients) and the control and single stressors
shared estimates the interaction for each order.

**Groups of units.** A stressor that takes a single level within every
group of a `groups` column was applied to whole groups (a whole-plot
factor). Each column gets one of four structures:

- `"blocks"`: no stressor is applied to whole groups; treatment
  combinations are compared within groups (complete or incomplete
  blocks).

- `"whole_units"`: some stressors are applied to whole groups, with at
  least two groups per combination of their levels. For effects that
  involve them the replication is the number of groups, not of units,
  and the random intercept for the groups is needed for valid inference.

- `"confounded"`: some combination of the stressors applied to whole
  groups was given to a single group, so effects that involve them
  cannot be separated from differences among groups (for example, one
  incubator per temperature, or one nitrate level per block). The sets
  affected are flagged in `sets$confounded`; a random intercept does not
  remove the confounding and is not added when it cannot be estimated.

- `"units"`: each group holds a single unit, so the column repeats
  `unit` and is not used.

## Assumptions

- The `control` level of each stressor represents its absence.

- With `background = "control"`, interactions are estimated with the
  stressors outside the set at control (interactions "in the absence of
  the other stressors"); an interaction can differ in the presence of
  another stressor, which `background = "all"` estimates.

- Each experimental unit receives a single treatment combination for the
  whole experiment. Treatments that change within a unit (sequential
  exposures) are reported as a problem, except stressors that start
  later than others, declared in `baseline` by stressor.

- Units and groups not given in `unit` and `groups` are assumed
  independent: the design cannot reveal a shared incubator or water
  supply that is not recorded in the data.

- Levels are categories. A gradient with many levels and few replicates
  is analysed level by level, not as a continuous predictor.

## Missing values and zeros

Missing values in the stressor, `unit`, `time`, `by` or `groups` columns
are an error: no treatment combination is guessed. The response is not
used here.

## References

Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B.
(2025). Hypothesis-driven research on multiple stressors: An analytical
framework for stressor interactions. *Ecology and Evolution*, 15,
e71959. [doi:10.1002/ece3.71959](https://doi.org/10.1002/ece3.71959)

## Examples

``` r
# Two stressors, two levels each, four replicate jars
jars <- expand.grid(
  salinity = c("ambient", "high"),
  warming = c("no", "yes"),
  replicate = 1:4,
  stringsAsFactors = FALSE
)
jars$jar <- seq_len(nrow(jars))
design <- describe_design(
  jars,
  stressors = c("salinity", "warming"),
  control = c(salinity = "ambient", warming = "no"),
  unit = "jar"
)
design
#> <must_design> 16 rows, 2 stressors, 4 observed treatment combinations
#> Stressors (control): salinity (ambient), warming (no)
#> Unit: jar | Time: none | By: none | Groups: none
#> Contrast sets: 1 candidate(s), 1 identifiable
#> Checks:
#>   info     Experimental unit: 'jar' (16 units, each in a single treatment
#>            combination).
#>   info     Units per treatment combination: 4.
#>   info     Identifiable contrast sets (net deviation): 1 with 2 stressors.
design$sets
#>   set_id context order          stressors                      levels
#> 1    S01             2 salinity x warming salinity=high | warming=yes
#>   background net_identifiable pure_identifiable missing_for_net
#> 1                        TRUE                NA                
#>   missing_for_pure min_units shared confounded
#> 1                          4  FALSE      FALSE

# A ray design: the mixture only at matching doses
ray <- data.frame(
  a = c(0, 2, 20, 0, 0, 2, 20),
  b = c(0, 0, 0, 2, 20, 2, 20)
)
ray <- ray[rep(seq_len(nrow(ray)), each = 3), ]
ray_design <- describe_design(ray, c("a", "b"), control = c(a = "0", b = "0"))
ray_design$sets
#>   set_id context order stressors      levels background net_identifiable
#> 1    S01             2     a x b   a=2 | b=2                        TRUE
#> 2    S02             2     a x b a=20 | b=20                        TRUE
#>   pure_identifiable missing_for_net missing_for_pure min_units shared
#> 1                NA                                          3   TRUE
#> 2                NA                                          3   TRUE
#>   confounded
#> 1      FALSE
#> 2      FALSE
ray_design$missing_cells
#>   context  a  b       cell n_stressed
#> 1          2 20 a=2 | b=20          2
#> 2         20  2 a=20 | b=2          2

# Temperature applied to whole rooms, one room per temperature: its
# effects cannot be separated from differences between rooms
rooms <- expand.grid(temperature = c("15", "20"), nutrients = c("no", "yes"),
                     jar = 1:3, stringsAsFactors = FALSE)
rooms$room <- paste0("R", rooms$temperature)
rooms$id <- seq_len(nrow(rooms))
room_design <- describe_design(rooms, c("temperature", "nutrients"),
  control = c(temperature = "15", nutrients = "no"), unit = "id",
  groups = "room")
room_design$groups
#>   group n_groups     between  structure groups_per_combination units_per_group
#> 1  room        2 temperature confounded                      1               6
#>   cells_per_group random
#> 1               2  FALSE
room_design$sets$confounded
#> [1] TRUE
```
