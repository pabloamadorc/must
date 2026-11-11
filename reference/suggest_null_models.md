# Suggest null models from what is known about the stressors

Applies the guidelines of Schäfer and Piggott (2018, Fig. 3) for the
choice of a null model, from what is known about the stressors and their
effects: whether the stressors act similarly, whether one of them acts
specifically on the physiology of the organisms, whether the effect is
mortality, whether the stressor-effect relationship or the distribution
of sensitivities is known, how the sensitivities of individuals to the
stressors are correlated, and how large the effects of the stressors
alone are relative to the effect limit (the largest possible effect,
such as 100% mortality or the loss of the whole response). A criterion
left unknown (`NA`) is tried with each of its values, and every null
model that one of them leads to is kept: the range of predictions that
the authors recommend when the information is lacking, which
[`compare_null_models()`](https://pabloamadorc.github.io/must/reference/compare_null_models.md)
then compares.

The result can be given as `null_models` to
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md),
which estimates the null models it suggests and keeps the criteria and
the reasons in `settings$null_model_choice`.

## Usage

``` r
suggest_null_models(
  x = NULL,
  similar_action = NA,
  specific_action = NA,
  mortality = NA,
  stressor_effect = NA,
  sensitivity_distribution = NA,
  correlation = NA,
  effects = NULL,
  near_limit = 0.9,
  level = 0.95,
  response = NULL,
  sensitivity = NULL,
  design = NULL,
  data = NULL,
  ...
)
```

## Arguments

- x:

  A model fitted with
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md),
  from which the effects of the stressors alone are estimated; a model
  fitted outside must, with its `design`, as in
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md);
  the output of
  [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md),
  to check which null models the design allows; or `NULL` (the default),
  to apply only the criteria given.

- similar_action:

  Do the stressors act similarly, with the same mode of action? Schäfer
  and Piggott (2018) distinguish stressors that act from outside on the
  whole organism (predation, floods), relatively nonspecifically on its
  physiology (competition, salinity, heat) and specifically on
  physiological processes (many toxicants). `TRUE`, `FALSE` or `NA`
  (unknown, the default).

- specific_action:

  Does one of the stressors act specifically on physiological processes,
  as many toxicants do? `TRUE`, `FALSE` or `NA`.

- mortality:

  Is the effect mortality? `TRUE`, `FALSE` or `NA`.

- stressor_effect:

  Is the stressor-effect relationship of the stressors known, from
  earlier studies or from enough levels of each stressor in the
  experiment (see Details)? `TRUE`, `FALSE` or `NA`.

- sensitivity_distribution:

  Is the distribution of the sensitivities (stress capacities) of the
  individuals known? `TRUE`, `FALSE` or `NA`.

- correlation:

  How the sensitivities of individuals to the stressors are correlated:
  `"positive"` (highly positively: the individuals sensitive to one
  stressor are sensitive to the other), `"negative"` (highly negatively:
  those sensitive to one are tolerant to the other), `"weak"` (neither)
  or `NA` (unknown, the default).

- effects:

  Effects of the stressors alone as proportions of the effect limit: a
  named numeric vector with one value per stressor, such as
  `c(salinity = 0.3, warming = 0.5)`. Values cannot exceed 1 (the effect
  limit); a negative value is a stressor that increases the response,
  and `NA` an unknown effect. `NULL` (the default) estimates them from
  `x` when it is a fitted model and leaves them unknown otherwise.

- near_limit:

  Proportion of the effect limit from which the effect of a stressor
  counts as close to it. Default 0.9 (see Details).

- level:

  Confidence level of the intervals of the effects estimated from `x`
  (default 0.95). A criterion that depends on the effects is unknown
  when an interval spans its threshold; with `NULL`, the point estimates
  decide.

- response:

  For a model fitted with
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md),
  the names of the responses whose effects are estimated; `NULL` (the
  default) for all. The other criteria apply to all of them: when they
  differ (mortality for one response, growth for another), call the
  function once per response.

- sensitivity:

  As in
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md):
  for a model fitted with
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md),
  the names of the responses for which an increase means a stronger
  effect (mortality, damage), whose proportions are then complemented
  (mortality becomes survival); for a model fitted outside must, `TRUE`
  or `FALSE`. Only used when the effects are estimated.

- design:

  For a model fitted outside must, the output of
  [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md)
  for its data.

- data, ...:

  For a model fitted outside must, as in
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md):
  the data it was fitted to, when emmeans cannot recover them, and
  further arguments to
  [`emmeans::ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.html).

## Value

An object of class `must_null_models`, a list with:

- `null_models`: the codes of the null models kept, for
  `estimate_interactions(null_models = )`.

- `models`: one row per null model of the guidelines (`"AD"`, `"MU"`,
  `"DO"`, `"SAM"`, `"CA"`) with its `status` (`"suggested"` when it is
  the only one reached, `"possible"`, or `"excluded"`), `reached` (in
  how many of the `evaluations`, the contrast sets and times of a fitted
  model, it is reached), `reason` (the steps that lead to it, or why it
  is not reached), `design` (whether the design allows it: `"yes"`,
  `"no"`, `"some sets"`, or `NA` without a design), `design_note` and
  `note`.

- `criteria`: the value of each criterion (`"yes"`, `"no"`, a
  correlation, `"below"`, `"above"` or `"equal"` for the sum of the
  effects, `"unknown"`, `"uncertain"` when an interval spans the
  threshold, or `"varies"` among sets and times) and its `source`
  (`"given"`, `"estimated"` from the fitted model, or `"unknown"`), with
  counts in `detail` when it varies.

- `decisive`: for each unknown criterion whose value would change the
  models kept, the models kept with each of its values.

- `effects`: one row per stressor (and response, contrast set and time,
  for a fitted model) with the effect (`estimate`), its interval
  (`lower`, `upper`, `level`, `df`), its `direction`, `status` and
  `reason`; for a fitted model, `scale` says whether the effect is on
  the response as analysed (`"response"`) or on its complement
  (`"complement"`, for the proportions described above, which have rows
  on both scales). `NULL` when no effects were given or estimated.

- `sets`: for a fitted model, one row per response, contrast set and
  time with the scale to which the guidelines were applied (`scale`),
  the `type` of the effects (on the response as analysed), the sum of
  the effects and its interval, whether one effect is close to the limit
  (`near_limit`), how the sum compares with the effect limit
  (`sum_vs_limit`: `"below"`, `"above"`, `"equal"`, or `NA` if unknown),
  whether a stressor alone increases the response (`increase`) and two
  have opposite effects (`opposing`), `similar_predictions` and the null
  models reached (`null_models`); `NULL` otherwise.

- `settings`: the thresholds, the source of the effects, the responses
  whose complement the guidelines were applied to (`complement`) and the
  `sensitivity` with which to estimate the interactions
  (`suggested_sensitivity`).

## Details

**The guidelines.** The steps of Fig. 3 of Schäfer and Piggott (2018),
in the order in which they are checked; the first whose conditions hold
gives the null model:

1.  The stressors act similarly and their stressor-effect relationship
    is known: concentration addition (`"CA"`). In must, also with
    effects of the stressors alone in the same direction (see below).

2.  One stressor acts specifically on physiology, the effect is
    mortality, and the stressor-effect relationship or the sensitivity
    distribution is known: the stress addition model (`"SAM"`). In must,
    also with no stressor alone that increases the response.

3.  The sensitivities are highly positively correlated, or the effect of
    one stressor is close to the effect limit: dominance (`"DO"`).

4.  The sensitivities are highly negatively correlated and the sum of
    the effects is below the effect limit: simple addition (`"AD"`).

5.  The effect is mortality: the multiplicative model (`"MU"`).

6.  The sum of the effects exceeds the effect limit: the multiplicative
    model (`"MU"`).

7.  Otherwise: simple addition (`"AD"`).

The multiplicative model is the effect addition of ecotoxicology (Bliss
independence): `"BL"` in
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
gives the same prediction, with the deviation relative to the reference
cell instead of on the log scale. When the mode of action or the
correlation of the sensitivities is unknown, the authors suggest
computing every null model that matches the remaining criteria; each
unknown criterion here takes all its values, and `decisive` lists those
whose value would change the models kept. They also note that the choice
matters less when the effect of one stressor is close to the effect
limit or about ten times larger than that of the other, since the null
models then predict similar joint effects (except perhaps the stress
addition model and concentration addition); `similar_predictions` in
`sets` flags those cases.

The guidelines were written for two stressors and for effects on
individuals (or populations, with caveats). For three or more stressors,
the same criteria are applied to all the stressors of a contrast set
(the effect of any of them close to the limit, the sum of all their
effects), and a larger effect counts as about ten times the others when
it is at least ten times their sum; dominance is defined in must for two
stressors only. Schäfer and Piggott (2018) give no value for an effect
close to the effect limit: `near_limit = 0.9` is a choice of this
package. They suggest at least five levels of a stressor, on an interval
or ratio scale, to establish its stressor-effect relationship in a
factorial experiment; with a design, `models` says how many levels of
each stressor were tested alone. Concentration addition in must fits its
curves to two or more doses, and the stress addition model assumes a
beta distribution of the stress capacities (`sam_shape` in
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)).

**Effects estimated from a fitted model.** For each contrast set and
sampling time, the effect of each stressor alone is 1 - y_i / y_0, from
the means of the model at the reference cell of the set (y_0) and at the
cell with that stressor alone (y_i): the reduction of the response as a
share of the reference, whose limit (the loss of the whole response)
is 1. The cells with two or more stressors are not used, so that the
choice does not depend on the interactions it serves to test. The
interval of each effect comes from the delta method on the log scale of
y_i / y_0, and that of their sum, the number of stressors minus the sum
of the ratios y_i / y_0, on the log scale of that sum, with the
covariance and the degrees of freedom that
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
uses by default (REML and the small-sample correction). One effect is
close to the limit if the lower bound of its interval reaches
`near_limit`, and none is if every upper bound is below it; otherwise
the criterion is unknown. The sum of the effects is below the effect
limit if the upper bound of its interval is, and above it if the lower
bound is. A cell whose mean is fixed at the boundary (only zeros, or
only successes) has no variance: the effect of a stressor alone in such
a cell gets no interval, and so can take any value up to the effect
limit, while the uncertainty of a reference cell at the boundary is left
out. Effects are not defined when the mean of the reference cell is
zero, when the response takes negative values, or when an increase means
a stronger effect in a response that is not a proportion
(`sensitivity`): such a response has no effect limit. The null models
kept are those reached in any set and time; `sets` gives them for each
one.

**Stressors that increase the response.** A stressor, in the sense of
Orr et al. (2020) and Madge Pimentel et al. (2025), is any driver of
change, whatever its direction: in the 170 studies reanalysed by Piggott
et al. (2015), 29% of the pairs of stressors had opposing effects and
22% two positive ones. Each effect of a stressor alone has a
`direction`: `"decrease"` or `"increase"` of the response when its
interval excludes zero (by its sign for effects given or point
estimates), `"none"` when it includes zero, and `"unknown"` without an
interval. Each set gets the interaction `type` of Piggott et al. (2015,
Table 1) from them: `"double negative"` (both decrease the response),
`"double positive"`, `"opposing"`, `"negative neutral"`,
`"positive neutral"` or `"double neutral"` (`"all"` instead of
`"double"` for three or more stressors). The guidelines take the effects
with their sign. The sum of the effects exceeds the effect limit exactly
when simple addition predicts a response below zero, y_0 (1 - sum of the
effects) \< 0, so an increase offsets the decreases; only a decrease can
be close to the effect limit. Concentration addition treats the
stressors as dilutions of one another, which is not possible with
opposite effects, and the stress addition model assumes stressors that
reduce the response (Liess et al. 2016): steps 1 and 2 also need, in
must, effects in the same direction and no stressor alone that increases
the response (`opposing` and `increase` in `sets`), as
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
already requires for those models. For a proportion fitted with `trials`
that the stressors alone increase and never decrease, the guidelines are
applied to its complement (1 - p, mortality for survival), where the
increases are reductions towards the limit 0: the mirror of the effect
limit for the upper limit 1, which simple addition can also cross, and
the scale on which the multiplicative model and dominance describe
independent and correlated effects (Madge Pimentel et al. 2025 discuss
the asymmetry between survival and mortality). The interactions must
then be estimated on that scale: `settings$suggested_sensitivity` gives
the `sensitivity` for
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md),
which checks it. For responses of whole communities, such as a total
density, some taxa often increase when others are harmed, through
interactions among them; the guidelines were written for individuals and
populations, and Schäfer and Piggott (2018) point instead to the
compositional null model of Thompson et al. (2018), which predicts the
community response from the null models of each taxon (not implemented
in must).

## References

Schäfer, R. B., & Piggott, J. J. (2018). Advancing understanding and
prediction in multiple stressor research through a mechanistic basis for
null models. *Global Change Biology*, 24(5), 1817–1826.
[doi:10.1111/gcb.14073](https://doi.org/10.1111/gcb.14073)

Liess, M., Foit, K., Knillmann, S., Schäfer, R. B., & Liess, H.-D.
(2016). Predicting the synergy of multiple stress effects. *Scientific
Reports*, 6, 32965.
[doi:10.1038/srep32965](https://doi.org/10.1038/srep32965)

Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B.
(2025). Hypothesis-driven research on multiple stressors: An analytical
framework for stressor interactions. *Ecology and Evolution*, 15,
e71959. [doi:10.1002/ece3.71959](https://doi.org/10.1002/ece3.71959)

Orr, J. A., Vinebrooke, R. D., Jackson, M. C., Kroeker, K. J., Kordas,
R. L., Mantyka-Pringle, C., Van den Brink, P. J., De Laender, F., Stoks,
R., Holmstrup, M., Matthaei, C. D., Monk, W. A., Penk, M. R., Leuzinger,
S., Schäfer, R. B., & Piggott, J. J. (2020). Towards a unified study of
multiple stressors: Divisions and common goals across research
disciplines. *Proceedings of the Royal Society B*, 287, 20200421.
[doi:10.1098/rspb.2020.0421](https://doi.org/10.1098/rspb.2020.0421)

Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015).
Reconceptualizing synergism and antagonism among multiple stressors.
*Ecology and Evolution*, 5(7), 1538–1547.
[doi:10.1002/ece3.1465](https://doi.org/10.1002/ece3.1465)

Thompson, P. L., MacLennan, M. M., & Vinebrooke, R. D. (2018). An
improved null model for assessing the net effects of multiple stressors
on communities. *Global Change Biology*, 24(1), 517–525.
[doi:10.1111/gcb.13852](https://doi.org/10.1111/gcb.13852)

## See also

[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md),
which takes the result as `null_models`, and
[`compare_null_models()`](https://pabloamadorc.github.io/must/reference/compare_null_models.md),
which compares the null models kept.

## Examples

``` r
# Two stressors with different modes of action, a nonlethal effect and
# effects of 30% and 50% of the effect limit; the correlation of the
# sensitivities is unknown
suggest_null_models(similar_action = FALSE, mortality = FALSE,
  stressor_effect = FALSE, effects = c(salinity = 0.3, warming = 0.5))
#> <must_null_models> after the guidelines of Schäfer and Piggott (2018), Fig. 3
#> Criteria:
#>   The stressors act similarly (same mode of action)                       no
#>   One stressor acts specifically on physiology                            unknown
#>   The effect is mortality                                                 no
#>   The stressor-effect relationship is known                               no
#>   The sensitivity distribution is known                                   unknown
#>   Correlation of the sensitivities                                        unknown
#>   The effect of one stressor is close to the effect limit (at least 90%)  no
#>   The sum of the effects compared with the effect limit                   below
#>   A stressor alone increases the response                                 no
#>   Stressors alone with opposite effects                                   no
#> Effects of the stressors alone (given): salinity = 0.3 (decrease),
#>   warming = 0.5 (decrease); type: double negative.
#> Null models:
#>   AD (simple addition): possible. Step 4: the sensitivities are highly
#>     negatively correlated and the sum of the effects is below the
#>     effect limit; step 7: none of the earlier steps applies.
#>   MU (multiplicative): excluded. The multiplicative model needs
#>     mortality as the effect (step 5) or the sum of the effects above
#>     the effect limit (step 6), but the effect is not mortality and the
#>     sum of the effects does not exceed the effect limit.
#>   DO (dominance): possible. Step 3: the sensitivities are highly
#>     positively correlated, or the effect of one stressor is close to
#>     the effect limit.
#>   SAM (stress addition): excluded. The stress addition model needs one
#>     stressor that acts specifically on physiology, mortality as the
#>     effect, stressors that do not increase the response and a known
#>     stressor-effect relationship or sensitivity distribution (step 2),
#>     but the effect is not mortality.
#>   CA (concentration addition): excluded. Concentration addition needs
#>     stressors that act similarly, with effects in the same direction,
#>     and a known stressor-effect relationship (step 1), but the
#>     stressors do not act similarly and their stressor-effect
#>     relationship is not known.
#> Unknown criteria that would change the choice:
#>   Correlation of the sensitivities (`correlation`): positive -> DO;
#>     negative -> AD; weak -> AD
#> For estimate_interactions(): null_models = c("AD", "DO")
# The same, with an effect close to the effect limit
suggest_null_models(similar_action = FALSE, mortality = FALSE,
  stressor_effect = FALSE, effects = c(salinity = 0.3, warming = 0.95))
#> <must_null_models> after the guidelines of Schäfer and Piggott (2018), Fig. 3
#> Criteria:
#>   The stressors act similarly (same mode of action)                       no
#>   One stressor acts specifically on physiology                            unknown
#>   The effect is mortality                                                 no
#>   The stressor-effect relationship is known                               no
#>   The sensitivity distribution is known                                   unknown
#>   Correlation of the sensitivities                                        unknown
#>   The effect of one stressor is close to the effect limit (at least 90%)  yes
#>   The sum of the effects compared with the effect limit                   above
#>   A stressor alone increases the response                                 no
#>   Stressors alone with opposite effects                                   no
#> Effects of the stressors alone (given): salinity = 0.30 (decrease),
#>   warming = 0.95 (decrease); type: double negative.
#> Null models:
#>   AD (simple addition): excluded. Not reached: an earlier step always
#>     applies (step 3, dominance).
#>   MU (multiplicative): excluded. Not reached: an earlier step always
#>     applies (step 3, dominance).
#>   DO (dominance): suggested. Step 3: the sensitivities are highly
#>     positively correlated, or the effect of one stressor is close to
#>     the effect limit.
#>   SAM (stress addition): excluded. The stress addition model needs one
#>     stressor that acts specifically on physiology, mortality as the
#>     effect, stressors that do not increase the response and a known
#>     stressor-effect relationship or sensitivity distribution (step 2),
#>     but the effect is not mortality.
#>   CA (concentration addition): excluded. Concentration addition needs
#>     stressors that act similarly, with effects in the same direction,
#>     and a known stressor-effect relationship (step 1), but the
#>     stressors do not act similarly and their stressor-effect
#>     relationship is not known.
#> For estimate_interactions(): null_models = "DO"

# Effects estimated from a fitted model
jars <- expand.grid(
  salinity = c("ambient", "high"),
  warming = c("no", "yes"),
  replicate = 1:5,
  stringsAsFactors = FALSE
)
jars$jar <- seq_len(nrow(jars))
mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
  ifelse(jars$warming == "yes", 0.4, 1)
jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
design <- describe_design(jars, c("salinity", "warming"),
  control = c(salinity = "ambient", warming = "no"), unit = "jar")
fits <- fit_responses(jars, "algae", design, residual_checks = FALSE,
  seed = 1, quiet = TRUE)
choice <- suggest_null_models(fits, similar_action = FALSE, mortality = FALSE,
  stressor_effect = FALSE)
choice
#> <must_null_models> after the guidelines of Schäfer and Piggott (2018), Fig. 3
#> Criteria:
#>   The stressors act similarly (same mode of action)                       no
#>   One stressor acts specifically on physiology                            unknown
#>   The effect is mortality                                                 no
#>   The stressor-effect relationship is known                               no
#>   The sensitivity distribution is known                                   unknown
#>   Correlation of the sensitivities                                        unknown
#>   The effect of one stressor is close to the effect limit (at least 90%)  no (estimated)
#>   The sum of the effects compared with the effect limit                   uncertain (estimated; an interval spans the threshold)
#>   A stressor alone increases the response                                 no (estimated)
#>   Stressors alone with opposite effects                                   no (estimated)
#> Effects of the stressors alone: estimated in 1 of 1 contrast sets and
#>   times, with 95% intervals.
#>   Directions of the effects (types of Piggott et al. 2015): double
#>   negative in 1.
#> Null models:
#>   AD (simple addition): possible. Step 4: the sensitivities are highly
#>     negatively correlated and the sum of the effects is below the
#>     effect limit; step 7: none of the earlier steps applies.
#>   MU (multiplicative): possible. Step 6: the sum of the effects exceeds
#>     the effect limit. The effect addition of ecotoxicology: "BL" (Bliss
#>     independence) gives the same prediction, with the deviation
#>     relative to the reference cell.
#>   DO (dominance): possible. Step 3: the sensitivities are highly
#>     positively correlated, or the effect of one stressor is close to
#>     the effect limit.
#>   SAM (stress addition): excluded. The stress addition model needs one
#>     stressor that acts specifically on physiology, mortality as the
#>     effect, stressors that do not increase the response and a known
#>     stressor-effect relationship or sensitivity distribution (step 2),
#>     but the effect is not mortality.
#>   CA (concentration addition): excluded. Concentration addition needs
#>     stressors that act similarly, with effects in the same direction,
#>     and a known stressor-effect relationship (step 1), but the
#>     stressors do not act similarly and their stressor-effect
#>     relationship is not known.
#> Unknown criteria that would change the choice:
#>   Correlation of the sensitivities (`correlation`): positive -> DO;
#>     negative -> AD, MU; weak -> AD, MU
#>   Sum of the effects vs the limit (`sum`): below -> AD, DO; above ->
#>     MU, DO
#> For estimate_interactions(): null_models = c("AD", "MU", "DO")
choice$effects[, c("stressor", "estimate", "lower", "upper")]
#>   stressor  estimate     lower     upper
#> 1 salinity 0.4974874 0.3610238 0.6048071
#> 2  warming 0.6030151 0.4848120 0.6940980
# The null models kept, with their justification in the settings
est <- estimate_interactions(fits, null_models = choice, n_sim = 2000, seed = 1)
est$settings$null_model_choice$criteria
#>                   criterion
#> 1            similar_action
#> 2           specific_action
#> 3                 mortality
#> 4           stressor_effect
#> 5  sensitivity_distribution
#> 6               correlation
#> 7                near_limit
#> 8                       sum
#> 9                  increase
#> 10                 opposing
#>                                                                  question
#> 1                       The stressors act similarly (same mode of action)
#> 2                            One stressor acts specifically on physiology
#> 3                                                 The effect is mortality
#> 4                               The stressor-effect relationship is known
#> 5                                   The sensitivity distribution is known
#> 6                                        Correlation of the sensitivities
#> 7  The effect of one stressor is close to the effect limit (at least 90%)
#> 8                   The sum of the effects compared with the effect limit
#> 9                                 A stressor alone increases the response
#> 10                                  Stressors alone with opposite effects
#>        value    source detail
#> 1         no     given       
#> 2    unknown   unknown       
#> 3         no     given       
#> 4         no     given       
#> 5    unknown   unknown       
#> 6    unknown   unknown       
#> 7         no estimated       
#> 8  uncertain estimated       
#> 9         no estimated       
#> 10        no estimated       
```
