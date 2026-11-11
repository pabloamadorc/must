# Agreement between null models on each interaction

Tells, for each interaction estimated under two or more null models,
whether all of them flag a deviation with the same sign, only some,
none, or whether they disagree on its direction. Whether a combined
effect is called an interaction depends on the null model (Madge
Pimentel et al. 2025): an interaction that every plausible null model
flags does not depend on that choice, whereas one flagged by only some
of them needs a mechanistic hypothesis to be interpreted (Schäfer &
Piggott 2018).

## Usage

``` r
compare_null_models(x, null_models = NULL, alpha = 0.05)
```

## Arguments

- x:

  A `must_interactions` object from
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md),
  estimated under two or more null models.

- null_models:

  Null models to compare, among those of `x`. Default `NULL`: all of
  them.

- alpha:

  Significance level applied to the adjusted p-values. Default 0.05.

## Value

A data frame with one row per response, `by` group, sampling time,
contrast set and estimand, and columns:

- `response`, `context` and the `by` columns, `time`, `set_id`, `order`,
  `stressors`, `levels`, `background` and `estimand`, as in
  `x$estimates`.

- One column per null model with the deviation it flags: `"+"` or `"-"`
  (adjusted p-value below `alpha`, with the sign of the estimate), `"0"`
  (tested, no deviation detected) or `NA` (not tested: not estimable,
  not applicable, not defined or not identifiable).

- `n_tested` and `n_deviating`: null models tested and null models that
  flag a deviation.

- `agreement`: `"robust"` (every null model tested flags a deviation,
  all with the same sign), `"model_dependent"` (some but not all of
  them, with the same sign), `"conflicting"` (deviations of opposite
  sign), `"none"` (no deviation detected under any null model) or `NA`
  (fewer than two null models tested).

- `direction`: `"positive"` or `"negative"` for robust and
  model-dependent deviations; `NA` otherwise.

The attributes `alpha` and `null_models` keep the significance level and
the null models compared, for
[`report_interactions()`](https://pabloamadorc.github.io/must/reference/report_interactions.md).

## Details

Deviations are decided as in
[`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)
and [`summary()`](https://rdrr.io/r/base/summary.html): the adjusted
p-value of each null model (adjusted within the null model, see
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md))
compared with `alpha`. With an equivalence `margin`, `"0"` covers both
equivalent and inconclusive interactions; the `verdict` column of
[`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)
tells them apart. The sign is that of the estimate on the scale of each
null model, negative when the combination falls below what the null
model predicts; whether that is a synergism or an antagonism depends on
the direction of the single effects (see
[`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)).
Net deviations and pure interactions are compared separately; null
models that do not define an estimand (the pure interaction under DO,
SAM and CA) are not tested for it.

"Robust" is relative to the null models compared, so choose those
plausible for the stressors at hand. In simulated experiments generated
without interaction by AD, MU, DO, SAM or CA (see "Default null model"
in
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)),
requiring the six null models to agree flagged false interactions in
0.6% of the experiments on average (2.3% at most), but detected only 15%
of the interactions that no mechanism explains. It is too strict a rule
for detecting interactions, but it identifies those that do not depend
on the null model. Model-dependent and conflicting interactions are
those where the mechanistic hypothesis decides; a conflict is common
when the data are compared with stress addition, which predicts larger
combined effects than the other null models.

## References

Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B.
(2025). Hypothesis-driven research on multiple stressors: An analytical
framework for stressor interactions. *Ecology and Evolution*, 15,
e71959. [doi:10.1002/ece3.71959](https://doi.org/10.1002/ece3.71959)

Schäfer, R. B., & Piggott, J. J. (2018). Advancing understanding and
prediction in multiple stressor research through a mechanistic basis for
null models. *Global Change Biology*, 24(5), 1817–1826.
[doi:10.1111/gcb.14073](https://doi.org/10.1111/gcb.14073)

## Examples

``` r
jars <- expand.grid(
  salinity = c("ambient", "high"),
  warming = c("no", "yes"),
  replicate = 1:5,
  stringsAsFactors = FALSE
)
jars$jar <- seq_len(nrow(jars))
mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
  ifelse(jars$warming == "yes", 0.5, 1)
jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
design <- describe_design(jars, c("salinity", "warming"),
  control = c(salinity = "ambient", warming = "no"), unit = "jar")
fits <- fit_responses(jars, "algae", design, residual_checks = FALSE,
  seed = 1, quiet = TRUE)
est <- estimate_interactions(fits, null_models = c("AD", "MU", "DO"), n_sim = 2000,
  seed = 1)
compare_null_models(est)
#>   response context time set_id order          stressors
#> 1    algae         <NA>    S01     2 salinity x warming
#>                        levels background estimand AD MU DO n_tested n_deviating
#> 1 salinity=high | warming=yes                 net  +  0  -        3           2
#>     agreement direction
#> 1 conflicting      <NA>
```
