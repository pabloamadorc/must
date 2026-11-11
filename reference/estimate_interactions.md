# Estimate single, combined and interaction effects under null models

Computes, for every identifiable contrast set of the design (see
[`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md)),
the effect of each stressor alone, the effect of their combination and
its deviation from what each null model predicts. This is
post-estimation inference: one model per response, one or more null
models tested from its predictions with the same criteria (Madge
Pimentel et al. 2025). The model is either the one selected by
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
or one fitted outside must, such as
[`stats::lm()`](https://rdrr.io/r/stats/lm.html),
[`stats::glm()`](https://rdrr.io/r/stats/glm.html),
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html),
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) or
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)
(see "Models fitted outside must"). By default the null model is simple
addition (see "Default null model").

## Usage

``` r
estimate_interactions(fit, ...)

# S3 method for class 'must_fit'
estimate_interactions(
  fit,
  null_models = "AD",
  inference = "simulation",
  correction = "small_sample",
  n_sim = 10000L,
  level = 0.95,
  p_adjust = "BH",
  overall = TRUE,
  sensitivity = NULL,
  seed = NULL,
  sam_shape = 3.2,
  margin = NULL,
  ...
)

# Default S3 method
estimate_interactions(
  fit,
  design,
  null_models = "AD",
  inference = "simulation",
  correction = "small_sample",
  n_sim = 10000L,
  level = 0.95,
  p_adjust = "BH",
  overall = TRUE,
  sensitivity = FALSE,
  seed = NULL,
  sam_shape = 3.2,
  margin = NULL,
  df = NULL,
  data = NULL,
  ...
)
```

## Arguments

- fit:

  A `must_fit` object from
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md),
  or a model fitted outside must whose predictors include the stressor
  columns of `design` (any model supported by the emmeans package).

- ...:

  For a model fitted outside must: further arguments passed to
  [`emmeans::ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.html),
  such as `at` (values of other numeric predictors), `offset`, `lmer.df`
  or `component`. Not used otherwise.

- null_models:

  Null models to evaluate: any of `"AD"` (simple addition), `"MU"`
  (multiplicative), `"DO"` (dominance), `"BL"` (Bliss independence),
  `"SAM"` (stress addition model) and `"CA"` (concentration addition, or
  Loewe additivity). Default `"AD"`, the null model that flagged the
  fewest false interactions across mechanisms in simulations (see
  "Default null model"); `c("AD", "MU", "DO")` gives the three null
  models of Madge Pimentel et al. (2025). It can also be the output of
  [`suggest_null_models()`](https://pabloamadorc.github.io/must/reference/suggest_null_models.md):
  its null models are estimated, and its criteria and reasons are kept
  in `settings$null_model_choice`.

- inference:

  `"simulation"` (default; Krinsky–Robb: draws of the fixed effects from
  their estimated distribution) or `"delta"` (delta method; not for
  `"DO"`).

- correction:

  `"small_sample"` (default) or `"none"`. For a `must_fit` with
  `"small_sample"`, the covariance of the fixed effects comes from the
  REML refit of the selected model and the sampling distribution is a t
  with the degrees of freedom of the experimental units (see Details);
  for a model fitted outside must, a t with the degrees of freedom given
  by emmeans (residual ones for a glm with an estimated dispersion) or
  by `df`. `"none"` uses the normal distribution (and, for a `must_fit`,
  the maximum-likelihood covariance), as in Madge Pimentel et al.
  (2025); with few replicates its tests reject a true null model too
  often.

- n_sim:

  Number of draws for `inference = "simulation"`, and for the tests of
  combinations with only zeros with either inference. Default 10000.

- level:

  Confidence level of the intervals. Default 0.95.

- p_adjust:

  Method of [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
  used across contrast sets and sampling times (see Details). Default
  `"BH"`.

- overall:

  If `TRUE` (default) and the design has two or more sampling times
  after any `baseline`, adds for each set the average of each effect
  over those times (`time = "overall"`).

- sensitivity:

  Whether an increase in the response means sensitivity (mortality,
  damage) rather than tolerance. For a `must_fit`: `NULL` (default) or
  the names of such responses; for a model fitted outside must: `TRUE`
  or `FALSE` (default). When the response is a proportion (`trials` in
  [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md);
  a binomial-type family or a logit, probit, cloglog or cauchit link in
  a model fitted outside must), every null model is then applied to its
  complement (1 - p, e.g. survival instead of mortality); for other
  responses the MU, DO and BL rows are `"not_applicable"`, because those
  null models assume that an increase means tolerance.

- seed:

  `NULL` (default) or a whole number. The simulation runs with this seed
  through
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html),
  which restores the random number state of the session afterwards. With
  `NULL`, a seed is drawn from the session's random number stream;
  either way it is stored in `settings$seed`, so the estimates can be
  repeated exactly.

- sam_shape:

  Shape of the symmetric beta distribution of the stress capacity of
  individuals in the stress addition model (both of its parameters).
  Default 3.2, the value fitted by Liess et al. (2016).

- margin:

  `NULL` (default) or a positive number: the equivalence margin for the
  relative deviation of each interaction (`relative`, a share of the
  response of the reference cell; 0.1 means 10%). With a margin, each
  interaction is also tested for equivalence (`p_equivalence`), and
  [`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)
  and [`summary()`](https://rdrr.io/r/base/summary.html) give its
  verdict (see "Equivalence"). It has no default: what deviation is
  negligible is a scientific decision, to be declared before the
  analysis.

- design:

  For a model fitted outside must: the output of
  [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md)
  for the data the model was fitted to. Its stressors, control levels,
  `time`, `by`, `baseline`, `max_order` and `background` define the
  contrast sets.

- df:

  For a model fitted outside must: `NULL` (default), for the degrees of
  freedom that emmeans gives for each contrast, or a single positive
  number used for every contrast (for example, experimental units minus
  treatment combinations). Only with `correction = "small_sample"`.
  emmeans gives infinite degrees of freedom to a
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) whose dispersion is
  estimated (quasi-Poisson, quasi-binomial and the other families except
  Poisson and binomial); with `NULL`, its residual degrees of freedom
  are used instead, as
  [summary.glm()](https://rdrr.io/r/stats/summary.glm.html) does.
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) keeps the
  normal distribution, as in its summary.

- data:

  For a model fitted outside must: the data it was fitted to, when
  emmeans cannot recover them from the model (for example, a
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) whose formula
  transforms a column, as in `factor(day)`). Default `NULL`.

## Value

An object of class `must_interactions`, a list with `estimates`,
`design`, `settings` and `diagnostics`. `estimates` is a data frame with
one row per response, contrast set, sampling time, null model and
quantity, and columns:

- `response`, `context` and the `by` columns, `time` (`NA` without time,
  `"overall"` for averages over times), `set_id`, `order`, `stressors`,
  `levels`, `background` (levels of the stressors held outside the set):
  what is estimated.

- `null_model` (`"AD"`, `"MU"`, `"DO"`, `"BL"`, `"SAM"`, `"CA"`) and
  `scale` (`"response"`, `"log(response)"` or `"response / control"`):
  what "no interaction" means and on which scale.

- `component`: `"single"` (one stressor alone vs the reference cell; its
  name in `stressor`), `"combined"` (all stressors of the set vs the
  reference) or `"interaction"`.

- `estimand` for interactions: `"net"` (deviation of the combination
  from the null-model combination of the single effects) or `"pure"`
  (highest-order interaction once all lower-order ones are removed; only
  for three or more stressors; the emergent interaction under BL).

- `estimate`, `se` (standard error; with `inference = "simulation"`, the
  standard deviation of the simulated estimates before they are scaled
  to a t distribution, comparable with the delta method), `lower`,
  `upper`, `level`, `df` (degrees of freedom of the t distribution;
  `Inf` for the normal), `p_value`, `p_adjusted`, `mc_se` (Monte Carlo
  standard error of a simulated p-value).

- `observed` and `expected` (interactions only, response scale): the
  fitted mean of the combination and the value the null model predicts.

- `relative`, `relative_lower` and `relative_upper` (interactions only):
  the relative deviation (observed - expected) / reference, the same
  quantity for every null model, and its interval at `level` (see
  "Equivalence").

- `p_equivalence` and `p_equivalence_adjusted` (interactions, with a
  `margin` only): p-value of the test of equivalence, raw and adjusted
  like `p_adjusted`.

- `rescaled` and `rescaled_class` (BL interactions only): the
  interaction relative to its extremes and its descriptive class (see
  "Bliss independence").

- `inference`, `status` (`"ok"`, `"boundary"` (tested from the
  probability of a combination with only zeros, see Details),
  `"not_estimable"`, `"not_identifiable"`, `"not_defined"` or
  `"not_applicable"`), `reason`, `note` (`"expected_below_zero"` when
  simple addition predicts an impossible negative value for a
  non-negative response; `"no_point_expectation"` when dominance gives
  no single expected value; `"reference_not_positive"` when the relative
  deviation cannot be computed because the reference mean (or one of its
  draws) is not positive; `"equivalence_not_tested"` when the closest
  means at the margin could not be found (see "Equivalence");
  `"confounded_with_groups"` for the sets flagged in `sets$confounded`
  by
  [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md);
  `"complement"` for responses analysed as 1 - p; several notes are
  separated by `"; "`).

`settings` holds the arguments used (with the seed and whether it was
given), the versions of must and R (and emmeans, for a model fitted
outside must), the criteria behind the null models when they came from
[`suggest_null_models()`](https://pabloamadorc.github.io/must/reference/suggest_null_models.md)
and, for a `must_fit`, a summary of the models (`fit`): what
[`report_interactions()`](https://pabloamadorc.github.io/must/reference/report_interactions.md)
describes.

For a model fitted outside must, `diagnostics$means` has the mean of
every treatment combination (and time) that was used: its estimate and
standard error on the link scale, the mean on the response scale, the
degrees of freedom given by emmeans and, if it could not be used, why
(`"nonestimable"`, `"no_se"` or `"boundary"`).

## Details

With y0 the mean of the reference cell (the control, or the background
combination of the other stressors), yi that of stressor i alone and yS
that of the combination of the stressors in S:

**Two stressors** (Madge Pimentel et al. 2025, eqs. 1–3 and 7–9):
\$\$AD:\\ y\_{AB} - y_A - y_B + y_0\$\$ \$\$MU:\\ \ln y\_{AB} - \ln
y_A - \ln y_B + \ln y_0\$\$ \$\$DO:\\ \min(y\_{AB}, y_0) - \min(y_A,
y_B)\$\$

**k stressors.** The net deviation compares the combination with the
null-model combination of the single effects only: \\y_K - \sum_i y_i +
(k-1) y_0\\ (AD) or the same on the log scale (MU). The pure k-th order
interaction is the inclusion–exclusion sum \\\sum\_{S \subseteq K}
(-1)^{k-\|S\|} f(y_S)\\, with f the identity (AD) or the logarithm (MU).
The net deviation equals the sum of the pure interactions of all subsets
of K with two or more stressors: a non-additive net deviation is not, by
itself, evidence of a k-way interaction. The dominance model is defined
for two stressors only.

**Single and combined effects** are differences from the reference cell
on the response scale (AD) or log ratios (MU). For DO only the
interaction is reported; its single effects are those of AD.

**Inference.** `"simulation"` follows the Krinsky–Robb method used by
Madge Pimentel et al. (2025): percentile intervals and two-sided
p-values \\2 \min(P(x \le 0), P(x \ge 0))\\ of the simulated estimates,
bounded below by 1/`n_sim`. `"delta"` uses Wald tests; it is not valid
for DO, which is not differentiable where the ranking of the means
changes. Uncertainty from the selection of the error distribution is not
included.

**Small-sample correction.** With few replicates, a normal distribution
with the maximum-likelihood covariance understates the uncertainty: in
simulated 2 x 2 Gaussian experiments with three replicates the AD test
then rejected a true null model in about 13% of the experiments instead
of 5%. With `correction = "small_sample"` the covariance comes from the
REML refit and the draws come from a multivariate t (normal draws scaled
by \\\sqrt{\nu / \chi^2\_\nu}\\), whose degrees of freedom follow the
strata of the design: observations minus estimated means without
repeated measures; experimental units minus the treatment combinations
they are nested in with repeated measures or subsamples (combination and
time when each unit is sampled at a single time); minus one per extra
block of each `groups` column with block structure in the model; and,
for groups of units to which some stressors are applied whole (whole
plots), the groups minus the combinations of those stressors form their
own stratum, which bounds the degrees of freedom of the sets that
involve them and is subtracted from those of the other sets. With
`dispersion = "cell"` in
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
the degrees of freedom of each contrast follow the Welch–Satterthwaite
approximation (Satterthwaite 1946). Families without a dispersion
parameter (Poisson, binomial) and without random effects keep the normal
distribution, as in [`stats::glm()`](https://rdrr.io/r/stats/glm.html).
For Gaussian models without random effects the corrected intervals of AD
contrasts are exact t intervals; with repeated measures, and for the
interactions between a stressor applied to whole groups and one applied
within them, the rule is conservative.

**Combinations with only zeros.** When the full combination of a set has
only zeros at a time (or, for binomial responses, only zeros or only
successes) and the model selected by
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
has a log or logit link, the null model still predicts a mean E for it
from the other cells. The p-value is twice the probability of observing
only zeros if the mean were E, averaged over the draws of E
(`status = "boundary"`; Poisson, NB2, Tweedie, binomial and
beta-binomial models). Random intercepts are integrated out:
observations that share a level of a random effect share one intercept
with the variance of all random terms, which is exact with a single
grouping (units or rooms, say) and conservative with nested or crossed
ones. The same test is used for a single stressor whose cell has only
zeros (against the reference cell). The estimate is computed with the
observed mean (zero), which makes MU estimates `-Inf`; intervals are not
given. If the reference cell, or a single-stressor cell that an
interaction needs, has only zeros, that contrast cannot be measured and
the row is `"not_estimable"`.

**Multiple testing.** p-values are adjusted within each response, null
model, component, stressor and estimand, across contrast sets, sampling
times and `by` groups; averages over times (`"overall"`) are adjusted
separately.

## Bliss independence

Under Bliss independence (`"BL"`; Bliss 1939) stressors act
independently, so the responses relative to the reference cell, \\w_S =
y_S / y_0\\ (relative fitness or survival), multiply. Single and
combined effects are \\w - 1\\, the net interaction is \\w_K - \prod_i
w_i\\ (the DA measure of Beppler et al. 2016, eqs. 2.1–2.2; Tekin et al.
2018) and, for three stressors, the pure interaction is the emergent
interaction E3 of Beppler et al. (2016, eq. 2.4): \$\$w\_{XYZ} - w_X
w\_{YZ} - w_Y w\_{XZ} - w_Z w\_{XY} + 2 w_X w_Y w_Z,\$\$ the part of the
three-way deviation that the pairwise interactions do not explain. It is
not defined here for four or more stressors (`status = "not_defined"`).
The net BL interaction tests the same hypothesis as the net MU
interaction on another scale (a ratio of the combined response to its
expectation under MU, a difference under BL); the emergent interaction
is a different hypothesis from the pure MU interaction, which removes
the pairwise interactions on the log scale.

`rescaled` divides each BL interaction by the absolute value of the same
expression with the combined response replaced by zero, when the
interaction is negative (lethal synergy), or by the smallest
single-stressor response, when it is positive (complete buffering), so
that -1 and 1 mark those two extremes (Segrè et al. 2005; Beppler et al.
2016, section 2.2.3). With a non-negative response it is at least -1;
when no stressor increases the response, it exceeds 1 if the combination
is less affected than its most affected single stressor.
`rescaled_class` labels values below -0.5 as `"synergy"`, above 0.5 as
`"antagonism"` and the rest as `"additive"`, the conservative cut-offs
of Beppler et al. (2016, section 2.2.4). Both describe the size of the
interaction only. Whether it deviates from Bliss independence is
decided, as for the other null models, by its interval and p-value,
because fixed cut-offs ignore the uncertainty of the estimate (Madge
Pimentel et al. 2025).

## Stress addition and concentration addition

**Stress addition (SAM).** Liess et al. (2016) assume that each
individual tolerates a general stress drawn from a beta distribution
with both shape parameters equal to `sam_shape`. The relative response
of each stressor alone, \\w_i\\ (survival relative to the reference
cell), becomes the general stress whose quantile in that distribution is
the mortality \\1 - w_i\\; the general stresses of the stressors add up,
and the predicted relative response of the combination is the share of
individuals that tolerate their sum, zero when it reaches 1 (Liess et
al. 2016, eqs. 2-5). The interaction is \\w_K\\ minus that prediction;
single and combined effects are relative, as under BL. A stressor that
increases the response exerts no stress; if one does so at the point
estimate, the interaction is `"not_applicable"`. The model predicts
larger combined effects than Bliss independence, so combinations that
Bliss independence calls additive can be antagonistic under it.

**Concentration addition (CA, Loewe additivity).** Each stressor acts as
a dilution of the others: a combination of doses \\x_i\\ has the
relative response w at which \\\sum_i x_i / EC_i(w) = 1\\, where
\\EC_i(w)\\ is the dose of stressor i alone that gives w (Loewe 1953;
Tang et al. 2015, eq. 2). Doses are the differences between the levels
and the control, so the stressor columns must hold numbers. \\EC_i\\
comes from a two-parameter log-logistic curve, \\w = 1 / (1 + (x /
EC\_{50})^b)\\, fitted by least squares to the relative responses at
every dose of the stressor tested alone, with the other stressors at the
reference. At least two such doses are needed, so CA is `"not_defined"`
in two-level factorials. The curves are refitted in each simulation
draw, so their uncertainty enters the intervals. A curve that does not
decrease with dose makes the interaction `"not_applicable"`; in single
draws it is taken as no dose equivalent. Doses tested alone with only
zeros keep their means at the boundary and still inform the curves.

Both models give only the net interaction (the pure interaction of three
or more stressors is `"not_defined"`), and their single and combined
effects are those of BL.

## Equivalence

A non-significant interaction is a failure to detect a deviation, not
evidence of additivity. To claim that a combination follows a null
model, declare before the analysis the largest deviation that would be
negligible, `margin`, as a share of the response of the reference cell.
Each interaction is then also tested against the hypotheses that its
relative deviation d = (observed - expected) / reference is at most
-`margin` and at least `margin`, with two one-sided tests (Schuirmann
1987; Lakens et al. 2018): `p_equivalence` is the larger of their
p-values.

Each one-sided test measures how far the estimated means of the cells
involved are from the closest means at which d equals the bound, in
standard errors (the Mahalanobis distance, signed), and compares that
distance with a t distribution with the degrees of freedom of the
interaction (normal if they are infinite). The distance is measured on a
scale on which the variance of a mean does not depend on the mean
(Anscombe 1948): the square root for Poisson, quasi-Poisson and other
counts whose variance grows in proportion to the mean, the inverse
hyperbolic sine of the square root for the negative binomial, the
arcsine of the square root for proportions, the logarithm for Gamma and
lognormal responses and the response itself for Gaussian ones (for other
families, and for models of a transformed response, the scale of the
model). When the bounds are linear in the means (AD with a constant
variance), the tests are t tests of observed - expected ± `margin` ×
reference, as in the t test of a ratio to a control (Hauschke et al.
1999). Tests built from the interval of d, or from its draws, declared
false equivalences too often with overdispersed counts, where the
standard error of a mean shrinks with the mean. The equivalence test
uses neither the draws nor the gradient, so it is the same with both
`inference` methods; near the margin, its verdict can differ from
whether the interval of d at 1 - 2 alpha lies within the margin. If the
closest means cannot be found, `p_equivalence` is `NA` with the note
`"equivalence_not_tested"`.

The relative deviation puts every null model on the same scale: under
BL, SAM and CA it is the interaction itself; under AD and DO, the
interaction divided by the reference; under MU, the observed minus the
multiplicative expectation, divided by the reference (the net BL
interaction). Its interval includes the uncertainty of the reference.
When the stressors raise the response far above the reference
(nutrients, say), relative deviations can be large. Combinations with
only zeros get a relative deviation but no interval and no test of
equivalence.

[`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) combine both
tests, with the adjusted p-values, into four verdicts: `"deviation"`
(deviation detected, not equivalent), `"negligible"` (a deviation
detected, but within the margin), `"equivalent"` (no deviation detected
and within the margin) and `"inconclusive"` (neither). Only
`"equivalent"` supports a claim that the stressors combine as the null
model predicts.

In simulated 2 x 2 experiments with four units per combination, the test
declared false equivalence (a true relative deviation exactly at the
margin of ±0.25) in 4–6% of the experiments, close to the nominal 5%,
for Gaussian responses under AD, survival under BL and overdispersed
counts under MU analysed with a quasi-Poisson
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) or with the negative
binomial of
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
(3–6% for counts with no, moderate or strong overdispersion). A Poisson
model fitted to the same overdispersed counts ignores their
overdispersion and gave 21%, which is why
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
uses it only when the negative binomial does not converge;
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html), whose
maximum-likelihood fit with the normal distribution underestimates the
variance with so few units, gave 10%. Equivalence needs replication:
with four units, a margin of ±0.1 was never reached; one of ±0.25 was
reached in about half (Gaussian) to three quarters (survival) of the
experiments that followed the null model exactly, but rarely (under 10%)
with overdispersed counts.

## Default null model

When there is a hypothesis about how the stressors act, name the null
model that represents it in `null_models` (Madge Pimentel et al. 2025);
without one, AD is used. It was chosen by simulation, with a rule fixed
before the results were seen: among the null models defined in two-level
factorials and with at least 80% of the highest average power against
interactions that no mechanism explains, the one with the fewest false
interactions on average over five mechanisms of no interaction (AD, MU
or BL, DO, SAM and CA). The experiments were 2 x 2 factorials and 3 x 3
dose factorials with four units per combination, survival out of 20
animals or overdispersed counts, and moderate or weak single effects,
analysed with a cell-means
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) (100 experiments per
scenario; p \< 0.05 without adjustment).

Each null model flagged interactions at about the nominal rate (3–5%)
when the data followed its own mechanism, but none was robust to every
mechanism: data generated by stress addition looked synergistic to all
the others (in 60% of the experiments under AD, 76% under MU and BL and
86% under DO), and data generated by the others looked antagonistic to
stress addition. On average, AD flagged false interactions in 20% of the
experiments, MU and BL in 24%, SAM in 27% and DO in 32%, with average
powers of 58%, 67%, 38% and 61%. In the dose design CA flagged fewer
(16%, against 23% under AD), with less power (51% against 63%). AD was
also chosen within each design, each size of the single effects and with
survival as the response. Without the stress addition mechanism, MU and
BL came within one percentage point of AD and, being more powerful,
would have been chosen; with counts only, the power of AD fell below 80%
of theirs, and they would have been chosen too. Requiring every null
model to agree (`agreement = "robust"` in
[`compare_null_models()`](https://pabloamadorc.github.io/must/reference/compare_null_models.md))
flagged under 3% false interactions, but only 15% of the interactions
that no mechanism explains.

The simulation is installed with the package:
`system.file("simulations", package = "must")` holds the script
(`null_model_robustness.R`), which reruns or extends it, the rates of
every scenario (`null_model_robustness.csv`) and its report
(`null_model_robustness.md`).

## Models fitted outside must

Any model supported by emmeans (Lenth & Piaskowski 2026) can be
analysed: [`stats::lm()`](https://rdrr.io/r/stats/lm.html),
[`stats::glm()`](https://rdrr.io/r/stats/glm.html),
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html),
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html),
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html),
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html) and
many more. The stressor, `time` and `by` columns of `design` must be
predictors of the model; each of their values is a level, also for a
numeric predictor, which is evaluated at each value observed. The mean
of each treatment combination (and time) is the estimated marginal mean
of emmeans: the prediction on the link scale with random effects at
zero, averaged with equal weights over the other factors of the model
and with other numeric predictors at their mean (unless set in `...`),
back-transformed with the inverse link or the transformation of the
response without bias adjustment (with `lm(log(y) ~ ...)` they are
medians of y, not means). The draws and the delta method use the
covariance of the coefficients of the model as it was fitted (no REML
refit).

With `correction = "small_sample"` the degrees of freedom are those
emmeans gives for the linear combination of coefficients of each
contrast (its gradient for non-linear contrasts): the residual degrees
of freedom for [`stats::lm()`](https://rdrr.io/r/stats/lm.html) and
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) with an estimated
dispersion and for
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) with any family
(observations minus effective degrees of freedom), Kenward–Roger or
Satterthwaite for
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) when pbkrtest
or lmerTest is installed, and infinite for Poisson and binomial models
fitted with [`stats::glm()`](https://rdrr.io/r/stats/glm.html) and for
most other mixed models. With few experimental units, infinite degrees
of freedom give tests that reject a true null model too often: `df` sets
them by hand, or
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
applies the correction described above.

The estimates inherit the structure of the model. A quantity the model
fixes at zero, such as the MU interaction of a log-link model without
the interaction term, is `"not_applicable"`. Means the model cannot
estimate (in rank-deficient models), or that lie at the boundary of the
link function (a combination with only zeros under a log link, or only
zeros or only successes under a logit link) make their contrasts
`"not_estimable"`: such models give no usable standard error for them,
whereas
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
tests them from the probability of that outcome.

## Assumptions

- Means are conditional predictions (random effects at zero). With a log
  link and a random intercept, marginal means differ by a constant
  factor: MU and BL estimates and all p-values are unchanged, AD and DO
  estimates change scale.

- MU, DO, BL, SAM and CA only represent their mechanistic assumptions
  when an increase in the response means tolerance (survival, abundance)
  rather than sensitivity (mortality); declare the latter in
  `sensitivity`. AD may predict impossible values below zero, in which
  case a deviation from it does not imply a mechanistic interaction
  (Madge Pimentel et al. 2025).

- A non-significant interaction is a failure to detect a deviation, not
  evidence of additivity; equivalence needs a `margin` (see
  "Equivalence").

## Missing values and zeros

MU needs positive means and BL a positive mean in the reference cell: if
a fitted or simulated mean is not positive (possible with an identity
link), the MU or BL rows of that set are `"not_estimable"`. Treatment
combinations without values of the response (at a time) make the sets
that need them `"not_estimable"` (at that time). Times declared as
`baseline` in
[`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md)
are not estimated.

## Reproducibility

With `inference = "simulation"` results depend on the random number
stream: the seed used is in `settings$seed`, and `mc_se` gives the Monte
Carlo error of each simulated p-value. The same seed may not give
identical draws in other versions of R or with other linear algebra
libraries.

## References

Beppler, C., Tekin, E., Mao, Z., White, C., McDiarmid, C., Vargas, E.,
Miller, J. H., Savage, V. M., & Yeh, P. J. (2016). Uncovering emergent
interactions in three-way combinations of stressors. *Journal of the
Royal Society Interface*, 13, 20160800.
[doi:10.1098/rsif.2016.0800](https://doi.org/10.1098/rsif.2016.0800)

Bliss, C. I. (1939). The toxicity of poisons applied jointly. *Annals of
Applied Biology*, 26(3), 585–615.
[doi:10.1111/j.1744-7348.1939.tb06990.x](https://doi.org/10.1111/j.1744-7348.1939.tb06990.x)

Krinsky, I., & Robb, A. L. (1986). On approximating the statistical
properties of elasticities. *The Review of Economics and Statistics*,
68(4), 715–719. [doi:10.2307/1924536](https://doi.org/10.2307/1924536)

Lenth, R. V., & Piaskowski, J. (2026). emmeans: Estimated marginal
means, aka least-squares means. R package version 2.0.4.
<https://CRAN.R-project.org/package=emmeans>

Liess, M., Foit, K., Knillmann, S., Schäfer, R. B., & Liess, H.-D.
(2016). Predicting the synergy of multiple stress effects. *Scientific
Reports*, 6, 32965.
[doi:10.1038/srep32965](https://doi.org/10.1038/srep32965)

Loewe, S. (1953). The problem of synergism and antagonism of combined
drugs. *Arzneimittelforschung*, 3, 285–290.

Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B.
(2025). Hypothesis-driven research on multiple stressors: An analytical
framework for stressor interactions. *Ecology and Evolution*, 15,
e71959. [doi:10.1002/ece3.71959](https://doi.org/10.1002/ece3.71959)

Anscombe, F. J. (1948). The transformation of Poisson, binomial and
negative-binomial data. *Biometrika*, 35(3–4), 246–254.
[doi:10.1093/biomet/35.3-4.246](https://doi.org/10.1093/biomet/35.3-4.246)

Hauschke, D., Kieser, M., & Hothorn, L. A. (1999). Proof of safety in
toxicology based on the ratio of two means for normally distributed
data. *Biometrical Journal*, 41(3), 295–304.
[doi:10.1002/(SICI)1521-4036(199906)41:3\<295::AID-BIMJ295\>3.0.CO;2-2](https://doi.org/10.1002/%28SICI%291521-4036%28199906%2941%3A3%3C295%3A%3AAID-BIMJ295%3E3.0.CO%3B2-2)

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
for psychological research: A tutorial. *Advances in Methods and
Practices in Psychological Science*, 1(2), 259–269.
[doi:10.1177/2515245918770963](https://doi.org/10.1177/2515245918770963)

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin*, 2(6), 110–114.
[doi:10.2307/3002019](https://doi.org/10.2307/3002019)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics*, 15(6), 657–680.
[doi:10.1007/BF01068419](https://doi.org/10.1007/BF01068419)

Segrè, D., DeLuna, A., Church, G. M., & Kishony, R. (2005). Modular
epistasis in yeast metabolism. *Nature Genetics*, 37, 77–83.
[doi:10.1038/ng1489](https://doi.org/10.1038/ng1489)

Tang, J., Wennerberg, K., & Aittokallio, T. (2015). What is synergy? The
Saariselkä agreement revisited. *Frontiers in Pharmacology*, 6, 181.
[doi:10.3389/fphar.2015.00181](https://doi.org/10.3389/fphar.2015.00181)

Tekin, E., Yeh, P. J., & Savage, V. M. (2018). General form for
interaction measures and framework for deriving higher-order emergent
effects. *Frontiers in Ecology and Evolution*, 6, 166.
[doi:10.3389/fevo.2018.00166](https://doi.org/10.3389/fevo.2018.00166)

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
# Simple addition (the default)
est <- estimate_interactions(fits, n_sim = 2000, seed = 1)
est
#> <must_interactions> 1 response(s); null models: AD
#> Inference: Krinsky-Robb simulation (2000 draws), small-sample correction (REML, t); 95% intervals; p-values adjusted with BH; seed 1
#> Scales: AD = differences (response scale)
#> Interactions (1 rows):
#>  response                 combination model estimand estimate lower upper  df
#>     algae salinity=high + warming=yes    AD      net      9.8  1.47  18.5 Inf
#>  p_adj status
#>  0.016     ok
est$estimates[est$estimates$component == "interaction",
  c("null_model", "estimate", "lower", "upper", "df", "p_adjusted")]
#>   null_model estimate    lower    upper  df p_adjusted
#> 4         AD      9.8 1.471868 18.47447 Inf      0.016
# The three null models and the procedure of Madge Pimentel et al. (2025),
# without the correction
estimate_interactions(fits, null_models = c("AD", "MU", "DO"), n_sim = 2000,
  seed = 1, correction = "none")
#> <must_interactions> 1 response(s); null models: AD, MU, DO
#> Inference: Krinsky-Robb simulation (2000 draws), no small-sample correction; 95% intervals; p-values adjusted with BH; seed 1
#> Scales: AD = differences (response scale); MU = log ratios; DO = dominance (response scale)
#> Interactions (3 rows):
#>  response                 combination model estimand  estimate   lower  upper
#>     algae salinity=high + warming=yes    AD      net   9.80000   1.470 18.500
#>     algae salinity=high + warming=yes    MU      net  -0.00501  -0.413  0.413
#>     algae salinity=high + warming=yes    DO      net -10.00000 -12.800 -4.730
#>   df  p_adj status
#>  Inf 0.0160     ok
#>  Inf 0.9700     ok
#>  Inf 0.0005     ok
# Additive, or only not detected? An equivalence margin of 10% of the control
eq <- estimate_interactions(fits, null_models = c("AD", "MU"), n_sim = 2000, seed = 1,
  margin = 0.1)
eq$estimates[eq$estimates$component == "interaction",
  c("null_model", "relative", "relative_lower", "relative_upper", "p_equivalence")]
#>   null_model    relative relative_lower relative_upper p_equivalence
#> 4         AD  0.24623114      0.0403576     0.42205268    0.92420415
#> 8         MU -0.00126261     -0.1144368     0.09664942    0.04634745
classify_interactions(eq)[, c("null_model", "class", "verdict")]
#>   null_model class    verdict
#> 1         AD    -A  deviation
#> 2         MU    AD equivalent
# Bliss independence, with its rescaled interaction
bl <- estimate_interactions(fits, null_models = "BL", n_sim = 2000, seed = 1)
bl$estimates[bl$estimates$component == "interaction",
  c("estimate", "lower", "upper", "p_value", "rescaled", "rescaled_class")]
#>      estimate      lower      upper p_value    rescaled rescaled_class
#> 4 -0.00126261 -0.1144368 0.09664942   0.967 -0.00500006       additive
# A model fitted outside must
model <- stats::glm(algae ~ salinity * warming, family = stats::poisson(), data = jars)
estimate_interactions(model, design, null_models = c("AD", "MU", "BL"), n_sim = 2000,
  seed = 1)
#> <must_interactions> 1 response(s); null models: AD, MU, BL
#> Model: glm, poisson (log); means from emmeans
#> Inference: Krinsky-Robb simulation (2000 draws), degrees of freedom from emmeans; 95% intervals; p-values adjusted with BH; seed 1
#> Scales: AD = differences (response scale); MU = log ratios; BL = relative to the reference cell (Bliss independence)
#> Interactions (3 rows):
#>  response                 combination model estimand estimate  lower   upper
#>     algae salinity=high + warming=yes    AD      net  9.80000  1.470 18.5000
#>     algae salinity=high + warming=yes    MU      net -0.00501 -0.413  0.4130
#>     algae salinity=high + warming=yes    BL      net -0.00126 -0.114  0.0966
#>   df p_adj rescaled status
#>  Inf 0.016       NA     ok
#>  Inf 0.970       NA     ok
#>  Inf 0.970   -0.005     ok
```
