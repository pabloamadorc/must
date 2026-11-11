# Fit a GLMM to each response, choosing its error distribution

Fits, for every response, one generalised linear mixed model per
candidate error distribution with
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html),
and keeps the one that fits the data best. The null models are not
involved here: following the post-estimation approach, the model is
chosen for the data and the null models are tested afterwards from its
predictions
([`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)).

## Usage

``` r
fit_responses(
  data,
  responses,
  design,
  families = NULL,
  response_types = NULL,
  trials = NULL,
  offset = NULL,
  time_structure = "interaction",
  dispersion = "common",
  random = NULL,
  control = NULL,
  residual_checks = TRUE,
  n_sim_residuals = 500L,
  seed = NULL,
  quiet = FALSE
)
```

## Arguments

- data:

  The data frame given to
  [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md)
  (same treatment combinations; other columns may be added).

- responses:

  Character vector with the names of one or more numeric response
  columns.

- design:

  A `must_design` object from
  [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md).
  Designs with checks of severity `"problem"` are refused.

- families:

  `NULL` (default) to choose the candidates from each response (see
  Details); a named list of family objects used for every response, e.g.
  `list(gaussian = gaussian(), gamma_log = Gamma("log"))`; or a named
  list with one such list per response.

- response_types:

  `NULL` (default) or a named character vector with `"count"` or
  `"continuous"` for some responses, to override the automatic detection
  used to choose the default candidates.

- trials:

  `NULL` (default) or a character vector named after responses that
  count successes out of a known number of trials (survivors out of the
  animals added, emerged out of the eggs), giving the column with the
  number of trials of each observation, e.g. `c(survivors = "n_start")`.
  These responses are modelled as binomial or beta-binomial and their
  means are probabilities.

- offset:

  `NULL` (default) or a character vector named after count responses
  giving the column with the sampling effort of each observation (volume
  filtered, area or time sampled), e.g. `c(daphnia = "volume_l")`. The
  model includes `offset(log(effort))`, only log-link candidates are
  used, and means are per unit of effort. Every value of the response
  needs an effort.

- time_structure:

  How treatment effects relate to sampling time when the design has
  `time`: `"interaction"` (default; a separate mean for each treatment
  combination and time, no assumption of constant effects) or
  `"additive"` (treatment effects constant over time on the scale of the
  link function).

- dispersion:

  Which observations share a dispersion (variance) parameter: `"common"`
  (default; one for all), `"by"` (one per `by` group) or `"cell"` (one
  per treatment combination and `by` group; the candidates with a
  dispersion parameter then need two or more units in every
  combination). Ignored for Poisson and binomial candidates, which have
  none.

- random:

  `NULL` (default) adds a random intercept for each `groups` column of
  the design that can be estimated (see
  [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md))
  and one per experimental unit when some unit has more than one row
  (repeated measures or subsamples); `"none"` removes random effects;
  any other string is used as the random part of the formula in lme4
  syntax, with column names of `data`, e.g.
  `"(1 | tank) + (1 | block)"`. Custom terms can also use the internal
  columns `.unit` (the units), `.time` (the sampling times, as a factor)
  and `.g_<name>` (the groups of each `groups` column, made unique
  within `by` groups).

- control:

  `NULL` (default) or the output of
  [`glmmTMB::glmmTMBControl()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMBControl.html),
  passed to every fit. By default the iteration and evaluation limits of
  the optimiser grow with the number of fixed effects (see Details); the
  rest of the settings are those of glmmTMB.

- residual_checks:

  If `TRUE` (default), simulation-based residual checks of the DHARMa
  package (uniformity, dispersion and zero inflation) decide which
  candidates may be selected. Needs DHARMa installed.

- n_sim_residuals:

  Number of simulations for the residual checks. Default 500.

- seed:

  `NULL` (default) or a whole number. The residual checks run with this
  seed through
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html),
  which restores the random number state of the session afterwards. With
  `NULL`, a seed is drawn from the session's random number stream;
  either way it is stored in `settings$seed`, so the selection can be
  repeated exactly.

- quiet:

  If `FALSE` (default), a message reports each response as it is fitted.

## Value

An object of class `must_fit`, a list with:

- `fits`: a list with one element per response, holding `response`,
  `status` (`"ok"` or `"failed"`), `reason`, `model` (the selected
  `glmmTMB` fit, by maximum likelihood), `model_reml` (the same model
  refitted by REML, used for the small-sample correction of
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md);
  `NULL` if the refit failed), `family` (name of the selected
  candidate), `link`, `response_type` and `type_source` (`"detected"`,
  `"declared"`, `"trials"` or `"families"`), `selection` (data frame
  with one row per candidate: `family`, `link`, `random` (random part
  used), `converged`, `aic`, `p_uniformity`, `p_dispersion`,
  `p_zero_inflation`, `n_tests_passed`, `residuals_ok`, `selected`,
  `message`), `checks_passed` (`TRUE`, `FALSE` if no candidate passed
  all the residual checks, or `NA` if they were not run), `formula` (of
  the selected candidate), `n_obs`, `n_missing` (rows without a response
  value) and `grid` (treatment combinations and times predicted later).

- `design`: the design.

- `settings`: the arguments used, with the seed (and whether it was
  given), the families of the candidates given in `families` and the
  versions of must, R, glmmTMB and, with the residual checks, DHARMa.

## Details

**Model.** Treatment combinations are coded as one factor (a cell-means
model), so any design is handled: full or partial factorials, ray
designs or a few combinations of many stressors. With
`time_structure = "interaction"` the fixed part is one mean per
treatment combination, time and `by` group; with `"additive"` it is one
mean per treatment combination and `by` group plus a common effect of
time.

**Default candidates.** Responses with `trials` are binomial: binomial
and beta-binomial, both with a logit link. Otherwise whole non-negative
numbers are treated as counts (the `type_source` says the type was
detected): Poisson and negative binomial (NB2), both with a log link,
but the Poisson is selected only if the negative binomial does not
converge, as happens when the counts show no overdispersion to estimate.
With overdispersed counts and four units per combination, the AIC chose
the Poisson in a fifth of simulated experiments (the residual checks
passed it too), and then the tests of
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
rejected true null models too often; fitting the negative binomial kept
them close to their nominal error, at a small cost in power when the
counts are truly Poisson. Continuous values: Gaussian, Gaussian with log
link, Gamma with log link and lognormal if all values are positive;
Gaussian, Gaussian with log link and Tweedie if there are zeros;
Gaussian only if there are negative values. With an `offset` only the
log-link candidates are kept. Use `response_types` when the detection is
wrong (a density rounded to whole numbers is not a count) and `families`
for other distributions (e.g.
[`glmmTMB::beta_family()`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html)
for proportions without a number of trials, such as cover).

**Optimiser.** A cell-means model has one parameter per treatment
combination (and time), which can be more than the default limits of the
optimiser of glmmTMB (300 iterations, 400 evaluations) allow. By default
they are raised to 20 iterations and 30 evaluations per fixed effect
(never below glmmTMB's defaults); a fit that reaches them still counts
as not converged.

**Combinations with only zeros.** On a log link, the mean of a treatment
combination (at a time) whose values are all zero is on the boundary of
the parameter space (minus infinity on the link scale); on a logit link
the same happens when all values are zero or all are successes. For
these candidates the mean of such a combination is fixed at that
boundary (`map` of glmmTMB), so the rest of the model is estimated
normally; their AIC counts the fixed parameters, so that it stays
comparable with the identity-link candidates. With
`time_structure = "additive"` only combinations with such values at
every time are fixed (otherwise their mean is not on the boundary).
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
tests the contrasts that involve these combinations from the probability
of the observed zeros.

**Random intercepts at zero.** With few units, the variance of a random
intercept is often estimated at zero, which often makes the Hessian
singular. If that happens to an automatic random intercept
(`random = NULL`), the candidate is refitted without it and the reduced
model is kept when it reaches the same maximum likelihood, whether or
not the full model converged; the `random` column of `selection` and the
`reason` say so. Terms given in `random` are never removed.

**Selection.** Candidates that fail to converge or have a non-positive
definite Hessian are discarded, and so is the default Poisson for counts
when the negative binomial is left. With `residual_checks = TRUE`, a
candidate passes when the three DHARMa tests give p \>= 0.05; if at
least one candidate passes, only those compete, and otherwise those
passing the most tests (`checks_passed = FALSE`: inspect that model
before using it). The candidate with the lowest AIC is selected.
Candidates whose AIC differs by less than 1e-4 fit equally well: in a
cell-means model without random effects the link does not change the
fitted means. The log link is then preferred, because its means stay
positive as the multiplicative and dominance null models need, unless
some combination has only zeros, in which case the identity link is
preferred; the `reason` says which.

**REML refit.** The selected model is refitted by restricted maximum
likelihood (REML), whose variance estimates account for the number of
fixed effects.
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
uses it for its default small-sample correction; the selection itself
uses the maximum-likelihood AIC.

## Assumptions

- A random intercept per unit implies the same correlation between any
  two samplings of a unit. Temporal autocorrelation is only modelled if
  a structure is given in `random`, e.g.
  `random = "ar1(.time + 0 | .unit)"` (first-order autoregressive; it
  treats the sampling times as equally spaced).

- With `dispersion = "common"`, all treatment combinations and `by`
  groups share one dispersion parameter on the scale of the family; with
  very different groups (species, populations), use `"by"` or `"cell"`.

- The selection among candidates is not accounted for in the uncertainty
  of later estimates.

## Missing values and zeros

Rows with a missing response are left out for that response only and
counted in `n_missing`. Zeros are kept as observed values and no
transformation is applied (see "Combinations with only zeros").

## References

Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg,
C. W., Nielsen, A., Skaug, H. J., Mächler, M., & Bolker, B. M. (2017).
glmmTMB balances speed and flexibility among packages for zero-inflated
generalized linear mixed modeling. *The R Journal*, 9(2), 378–400.
[doi:10.32614/RJ-2017-066](https://doi.org/10.32614/RJ-2017-066)

Madge Pimentel, I., et al. (2025). Hypothesis-driven research on
multiple stressors: An analytical framework for stressor interactions.
*Ecology and Evolution*, 15, e71959.
[doi:10.1002/ece3.71959](https://doi.org/10.1002/ece3.71959)

## Examples

``` r
jars <- expand.grid(
  salinity = c("ambient", "high"),
  warming = c("no", "yes"),
  replicate = 1:5,
  stringsAsFactors = FALSE
)
jars$jar <- seq_len(nrow(jars))
# Algal counts: each stressor halves the expected count
mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
  ifelse(jars$warming == "yes", 0.5, 1)
jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
# Survivors out of 20 animals per jar
jars$n_start <- 20
jars$survivors <- c(18, 15, 16, 9)[match(paste(jars$salinity, jars$warming),
  c("ambient no", "high no", "ambient yes", "high yes"))] +
  c(-1, 0, 1, 0, 0)[jars$replicate]
design <- describe_design(jars, c("salinity", "warming"),
  control = c(salinity = "ambient", warming = "no"), unit = "jar")
fits <- fit_responses(jars, c("algae", "survivors"), design,
  trials = c(survivors = "n_start"), residual_checks = FALSE, seed = 1,
  quiet = TRUE)
fits
#> <must_fit> 2 response(s); time structure: interaction; dispersion: common; random effects: none
#> Residual checks: not run; seed: 1
#>   algae: nbinom2 (log link; count, detected), AIC 114.2, n = 20, missing = 0, residual checks: not run
#>   survivors: binomial (logit link; binomial, trials), AIC 71.4, n = 20, missing = 0, residual checks: not run
fits$fits$algae$selection
#>    family link random converged      aic p_uniformity p_dispersion
#> 1 poisson  log   none      TRUE 112.2431           NA           NA
#> 2 nbinom2  log   none      TRUE 114.2431           NA           NA
#>   p_zero_inflation n_tests_passed residuals_ok selected message
#> 1               NA             NA           NA    FALSE        
#> 2               NA             NA           NA     TRUE        
```
