# must 0.3.0

* `plot()` of `must_interactions` objects has been redesigned and draws four
  kinds of figure, chosen with `type`:
  * `"effects"` (the new default): for each interaction, the effect of each
    stressor alone, the combined effect expected under the null model, the
    combined effect observed and the interaction, with its class at the
    right: how the interaction arises.
  * `"interactions"`: the interactions alone, one row each, on the scale of
    the null model or as relative deviations (`scale`), optionally under
    headers by set, time or group (`group_by`).
  * `"models"`: the interaction under each null model as a deviation
    relative to the reference, with its class: whether the conclusion
    depends on the null model.
  * `"time"`: one figure per set with the single and combined effects and
    the interaction over the sampling times.

  Sets, times and groups are chosen with `set_id`, `time` and `context`;
  labels with `labels`, `stressor_labels`, `label_style` and `label_width`.
  Classes are written as codes (`S+`, `S-`, `A+`, `A-`, `n.d.`) or names
  (`classes`). Long figures are split into pages (`per_page`). Infinite
  estimates, estimates outside given limits and intervals cut by the limits
  end in arrows (`type = "models"` cuts intervals that reach far beyond the
  estimates). `plot()` returns the rows drawn invisibly, and arguments it
  does not use are now an error instead of being passed to
  `plot.default()`.

* `estimate_interactions()` no longer fails with glmmTMB 1.1.11 to 1.1.15
  when two or more treatment combinations have only zeros (or only
  successes): `vcov()` of those versions fails for a model with two or more
  parameters fixed at the boundary, and the covariance of the fixed effects
  is now built from that of the estimated parameters.

* Andreu Rico is now an author of the package: `citation("must")` and the
  methods written by `report_interactions()` cite both authors.

# must 0.2.0

* The package has a website, <https://pabloamadorc.github.io/must/>, with
  the help pages, the vignettes and the case studies.

* Data: `polazzo_zooplankton` now has every taxon of the source (its eight
  columns are unchanged), and two new tables hold the other communities of
  that experiment (`polazzo_phytoplankton`, `polazzo_macroinvertebrates`).
  New data of Hermann et al. (2024), an outdoor mesocosm experiment with an
  insecticide (imidacloprid) and three temperature regimes (ambient,
  elevated and heatwaves): zooplankton, phytoplankton, macroinvertebrates,
  water quality, nutrients, leaf litter decomposition and imidacloprid
  concentrations (`hermann_zooplankton` and six more tables). Community
  tables have one column per taxon, with its name in the source file in
  the attribute `"taxa"`.

* New data of six experiments with different designs, from two to five
  stressors: drought and plant invasion on field plots in blocks
  (`alba_arthropods`), salinisation and nutrient enrichment in two orders
  (`fuggle_biofilms`), temperature, photoperiod and microbiome removal in
  two diatoms (`giesler_diatoms`), temperature, salinity and oxygen in
  horseshoe crab embryos (`vasquez_embryos`), four stressors on stream
  mesocosms fed by header tanks (`hunn_invertebrates`) and five global
  changes on split plots (`niboyet_denitrification`). The article "Case
  studies" of the package website reanalyses each of them.

* New vignette, "Two mesocosm experiments" (`vignette("mesocosms")`): the
  whole workflow, from the design to the methods paragraph, on the
  phytoplankton and macroinvertebrates of Polazzo & Rico (2021) and the
  zooplankton, chironomids and leaf litter decomposition of Hermann et al.
  (2024), with baselines common to all stressors or per stressor, a density
  declared continuous, a taxon eliminated by a stressor, the null models
  suggested for each response and a test of equivalence.

* New `report_interactions()`: writes the methods of an analysis from the
  settings of the output of `estimate_interactions()` and, if given, of
  `classify_interactions()` and `compare_null_models()`: the design, the
  model of each response and how its error distribution was chosen, the null
  models and why (the criteria of `suggest_null_models()` or a
  `justification`), how the interactions were estimated and tested, and how
  they were classified, with the works cited. A table lists every decision
  with the function and argument that control it, the value used, its
  default and its source (default, set by the user, chosen by must,
  estimated, suggested or recorded). `format()` gives plain text or markdown
  for R Markdown and Quarto documents. To support it, `fit_responses()` and
  `estimate_interactions()` record the versions of must, R and the packages
  used and whether the seed was given, `estimate_interactions()` keeps a
  summary of the models in `settings$fit`, and the outputs of
  `classify_interactions()` and `compare_null_models()` keep `alpha` (and
  the null models compared) as attributes.

* New logo: the red hexagon with the name of the package, on a white
  background.

* New `suggest_null_models()`: the guidelines of Schäfer & Piggott (2018,
  Fig. 3) for the choice of a null model, from whether the stressors act
  similarly, whether one acts specifically on physiology, whether the effect
  is mortality, whether the stressor-effect relationship or the distribution
  of sensitivities is known, how the sensitivities are correlated and how
  large the effects of the stressors alone are relative to the effect limit.
  Unknown criteria take every value, and all the null models they lead to
  are kept, with the steps that lead to each, why the others are excluded
  and which unknown criteria would change the choice. The effects can be
  given or estimated from a fitted model (also one fitted outside must) in
  every contrast set and time, from the cells with one stressor alone, with
  intervals: a criterion whose interval spans its threshold stays unknown.
  Stressors may increase the response: each effect has a direction and each
  set the interaction type of Piggott et al. (2015); concentration addition
  needs effects in the same direction and the stress addition model no
  increase; and for a proportion that the stressors increase, the guidelines
  are applied to its complement, on which the interactions are then
  estimated. `estimate_interactions()` takes the result as `null_models`,
  checks its `sensitivity` and keeps its criteria and reasons in
  `settings$null_model_choice`.

* `describe_design()` gains `shared`, for rows that belong to every group of
  a `by` column: a control common to two exposure regimes, to two orders of
  application of the stressors, or to upstream and downstream ponds. Those
  rows enter the analysis once and supply the cells that each group lacks,
  instead of being copied into every group. `fit_responses()` and
  `estimate_interactions()` (also for models fitted outside must) use them
  in the sets of each group; `set_cells` gives the group in which each cell
  was observed.

* `baseline` in `describe_design()` can be a named list with the times
  before each stressor was applied, for stressors that started later than
  others (drying started after the warming and shading treatments): an
  interaction is not estimated before all of its stressors were applied,
  and the interactions among the other stressors are. With
  `time_structure = "additive"`, those samplings are coded with the
  stressors not yet applied at their control.

* `describe_design()` lists as missing, within each `by` group, only the
  combinations of the levels used in that group, and checks the labels of
  the treatments: levels that differ only in letter case or spaces
  ("control" and "Control") and labels that are not valid UTF-8 text (a file
  read with the wrong encoding) give a warning, and levels of a stressor
  observed in only some `by` groups, which may be one treatment under two
  labels, are pointed out.

* `describe_design()` has been checked on 36 published experiments with very
  different designs, which motivated the changes above.

* Equivalence: `estimate_interactions()` gains `margin` (no default), the
  largest negligible deviation as a share of the response of the reference
  cell. Every interaction now has its relative deviation, (observed −
  expected) / reference, the same quantity for every null model
  (`relative`, `relative_lower`, `relative_upper`); with a margin it is also
  tested for equivalence with two one-sided tests (`p_equivalence`,
  `p_equivalence_adjusted`), each from the distance to the margin on a scale
  on which the variance of a mean does not depend on the mean (the same
  with both `inference` methods). `classify_interactions()` (`verdict`) and
  `summary()` combine both tests into four verdicts: deviation, negligible
  deviation, equivalent and inconclusive, so that "additive" can be a
  claim rather than a failure to detect a deviation.

* `fit_responses()` no longer chooses between the Poisson and the negative
  binomial by AIC for counts with its default candidates: it fits the
  negative binomial and uses the Poisson only if that does not converge, as
  happens when the counts show no overdispersion (`reason` says so). With
  overdispersed counts and four units per combination, the AIC chose the
  Poisson in a fifth of simulated experiments, and the residual checks did
  not exclude it, so the tests rejected true null models too often (7% at a
  nominal 5%); with the new rule they stayed at 5–6%. Candidates given in
  `families` are still chosen by AIC.

* For a `glm()` whose dispersion is estimated (quasi-Poisson,
  quasi-binomial and the other families except Poisson and binomial),
  `estimate_interactions()` now uses a t distribution with its residual
  degrees of freedom, as `summary.glm()` does; emmeans gives such models
  infinite degrees of freedom, which made their tests reject a true null
  model too often with few units. Results for these models change.

* New `compare_null_models()`: for each interaction estimated under two or
  more null models, whether all of them flag a deviation with the same sign
  (`"robust"`), only some (`"model_dependent"`), none (`"none"`) or some in
  opposite directions (`"conflicting"`), with what each null model
  flags. `summary()` counts these cases when there are two or more null
  models.

* The default of `null_models` in `estimate_interactions()` is now `"AD"`
  (it was `c("AD", "MU", "DO")`): in simulated experiments generated without
  interaction by five mechanisms (AD, MU or BL, DO, SAM and CA), AD flagged
  the fewest false interactions on average among the null models that
  two-level factorials can test, under a rule fixed before the simulations
  (see "Default null model" in `?estimate_interactions` and
  `vignette("must")`). `null_models = c("AD", "MU", "DO")` gives the
  previous results. The simulation script, the rates of every scenario and
  the report are installed in `system.file("simulations", package = "must")`.
* `estimate_interactions()` gains the Bliss independence null model
  (`null_models = "BL"`; Bliss 1939): single and combined effects relative
  to the reference cell, the net interaction wK - prod(wi) and, for three
  stressors, the emergent interaction of Beppler et al. (2016), with its
  rescaled value (`rescaled`) and the descriptive class of Beppler et al.
  (`rescaled_class`). As for the other null models, a deviation is judged
  from the interval and the p-value; the cut-offs only describe its size.
  The emergent interaction is `not_defined` for four or more stressors.
* `estimate_interactions()` is now generic and also analyses models fitted
  outside must (`lm()`, `glm()`, `mgcv::gam()`, `lme4::glmer()`,
  `glmmTMB::glmmTMB()` and the other models supported by emmeans), given the
  output of `describe_design()` for their data. The means of the treatment
  combinations are those of emmeans, the degrees of freedom those that
  emmeans gives or those set in the new `df` argument, and
  `diagnostics$means` lists the means used. Quantities that the model fixes
  at zero are `not_applicable`, and means that the model cannot estimate or
  that lie at the boundary of its link make their contrasts `not_estimable`.
  emmeans is a suggested package.
* Two more null models from ecotoxicology: the stress addition model
  (`null_models = "SAM"`; Liess et al. 2016), whose beta distribution of
  stress capacity is set with `sam_shape`, and concentration addition
  (`"CA"`, Loewe additivity), which fits a log-logistic dose-response curve
  to the doses of each stressor tested alone and needs at least two of them.
  Both give the net interaction on the scale relative to the reference cell.
* `classify_interactions()` also classifies the BL, SAM and CA interactions.

# must 0.1.0

First public release.

* `describe_design()` reports, before any model is fitted, the observed and
  missing treatment combinations, the contrast sets the design identifies
  (net deviation and, for three or more stressors, the pure highest-order
  interaction), cells shared by several sets, the replication unit,
  subsamples, repeated measures and unbalanced sampling. With `groups` it
  reports the stressors applied to whole groups of units (blocks, whole
  plots, rooms) and flags the contrast sets confounded with differences
  among groups. `baseline` declares pre-treatment sampling times, and
  `background = "all"` adds contrast sets with the other stressors at their
  non-control levels.
* `fit_responses()` fits a cell-means GLMM with glmmTMB to one or several
  responses and selects the error distribution by convergence, DHARMa
  residual checks and AIC. It handles numbers of successes out of a number
  of trials (`trials`), counts per unit of sampling effort (`offset`), one
  dispersion per `by` group or per treatment combination (`dispersion`) and
  treatment combinations with only zeros (or only successes), whose means
  are fixed at the boundary of the link. Random intercepts for repeated
  units and estimable groups are added automatically, and removed when
  their variance is estimated at zero.
* `estimate_interactions()` computes single, combined and interaction
  effects under the simple addition (AD), multiplicative (MU) and dominance
  (DO) null models, with Krinsky–Robb simulation or the delta method and
  adjusted p-values. By default it corrects for small samples (REML
  covariance and t distributions with the degrees of freedom of the
  experimental units); `correction = "none"` gives the asymptotic procedure
  of Madge Pimentel et al. (2025). Contrasts that involve a combination with
  only zeros are tested from the probability of that outcome under the null
  model, and designs with several sampling times also get the average of
  each effect over the times.
* `classify_interactions()` applies the directional classification of
  Piggott et al. (2015) under AD and MU, extended to sets of three or more
  stressors.
* `summary()` and `plot()` methods for `must_interactions` objects, the
  example data set `polazzo_zooplankton` and `citation("must")`.
