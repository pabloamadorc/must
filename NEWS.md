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
