#' Estimate single, combined and interaction effects under null models
#'
#' Computes, for every identifiable contrast set of the design (see
#' [describe_design()]), the effect of each stressor alone, the effect of
#' their combination and its deviation from what each null model predicts.
#' This is post-estimation inference: one model per response, one or more
#' null models tested from its predictions with the same criteria (Madge
#' Pimentel et al. 2025). The model is either the one selected by
#' [fit_responses()] or one fitted outside must, such as [stats::lm()],
#' [stats::glm()], `mgcv::gam()`, `lme4::glmer()` or [glmmTMB::glmmTMB()]
#' (see "Models fitted outside must"). By default the null model is simple
#' addition (see "Default null model").
#'
#' @param fit A `must_fit` object from [fit_responses()], or a model fitted
#'   outside must whose predictors include the stressor columns of `design`
#'   (any model supported by the emmeans package).
#' @param design For a model fitted outside must: the output of
#'   [describe_design()] for the data the model was fitted to. Its stressors,
#'   control levels, `time`, `by`, `baseline`, `max_order` and `background`
#'   define the contrast sets.
#' @param null_models Null models to evaluate: any of `"AD"` (simple
#'   addition), `"MU"` (multiplicative), `"DO"` (dominance), `"BL"` (Bliss
#'   independence), `"SAM"` (stress addition model) and `"CA"`
#'   (concentration addition, or Loewe additivity). Default `"AD"`, the null
#'   model that flagged the fewest false interactions across mechanisms in
#'   simulations (see "Default null model"); `c("AD", "MU", "DO")` gives the
#'   three null models of Madge Pimentel et al. (2025). It can also be the
#'   output of [suggest_null_models()]: its null models are estimated, and its
#'   criteria and reasons are kept in `settings$null_model_choice`.
#' @param inference `"simulation"` (default; Krinsky–Robb: draws of the fixed
#'   effects from their estimated distribution) or `"delta"` (delta method;
#'   not for `"DO"`).
#' @param correction `"small_sample"` (default) or `"none"`. For a `must_fit`
#'   with `"small_sample"`, the covariance of the fixed effects comes from the
#'   REML refit of the selected model and the sampling distribution is a t
#'   with the degrees of freedom of the experimental units (see Details); for
#'   a model fitted outside must, a t with the degrees of freedom given by
#'   emmeans (residual ones for a glm with an estimated dispersion) or by
#'   `df`. `"none"` uses the normal distribution (and, for a
#'   `must_fit`, the maximum-likelihood covariance), as in Madge Pimentel et
#'   al. (2025); with few replicates its tests reject a true null model too
#'   often.
#' @param n_sim Number of draws for `inference = "simulation"`, and for the
#'   tests of combinations with only zeros with either inference. Default
#'   10000.
#' @param level Confidence level of the intervals. Default 0.95.
#' @param p_adjust Method of [stats::p.adjust()] used across contrast sets
#'   and sampling times (see Details). Default `"BH"`.
#' @param overall If `TRUE` (default) and the design has two or more sampling
#'   times after any `baseline`, adds for each set the average of each
#'   effect over those times (`time = "overall"`).
#' @param sensitivity Whether an increase in the response means sensitivity
#'   (mortality, damage) rather than tolerance. For a `must_fit`: `NULL`
#'   (default) or the names of such responses; for a model fitted outside
#'   must: `TRUE` or `FALSE` (default). When the response is a proportion
#'   (`trials` in [fit_responses()]; a binomial-type family or a logit,
#'   probit, cloglog or cauchit link in a model fitted outside must), every
#'   null model is then applied to its complement (1 - p, e.g. survival
#'   instead of mortality); for other responses the MU, DO and BL rows are
#'   `"not_applicable"`, because those null models assume that an increase
#'   means tolerance.
#' @param seed `NULL` (default) or a whole number. The simulation runs with
#'   this seed through [withr::with_seed()], which restores the random
#'   number state of the session afterwards. With `NULL`, a seed is drawn
#'   from the session's random number stream; either way it is stored in
#'   `settings$seed`, so the estimates can be repeated exactly.
#' @param sam_shape Shape of the symmetric beta distribution of the stress
#'   capacity of individuals in the stress addition model (both of its
#'   parameters). Default 3.2, the value fitted by Liess et al. (2016).
#' @param margin `NULL` (default) or a positive number: the equivalence
#'   margin for the relative deviation of each interaction (`relative`, a
#'   share of the response of the reference cell; 0.1 means 10%). With a
#'   margin, each interaction is also tested for equivalence
#'   (`p_equivalence`), and [classify_interactions()] and `summary()` give
#'   its verdict (see "Equivalence"). It has no default: what deviation is
#'   negligible is a scientific decision, to be declared before the analysis.
#' @param df For a model fitted outside must: `NULL` (default), for the
#'   degrees of freedom that emmeans gives for each contrast, or a single
#'   positive number used for every contrast (for example, experimental units
#'   minus treatment combinations). Only with `correction = "small_sample"`.
#'   emmeans gives infinite degrees of freedom to a [stats::glm()] whose
#'   dispersion is estimated (quasi-Poisson, quasi-binomial and the other
#'   families except Poisson and binomial); with `NULL`, its residual degrees
#'   of freedom are used instead, as [summary.glm()][stats::summary.glm] does.
#'   `MASS::glm.nb()` keeps the normal distribution, as in its summary.
#' @param data For a model fitted outside must: the data it was fitted to,
#'   when emmeans cannot recover them from the model (for example, a
#'   `mgcv::gam()` whose formula transforms a column, as in `factor(day)`).
#'   Default `NULL`.
#' @param ... For a model fitted outside must: further arguments passed to
#'   [emmeans::ref_grid()], such as `at` (values of other numeric
#'   predictors), `offset`, `lmer.df` or `component`. Not used otherwise.
#'
#' @return An object of class `must_interactions`, a list with
#'   `estimates`, `design`, `settings` and `diagnostics`. `estimates` is a
#'   data frame with one row per response, contrast set, sampling time, null
#'   model and quantity, and columns:
#'   * `response`, `context` and the `by` columns, `time` (`NA` without
#'     time, `"overall"` for averages over times), `set_id`, `order`,
#'     `stressors`, `levels`, `background` (levels of the stressors held
#'     outside the set): what is estimated.
#'   * `null_model` (`"AD"`, `"MU"`, `"DO"`, `"BL"`, `"SAM"`, `"CA"`) and
#'     `scale` (`"response"`, `"log(response)"` or `"response / control"`):
#'     what "no interaction" means and on which scale.
#'   * `component`: `"single"` (one stressor alone vs the reference cell;
#'     its name in `stressor`), `"combined"` (all stressors of the set vs the
#'     reference) or `"interaction"`.
#'   * `estimand` for interactions: `"net"` (deviation of the combination
#'     from the null-model combination of the single effects) or `"pure"`
#'     (highest-order interaction once all lower-order ones are removed;
#'     only for three or more stressors; the emergent interaction under BL).
#'   * `estimate`, `se` (standard error; with `inference = "simulation"`,
#'     the standard deviation of the simulated estimates before they are
#'     scaled to a t distribution, comparable with the delta method),
#'     `lower`, `upper`, `level`, `df` (degrees of freedom of the t
#'     distribution; `Inf` for the normal), `p_value`, `p_adjusted`, `mc_se`
#'     (Monte Carlo standard error of a simulated p-value).
#'   * `observed` and `expected` (interactions only, response scale): the
#'     fitted mean of the combination and the value the null model predicts.
#'   * `relative`, `relative_lower` and `relative_upper` (interactions only):
#'     the relative deviation (observed - expected) / reference, the same
#'     quantity for every null model, and its interval at `level` (see
#'     "Equivalence").
#'   * `p_equivalence` and `p_equivalence_adjusted` (interactions, with a
#'     `margin` only): p-value of the test of equivalence, raw and adjusted
#'     like `p_adjusted`.
#'   * `rescaled` and `rescaled_class` (BL interactions only): the
#'     interaction relative to its extremes and its descriptive class (see
#'     "Bliss independence").
#'   * `inference`, `status` (`"ok"`, `"boundary"` (tested from the
#'     probability of a combination with only zeros, see Details),
#'     `"not_estimable"`, `"not_identifiable"`, `"not_defined"` or
#'     `"not_applicable"`), `reason`, `note` (`"expected_below_zero"` when
#'     simple addition predicts an impossible negative value for a
#'     non-negative response; `"no_point_expectation"` when dominance gives no
#'     single expected value; `"reference_not_positive"` when the relative
#'     deviation cannot be computed because the reference mean (or one of its
#'     draws) is not positive; `"equivalence_not_tested"` when the closest
#'     means at the margin could not be found (see "Equivalence");
#'     `"confounded_with_groups"` for the sets flagged
#'     in `sets$confounded` by [describe_design()]; `"complement"` for
#'     responses analysed as 1 - p; several notes are separated by `"; "`).
#'
#'   `settings` holds the arguments used (with the seed and whether it was
#'   given), the versions of must and R (and emmeans, for a model fitted
#'   outside must), the criteria behind the null models when they came from
#'   [suggest_null_models()] and, for a `must_fit`, a summary of the models
#'   (`fit`): what [report_interactions()] describes.
#'
#'   For a model fitted outside must, `diagnostics$means` has the mean of
#'   every treatment combination (and time) that was used: its estimate and
#'   standard error on the link scale, the mean on the response scale, the
#'   degrees of freedom given by emmeans and, if it could not be used, why
#'   (`"nonestimable"`, `"no_se"` or `"boundary"`).
#'
#' @details
#' With y0 the mean of the reference cell (the control, or the background
#' combination of the other stressors), yi that of stressor i alone and yS
#' that of the combination of the stressors in S:
#'
#' **Two stressors** (Madge Pimentel et al. 2025, eqs. 1–3 and 7–9):
#' \deqn{AD:\ y_{AB} - y_A - y_B + y_0}{AD: yAB - yA - yB + y0}
#' \deqn{MU:\ \ln y_{AB} - \ln y_A - \ln y_B + \ln y_0}{MU: ln yAB - ln yA - ln yB + ln y0}
#' \deqn{DO:\ \min(y_{AB}, y_0) - \min(y_A, y_B)}{DO: min(yAB, y0) - min(yA, yB)}
#'
#' **k stressors.** The net deviation compares the combination with the
#' null-model combination of the single effects only:
#' \eqn{y_K - \sum_i y_i + (k-1) y_0}{yK - sum(yi) + (k - 1) y0} (AD) or the
#' same on the log scale (MU). The pure k-th order interaction is the
#' inclusion–exclusion sum \eqn{\sum_{S \subseteq K} (-1)^{k-|S|} f(y_S)}{sum
#' over subsets S of K of (-1)^(k - |S|) f(yS)}, with f the identity (AD) or
#' the logarithm (MU). The net deviation equals the sum of the pure
#' interactions of all subsets of K with two or more stressors: a
#' non-additive net deviation is not, by itself, evidence of a k-way
#' interaction. The dominance model is defined for two stressors only.
#'
#' **Single and combined effects** are differences from the reference cell
#' on the response scale (AD) or log ratios (MU). For DO only the
#' interaction is reported; its single effects are those of AD.
#'
#' **Inference.** `"simulation"` follows the Krinsky–Robb method used by
#' Madge Pimentel et al. (2025): percentile intervals and two-sided p-values
#' \eqn{2 \min(P(x \le 0), P(x \ge 0))}{2 min(P(x <= 0), P(x >= 0))} of the
#' simulated estimates, bounded below by 1/`n_sim`. `"delta"` uses Wald
#' tests; it is not valid for DO, which is not differentiable where the
#' ranking of the means changes. Uncertainty from the selection of the error
#' distribution is not included.
#'
#' **Small-sample correction.** With few replicates, a normal distribution
#' with the maximum-likelihood covariance understates the uncertainty: in
#' simulated 2 x 2 Gaussian experiments with three replicates the AD test then
#' rejected a true null model in about 13% of the experiments instead of 5%.
#' With
#' `correction = "small_sample"` the covariance comes from the REML refit and
#' the draws come from a multivariate t (normal draws scaled by
#' \eqn{\sqrt{\nu / \chi^2_\nu}}{sqrt(df / chisq(df))}), whose degrees of
#' freedom follow the strata of the design: observations minus estimated
#' means without repeated measures; experimental units minus the treatment
#' combinations they are nested in with repeated measures or subsamples
#' (combination and time when each unit is sampled at a single time); minus
#' one per extra block of each `groups` column with block structure in the
#' model; and, for groups of units to which some stressors are applied whole
#' (whole plots), the groups minus the combinations of those stressors form
#' their own stratum, which bounds the degrees of freedom of the sets that
#' involve them and is subtracted from those of the other sets. With
#' `dispersion = "cell"` in [fit_responses()] the degrees of freedom of each
#' contrast follow the Welch–Satterthwaite approximation (Satterthwaite
#' 1946). Families without a dispersion parameter (Poisson, binomial) and
#' without random effects keep the normal distribution, as in [stats::glm()].
#' For Gaussian models without random effects the corrected intervals of AD
#' contrasts are exact t intervals; with repeated measures, and for the
#' interactions between a stressor applied to whole groups and one applied
#' within them, the rule is conservative.
#'
#' **Combinations with only zeros.** When the full combination of a set has
#' only zeros at a time (or, for binomial responses, only zeros or only
#' successes) and the model selected by [fit_responses()] has a log or logit
#' link, the null model still predicts a mean E for it from the other cells.
#' The p-value is twice the probability of observing only zeros if the mean
#' were E, averaged over the draws of E (`status = "boundary"`; Poisson, NB2,
#' Tweedie, binomial and beta-binomial models). Random intercepts are
#' integrated out: observations that share a level of a random effect share
#' one intercept with the variance of all random terms, which is exact with a
#' single grouping (units or rooms, say) and conservative with nested or
#' crossed ones. The same test is used for a single stressor whose cell has
#' only zeros (against the reference cell). The estimate is computed with
#' the observed mean (zero), which makes MU estimates `-Inf`; intervals are
#' not given. If the reference cell, or a single-stressor cell that an
#' interaction needs, has only zeros, that contrast cannot be measured and
#' the row is `"not_estimable"`.
#'
#' **Multiple testing.** p-values are adjusted within each response, null
#' model, component, stressor and estimand, across contrast sets, sampling
#' times and `by` groups; averages over times (`"overall"`) are adjusted
#' separately.
#'
#' @section Bliss independence:
#' Under Bliss independence (`"BL"`; Bliss 1939) stressors act
#' independently, so the responses relative to the reference cell,
#' \eqn{w_S = y_S / y_0}{wS = yS / y0} (relative fitness or survival),
#' multiply. Single and combined effects are \eqn{w - 1}{w - 1}, the net
#' interaction is \eqn{w_K - \prod_i w_i}{wK - prod(wi)} (the DA measure of
#' Beppler et al. 2016, eqs. 2.1–2.2; Tekin et al. 2018) and, for three
#' stressors, the pure interaction is the emergent interaction E3 of Beppler
#' et al. (2016, eq. 2.4):
#' \deqn{w_{XYZ} - w_X w_{YZ} - w_Y w_{XZ} - w_Z w_{XY} + 2 w_X w_Y w_Z,}{wXYZ
#' - wX wYZ - wY wXZ - wZ wXY + 2 wX wY wZ,}
#' the part of the three-way deviation that the pairwise interactions do not
#' explain. It is not defined here for four or more stressors
#' (`status = "not_defined"`). The net BL interaction tests the same
#' hypothesis as the net MU interaction on another scale (a ratio of the
#' combined response to its expectation under MU, a difference under BL);
#' the emergent interaction is a different hypothesis from the pure MU
#' interaction, which removes the pairwise interactions on the log scale.
#'
#' `rescaled` divides each BL interaction by the absolute value of the same
#' expression with the combined response replaced by zero, when the
#' interaction is negative (lethal synergy), or by the smallest
#' single-stressor response, when it is positive (complete buffering), so
#' that -1 and 1 mark those two extremes (Segrè et al. 2005; Beppler et al.
#' 2016, section 2.2.3). With a non-negative response it is at least -1; when
#' no stressor increases the response, it exceeds 1 if the combination is
#' less affected than its most affected single stressor. `rescaled_class`
#' labels values below -0.5 as `"synergy"`, above 0.5 as `"antagonism"` and
#' the rest as `"additive"`, the conservative cut-offs of Beppler et al.
#' (2016, section 2.2.4). Both
#' describe the size of the interaction only. Whether it deviates from
#' Bliss independence is decided, as for the other null models, by its
#' interval and p-value, because fixed cut-offs ignore the uncertainty of the
#' estimate (Madge Pimentel et al. 2025).
#'
#' @section Stress addition and concentration addition:
#' **Stress addition (SAM).** Liess et al. (2016) assume that each individual
#' tolerates a general stress drawn from a beta distribution with both shape
#' parameters equal to `sam_shape`. The relative response of each stressor
#' alone, \eqn{w_i}{wi} (survival relative to the reference cell), becomes
#' the general stress whose quantile in that distribution is the mortality
#' \eqn{1 - w_i}{1 - wi}; the general stresses of the stressors add up, and
#' the predicted relative response of the combination is the share of
#' individuals that tolerate their sum, zero when it reaches 1 (Liess et al.
#' 2016, eqs. 2-5). The interaction is \eqn{w_K}{wK} minus that prediction;
#' single and combined effects are relative, as under BL. A stressor that
#' increases the response exerts no stress; if one does so at the point
#' estimate, the interaction is `"not_applicable"`. The model predicts larger
#' combined effects than Bliss independence, so combinations that Bliss
#' independence calls additive can be antagonistic under it.
#'
#' **Concentration addition (CA, Loewe additivity).** Each stressor acts as a
#' dilution of the others: a combination of doses \eqn{x_i}{xi} has the
#' relative response w at which
#' \eqn{\sum_i x_i / EC_i(w) = 1}{sum(xi / ECi(w)) = 1}, where
#' \eqn{EC_i(w)}{ECi(w)} is the dose of stressor i alone that gives w (Loewe
#' 1953; Tang et al. 2015, eq. 2). Doses are the differences between the
#' levels and the control, so the stressor columns must hold numbers.
#' \eqn{EC_i}{ECi} comes from a two-parameter log-logistic curve,
#' \eqn{w = 1 / (1 + (x / EC_{50})^b)}{w = 1 / (1 + (x / EC50)^b)}, fitted by
#' least squares to the relative responses at every dose of the stressor
#' tested alone, with the other stressors at the reference. At least two such
#' doses are needed, so CA is `"not_defined"` in two-level factorials. The
#' curves are refitted in each simulation draw, so their uncertainty enters
#' the intervals. A curve that does not decrease with dose makes the
#' interaction `"not_applicable"`; in single draws it is taken as no dose
#' equivalent. Doses tested alone with only zeros keep their means at the
#' boundary and still inform the curves.
#'
#' Both models give only the net interaction (the pure interaction of three or
#' more stressors is `"not_defined"`), and their single and combined effects
#' are those of BL.
#'
#' @section Equivalence:
#' A non-significant interaction is a failure to detect a deviation, not
#' evidence of additivity. To claim that a combination follows a null model,
#' declare before the analysis the largest deviation that would be
#' negligible, `margin`, as a share of the response of the reference cell.
#' Each interaction is then also tested against the hypotheses that its
#' relative deviation d = (observed - expected) / reference is at most
#' -`margin` and at least `margin`, with two one-sided tests (Schuirmann
#' 1987; Lakens et al. 2018): `p_equivalence` is the larger of their
#' p-values.
#'
#' Each one-sided test measures how far the estimated means of the cells
#' involved are from the closest means at which d equals the bound, in
#' standard errors (the Mahalanobis distance, signed), and compares that
#' distance with a t distribution with the degrees of freedom of the
#' interaction (normal if they are infinite). The distance is measured on a
#' scale on which the variance of a mean does not depend on the mean
#' (Anscombe 1948): the square root for Poisson, quasi-Poisson and other
#' counts whose variance grows in proportion to the mean, the inverse
#' hyperbolic sine of the square root for the negative binomial, the arcsine
#' of the square root for proportions, the logarithm for Gamma and lognormal
#' responses and the response itself for Gaussian ones (for other families,
#' and for models of a transformed response, the scale of the model). When
#' the bounds are linear in the means (AD with a constant variance), the
#' tests are t tests of observed - expected ± `margin` × reference, as in
#' the t test of a ratio to a control (Hauschke et al. 1999). Tests built
#' from the interval of d, or from its draws, declared false equivalences too
#' often with overdispersed counts, where the standard error of a mean
#' shrinks with the mean. The equivalence test uses neither the draws nor
#' the gradient, so it is the same with both `inference` methods; near the
#' margin, its verdict can differ from whether the interval of d at 1 - 2
#' alpha lies within the margin. If the closest means cannot be found,
#' `p_equivalence` is `NA` with the note `"equivalence_not_tested"`.
#'
#' The relative deviation puts every null model on the same scale: under BL,
#' SAM and CA it is the interaction itself; under AD and DO, the interaction
#' divided by the reference; under MU, the observed minus the multiplicative
#' expectation, divided by the reference (the net BL interaction). Its
#' interval includes the uncertainty of the reference. When the stressors
#' raise the response far above the reference (nutrients, say), relative
#' deviations can be large. Combinations with only zeros get a relative
#' deviation but no interval and no test of equivalence.
#'
#' [classify_interactions()] and `summary()` combine both tests, with the
#' adjusted p-values, into four verdicts: `"deviation"` (deviation detected,
#' not equivalent), `"negligible"` (a deviation detected, but within the
#' margin), `"equivalent"` (no deviation detected and within the margin) and
#' `"inconclusive"` (neither). Only `"equivalent"` supports a claim that the
#' stressors combine as the null model predicts.
#'
#' In simulated 2 x 2 experiments with four units per combination, the test
#' declared false equivalence (a true relative deviation exactly at the
#' margin of ±0.25) in 4–6% of the experiments, close to the nominal 5%, for
#' Gaussian responses under AD, survival under BL and overdispersed counts
#' under MU analysed with a quasi-Poisson [stats::glm()] or with the
#' negative binomial of [fit_responses()] (3–6% for counts with no, moderate
#' or strong overdispersion). A Poisson model fitted to the same overdispersed counts
#' ignores their overdispersion and gave 21%, which is why
#' [fit_responses()] uses it only when the negative binomial does not
#' converge; `MASS::glm.nb()`, whose maximum-likelihood fit with the normal
#' distribution underestimates the variance with so few units, gave 10%.
#' Equivalence needs replication: with four units, a margin of
#' ±0.1 was never reached; one of ±0.25 was reached in about half (Gaussian)
#' to three quarters (survival) of the experiments that followed the null
#' model exactly, but rarely (under 10%) with overdispersed counts.
#'
#' @section Default null model:
#' When there is a hypothesis about how the stressors act, name the null model
#' that represents it in `null_models` (Madge Pimentel et al. 2025); without
#' one, AD is used. It was chosen by simulation, with a rule fixed before the
#' results were seen: among the null models defined in two-level factorials
#' and with at least 80% of the highest average power against interactions
#' that no mechanism explains, the one with the fewest false interactions on
#' average over five mechanisms of no interaction (AD, MU or BL, DO, SAM and
#' CA). The experiments were 2 x 2 factorials and 3 x 3 dose factorials with
#' four units per combination, survival out of 20 animals or overdispersed
#' counts, and moderate or weak single effects, analysed with a cell-means
#' [stats::glm()] (100 experiments per scenario; p < 0.05 without
#' adjustment).
#'
#' Each null model flagged interactions at about the nominal rate (3–5%)
#' when the data followed its own mechanism, but none was robust to every
#' mechanism: data generated by stress addition looked synergistic to
#' all the others (in 60% of the experiments under AD, 76% under MU and BL
#' and 86% under DO), and data generated by the others looked antagonistic to
#' stress addition. On average, AD flagged false interactions in 20% of the
#' experiments, MU and BL in 24%, SAM in 27% and DO in 32%, with average
#' powers of 58%, 67%, 38% and 61%. In the dose design CA flagged fewer (16%,
#' against 23% under AD), with less power (51% against 63%). AD was also
#' chosen within each design, each size of the single effects and with
#' survival as the response. Without the stress addition mechanism, MU and BL
#' came within one percentage point of AD and, being more powerful, would
#' have been chosen; with counts only, the power of AD fell below 80% of
#' theirs, and they would have been chosen too. Requiring every null model
#' to agree (`agreement = "robust"` in [compare_null_models()]) flagged under
#' 3% false interactions, but only 15% of the interactions that no mechanism
#' explains.
#'
#' The simulation is installed with the package:
#' `system.file("simulations", package = "must")` holds the script
#' (`null_model_robustness.R`), which reruns or extends it, the rates of every
#' scenario (`null_model_robustness.csv`) and its report
#' (`null_model_robustness.md`).
#'
#' @section Models fitted outside must:
#' Any model supported by emmeans (Lenth & Piaskowski 2026) can be analysed:
#' [stats::lm()], [stats::glm()], `mgcv::gam()`, `lme4::lmer()`,
#' `lme4::glmer()`, [glmmTMB::glmmTMB()] and many more. The stressor, `time`
#' and `by` columns of `design` must be predictors of the model; each of their
#' values is a level, also for a numeric predictor, which is evaluated at
#' each value observed. The mean of each treatment combination (and time)
#' is the estimated marginal mean of emmeans: the prediction on the link
#' scale with random effects at zero, averaged with equal weights over the
#' other factors of the model and with other numeric predictors at their
#' mean (unless set in `...`), back-transformed with the inverse link or
#' the transformation of the response without bias adjustment (with
#' `lm(log(y) ~ ...)` they are medians of y, not means). The draws and the
#' delta method use the covariance of the coefficients of the model as it
#' was fitted (no REML refit).
#'
#' With `correction = "small_sample"` the degrees of freedom are those
#' emmeans gives for the linear combination of coefficients of each contrast
#' (its gradient for non-linear contrasts): the residual degrees of freedom
#' for [stats::lm()] and [stats::glm()] with an estimated dispersion and for
#' `mgcv::gam()` with any family (observations minus effective degrees of
#' freedom), Kenward–Roger or Satterthwaite for `lme4::lmer()` when pbkrtest
#' or lmerTest is installed, and infinite for Poisson and binomial models
#' fitted with [stats::glm()] and for most other mixed models. With few
#' experimental units, infinite degrees of freedom give tests that reject a
#' true null model too often: `df` sets them by hand, or [fit_responses()]
#' applies the correction described above.
#'
#' The estimates inherit the structure of the model. A quantity the model
#' fixes at zero, such as the MU interaction of a log-link model without
#' the interaction term, is `"not_applicable"`. Means the model cannot
#' estimate (in rank-deficient models), or that lie at the boundary of the
#' link function (a combination with only zeros under a log link, or only
#' zeros or only successes under a logit link) make their contrasts
#' `"not_estimable"`: such models give no usable standard error for them,
#' whereas [fit_responses()] tests them from the probability of that
#' outcome.
#'
#' @section Assumptions:
#' * Means are conditional predictions (random effects at zero). With a log
#'   link and a random intercept, marginal means differ by a constant
#'   factor: MU and BL estimates and all p-values are unchanged, AD and DO
#'   estimates change scale.
#' * MU, DO, BL, SAM and CA only represent their mechanistic assumptions when
#'   an increase in the response means tolerance (survival, abundance) rather
#'   than sensitivity (mortality); declare the latter in `sensitivity`. AD
#'   may predict impossible values below zero, in which case a deviation
#'   from it does not imply a mechanistic interaction (Madge Pimentel et al.
#'   2025).
#' * A non-significant interaction is a failure to detect a deviation, not
#'   evidence of additivity; equivalence needs a `margin` (see
#'   "Equivalence").
#'
#' @section Missing values and zeros:
#' MU needs positive means and BL a positive mean in the reference cell: if a
#' fitted or simulated mean is not positive (possible with an identity link),
#' the MU or BL rows of that set are `"not_estimable"`. Treatment
#' combinations without values of the response (at a time) make the sets
#' that need them `"not_estimable"` (at that time). Times declared as
#' `baseline` in [describe_design()] are not estimated.
#'
#' @section Reproducibility:
#' With `inference = "simulation"` results depend on the random number
#' stream: the seed used is in `settings$seed`, and `mc_se` gives the Monte
#' Carlo error of each simulated p-value. The same seed may not give
#' identical draws in other versions of R or with other linear algebra
#' libraries.
#'
#' @references
#' Beppler, C., Tekin, E., Mao, Z., White, C., McDiarmid, C., Vargas, E.,
#' Miller, J. H., Savage, V. M., & Yeh, P. J. (2016). Uncovering emergent
#' interactions in three-way combinations of stressors. *Journal of the
#' Royal Society Interface*, 13, 20160800. \doi{10.1098/rsif.2016.0800}
#'
#' Bliss, C. I. (1939). The toxicity of poisons applied jointly. *Annals of
#' Applied Biology*, 26(3), 585–615.
#' \doi{10.1111/j.1744-7348.1939.tb06990.x}
#'
#' Krinsky, I., & Robb, A. L. (1986). On approximating the statistical
#' properties of elasticities. *The Review of Economics and Statistics*,
#' 68(4), 715–719. \doi{10.2307/1924536}
#'
#' Lenth, R. V., & Piaskowski, J. (2026). emmeans: Estimated marginal means,
#' aka least-squares means. R package version 2.0.4.
#' \url{https://CRAN.R-project.org/package=emmeans}
#'
#' Liess, M., Foit, K., Knillmann, S., Schäfer, R. B., & Liess, H.-D. (2016).
#' Predicting the synergy of multiple stress effects. *Scientific Reports*,
#' 6, 32965. \doi{10.1038/srep32965}
#'
#' Loewe, S. (1953). The problem of synergism and antagonism of combined
#' drugs. *Arzneimittelforschung*, 3, 285–290.
#'
#' Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
#' J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025).
#' Hypothesis-driven research on multiple stressors: An analytical framework
#' for stressor interactions. *Ecology and Evolution*, 15, e71959.
#' \doi{10.1002/ece3.71959}
#'
#' Anscombe, F. J. (1948). The transformation of Poisson, binomial and
#' negative-binomial data. *Biometrika*, 35(3–4), 246–254.
#' \doi{10.1093/biomet/35.3-4.246}
#'
#' Hauschke, D., Kieser, M., & Hothorn, L. A. (1999). Proof of safety in
#' toxicology based on the ratio of two means for normally distributed data.
#' *Biometrical Journal*, 41(3), 295–304.
#' \doi{10.1002/(SICI)1521-4036(199906)41:3<295::AID-BIMJ295>3.0.CO;2-2}
#'
#' Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
#' for psychological research: A tutorial. *Advances in Methods and
#' Practices in Psychological Science*, 1(2), 259–269.
#' \doi{10.1177/2515245918770963}
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin*, 2(6), 110–114.
#' \doi{10.2307/3002019}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#' procedure and the power approach for assessing the equivalence of average
#' bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics*,
#' 15(6), 657–680. \doi{10.1007/BF01068419}
#'
#' Segrè, D., DeLuna, A., Church, G. M., & Kishony, R. (2005). Modular
#' epistasis in yeast metabolism. *Nature Genetics*, 37, 77–83.
#' \doi{10.1038/ng1489}
#'
#' Tang, J., Wennerberg, K., & Aittokallio, T. (2015). What is synergy? The
#' Saariselkä agreement revisited. *Frontiers in Pharmacology*, 6, 181.
#' \doi{10.3389/fphar.2015.00181}
#'
#' Tekin, E., Yeh, P. J., & Savage, V. M. (2018). General form for
#' interaction measures and framework for deriving higher-order emergent
#' effects. *Frontiers in Ecology and Evolution*, 6, 166.
#' \doi{10.3389/fevo.2018.00166}
#'
#' @examples
#' jars <- expand.grid(
#'   salinity = c("ambient", "high"),
#'   warming = c("no", "yes"),
#'   replicate = 1:5,
#'   stringsAsFactors = FALSE
#' )
#' jars$jar <- seq_len(nrow(jars))
#' mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
#'   ifelse(jars$warming == "yes", 0.5, 1)
#' jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
#' design <- describe_design(jars, c("salinity", "warming"),
#'   control = c(salinity = "ambient", warming = "no"), unit = "jar")
#' fits <- fit_responses(jars, "algae", design, residual_checks = FALSE,
#'   seed = 1, quiet = TRUE)
#' # Simple addition (the default)
#' est <- estimate_interactions(fits, n_sim = 2000, seed = 1)
#' est
#' est$estimates[est$estimates$component == "interaction",
#'   c("null_model", "estimate", "lower", "upper", "df", "p_adjusted")]
#' # The three null models and the procedure of Madge Pimentel et al. (2025),
#' # without the correction
#' estimate_interactions(fits, null_models = c("AD", "MU", "DO"), n_sim = 2000,
#'   seed = 1, correction = "none")
#' # Additive, or only not detected? An equivalence margin of 10% of the control
#' eq <- estimate_interactions(fits, null_models = c("AD", "MU"), n_sim = 2000, seed = 1,
#'   margin = 0.1)
#' eq$estimates[eq$estimates$component == "interaction",
#'   c("null_model", "relative", "relative_lower", "relative_upper", "p_equivalence")]
#' classify_interactions(eq)[, c("null_model", "class", "verdict")]
#' # Bliss independence, with its rescaled interaction
#' bl <- estimate_interactions(fits, null_models = "BL", n_sim = 2000, seed = 1)
#' bl$estimates[bl$estimates$component == "interaction",
#'   c("estimate", "lower", "upper", "p_value", "rescaled", "rescaled_class")]
#' @examplesIf requireNamespace("emmeans", quietly = TRUE)
#' # A model fitted outside must
#' model <- stats::glm(algae ~ salinity * warming, family = stats::poisson(), data = jars)
#' estimate_interactions(model, design, null_models = c("AD", "MU", "BL"), n_sim = 2000,
#'   seed = 1)
#' @export
estimate_interactions <- function(fit, ...) {
  UseMethod("estimate_interactions")
}

#' @rdname estimate_interactions
#' @export
estimate_interactions.must_fit <- function(fit, null_models = "AD",
                                           inference = "simulation",
                                           correction = "small_sample", n_sim = 10000L,
                                           level = 0.95, p_adjust = "BH", overall = TRUE,
                                           sensitivity = NULL, seed = NULL, sam_shape = 3.2,
                                           margin = NULL, ...) {
  call <- generic_call("estimate_interactions")
  check_dots_empty(..., call = call)
  opts <- check_estimate_options(null_models, inference, correction, n_sim, level, p_adjust,
                                 overall, seed, sam_shape, margin, call = call)
  sensitivity <- check_sensitivity(sensitivity, names(fit$fits), call = call)
  check_choice_sensitivity(opts$null_model_choice, sensitivity, call = call)
  seed_given <- !is.null(seed)
  opts$seed <- resolve_seed(seed)
  run <- function() {
    lapply(fit$fits, function(f) {
      estimate_one_response(f, design = fit$design, opts = opts,
                            sensitive = f$response %in% sensitivity)
    })
  }
  parts <- withr::with_seed(opts$seed, run())
  failed <- vapply(fit$fits, function(f) !identical(f$status, "ok"), logical(1))
  new_interactions(
    parts, fit$design, opts, sensitivity,
    settings = list(
      prediction = "conditional (random effects at zero)",
      time_structure = fit$settings$time_structure,
      families = vapply(fit$fits, function(f) {
        if (identical(f$status, "ok")) paste0(f$family, " (", f$link, ")") else NA_character_
      }, character(1)),
      fit = fit_summary(fit),
      versions = analysis_versions(),
      seed_given = seed_given
    ),
    diagnostics = list(
      failed_responses = data.frame(
        response = names(fit$fits)[failed],
        reason = vapply(fit$fits[failed], `[[`, character(1), "reason"),
        stringsAsFactors = FALSE
      ),
      reml = vapply(fit$fits, function(f) !is.null(f$model_reml), logical(1))
    )
  )
}

# Checks of the options shared by every method. The seed is checked here and
# drawn by the caller after its own checks, so that an input error does not
# advance the random number stream.
check_estimate_options <- function(null_models, inference, correction, n_sim, level,
                                   p_adjust, overall, seed, sam_shape, margin = NULL,
                                   call = sys.call(-1)) {
  # The output of suggest_null_models(): its null models, and its criteria and
  # reasons for the settings.
  choice <- NULL
  if (inherits(null_models, "must_null_models")) {
    choice <- null_model_choice(null_models)
    null_models <- null_models$null_models
  }
  null_models <- check_null_models(null_models, call = call)
  check_choice(inference, c("simulation", "delta"), call = call)
  check_choice(correction, c("small_sample", "none"), call = call)
  if (identical(inference, "delta") && "DO" %in% null_models) {
    abort_input(paste(
      "The delta method is not valid for the dominance null model, which is",
      "not differentiable where the ranking of the means changes (Madge",
      "Pimentel et al. 2025). Use inference = \"simulation\" or drop \"DO\"",
      if (is.null(choice)) {
        "from `null_models`."
      } else {
        "from the null models suggested, as in setdiff(<suggestion>$null_models, \"DO\")."
      }
    ), call = call)
  }
  check_number(n_sim, lower = 100, call = call)
  if (n_sim != round(n_sim)) {
    abort_input("`n_sim` must be a whole number.", call = call)
  }
  check_number(level, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE,
               call = call)
  check_choice(p_adjust, stats::p.adjust.methods, call = call)
  check_flag(overall, call = call)
  check_seed(seed, call = call)
  check_number(sam_shape, lower = 0, include_lower = FALSE, call = call)
  if (!is.null(margin)) {
    check_number(margin, lower = 0, include_lower = FALSE, call = call)
    if (!is.finite(margin)) {
      abort_input("`margin` must be a finite number.", call = call)
    }
  }
  list(null_models = null_models, inference = inference, correction = correction,
       n_sim = as.integer(n_sim), level = level, p_adjust = p_adjust, overall = overall,
       seed = NULL, sam_shape = sam_shape, margin = margin, null_model_choice = choice)
}

# The output object, from the rows estimated for each response.
new_interactions <- function(parts, design, opts, sensitivity, settings, diagnostics) {
  estimates <- bind_rows_base(parts, template = empty_estimates(design))
  estimates <- adjust_p_values(estimates, opts$p_adjust)
  structure(
    list(
      estimates = estimates,
      design = design,
      settings = c(
        opts[c("null_models", "inference", "correction", "n_sim", "level", "p_adjust",
               "overall")],
        list(sensitivity = sensitivity, seed = opts$seed, sam_shape = opts$sam_shape,
             margin = opts$margin, null_model_choice = opts$null_model_choice),
        settings
      ),
      diagnostics = diagnostics
    ),
    class = "must_interactions"
  )
}

check_null_models <- function(null_models, call = sys.call(-1)) {
  valid <- c("AD", "MU", "DO", "BL", "SAM", "CA")
  if (!is.character(null_models) || length(null_models) == 0L || anyNA(null_models) ||
      !all(null_models %in% valid) || anyDuplicated(null_models) > 0L) {
    abort_input(
      sprintf(
        "`null_models` must contain one or more of %s, without repetitions.",
        format_items(quote_items(valid), max = Inf)
      ),
      call = call
    )
  }
  valid[valid %in% null_models]
}

check_sensitivity <- function(sensitivity, responses, call = sys.call(-1)) {
  if (is.null(sensitivity)) {
    return(character())
  }
  if (!is.character(sensitivity) || anyNA(sensitivity) || !all(sensitivity %in% responses)) {
    abort_input(
      sprintf(
        "`sensitivity` must be NULL or names of fitted responses (%s).",
        format_items(quote_items(responses), max = Inf)
      ),
      call = call
    )
  }
  unique(sensitivity)
}

empty_estimates <- function(design) {
  out <- data.frame(
    response = character(), context = character(), stringsAsFactors = FALSE
  )
  for (b in design$spec$by) {
    out[[b]] <- character()
  }
  cbind(out, data.frame(
    time = character(), set_id = character(), order = integer(),
    stressors = character(), levels = character(), background = character(),
    null_model = character(), scale = character(), component = character(),
    stressor = character(), estimand = character(), estimate = numeric(),
    se = numeric(), lower = numeric(), upper = numeric(), level = numeric(),
    df = numeric(), p_value = numeric(), mc_se = numeric(), observed = numeric(),
    expected = numeric(), relative = numeric(), relative_lower = numeric(),
    relative_upper = numeric(), p_equivalence = numeric(), rescaled = numeric(),
    rescaled_class = character(), inference = character(), status = character(),
    reason = character(), note = character(), stringsAsFactors = FALSE
  ))
}

estimate_one_response <- function(f, design, opts, sensitive) {
  if (!identical(f$status, "ok")) {
    return(NULL)
  }
  model <- if (identical(opts$correction, "small_sample") && !is.null(f$model_reml)) {
    f$model_reml
  } else {
    f$model
  }
  complement <- sensitive && identical(f$kind, "binomial")
  pred <- prepare_predictions(f, model, opts$n_sim, complement = complement)
  estimate_sets(f, pred, design, opts, sensitive, baseline = baseline_times(design$spec),
                set_df = function(s) {
                  set_df(f$df_info, s, opts$correction,
                         contexts = unique(set_cell_contexts(design$set_cells[[s$set_id]], s)))
                })
}

# Estimates of every identifiable contrast set of one response, from the
# predictions `pred` at the cells of f$grid. `set_df(s)` gives the degrees of
# freedom of a set, or NA when each contrast gets its own (from
# pred$contrast_df()); `baseline` gives, for each stressor, the times before
# it was applied, at which the sets that involve it are skipped.
estimate_sets <- function(f, pred, design, opts, sensitive, baseline, set_df) {
  key <- paste_key(f$grid$.context, f$grid$.cell, f$grid$.time)
  sets <- design$sets[design$sets$net_identifiable, , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(sets))) {
    s <- sets[i, , drop = FALSE]
    sc <- design$set_cells[[s$set_id]]
    sc_context <- set_cell_contexts(sc, s)
    times <- unique(as.character(f$grid$.time[f$grid$.context == s$context]))
    times <- times[!times %in% set_baseline(baseline, design, s)]
    df <- set_df(s)
    per_time <- lapply(times, function(tt) {
      idx <- match(paste_key(sc_context, sc$cell, tt), key)
      set_estimates(s, sc, idx, tt, pred, opts, f, design, df, sensitive, key)
    })
    rows <- c(rows, lapply(per_time, `[[`, "rows"))
    if (opts$overall && length(times) >= 2L) {
      rows[[length(rows) + 1L]] <- overall_estimates(s, per_time, pred, opts$inference,
                                                     opts$level, f, design,
                                                     margin = opts$margin)
    }
  }
  bind_rows_base(rows)
}

# Context in which each cell of a set was observed: its own `by` group or the
# group of rows shared with it (designs made before `shared` have no column).
set_cell_contexts <- function(sc, s) {
  if (is.null(sc$context)) rep(s$context, nrow(sc)) else sc$context
}

# Times at which a set is not estimated: the baseline times of the stressors
# away from control in its combination, those of the set and those of its
# background. `baseline` is a list by stressor (or a vector for all).
set_baseline <- function(baseline, design, s) {
  if (is.null(baseline)) {
    return(character())
  }
  if (!is.list(baseline)) {
    return(as.character(baseline))
  }
  unique(unlist(baseline[set_stressors_away(design, s)], use.names = FALSE))
}

# Stressors away from control in the combination of a set.
set_stressors_away <- function(design, s) {
  spec <- design$spec
  row <- design$cells[design$cells$context == s$context & design$cells$cell == s$levels, ,
                      drop = FALSE]
  if (nrow(row) == 0L) {
    return(strsplit(s$stressors, " x ", fixed = TRUE)[[1L]])
  }
  spec$stressors[vapply(spec$stressors, function(x) {
    !identical(row[[x]][[1L]], spec$control[[x]])
  }, logical(1))]
}

estimate_row <- function(s, f, design, time_value, null_model, sp, values) {
  row <- data.frame(response = f$response, context = s$context, stringsAsFactors = FALSE)
  for (b in design$spec$by) {
    row[[b]] <- s[[b]]
  }
  cbind(row, data.frame(
    time = as.character(time_value), set_id = s$set_id, order = s$order,
    stressors = s$stressors, levels = s$levels, background = s$background,
    null_model = null_model, scale = null_model_scale(null_model),
    component = sp$component, stressor = sp$stressor, estimand = sp$estimand,
    estimate = values$estimate, se = values$se, lower = values$lower,
    upper = values$upper, level = values$level, df = values$df,
    p_value = values$p_value, mc_se = values$mc_se, observed = values$observed,
    expected = values$expected, relative = values$relative,
    relative_lower = values$relative_lower, relative_upper = values$relative_upper,
    p_equivalence = values$p_equivalence, rescaled = values$rescaled,
    rescaled_class = values$rescaled_class, inference = values$inference, status = values$status,
    reason = values$reason, note = values$note, stringsAsFactors = FALSE
  ))
}

empty_values <- function(status, reason, level, inference) {
  list(estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
       level = level, df = NA_real_, p_value = NA_real_, mc_se = NA_real_,
       observed = NA_real_, expected = NA_real_, relative = NA_real_,
       relative_lower = NA_real_, relative_upper = NA_real_, p_equivalence = NA_real_,
       rescaled = NA_real_, rescaled_class = NA_character_, inference = inference,
       status = status, reason = reason, note = "", draws = NULL, grad = NULL)
}

set_estimates <- function(s, sc, idx, time_value, pred, opts, f, design, df, sensitive, key) {
  inference <- opts$inference
  level <- opts$level
  m <- s$order
  in_set <- strsplit(s$stressors, " x ", fixed = TRUE)[[1L]]
  w <- contrast_weights(m)
  missing_now <- sc$cell[sc$for_net & is.na(idx)]
  has_time <- !is.na(time_value)
  control <- design$spec$control
  ca <- if ("CA" %in% opts$null_models) {
    ca_layout(s, sc, idx, time_value, design, key, in_set)
  } else {
    NULL
  }
  rows <- list()
  parts <- list()
  for (nm in opts$null_models) {
    for (sp in component_specs(nm, m, in_set, w)) {
      if (identical(nm, "SAM")) {
        sp$shape <- opts$sam_shape
      }
      ca_interaction <- identical(nm, "CA") && identical(sp$estimand, "net")
      if (ca_interaction && isTRUE(ca$ok)) {
        sp$cols <- ca$cols
        sp$labels <- ca$labels
        sp$ca <- ca$ca
      }
      values <- if (identical(nm, "DO") && m > 2L) {
        empty_values("not_defined", "The dominance null model is defined for two stressors only.",
                     level, inference)
      } else if (identical(sp$estimand, "pure") && nm %in% c("SAM", "CA")) {
        empty_values("not_defined", sprintf(
          paste(
            "%s gives no pure (emergent) interaction of three or more stressors:",
            "only the net deviation is estimated."
          ),
          if (identical(nm, "SAM")) "The stress addition model" else "Concentration addition"
        ), level, inference)
      } else if (identical(sp$estimand, "pure") && !isTRUE(s$pure_identifiable)) {
        lacking <- strsplit(s$missing_for_pure, "; ", fixed = TRUE)[[1L]]
        empty_values("not_identifiable", sprintf(
          paste(
            "The design lacks lower-order combinations (%s): the pure",
            "interaction cannot be separated from them."
          ),
          format_items(quote_items(short_cells(lacking, control)), max = 3L)
        ), level, inference)
      } else if (identical(nm, "BL") && identical(sp$estimand, "pure") && m > 3L) {
        empty_values("not_defined", paste(
          "The emergent interaction under Bliss independence is defined for up",
          "to three stressors (Beppler et al. 2016; Tekin et al. 2018)."
        ), level, inference)
      } else if (sensitive && !identical(f$kind, "binomial") &&
                 nm %in% c("MU", "DO", relative_null_models())) {
        empty_values("not_applicable", paste(
          "An increase in this response means sensitivity: the multiplicative,",
          "dominance, Bliss, stress addition and concentration addition null",
          "models assume that an increase means tolerance (Madge Pimentel et al.",
          "2025). Analyse a response where higher means better, such as survival."
        ), level, inference)
      } else if (length(missing_now) > 0L) {
        empty_values("not_estimable", sprintf(
          if (has_time) "No values of the response at this time in %s." else "No values of the response in %s.",
          format_items(quote_items(short_cells(missing_now, control)), max = 3L)
        ), level, inference)
      } else if (identical(sp$estimand, "pure") && anyNA(idx)) {
        empty_values("not_estimable", paste0(
          "Some lower-order combinations have no values of the response",
          if (has_time) " at this time." else "."
        ), level, inference)
      } else if (ca_interaction && !isTRUE(ca$ok)) {
        empty_values(ca$status, ca$reason, level, inference)
      } else if (identical(inference, "simulation") && nzchar(pred$sim_problem())) {
        empty_values("not_estimable", pred$sim_problem(), level, inference)
      } else {
        compute_component(nm, sp, idx, m, pred, inference, level, f, df, sc, design,
                          margin = opts$margin)
      }
      if (isTRUE(s$confounded)) {
        values$note <- add_note(values$note, "confounded_with_groups")
      }
      if (pred$complement) {
        values$note <- add_note(values$note, "complement")
      }
      rows[[length(rows) + 1L]] <- estimate_row(s, f, design, time_value, nm, sp, values)
      parts[[component_key(nm, sp)]] <- values
    }
  }
  list(rows = bind_rows_base(rows), parts = parts)
}

component_key <- function(null_model, sp) {
  paste(null_model, sp$component, sp$stressor, sp$estimand, sep = "|")
}

add_note <- function(note, extra) {
  if (nzchar(note)) paste(note, extra, sep = "; ") else extra
}

# Which quantities are reported for each null model.
component_specs <- function(null_model, m, in_set, w) {
  interaction <- list(list(component = "interaction", stressor = NA_character_,
                           estimand = "net", weights = w$net))
  if (identical(null_model, "DO")) {
    return(interaction)
  }
  singles <- lapply(seq_len(m), function(i) {
    list(component = "single", stressor = in_set[[i]], estimand = NA_character_,
         weights = w$single[[i]])
  })
  combined <- list(list(component = "combined", stressor = NA_character_,
                        estimand = NA_character_, weights = w$combined))
  pure <- if (m >= 3L) {
    list(list(component = "interaction", stressor = NA_character_, estimand = "pure",
              weights = w$pure))
  } else {
    list()
  }
  c(singles, combined, interaction, pure)
}

# Cells that concentration addition needs for a set at one time: the
# reference cell, every dose of each stressor of the set tested alone (the
# other stressors at the reference) and the full combination. Doses are the
# differences between the levels and the control. Returns the grid columns
# (`cols`), their labels and, in `ca`, the positions of each stressor's doses
# in `cols`, their logarithms and the doses of the combination; or a status
# and reason if concentration addition cannot be computed.
ca_layout <- function(s, sc, idx, time_value, design, key, in_set) {
  spec <- design$spec
  fail <- function(status, reason) list(ok = FALSE, status = status, reason = reason)
  # The cells of the set's `by` group and of the rows shared with it, its own
  # first.
  pool <- pool_of(context_pools(design$cells, spec), s$context)
  cells <- design$cells[design$cells$context %in% pool, , drop = FALSE]
  cells <- cells[order(match(cells$context, pool), method = "radix"), , drop = FALSE]
  reference <- cells[match(sc$cell[[1L]], cells$cell), spec$stressors, drop = FALSE]
  full <- cells[match(sc$cell[[nrow(sc)]], cells$cell), spec$stressors, drop = FALSE]
  if (nrow(reference) != 1L || anyNA(reference) || nrow(full) != 1L || anyNA(full)) {
    return(fail("not_estimable", "The reference cell or the combination is not in the design."))
  }
  in_background <- nzchar(s$background)
  cols <- idx[[1L]]
  labels <- sc$cell[[1L]]
  curve_cols <- list()
  u <- list()
  x <- numeric(length(in_set))
  for (i in seq_along(in_set)) {
    st <- in_set[[i]]
    lv <- spec$levels[[st]]
    num <- suppressWarnings(as.numeric(lv))
    ctrl <- suppressWarnings(as.numeric(spec$control[[st]]))
    if (anyNA(num) || is.na(ctrl)) {
      return(fail("not_defined", sprintf(
        "Concentration addition needs numeric doses: the levels of '%s' are not numbers.", st
      )))
    }
    dose <- num - ctrl
    x[[i]] <- as.numeric(full[[st]]) - ctrl
    if (!isTRUE(x[[i]] > 0)) {
      return(fail("not_applicable", sprintf(
        "Concentration addition needs doses above the control; '%s' is below it here.", st
      )))
    }
    alone <- lv[dose > 0]
    lab <- vapply(alone, function(l) {
      r <- reference
      r[[st]] <- l
      make_cell_label(r, spec$stressors)
    }, character(1), USE.NAMES = FALSE)
    n_design <- sum(lab %in% cells$cell)
    if (n_design < 2L) {
      return(fail("not_defined", sprintf(
        paste(
          "Concentration addition needs a dose-response curve for each stressor:",
          "at least two doses of '%s' alone%s; the design has %d."
        ),
        st, if (in_background) " in this background" else "", n_design
      )))
    }
    pos <- match(paste_key(cells$context[match(lab, cells$cell)], lab, time_value), key)
    usable <- !is.na(pos)
    if (sum(usable) < 2L) {
      return(fail("not_estimable", sprintf(
        "Fewer than two doses of '%s' alone have values of the response%s.",
        st, if (is.na(time_value)) "" else " at this time"
      )))
    }
    curve_cols[[i]] <- length(cols) + seq_len(sum(usable))
    cols <- c(cols, pos[usable])
    labels <- c(labels, lab[usable])
    u[[i]] <- log(dose[dose > 0][usable])
  }
  list(
    ok = TRUE,
    cols = c(cols, idx[[length(idx)]]),
    labels = c(labels, sc$cell[[length(idx)]]),
    ca = list(curve_cols = curve_cols, u = u, x = x, stressors = in_set)
  )
}

compute_component <- function(null_model, sp, idx, m, pred, inference, level, f, df, sc,
                              design, margin = NULL) {
  use <- if (identical(null_model, "DO")) 1:4 else which(sp$weights != 0)
  # The concentration-addition interaction also uses the doses of each
  # stressor tested alone (sp$cols, see ca_layout()).
  cols <- if (is.null(sp$cols)) idx[use] else sp$cols
  labels <- if (is.null(sp$labels)) sc$cell[use] else sp$labels
  problem <- pred$cell_problem[cols]
  if (any(nzchar(problem))) {
    return(empty_values("not_estimable", cell_problem_reason(problem, labels, design),
                        level, inference))
  }
  at_boundary <- pred$at_boundary[cols]
  if (!is.null(sp$ca)) {
    # Doses alone with only zeros keep their means fixed at the boundary and
    # still inform the curves; the reference and the combination do not.
    at_boundary[-c(1L, length(cols))] <- FALSE
  }
  if (any(at_boundary)) {
    return(boundary_component(null_model, sp, idx, use, m, pred, level, f, df, sc, design))
  }
  mu <- pred$mu[cols]
  fun <- contrast_function(null_model, sp, m, use)
  weights <- if (identical(null_model, "DO")) dominance_weights(mu) else sp$weights[use]
  # Gradient of the contrast with respect to the linear predictor of its cells.
  gradient <- function() {
    if (null_model %in% relative_null_models()) {
      numeric_gradient(fun, mu) * pred$dmu[cols]
    } else {
      delta_gradient(weights, mu, pred$dmu[cols], if (identical(null_model, "MU")) "MU" else "AD")
    }
  }
  g0 <- gradient()
  if (fixed_by_model(g0, pred$V_eta[cols, cols, drop = FALSE])) {
    return(empty_values("not_applicable", paste(
      "The model fixes this quantity at zero (it has no term for it, such as an",
      "interaction term on the scale of the link function), so it is not tested.",
      "Fit a model that includes it."
    ), level, inference))
  }
  if (is.na(df)) {
    df <- pred$contrast_df(g0, cols)
  }
  if (is.na(df) || df < 1) {
    return(empty_values(
      "not_estimable",
      paste(
        "No degrees of freedom are left to estimate the variance of this",
        "contrast (too few replicated units)."
      ),
      level, inference
    ))
  }
  draws <- if (identical(inference, "simulation")) pred$mu_draws(cols, df) else NULL
  # The standard error is taken from the normal draws (before they are scaled
  # to a t distribution, whose standard deviation is larger than its scale and
  # infinite with two or fewer degrees of freedom).
  normal <- if (!is.null(draws) && is.finite(df)) pred$mu_draws(cols, Inf) else draws
  positive <- switch(null_model, MU = seq_along(cols), BL = , SAM = , CA = 1L, integer())
  if (length(positive) > 0L) {
    bad_point <- any(!is.finite(mu[positive]) | mu[positive] <= 0)
    bad_sims <- !is.null(draws) && any(draws[, positive] <= 0)
    if (bad_point || bad_sims) {
      return(empty_values(
        "not_estimable",
        if (identical(null_model, "MU")) {
          "The multiplicative model needs positive means; some fitted or simulated means are not positive."
        } else {
          sprintf(
            paste(
              "The %s needs a positive mean in the reference cell; some fitted or",
              "simulated values are not positive."
            ),
            null_model_name(null_model)
          )
        },
        level, inference
      ))
    }
  }
  unsuited <- relative_assumption(null_model, sp, mu, m, labels, design)
  if (!is.null(unsuited)) {
    return(empty_values("not_applicable", unsuited, level, inference))
  }
  estimate <- fun(matrix(mu, nrow = 1L))
  grad <- NULL
  component_draws <- NULL
  normal_draws <- NULL
  if (identical(inference, "simulation")) {
    component_draws <- fun(draws)
    normal_draws <- fun(normal)
    stats_list <- summarise_draws(component_draws, level)
    stats_list$se <- stats::sd(normal_draws)
  } else {
    g <- gradient()
    se <- sqrt(max(drop(t(g) %*% pred$V_eta[cols, cols, drop = FALSE] %*% g), 0))
    stats_list <- summarise_wald(estimate, se, df, level)
    grad <- numeric(length(pred$mu))
    grad[cols] <- g
  }
  observed <- NA_real_
  expected <- NA_real_
  rescaled <- NA_real_
  note <- ""
  if (identical(sp$component, "interaction")) {
    observed <- pred$mu[[idx[[length(idx)]]]]
    y0 <- pred$mu[[idx[[1L]]]]
    if (null_model %in% relative_null_models()) {
      # Interaction = observed / y0 - expected / y0.
      expected <- observed - y0 * estimate
    } else if (identical(sp$estimand, "pure")) {
      expected <- if (identical(null_model, "MU")) observed / exp(estimate) else observed - estimate
    } else {
      expected <- expected_combined(pred$mu[idx[seq_len(m + 1L)]], m, null_model)
      if (identical(null_model, "AD") && isTRUE(f$non_negative) && is.finite(expected) &&
          expected < 0) {
        note <- "expected_below_zero"
      }
      if (identical(null_model, "DO") && is.na(expected)) {
        note <- "no_point_expectation"
      }
    }
    if (identical(null_model, "BL")) {
      rescaled <- bliss_rescaled(estimate, full_means(mu, use, length(idx)), m)
    }
  }
  rel <- if (identical(sp$component, "interaction")) {
    relative_deviation(null_model, sp, m, use, fun, mu, draws, inference, level, df, pred,
                       cols, margin)
  } else {
    list(relative = NA_real_, relative_lower = NA_real_, relative_upper = NA_real_,
         p_equivalence = NA_real_, note = "")
  }
  if (nzchar(rel$note)) {
    note <- add_note(note, rel$note)
  }
  list(
    estimate = estimate, se = stats_list$se, lower = stats_list$lower,
    upper = stats_list$upper, level = level, df = df, p_value = stats_list$p_value,
    mc_se = if (identical(inference, "simulation")) mc_error(stats_list$p_value, pred$n_sim) else NA_real_,
    observed = observed, expected = expected, relative = rel$relative,
    relative_lower = rel$relative_lower, relative_upper = rel$relative_upper,
    p_equivalence = rel$p_equivalence, rescaled = rescaled,
    rescaled_class = bliss_class(rescaled), inference = inference, status = "ok",
    reason = "", note = note, draws = component_draws, normal_draws = normal_draws,
    grad = grad, cols = cols, fun = fun, null_model = null_model,
    rel_draws = rel$rel_draws, rel_grad = rel$rel_grad, rel_fun = rel$rel_fun
  )
}

# Relative deviation of an interaction, (observed - expected) / reference: the
# deviation from the null model as a share of the response of the reference
# cell, on the same scale for every null model. Under BL, SAM and CA it is
# the interaction itself; under AD and DO, the interaction divided by the
# reference; under MU, the observed minus the multiplicative expectation,
# divided by the reference. Returns it with its interval (from the draws or
# the gradient, like the interaction) and, with a `margin`, the p-value of
# the equivalence test (see equivalence_p_signed()), and what is needed to
# average it over times.
relative_deviation <- function(null_model, sp, m, use, fun, mu, draws, inference, level, df,
                               pred, cols, margin) {
  out <- list(relative = NA_real_, relative_lower = NA_real_, relative_upper = NA_real_,
              p_equivalence = NA_real_, rel_draws = NULL, rel_grad = NULL, rel_fun = NULL,
              note = "")
  # The reference is the first cell of every contrast (see compute_component()).
  if (!is.finite(mu[[1L]]) || mu[[1L]] <= 0 || (!is.null(draws) && any(draws[, 1L] <= 0))) {
    out$note <- "reference_not_positive"
    return(out)
  }
  rel_fun <- relative_function(null_model, sp, m, use, fun)
  out$rel_fun <- rel_fun
  out$relative <- rel_fun(matrix(mu, nrow = 1L))
  if (identical(inference, "simulation")) {
    d <- rel_fun(draws)
    s <- summarise_draws(d, level)
    out$rel_draws <- d
  } else {
    g <- numeric_gradient(rel_fun, mu) * pred$dmu[cols]
    se <- sqrt(max(drop(t(g) %*% pred$V_eta[cols, cols, drop = FALSE] %*% g), 0))
    s <- summarise_wald(out$relative, se, df, level)
    grad <- numeric(length(pred$mu))
    grad[cols] <- g
    out$rel_grad <- grad
  }
  if (!is.null(margin)) {
    out$p_equivalence <- equivalence_p_signed(rel_fun, pred, cols, df, margin)
    if (is.na(out$p_equivalence)) {
      out$note <- "equivalence_not_tested"
    }
  }
  out$relative_lower <- s$lower
  out$relative_upper <- s$upper
  out
}

# Function from the means of the cells of a contrast (one row per evaluation)
# to its relative deviation.
relative_function <- function(null_model, sp, m, use, fun) {
  if (null_model %in% relative_null_models()) {
    return(fun)
  }
  n <- length(sp$weights)
  function(D) {
    D <- matrix(D, ncol = length(use))
    if (identical(null_model, "MU")) {
      M <- full_means(D, use, n)
      return((M[, n] - expected_full(M, "interaction", sp$estimand, "MU", m, sp$weights)) / M[, 1L])
    }
    fun(D) / D[, 1L]
  }
}

# Function that maps a matrix of means of the cells `use` of a contrast (one
# row per evaluation: the point estimate or a simulation draw) to the values
# of the contrast.
contrast_function <- function(null_model, sp, m, use) {
  if (!is.null(sp$ca)) {
    return(function(D) ca_contrast(D, sp))
  }
  switch(
    null_model,
    DO = function(D) dominance_interaction(D),
    BL = ,
    CA = function(D) bliss_contrast(full_means(D, use, length(sp$weights)), m, sp),
    SAM = function(D) sam_contrast(full_means(D, use, length(sp$weights)), m, sp),
    function(D) drop(transform_means(D, null_model) %*% sp$weights[use])
  )
}

# Name of a null model in messages.
null_model_name <- function(null_model) {
  switch(null_model, AD = "simple addition model", MU = "multiplicative model",
         DO = "dominance model", BL = "Bliss independence model",
         SAM = "stress addition model", CA = "concentration addition model")
}

# NULL if the point estimates suit the assumption of a null model on the
# relative scale, or the reason why they do not: the stress addition model
# needs stressors that reduce the response (Liess et al. 2016), and
# concentration addition, single-stressor curves that decrease with dose.
relative_assumption <- function(null_model, sp, mu, m, labels, design) {
  if (!identical(sp$component, "interaction")) {
    return(NULL)
  }
  if (identical(null_model, "SAM")) {
    W <- mu[1L + seq_len(m)] / mu[[1L]]
    raising <- labels[1L + seq_len(m)][W > 1]
    if (length(raising) > 0L) {
      return(sprintf(
        paste(
          "The stress addition model needs stressors that reduce the response: %s",
          "increase%s it, so it exerts no general stress (Liess et al. 2016)."
        ),
        format_items(quote_items(short_cells(raising, design$spec$control)), max = 3L),
        if (length(raising) == 1L) "s" else ""
      ))
    }
  }
  if (!is.null(sp$ca)) {
    curves <- ca_curves(matrix(mu / mu[[1L]], nrow = 1L), sp)
    flat <- vapply(curves, function(cv) !(cv[1L, 2L] < 0), logical(1))
    if (any(flat)) {
      stressors <- sp$ca$stressors
      return(sprintf(
        paste(
          "Concentration addition needs dose-response curves that decrease with dose;",
          "the one fitted to %s alone does not."
        ),
        format_items(quote_items(stressors[flat]))
      ))
    }
  }
  NULL
}

# Matrix of means with one column per subset of the set (NA for the subsets
# that a contrast does not use) from the means of the columns `use`.
full_means <- function(D, use, n_subsets) {
  D <- matrix(D, ncol = length(use))
  M <- matrix(NA_real_, nrow = nrow(D), ncol = n_subsets)
  M[, use] <- D
  M
}

# Central-difference gradient of fun() (applied to a one-row matrix) at `x`.
numeric_gradient <- function(fun, x) {
  vapply(seq_along(x), function(j) {
    h <- 1e-6 * max(abs(x[[j]]), 1e-8)
    up <- x
    down <- x
    up[[j]] <- up[[j]] + h
    down[[j]] <- down[[j]] - h
    (fun(matrix(up, nrow = 1L)) - fun(matrix(down, nrow = 1L))) / (2 * h)
  }, numeric(1))
}

# Local weights of the dominance interaction min(yAB, y0) - min(yA, yB) at the
# point estimates: which cells it depends on.
dominance_weights <- function(mu) {
  c(as.numeric(mu[[4L]] >= mu[[1L]]), -as.numeric(mu[[2L]] <= mu[[3L]]),
    -as.numeric(mu[[2L]] > mu[[3L]]), as.numeric(mu[[4L]] < mu[[1L]]))
}

# TRUE when a contrast with gradient `g` (with respect to the linear predictor
# of its cells, whose covariance is `V`) has no variance, relative to the
# variance contributions of its cells: the model constrains it, as a model
# without an interaction term does with the interaction on the link scale.
fixed_by_model <- function(g, V) {
  if (!all(is.finite(g)) || !all(is.finite(V))) {
    return(FALSE)
  }
  scale <- sum(g^2 * diag(V))
  if (!isTRUE(scale > 0)) {
    return(FALSE)
  }
  drop(t(g) %*% V %*% g) <= 1e-10 * scale
}

# Reason for a contrast that needs a mean the model cannot provide (see
# external_problems()); `cells` are the labels of the cells of the contrast.
cell_problem_reason <- function(problem, cells, design) {
  first <- problem[nzchar(problem)][[1L]]
  where <- format_items(quote_items(short_cells(cells[problem == first],
                                                design$spec$control)), max = 3L)
  switch(
    first,
    nonestimable = sprintf(
      "The model cannot estimate the mean of %s from its terms (it is not estimable).", where
    ),
    no_se = sprintf("The model gives no standard error for the mean of %s.", where),
    boundary = sprintf(
      paste(
        "The mean of %s is at the boundary of the link function (only zeros, or only",
        "successes) and the model gives no usable standard error for it.",
        "fit_responses() tests such combinations from the probability of that outcome."
      ),
      where
    ),
    abort_internal(sprintf("unknown problem '%s' with a mean.", first))
  )
}

# Monte Carlo standard error of a two-sided simulated p-value p = 2q, where q
# is a proportion of n_sim draws.
mc_error <- function(p, n_sim) {
  if (is.na(p)) NA_real_ else sqrt(p * (2 - p) / n_sim)
}

# Contrast that uses a combination at the boundary (only zeros, or only
# successes): see R/boundary.R. The cell tested is the one away from the
# reference: the single stressor for single effects, the full combination
# otherwise.
boundary_component <- function(null_model, sp, idx, use, m, pred, level, f, df, sc, design) {
  inference <- "boundary test"
  # The cells of the contrast: those of the set or, for the concentration-
  # addition interaction, the reference, the doses alone and the combination.
  ca <- !is.null(sp$ca)
  cells <- if (ca) sp$cols else idx
  cell_labels <- if (ca) sp$labels else sc$cell
  used <- if (ca) seq_along(cells) else use
  target <- if (identical(sp$component, "single")) which(sp$weights == 1) else length(cells)
  boundary_now <- pred$at_boundary[cells[used]]
  if (ca) {
    boundary_now[-c(1L, length(used))] <- FALSE
  }
  at <- used[boundary_now]
  control <- design$spec$control
  if (!identical(as.integer(at), as.integer(target))) {
    lacking <- cell_labels[setdiff(at, target)]
    return(empty_values("not_estimable", sprintf(
      paste(
        "%s had only zeros: with the response eliminated in the reference cell",
        "or by a stressor alone, this contrast cannot be measured."
      ),
      format_items(quote_items(short_cells(lacking, control)), max = 3L)
    ), level, inference))
  }
  target_col <- cells[[target]]
  rows <- f$boundary_rows[as.character(f$boundary_rows$.cell_time) ==
                            as.character(pred$grid$.cell_time[[target_col]]), , drop = FALSE]
  params <- boundary_parameters(pred$model, rows)
  if (!params$family %in% boundary_families()) {
    return(empty_values("not_estimable", sprintf(
      paste(
        "A combination had only zeros, to which the %s family gives no",
        "probability: its contrasts cannot be tested."
      ),
      params$family
    ), level, inference))
  }
  type <- pred$boundary_type[[target_col]]
  observed_raw <- if (identical(type, "zero")) 0 else 1
  observed <- if (pred$complement) 1 - observed_raw else observed_raw
  others <- setdiff(seq_along(cells), target)
  other_weights <- if (identical(null_model, "DO") || ca) {
    rep(1, length(others))
  } else {
    sp$weights[others]
  }
  df_used <- df
  if (is.na(df_used)) {
    scale <- if (null_model %in% c("MU", relative_null_models())) "MU" else "AD"
    g <- delta_gradient(other_weights, pred$mu[cells[others]], pred$dmu[cells[others]], scale)
    df_used <- pred$contrast_df(g, cells[others])
  }
  if (is.na(df_used) || df_used < 1) {
    return(empty_values(
      "not_estimable",
      paste(
        "No degrees of freedom are left to estimate the variance of the expected",
        "mean of this contrast (too few replicated units)."
      ),
      level, inference
    ))
  }
  draws_other <- pred$mu_draws(cells[others], df_used)
  if (is.null(draws_other)) {
    return(empty_values("not_estimable", pred$sim_problem(), level, inference))
  }
  expected_of <- function(M) {
    if (ca) {
      M[, 1L] * ca_expected(ca_curves(M / M[, 1L], sp), sp$ca$x)
    } else {
      expected_full(M, sp$component, sp$estimand, null_model, m, sp$weights, shape = sp$shape)
    }
  }
  M <- matrix(NA_real_, nrow = nrow(draws_other), ncol = length(cells))
  M[, others] <- draws_other
  expected_draws <- expected_of(M)
  mu_point <- pred$mu[cells]
  mu_point[[target]] <- NA_real_
  expected <- expected_of(matrix(mu_point, nrow = 1L))
  # The probability is computed on the scale of the data (not the complement).
  mean_data <- if (pred$complement) 1 - expected_draws else expected_draws
  prob <- prob_all_boundary(mean_data, rows, params$family, params$link, params$disp,
                            params$power, params$sd_re, type, factors = params$factors)
  ok <- is.finite(prob)
  p_value <- min(1, 2 * mean(prob[ok]))
  mu_obs <- pred$mu[cells]
  mu_obs[[target]] <- observed
  estimate <- if (ca) {
    ca_contrast(matrix(mu_obs, nrow = 1L), sp)
  } else {
    contrast_function(null_model, sp, m, use)(matrix(mu_obs[use], nrow = 1L))
  }
  is_interaction <- identical(sp$component, "interaction")
  rescaled <- if (is_interaction && identical(null_model, "BL")) {
    bliss_rescaled(estimate, matrix(mu_obs, nrow = 1L), m)
  } else {
    NA_real_
  }
  # Relative deviation from the observed boundary value; it has no interval
  # and is not tested for equivalence.
  relative <- NA_real_
  reference <- pred$mu[[cells[[1L]]]]
  if (is_interaction && is.finite(reference) && reference > 0) {
    relative <- if (null_model %in% relative_null_models()) {
      estimate
    } else if (identical(null_model, "MU")) {
      (observed - expected[[1L]]) / reference
    } else {
      estimate / reference
    }
  }
  list(
    estimate = estimate, se = NA_real_, lower = NA_real_, upper = NA_real_, level = level,
    df = df_used, p_value = p_value,
    mc_se = if (p_value < 1) 2 * stats::sd(prob[ok]) / sqrt(sum(ok)) else NA_real_,
    observed = if (is_interaction) observed else NA_real_,
    expected = if (is_interaction) expected[[1L]] else NA_real_,
    relative = relative, relative_lower = NA_real_, relative_upper = NA_real_,
    p_equivalence = NA_real_,
    rescaled = rescaled, rescaled_class = bliss_class(rescaled),
    inference = inference, status = "boundary",
    reason = sprintf(
      paste(
        "%s had only %s: the p-value is twice the probability of that outcome if",
        "the null model held, averaged over the uncertainty of the expected mean;",
        "no interval is given."
      ),
      quote_items(short_cells(cell_labels[[target]], control)),
      if (identical(type, "zero")) "zeros" else "successes"
    ),
    note = "", draws = NULL, grad = NULL
  )
}

# Averages of each effect over the sampling times of a set.
overall_estimates <- function(s, per_time, pred, inference, level, f, design,
                              margin = NULL) {
  keys <- names(per_time[[1L]]$parts)
  first_rows <- per_time[[1L]]$rows
  out <- list()
  for (k in seq_along(keys)) {
    parts <- lapply(per_time, function(p) p$parts[[keys[[k]]]])
    statuses <- vapply(parts, `[[`, character(1), "status")
    sp <- list(
      component = first_rows$component[[k]],
      stressor = first_rows$stressor[[k]],
      estimand = first_rows$estimand[[k]]
    )
    null_model <- first_rows$null_model[[k]]
    values <- if (all(statuses == "ok")) {
      combine_over_times(parts, pred, inference, level, margin)
    } else if (length(unique(statuses)) == 1L &&
               statuses[[1L]] %in% c("not_defined", "not_identifiable", "not_applicable")) {
      empty_values(statuses[[1L]], parts[[1L]]$reason, level, inference)
    } else {
      empty_values("not_estimable", "Not estimable (with an interval) at every sampling time.",
                   level, inference)
    }
    if (isTRUE(s$confounded)) {
      values$note <- add_note(values$note, "confounded_with_groups")
    }
    if (pred$complement) {
      values$note <- add_note(values$note, "complement")
    }
    out[[k]] <- estimate_row(s, f, design, "overall", null_model, sp, values)
  }
  bind_rows_base(out)
}

combine_over_times <- function(parts, pred, inference, level, margin = NULL) {
  estimate <- mean(vapply(parts, `[[`, numeric(1), "estimate"))
  dfs <- vapply(parts, `[[`, numeric(1), "df")
  df <- min(dfs)
  if (identical(inference, "simulation")) {
    by_time <- if (length(unique(dfs)) == 1L) {
      lapply(parts, `[[`, "draws")
    } else {
      # Times with their own degrees of freedom (one dispersion per cell):
      # all times are redrawn with the smallest, so that the average follows
      # the t distribution it reports.
      lapply(parts, function(p) p$fun(pred$mu_draws(p$cols, df)))
    }
    draws <- rowMeans(do.call(cbind, by_time))
    stats_list <- summarise_draws(draws, level)
    stats_list$se <- stats::sd(rowMeans(do.call(cbind, lapply(parts, `[[`, "normal_draws"))))
    mc <- mc_error(stats_list$p_value, pred$n_sim)
  } else {
    g <- Reduce(`+`, lapply(parts, `[[`, "grad")) / length(parts)
    se <- sqrt(max(drop(t(g) %*% pred$V_eta %*% g), 0))
    stats_list <- summarise_wald(estimate, se, df, level)
    mc <- NA_real_
  }
  rel <- combine_relative_over_times(parts, pred, inference, level, df, dfs, margin)
  list(estimate = estimate, se = stats_list$se, lower = stats_list$lower,
       upper = stats_list$upper, level = level, df = df, p_value = stats_list$p_value,
       mc_se = mc, observed = NA_real_, expected = NA_real_, relative = rel$relative,
       relative_lower = rel$relative_lower, relative_upper = rel$relative_upper,
       p_equivalence = rel$p_equivalence, rescaled = NA_real_,
       rescaled_class = NA_character_, inference = inference,
       status = "ok", reason = "", note = rel$note, draws = NULL, grad = NULL)
}

# Average of the relative deviations of an interaction over times, with its
# interval from the same draws (or gradients) as the interaction itself, and
# its equivalence test on the means of the cells of every time.
combine_relative_over_times <- function(parts, pred, inference, level, df, dfs, margin) {
  out <- list(relative = NA_real_, relative_lower = NA_real_, relative_upper = NA_real_,
              p_equivalence = NA_real_, note = "")
  funs <- lapply(parts, `[[`, "rel_fun")
  if (any(vapply(funs, is.null, logical(1)))) {
    if (any(vapply(parts, function(p) isTRUE(grepl("reference_not_positive", p$note,
                                                   fixed = TRUE)), logical(1)))) {
      out$note <- "reference_not_positive"
    }
    return(out)
  }
  out$relative <- mean(vapply(parts, `[[`, numeric(1), "relative"))
  if (identical(inference, "simulation")) {
    by_time <- if (length(unique(dfs)) == 1L) {
      lapply(parts, `[[`, "rel_draws")
    } else {
      lapply(parts, function(p) p$rel_fun(pred$mu_draws(p$cols, df)))
    }
    d <- rowMeans(do.call(cbind, by_time))
    s <- summarise_draws(d, level)
  } else {
    g <- Reduce(`+`, lapply(parts, `[[`, "rel_grad")) / length(parts)
    se <- sqrt(max(drop(t(g) %*% pred$V_eta %*% g), 0))
    s <- summarise_wald(out$relative, se, df, level)
  }
  if (!is.null(margin)) {
    cols <- sort(unique(unlist(lapply(parts, `[[`, "cols"))))
    pos <- lapply(parts, function(p) match(p$cols, cols))
    average <- function(M) {
      M <- matrix(M, ncol = length(cols))
      Reduce(`+`, lapply(seq_along(parts), function(k) {
        funs[[k]](M[, pos[[k]], drop = FALSE])
      })) / length(parts)
    }
    out$p_equivalence <- equivalence_p_signed(average, pred, cols, df, margin)
    if (is.na(out$p_equivalence)) {
      out$note <- "equivalence_not_tested"
    }
  }
  out$relative_lower <- s$lower
  out$relative_upper <- s$upper
  out
}

# Adjusted p-values of the tests of deviation (`p_adjusted`) and of
# equivalence (`p_equivalence_adjusted`), within the same families, each
# placed after its raw p-value.
adjust_p_values <- function(estimates, method) {
  estimates$p_adjusted <- rep(NA_real_, nrow(estimates))
  estimates$p_equivalence_adjusted <- rep(NA_real_, nrow(estimates))
  if (nrow(estimates) > 0L) {
    group <- paste_key(estimates$response, estimates$null_model, estimates$component,
                       estimates$stressor, estimates$estimand,
                       estimates$time %in% "overall")
    tested <- estimates$status %in% c("ok", "boundary")
    for (pair in list(c("p_value", "p_adjusted"),
                      c("p_equivalence", "p_equivalence_adjusted"))) {
      ok <- tested & is.finite(estimates[[pair[[1L]]]])
      for (g in unique(group[ok])) {
        sel <- ok & group == g
        estimates[[pair[[2L]]]][sel] <- stats::p.adjust(estimates[[pair[[1L]]]][sel],
                                                        method = method)
      }
    }
  }
  place_after <- function(df, column, after) {
    rest <- setdiff(names(df), column)
    at <- match(after, rest)
    df[, c(rest[seq_len(at)], column, rest[-seq_len(at)]), drop = FALSE]
  }
  estimates <- place_after(estimates, "p_adjusted", "p_value")
  place_after(estimates, "p_equivalence_adjusted", "p_equivalence")
}

# Distribution of the tests, for print methods: the small-sample correction
# of fit_responses() models, or the degrees of freedom of a model fitted
# outside must (`settings$model`).
correction_label <- function(s, detail = FALSE) {
  if (!identical(s$correction, "small_sample")) {
    return("no small-sample correction")
  }
  if (is.null(s$model)) {
    return(if (detail) "small-sample correction (REML, t)" else "small-sample correction")
  }
  if (is.numeric(s$df)) {
    sprintf("t with %s degrees of freedom (`df`)", format(s$df))
  } else if (identical(s$df, "residual")) {
    "t with the residual degrees of freedom of the glm (estimated dispersion)"
  } else {
    "degrees of freedom from emmeans"
  }
}

#' @export
print.must_interactions <- function(x, n = 12L, ...) {
  s <- x$settings
  cat(sprintf(
    "<must_interactions> %d response(s); null models: %s%s\n",
    length(s$families), paste(s$null_models, collapse = ", "),
    if (is.null(s$null_model_choice)) "" else " (from suggest_null_models())"
  ))
  # Models fitted outside must record their class in `settings$model`.
  if (!is.null(s$model)) {
    cat(sprintf("Model: %s, %s; means from emmeans\n", s$model, s$families[[1L]]))
  }
  cat(sprintf(
    "Inference: %s, %s; %s%% intervals; p-values adjusted with %s; seed %s\n",
    if (identical(s$inference, "simulation")) sprintf("Krinsky-Robb simulation (%d draws)", s$n_sim) else "delta method",
    correction_label(s, detail = TRUE), format(100 * s$level), s$p_adjust, format(s$seed)
  ))
  scales <- c(AD = "AD = differences (response scale)", MU = "MU = log ratios",
              DO = "DO = dominance (response scale)",
              BL = "BL = relative to the reference cell (Bliss independence)",
              SAM = "SAM = relative to the reference cell (stress addition)",
              CA = "CA = relative to the reference cell (concentration addition)")
  cat("Scales: ", paste(scales[s$null_models], collapse = "; "), "\n", sep = "")
  if (!is.null(s$margin)) {
    cat(sprintf(paste("Equivalence margin: +/-%s of the reference cell for the relative",
                      "deviations (`relative`)\n"), format(s$margin)))
  }
  est <- x$estimates
  int <- est[est$component == "interaction", , drop = FALSE]
  if (nrow(int) == 0L) {
    cat("No interaction could be estimated.\n")
  } else {
    show <- data.frame(
      response = int$response,
      time = int$time,
      combination = set_combination(int$levels, int$stressors),
      background = compact_levels(int$background, x$design$spec$control),
      model = int$null_model,
      estimand = int$estimand,
      estimate = signif(int$estimate, 3),
      lower = signif(int$lower, 3),
      upper = signif(int$upper, 3),
      df = signif(int$df, 3),
      p_adj = signif(int$p_adjusted, 2),
      relative = signif(int$relative, 2),
      p_eq_adj = signif(int$p_equivalence_adjusted, 2),
      rescaled = signif(int$rescaled, 2),
      status = int$status,
      stringsAsFactors = FALSE
    )
    if (is.null(s$margin)) {
      show$relative <- NULL
      show$p_eq_adj <- NULL
    }
    if (!"BL" %in% s$null_models) {
      show$rescaled <- NULL
    }
    if (all(is.na(show$time))) {
      show$time <- NULL
    }
    if (!any(nzchar(show$background))) {
      show$background <- NULL
    }
    cat(sprintf("Interactions (%d rows%s):\n", nrow(show),
                if (nrow(show) > n) sprintf(", first %d shown", n) else ""))
    print(utils::head(show, n), row.names = FALSE)
  }
  warn <- x$design$checks[x$design$checks$severity != "info", , drop = FALSE]
  shared <- x$design$checks[x$design$checks$check == "shared_cells", , drop = FALSE]
  notes <- rbind(warn, shared)
  if (nrow(notes) > 0L) {
    cat("Design notes:\n")
    print_checks(notes)
  }
  if (nrow(x$diagnostics$failed_responses) > 0L) {
    cat("Responses without a fitted model: ",
        paste(x$diagnostics$failed_responses$response, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
