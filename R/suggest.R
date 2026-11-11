# Null models from what is known about the stressors ---------------------------
#
# The guidelines of Schäfer & Piggott (2018, Fig. 3) as seven steps checked in
# order; the first one whose conditions hold gives the null model. Criteria
# left unknown (NA) take each of their values, and every null model reached
# with some of them is kept: the range of predictions that the authors
# suggest when information is lacking (their section 4).

#' Suggest null models from what is known about the stressors
#'
#' @description
#' Applies the guidelines of Schäfer and Piggott (2018, Fig. 3) for the choice
#' of a null model, from what is known about the stressors and their effects:
#' whether the stressors act similarly, whether one of them acts specifically
#' on the physiology of the organisms, whether the effect is mortality,
#' whether the stressor-effect relationship or the distribution of
#' sensitivities is known, how the sensitivities of individuals to the
#' stressors are correlated, and how large the effects of the stressors alone
#' are relative to the effect limit (the largest possible effect, such as
#' 100% mortality or the loss of the whole response). A criterion left
#' unknown (`NA`) is tried with each of its values, and every null model that
#' one of them leads to is kept: the range of predictions that the authors
#' recommend when the information is lacking, which [compare_null_models()]
#' then compares.
#'
#' The result can be given as `null_models` to [estimate_interactions()],
#' which estimates the null models it suggests and keeps the criteria and the
#' reasons in `settings$null_model_choice`.
#'
#' @param x A model fitted with [fit_responses()], from which the effects of
#'   the stressors alone are estimated; a model fitted outside must, with its
#'   `design`, as in [estimate_interactions()]; the output of
#'   [describe_design()], to check which null models the design allows; or
#'   `NULL` (the default), to apply only the criteria given.
#' @param similar_action Do the stressors act similarly, with the same mode of
#'   action? Schäfer and Piggott (2018) distinguish stressors that act from
#'   outside on the whole organism (predation, floods), relatively
#'   nonspecifically on its physiology (competition, salinity, heat) and
#'   specifically on physiological processes (many toxicants). `TRUE`,
#'   `FALSE` or `NA` (unknown, the default).
#' @param specific_action Does one of the stressors act specifically on
#'   physiological processes, as many toxicants do? `TRUE`, `FALSE` or `NA`.
#' @param mortality Is the effect mortality? `TRUE`, `FALSE` or `NA`.
#' @param stressor_effect Is the stressor-effect relationship of the stressors
#'   known, from earlier studies or from enough levels of each stressor in the
#'   experiment (see Details)? `TRUE`, `FALSE` or `NA`.
#' @param sensitivity_distribution Is the distribution of the sensitivities
#'   (stress capacities) of the individuals known? `TRUE`, `FALSE` or `NA`.
#' @param correlation How the sensitivities of individuals to the stressors
#'   are correlated: `"positive"` (highly positively: the individuals
#'   sensitive to one stressor are sensitive to the other), `"negative"`
#'   (highly negatively: those sensitive to one are tolerant to the other),
#'   `"weak"` (neither) or `NA` (unknown, the default).
#' @param effects Effects of the stressors alone as proportions of the effect
#'   limit: a named numeric vector with one value per stressor, such as
#'   `c(salinity = 0.3, warming = 0.5)`. Values cannot exceed 1 (the effect
#'   limit); a negative value is a stressor that increases the response, and
#'   `NA` an unknown effect. `NULL` (the default) estimates them from `x` when
#'   it is a fitted model and leaves them unknown otherwise.
#' @param near_limit Proportion of the effect limit from which the effect of a
#'   stressor counts as close to it. Default 0.9 (see Details).
#' @param level Confidence level of the intervals of the effects estimated
#'   from `x` (default 0.95). A criterion that depends on the effects is
#'   unknown when an interval spans its threshold; with `NULL`, the point
#'   estimates decide.
#' @param response For a model fitted with [fit_responses()], the names of
#'   the responses whose effects are estimated; `NULL` (the default) for all.
#'   The other criteria apply to all of them: when they differ (mortality for
#'   one response, growth for another), call the function once per response.
#' @param sensitivity As in [estimate_interactions()]: for a model fitted
#'   with [fit_responses()], the names of the responses for which an increase
#'   means a stronger effect (mortality, damage), whose proportions are then
#'   complemented (mortality becomes survival); for a model fitted outside
#'   must, `TRUE` or `FALSE`. Only used when the effects are estimated.
#' @param design For a model fitted outside must, the output of
#'   [describe_design()] for its data.
#' @param data,... For a model fitted outside must, as in
#'   [estimate_interactions()]: the data it was fitted to, when emmeans cannot
#'   recover them, and further arguments to [emmeans::ref_grid()].
#'
#' @details
#' **The guidelines.** The steps of Fig. 3 of Schäfer and Piggott (2018), in
#' the order in which they are checked; the first whose conditions hold gives
#' the null model:
#'
#' 1. The stressors act similarly and their stressor-effect relationship is
#'    known: concentration addition (`"CA"`). In must, also with effects of
#'    the stressors alone in the same direction (see below).
#' 2. One stressor acts specifically on physiology, the effect is mortality,
#'    and the stressor-effect relationship or the sensitivity distribution is
#'    known: the stress addition model (`"SAM"`). In must, also with no
#'    stressor alone that increases the response.
#' 3. The sensitivities are highly positively correlated, or the effect of
#'    one stressor is close to the effect limit: dominance (`"DO"`).
#' 4. The sensitivities are highly negatively correlated and the sum of the
#'    effects is below the effect limit: simple addition (`"AD"`).
#' 5. The effect is mortality: the multiplicative model (`"MU"`).
#' 6. The sum of the effects exceeds the effect limit: the multiplicative
#'    model (`"MU"`).
#' 7. Otherwise: simple addition (`"AD"`).
#'
#' The multiplicative model is the effect addition of ecotoxicology (Bliss
#' independence): `"BL"` in [estimate_interactions()] gives the same
#' prediction, with the deviation relative to the reference cell instead of
#' on the log scale. When the mode of action or the correlation of the
#' sensitivities is unknown, the authors suggest computing every null model
#' that matches the remaining criteria; each unknown criterion here takes all
#' its values, and `decisive` lists those whose value would change the
#' models kept. They also note that the choice matters less when the effect
#' of one stressor is close to the effect limit or about ten times larger
#' than that of the other, since the null models then predict similar joint
#' effects (except perhaps the stress addition model and concentration
#' addition); `similar_predictions` in `sets` flags those cases.
#'
#' The guidelines were written for two stressors and for effects on
#' individuals (or populations, with caveats). For three or more stressors,
#' the same criteria are applied to all the stressors of a contrast set (the
#' effect of any of them close to the limit, the sum of all their effects),
#' and a larger effect counts as about ten times the others when it is at
#' least ten times their sum; dominance is defined in must for two stressors
#' only. Schäfer and Piggott (2018) give no value for an effect close to the
#' effect limit: `near_limit = 0.9` is a choice of this package. They suggest
#' at least five levels of a stressor, on an interval or ratio scale, to
#' establish its stressor-effect relationship in a factorial experiment;
#' with a design, `models` says how many levels of each stressor were tested
#' alone. Concentration addition in must fits its curves to two or more
#' doses, and the stress addition model assumes a beta distribution of the
#' stress capacities (`sam_shape` in [estimate_interactions()]).
#'
#' **Effects estimated from a fitted model.** For each contrast set and
#' sampling time, the effect of each stressor alone is 1 - y_i / y_0, from
#' the means of the model at the reference cell of the set (y_0) and at the
#' cell with that stressor alone (y_i): the reduction of the response as a
#' share of the reference, whose limit (the loss of the whole response) is 1.
#' The cells with two or more stressors are not used, so that the choice does
#' not depend on the interactions it serves to test. The interval of each
#' effect comes from the delta method on the log scale of y_i / y_0, and that
#' of their sum, the number of stressors minus the sum of the ratios
#' y_i / y_0, on the log scale of that sum, with the covariance and the
#' degrees of freedom that [estimate_interactions()] uses by default (REML
#' and the small-sample correction). One effect is close to the limit if the
#' lower bound of its interval reaches `near_limit`, and none is if every
#' upper bound is below it; otherwise the criterion is unknown. The sum of the
#' effects is below the effect limit if the upper bound of its interval is,
#' and above it if the lower bound is. A cell whose mean is fixed
#' at the boundary (only zeros, or only successes) has no variance: the effect
#' of a stressor alone in such a cell gets no interval, and so can take any
#' value up to the effect limit, while the uncertainty of a reference cell at
#' the boundary is left out. Effects are not defined when the mean of the
#' reference cell is zero, when the response takes negative values, or when
#' an increase means a stronger effect in a response that is not a
#' proportion (`sensitivity`): such a response has no effect limit. The null
#' models kept are those reached in any set and time; `sets` gives them for
#' each one.
#'
#' **Stressors that increase the response.** A stressor, in the sense of Orr et
#' al. (2020) and Madge Pimentel et al. (2025), is any driver of change,
#' whatever its direction: in the 170 studies reanalysed by Piggott et al.
#' (2015), 29% of the pairs of stressors had opposing effects and 22% two
#' positive ones. Each effect of a stressor alone has a `direction`: `"decrease"` or
#' `"increase"` of the response when its interval excludes zero (by its sign
#' for effects given or point estimates), `"none"` when it includes zero, and
#' `"unknown"` without an interval. Each set gets the interaction `type` of
#' Piggott et al. (2015, Table 1) from them: `"double negative"` (both
#' decrease the response), `"double positive"`, `"opposing"`, `"negative
#' neutral"`, `"positive neutral"` or `"double neutral"` (`"all"` instead of
#' `"double"` for three or more stressors). The guidelines take the effects
#' with their sign. The sum of the effects exceeds the effect limit exactly
#' when simple addition predicts a response below zero, y_0 (1 - sum of the
#' effects) < 0, so an increase offsets the decreases; only a decrease can be
#' close to the effect limit. Concentration addition treats the stressors as
#' dilutions of one another, which is not possible with opposite effects, and
#' the stress addition model assumes stressors that reduce the response
#' (Liess et al. 2016): steps 1 and 2 also need, in must, effects in the
#' same direction and no stressor alone that increases the response
#' (`opposing` and `increase` in `sets`), as [estimate_interactions()]
#' already requires for those models. For a proportion fitted with `trials`
#' that the stressors alone increase and never decrease, the guidelines are
#' applied to its complement (1 - p, mortality for survival), where the
#' increases are reductions towards the limit 0: the mirror of the effect
#' limit for the upper limit 1, which simple addition can also cross, and
#' the scale on which the multiplicative model and dominance describe
#' independent and correlated effects (Madge Pimentel et al. 2025 discuss
#' the asymmetry between survival and mortality). The interactions must
#' then be estimated on that scale: `settings$suggested_sensitivity` gives
#' the `sensitivity` for [estimate_interactions()], which checks it. For
#' responses of whole communities, such as a total density, some taxa often
#' increase when others are harmed, through interactions among them; the
#' guidelines were written for individuals and populations, and Schäfer and
#' Piggott (2018) point instead to the compositional null model of Thompson
#' et al. (2018), which predicts the community response from the null models
#' of each taxon (not implemented in must).
#'
#' @return An object of class `must_null_models`, a list with:
#' * `null_models`: the codes of the null models kept, for
#'   `estimate_interactions(null_models = )`.
#' * `models`: one row per null model of the guidelines (`"AD"`, `"MU"`,
#'   `"DO"`, `"SAM"`, `"CA"`) with its `status` (`"suggested"` when it is
#'   the only one reached, `"possible"`, or `"excluded"`), `reached` (in how
#'   many of the `evaluations`, the contrast sets and times of a fitted
#'   model, it is reached), `reason` (the steps that lead to it, or why it is
#'   not reached), `design` (whether the design allows it: `"yes"`, `"no"`,
#'   `"some sets"`, or `NA` without a design), `design_note` and `note`.
#' * `criteria`: the value of each criterion (`"yes"`, `"no"`, a
#'   correlation, `"below"`, `"above"` or `"equal"` for the sum of the
#'   effects, `"unknown"`, `"uncertain"` when an interval spans the
#'   threshold, or `"varies"` among sets and times) and its `source`
#'   (`"given"`, `"estimated"` from the fitted model, or `"unknown"`), with
#'   counts in `detail` when it varies.
#' * `decisive`: for each unknown criterion whose value would change the
#'   models kept, the models kept with each of its values.
#' * `effects`: one row per stressor (and response, contrast set and time,
#'   for a fitted model) with the effect (`estimate`), its interval
#'   (`lower`, `upper`, `level`, `df`), its `direction`, `status` and
#'   `reason`; for a fitted model, `scale` says whether the effect is on the
#'   response as analysed (`"response"`) or on its complement
#'   (`"complement"`, for the proportions described above, which have rows on
#'   both scales). `NULL` when no effects were given or estimated.
#' * `sets`: for a fitted model, one row per response, contrast set and time
#'   with the scale to which the guidelines were applied (`scale`), the
#'   `type` of the effects (on the response as analysed), the sum of the
#'   effects and its interval, whether one effect is close to the limit
#'   (`near_limit`), how the sum compares with the effect limit
#'   (`sum_vs_limit`: `"below"`, `"above"`, `"equal"`, or `NA` if unknown),
#'   whether a stressor alone increases the response (`increase`) and two
#'   have opposite effects (`opposing`), `similar_predictions` and the null
#'   models reached (`null_models`); `NULL` otherwise.
#' * `settings`: the thresholds, the source of the effects, the responses
#'   whose complement the guidelines were applied to (`complement`) and the
#'   `sensitivity` with which to estimate the interactions
#'   (`suggested_sensitivity`).
#'
#' @references
#' Schäfer, R. B., & Piggott, J. J. (2018). Advancing understanding and
#' prediction in multiple stressor research through a mechanistic basis for
#' null models. *Global Change Biology*, 24(5), 1817–1826.
#' \doi{10.1111/gcb.14073}
#'
#' Liess, M., Foit, K., Knillmann, S., Schäfer, R. B., & Liess, H.-D. (2016).
#' Predicting the synergy of multiple stress effects. *Scientific Reports*,
#' 6, 32965. \doi{10.1038/srep32965}
#'
#' Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
#' J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025).
#' Hypothesis-driven research on multiple stressors: An analytical framework
#' for stressor interactions. *Ecology and Evolution*, 15, e71959.
#' \doi{10.1002/ece3.71959}
#'
#' Orr, J. A., Vinebrooke, R. D., Jackson, M. C., Kroeker, K. J., Kordas, R.
#' L., Mantyka-Pringle, C., Van den Brink, P. J., De Laender, F., Stoks, R.,
#' Holmstrup, M., Matthaei, C. D., Monk, W. A., Penk, M. R., Leuzinger, S.,
#' Schäfer, R. B., & Piggott, J. J. (2020). Towards a unified study of
#' multiple stressors: Divisions and common goals across research
#' disciplines. *Proceedings of the Royal Society B*, 287, 20200421.
#' \doi{10.1098/rspb.2020.0421}
#'
#' Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015).
#' Reconceptualizing synergism and antagonism among multiple stressors.
#' *Ecology and Evolution*, 5(7), 1538–1547. \doi{10.1002/ece3.1465}
#'
#' Thompson, P. L., MacLennan, M. M., & Vinebrooke, R. D. (2018). An
#' improved null model for assessing the net effects of multiple stressors
#' on communities. *Global Change Biology*, 24(1), 517–525.
#' \doi{10.1111/gcb.13852}
#'
#' @seealso [estimate_interactions()], which takes the result as
#'   `null_models`, and [compare_null_models()], which compares the null
#'   models kept.
#'
#' @examples
#' # Two stressors with different modes of action, a nonlethal effect and
#' # effects of 30% and 50% of the effect limit; the correlation of the
#' # sensitivities is unknown
#' suggest_null_models(similar_action = FALSE, mortality = FALSE,
#'   stressor_effect = FALSE, effects = c(salinity = 0.3, warming = 0.5))
#' # The same, with an effect close to the effect limit
#' suggest_null_models(similar_action = FALSE, mortality = FALSE,
#'   stressor_effect = FALSE, effects = c(salinity = 0.3, warming = 0.95))
#'
#' # Effects estimated from a fitted model
#' jars <- expand.grid(
#'   salinity = c("ambient", "high"),
#'   warming = c("no", "yes"),
#'   replicate = 1:5,
#'   stringsAsFactors = FALSE
#' )
#' jars$jar <- seq_len(nrow(jars))
#' mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
#'   ifelse(jars$warming == "yes", 0.4, 1)
#' jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
#' design <- describe_design(jars, c("salinity", "warming"),
#'   control = c(salinity = "ambient", warming = "no"), unit = "jar")
#' fits <- fit_responses(jars, "algae", design, residual_checks = FALSE,
#'   seed = 1, quiet = TRUE)
#' choice <- suggest_null_models(fits, similar_action = FALSE, mortality = FALSE,
#'   stressor_effect = FALSE)
#' choice
#' choice$effects[, c("stressor", "estimate", "lower", "upper")]
#' # The null models kept, with their justification in the settings
#' est <- estimate_interactions(fits, null_models = choice, n_sim = 2000, seed = 1)
#' est$settings$null_model_choice$criteria
#' @export
suggest_null_models <- function(x = NULL, similar_action = NA, specific_action = NA,
                                mortality = NA, stressor_effect = NA,
                                sensitivity_distribution = NA, correlation = NA,
                                effects = NULL, near_limit = 0.9, level = 0.95,
                                response = NULL, sensitivity = NULL, design = NULL,
                                data = NULL, ...) {
  call <- sys.call()
  if (inherits(x, "must_interactions")) {
    abort_input(paste(
      "`x` is the output of estimate_interactions(). Give the fitted model",
      "(the output of fit_responses()) or the design instead."
    ), call = call)
  }
  kind <- if (is.null(x)) {
    "none"
  } else if (inherits(x, "must_fit")) {
    "fit"
  } else if (inherits(x, "must_design")) {
    "design"
  } else {
    "external"
  }
  if (identical(kind, "external")) {
    if (is.null(design) || !inherits(design, "must_design")) {
      abort_input(sprintf(
        paste(
          "`x` is %s. Give a model fitted with fit_responses(), the output of",
          "describe_design(), NULL, or a model fitted outside must with the",
          "output of describe_design() for its data as `design`."
        ),
        describe_value(x)
      ), call = call)
    }
    if (!is.null(data)) {
      check_data_frame(data, call = call)
    }
  } else {
    check_dots_empty(..., call = call)
    if (!is.null(design)) {
      abort_input(paste(
        "`design` is only used with a model fitted outside must; the design of",
        "a model fitted with fit_responses() is part of it."
      ), call = call)
    }
    if (!is.null(data)) {
      abort_input("`data` is only used with a model fitted outside must.", call = call)
    }
  }
  the_design <- switch(kind, none = NULL, fit = x$design, design = x, external = design)
  criteria <- list(
    similar_action = check_criterion(similar_action, "similar_action", call),
    specific_action = check_criterion(specific_action, "specific_action", call),
    mortality = check_criterion(mortality, "mortality", call),
    stressor_effect = check_criterion(stressor_effect, "stressor_effect", call),
    sensitivity_distribution = check_criterion(sensitivity_distribution,
                                               "sensitivity_distribution", call),
    correlation = check_correlation(correlation, call)
  )
  check_number(near_limit, lower = 0, upper = 1, include_lower = FALSE, call = call)
  if (!is.null(level)) {
    check_number(level, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE,
                 call = call)
  }
  effects <- check_effects(effects, the_design, call)
  if (!is.null(effects) && kind %in% c("fit", "external")) {
    abort_input(paste(
      "Give `effects` or a fitted model in `x`, not both. To use your own",
      "effects with the design of a fitted model, give the design (the output",
      "of describe_design()) as `x`."
    ), call = call)
  }
  responses <- NULL
  if (identical(kind, "fit")) {
    responses <- check_response_names(response, names(x$fits), call)
    sensitivity <- check_sensitivity(sensitivity, names(x$fits), call = call)
  } else {
    if (!is.null(response)) {
      abort_input("`response` is only used with a model fitted with fit_responses().",
                  call = call)
    }
    if (identical(kind, "external")) {
      if (is.null(sensitivity)) {
        sensitivity <- FALSE
      }
      check_flag(sensitivity, call = call)
    } else if (!is.null(sensitivity)) {
      abort_input(paste(
        "`sensitivity` is only used when the effects are estimated from a fitted",
        "model in `x`."
      ), call = call)
    }
  }
  estimated <- NULL
  if (identical(kind, "fit")) {
    estimated <- fit_effects(x, responses, sensitivity, level, near_limit)
    fitted <- vapply(x$fits[responses], function(f) identical(f$status, "ok"), logical(1))
    responses <- responses[fitted]
  } else if (identical(kind, "external")) {
    check_installed("emmeans", "suggest null models from a model fitted outside must",
                    call = call)
    ext <- external_predictions(x, design, data, 100L, sensitivity, call, ...)
    responses <- ext$f$response
    dots <- list(...)
    flipped_pred <- function() {
      do.call(external_predictions,
              c(list(x, design, data, 100L, !sensitivity, call), dots))$pred
    }
    estimated <- response_effects(
      ext$f, function(flip) if (flip) flipped_pred() else ext$pred,
      design, ext$baseline, set_df = function(s) NA_real_, level = level,
      near_limit = near_limit, sensitive = sensitivity
    )
    estimated$complement <- if (isTRUE(estimated$complement)) responses else character()
  }
  sets <- NULL
  if (!is.null(estimated)) {
    sets <- estimated$sets
    ok <- sets$status == "ok"
    evals <- if (any(ok)) {
      sets[ok, c("near_limit", "sum_vs_limit", "increase", "opposing"), drop = FALSE]
    } else {
      data.frame(near_limit = NA, sum_vs_limit = NA_character_, increase = NA, opposing = NA,
                 stringsAsFactors = FALSE)
    }
    names(evals)[names(evals) == "sum_vs_limit"] <- "sum"
    rownames(evals) <- NULL
    origin <- if (any(ok)) "estimated" else "unknown"
  } else {
    ev <- effect_criteria(effects, lower = NULL, upper = NULL, sum = NA_real_,
                          sum_lower = NA_real_, sum_upper = NA_real_,
                          near_limit = near_limit, use_points = TRUE)
    evals <- data.frame(near_limit = ev$near_limit, sum = ev$sum, increase = ev$increase,
                        opposing = ev$opposing, stringsAsFactors = FALSE)
    origin <- if (is.null(effects)) "unknown" else "given"
  }
  reach <- lapply(seq_len(nrow(evals)), function(i) sp_scenarios(sp_evaluation(criteria, evals, i)))
  if (!is.null(sets)) {
    sets$null_models <- NA_character_
    sets$null_models[sets$status == "ok"] <- vapply(reach, function(g) {
      paste(sp_order(g$null_model), collapse = ", ")
    }, character(1))[seq_len(sum(sets$status == "ok"))]
  }
  models <- sp_models_table(reach, criteria, evals, the_design)
  complement <- if (is.null(estimated)) character() else estimated$complement
  structure(
    list(
      null_models = models$null_model[models$status != "excluded"],
      models = models,
      criteria = sp_criteria_table(criteria, evals, origin, near_limit, complement),
      decisive = sp_decisive(criteria, evals),
      effects = if (!is.null(estimated)) estimated$effects else user_effects_table(effects),
      sets = sets,
      settings = list(
        guide = "Sch\u00e4fer and Piggott (2018), Fig. 3",
        near_limit = near_limit,
        level = level,
        effects = origin,
        responses = responses,
        sensitivity = if (kind %in% c("fit", "external")) sensitivity else NULL,
        complement = complement,
        suggested_sensitivity = suggested_sensitivity(kind, sensitivity, complement),
        design = !is.null(the_design)
      )
    ),
    class = "must_null_models"
  )
}

# The `sensitivity` with which the interactions should be estimated: the one
# used for the effects, with the responses whose complement the guidelines
# were applied to switched (a proportion analysed as p is then analysed as
# 1 - p, and the other way round).
suggested_sensitivity <- function(kind, sensitivity, complement) {
  if (identical(kind, "fit")) {
    sort(c(setdiff(sensitivity, complement), setdiff(complement, sensitivity)))
  } else if (identical(kind, "external")) {
    xor(isTRUE(sensitivity), length(complement) > 0L)
  } else {
    NULL
  }
}

# Checks of the criteria ------------------------------------------------------

check_criterion <- function(x, arg, call) {
  if (length(x) != 1L || !is.atomic(x) || !(is.na(x) || isTRUE(x) || isFALSE(x))) {
    abort_input(sprintf("`%s` must be TRUE, FALSE or NA (unknown); it is %s.", arg,
                        describe_value(x)), call = call)
  }
  if (is.na(x)) NA else x
}

check_correlation <- function(x, call) {
  choices <- c("positive", "negative", "weak")
  if (length(x) == 1L && is.atomic(x) && is.na(x)) {
    return(NA_character_)
  }
  if (!is.character(x) || length(x) != 1L || !x %in% choices) {
    abort_input(sprintf(
      "`correlation` must be %s or NA (unknown); it is %s.",
      format_items(paste0("\"", choices, "\""), max = Inf), describe_value(x)
    ), call = call)
  }
  x
}

check_effects <- function(effects, design, call) {
  if (is.null(effects)) {
    return(NULL)
  }
  nm <- names(effects)
  if (!(is.numeric(effects) || (is.logical(effects) && all(is.na(effects)))) ||
      length(effects) < 2L || is.null(nm) || anyNA(nm) || any(!nzchar(nm)) ||
      anyDuplicated(nm) > 0L) {
    abort_input(paste(
      "`effects` must be a named numeric vector with one effect per stressor (at",
      "least two), such as c(salinity = 0.3, warming = 0.5)."
    ), call = call)
  }
  values <- as.numeric(effects)
  bad <- !is.na(values) & (!is.finite(values) | values > 1)
  if (any(bad)) {
    abort_input(sprintf(
      paste(
        "`effects` are proportions of the effect limit and cannot exceed 1 (the",
        "effect limit): %s."
      ),
      format_items(paste0(nm[bad], " = ", format(values[bad])), max = 3L)
    ), call = call)
  }
  if (!is.null(design)) {
    unknown <- setdiff(nm, design$spec$stressors)
    if (length(unknown) > 0L) {
      abort_input(sprintf(
        "`effects` names %s, which %s not a stressor of the design (%s).",
        format_items(quote_items(unknown), max = 3L),
        if (length(unknown) == 1L) "is" else "are",
        format_items(quote_items(design$spec$stressors), max = Inf)
      ), call = call)
    }
  }
  stats::setNames(values, nm)
}

check_response_names <- function(response, responses, call) {
  if (is.null(response)) {
    return(responses)
  }
  if (!is.character(response) || length(response) == 0L || anyNA(response) ||
      !all(response %in% responses)) {
    abort_input(sprintf(
      "`response` must be NULL or names of fitted responses (%s).",
      format_items(quote_items(responses), max = Inf)
    ), call = call)
  }
  unique(response)
}

# The guidelines (Schäfer & Piggott 2018, Fig. 3) -------------------------------

sp_codes <- function() {
  c("AD", "MU", "DO", "SAM", "CA")
}

# The order of the null models in must.
sp_order <- function(x) {
  codes <- c("AD", "MU", "DO", "BL", "SAM", "CA")
  codes[codes %in% x]
}

sp_names <- function() {
  c(AD = "simple addition", MU = "multiplicative", DO = "dominance",
    SAM = "stress addition", CA = "concentration addition")
}

# Values that an unknown criterion can take. The sum of the effects can also
# equal the effect limit (only with effects given): no step chooses it, and
# it leads to models also reached with the other values.
sp_domains <- function() {
  list(similar_action = c(TRUE, FALSE), specific_action = c(TRUE, FALSE),
       mortality = c(TRUE, FALSE), stressor_effect = c(TRUE, FALSE),
       sensitivity_distribution = c(TRUE, FALSE),
       correlation = c("positive", "negative", "weak"),
       near_limit = c(TRUE, FALSE), sum = c("below", "above"),
       increase = c(TRUE, FALSE), opposing = c(TRUE, FALSE))
}

# The criteria of evaluation `i`: those given and those that depend on the
# effects (one row of `evals`).
sp_evaluation <- function(criteria, evals, i) {
  c(criteria, list(near_limit = evals$near_limit[[i]], sum = evals$sum[[i]],
                   increase = evals$increase[[i]], opposing = evals$opposing[[i]]))
}

# Null model of each step.
sp_step_models <- function() {
  c("CA", "SAM", "DO", "AD", "MU", "MU", "AD")
}

sp_step_texts <- function() {
  c(
    paste(
      "the stressors act similarly, with effects in the same direction, and their",
      "stressor-effect relationship is known"
    ),
    paste(
      "one stressor acts specifically on physiology, the effect is mortality, no",
      "stressor alone increases the response, and the stressor-effect relationship",
      "or the sensitivity distribution is known"
    ),
    paste(
      "the sensitivities are highly positively correlated, or the effect of one",
      "stressor is close to the effect limit"
    ),
    paste(
      "the sensitivities are highly negatively correlated and the sum of the",
      "effects is below the effect limit"
    ),
    "the effect is mortality",
    "the sum of the effects exceeds the effect limit",
    "none of the earlier steps applies"
  )
}

# The step at which the guidelines stop for a complete set of criteria.
sp_step <- function(cr) {
  if (cr$similar_action && cr$stressor_effect && !cr$opposing) {
    return(1L)
  }
  if (cr$specific_action && cr$mortality && !cr$increase &&
      (cr$stressor_effect || cr$sensitivity_distribution)) {
    return(2L)
  }
  if (identical(cr$correlation, "positive") || cr$near_limit) {
    return(3L)
  }
  if (identical(cr$correlation, "negative") && identical(cr$sum, "below")) {
    return(4L)
  }
  if (cr$mortality) {
    return(5L)
  }
  if (identical(cr$sum, "above")) {
    return(6L)
  }
  7L
}

# Every combination of values of the unknown criteria, with its step and null
# model.
sp_scenarios <- function(cr) {
  domains <- sp_domains()
  values <- lapply(names(domains), function(n) {
    if (is.na(cr[[n]])) domains[[n]] else cr[[n]]
  })
  names(values) <- names(domains)
  grid <- expand.grid(values, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  # Opposite effects need one that increases the response.
  grid <- grid[!(grid$opposing & !grid$increase), , drop = FALSE]
  grid$step <- vapply(seq_len(nrow(grid)), function(i) {
    sp_step(as.list(grid[i, , drop = FALSE]))
  }, integer(1))
  grid$null_model <- sp_step_models()[grid$step]
  grid
}

# Known criteria that rule out the conditions of a step (none if they can
# hold).
sp_unmet <- function(step, cr) {
  out <- switch(
    step,
    c(if (isFALSE(cr$similar_action)) "the stressors do not act similarly",
      if (isTRUE(cr$opposing)) "the stressors alone have opposite effects",
      if (isFALSE(cr$stressor_effect)) "their stressor-effect relationship is not known"),
    c(if (isFALSE(cr$specific_action)) "no stressor acts specifically on physiology",
      if (isFALSE(cr$mortality)) "the effect is not mortality",
      if (isTRUE(cr$increase)) "a stressor alone increases the response",
      if (isFALSE(cr$stressor_effect) && isFALSE(cr$sensitivity_distribution)) {
        "neither the stressor-effect relationship nor the sensitivity distribution is known"
      }),
    if (!is.na(cr$correlation) && !identical(cr$correlation, "positive") &&
        isFALSE(cr$near_limit)) {
      c("the sensitivities are not highly positively correlated",
        "no effect is close to the effect limit")
    },
    c(if (!is.na(cr$correlation) && !identical(cr$correlation, "negative")) {
      "the sensitivities are not highly negatively correlated"
    },
    if (!is.na(cr$sum) && !identical(cr$sum, "below")) {
      "the sum of the effects is not below the effect limit"
    }),
    if (isFALSE(cr$mortality)) "the effect is not mortality",
    if (!is.na(cr$sum) && !identical(cr$sum, "above")) {
      "the sum of the effects does not exceed the effect limit"
    },
    NULL
  )
  if (is.null(out)) character() else out
}

sp_requirements <- function() {
  c(
    CA = paste(
      "Concentration addition needs stressors that act similarly, with effects in",
      "the same direction, and a known stressor-effect relationship (step 1)"
    ),
    SAM = paste(
      "The stress addition model needs one stressor that acts specifically on",
      "physiology, mortality as the effect, stressors that do not increase the",
      "response and a known stressor-effect relationship or sensitivity",
      "distribution (step 2)"
    ),
    DO = paste(
      "Dominance needs highly positively correlated sensitivities or the effect",
      "of one stressor close to the effect limit (step 3)"
    ),
    AD = paste(
      "Simple addition needs highly negatively correlated sensitivities with the",
      "sum of the effects below the effect limit (step 4), or that no earlier",
      "step applies (step 7)"
    ),
    MU = paste(
      "The multiplicative model needs mortality as the effect (step 5) or the",
      "sum of the effects above the effect limit (step 6)"
    )
  )
}

sp_step_label <- function(steps) {
  sprintf("step %d, %s", steps, sp_names()[sp_step_models()[steps]])
}

# "a", "a and b", "a, b and c"
and_list <- function(x) {
  n <- length(x)
  if (n <= 1L) {
    return(paste(x, collapse = ""))
  }
  paste(paste(x[-n], collapse = ", "), "and", x[[n]])
}

# Why a null model is reached (the steps that lead to it) or not.
sp_reason <- function(model, steps_reached, cr) {
  steps <- which(sp_step_models() == model)
  hit <- intersect(steps, steps_reached)
  if (length(hit) > 0L) {
    texts <- sprintf("step %d: %s", hit, sp_step_texts()[hit])
    out <- paste(texts, collapse = "; ")
    return(paste0(toupper(substr(out, 1L, 1L)), substr(out, 2L, nchar(out)), "."))
  }
  unmet <- lapply(steps, sp_unmet, cr = cr)
  if (all(lengths(unmet) > 0L)) {
    return(sprintf("%s, but %s.", sp_requirements()[[model]], and_list(unique(unlist(unmet)))))
  }
  open <- max(steps[lengths(unmet) == 0L])
  earlier <- steps_reached[steps_reached < open]
  sprintf("Not reached: an earlier step always applies (%s).",
          paste(sp_step_label(earlier), collapse = "; "))
}

# Criteria that depend on the effects ------------------------------------------

# Whether one effect is close to the effect limit and how their sum compares
# with it, from intervals (or from the point estimates). Unknown effects can
# take any value up to the effect limit.
effect_criteria <- function(estimate, lower, upper, sum, sum_lower, sum_upper, near_limit,
                            use_points) {
  if (length(estimate) == 0L) {
    return(list(near_limit = NA, sum = NA_character_, similar = NA, increase = NA,
                opposing = NA, type = NA_character_))
  }
  directions <- effect_directions(estimate, lower, upper, use_points)
  lo <- if (use_points) estimate else lower
  hi <- if (use_points) estimate else upper
  lo[is.na(lo)] <- -Inf
  hi[is.na(hi)] <- 1
  near <- if (any(lo >= near_limit)) TRUE else if (all(hi < near_limit)) FALSE else NA
  if (use_points && !anyNA(estimate)) {
    total <- sum(estimate)
    sum_class <- if (abs(total - 1) < 1e-8) "equal" else if (total < 1) "below" else "above"
  } else {
    s_lo <- if (use_points || is.na(sum_lower)) sum(lo) else sum_lower
    s_hi <- if (use_points || is.na(sum_upper)) sum(hi) else sum_upper
    sum_class <- if (s_hi < 1) "below" else if (s_lo > 1) "above" else NA_character_
  }
  c(list(near_limit = near, sum = sum_class, similar = similar_predictions(estimate, near_limit)),
    direction_criteria(directions))
}

# Direction of each effect of a stressor alone: "decrease" or "increase" of the
# response when its interval excludes zero (or by its sign, for point
# estimates), "none" when the interval includes zero, "unknown" without an
# interval or an estimate.
effect_directions <- function(estimate, lower, upper, use_points) {
  if (use_points) {
    out <- ifelse(is.na(estimate), "unknown",
                  ifelse(estimate > 0, "decrease", ifelse(estimate < 0, "increase", "none")))
  } else {
    out <- ifelse(is.na(lower) | is.na(upper), "unknown",
                  ifelse(lower > 0, "decrease", ifelse(upper < 0, "increase", "none")))
  }
  unname(as.character(out))
}

# Whether some stressor alone increases the response, whether two have
# opposite effects (TRUE, FALSE or NA when unknown directions leave it open),
# and the type of the set after Piggott et al. (2015).
direction_criteria <- function(directions) {
  n_inc <- sum(directions == "increase")
  n_dec <- sum(directions == "decrease")
  n_unk <- sum(directions == "unknown")
  increase <- if (n_inc > 0L) TRUE else if (n_unk == 0L) FALSE else NA
  possible <- (n_inc > 0L && n_unk > 0L) || (n_dec > 0L && n_unk > 0L) || n_unk >= 2L
  opposing <- if (n_inc > 0L && n_dec > 0L) TRUE else if (possible) NA else FALSE
  list(increase = increase, opposing = opposing, type = effect_type(directions))
}

# Interaction types of Piggott et al. (2015, Table 1) from the directions of
# the single effects ("negative" is a decrease of the response); for three or
# more stressors, "all" instead of "double".
effect_type <- function(directions) {
  if (length(directions) == 0L) {
    return(NA_character_)
  }
  has <- function(d) any(directions == d)
  if (has("increase") && has("decrease")) {
    return("opposing")
  }
  if (has("unknown")) {
    return(NA_character_)
  }
  all_of <- if (length(directions) == 2L) "double" else "all"
  if (all(directions == "decrease")) {
    paste(all_of, "negative")
  } else if (all(directions == "increase")) {
    paste(all_of, "positive")
  } else if (has("decrease")) {
    "negative neutral"
  } else if (has("increase")) {
    "positive neutral"
  } else {
    paste(all_of, "neutral")
  }
}

# Whether the null models predict similar joint effects (Schäfer & Piggott
# 2018, section 4): one effect close to the limit, or one at least ten times
# the sum of the others (all effects not negative).
similar_predictions <- function(estimate, near_limit) {
  if (length(estimate) < 2L || anyNA(estimate)) {
    return(NA)
  }
  if (any(estimate >= near_limit)) {
    return(TRUE)
  }
  if (any(estimate < 0)) {
    return(FALSE)
  }
  largest <- max(estimate)
  largest > 0 && largest >= 10 * (sum(estimate) - largest)
}

# Effects estimated from a fitted model ------------------------------------------

fit_effects <- function(fit, responses, sensitivity, level, near_limit) {
  parts <- lapply(responses, function(r) {
    f <- fit$fits[[r]]
    if (!identical(f$status, "ok")) {
      return(NULL)
    }
    model <- if (!is.null(f$model_reml)) f$model_reml else f$model
    sensitive <- r %in% sensitivity
    analysed <- sensitive && identical(f$kind, "binomial")
    part <- response_effects(
      f, function(flip) prepare_predictions(f, model, 100L, complement = xor(analysed, flip)),
      fit$design, baseline_times(fit$design$spec), level = level,
      near_limit = near_limit, sensitive = sensitive,
      set_df = function(s) {
        set_df(f$df_info, s, "small_sample",
               contexts = unique(set_cell_contexts(fit$design$set_cells[[s$set_id]], s)))
      }
    )
    part$response <- r
    part
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  list(
    effects = bind_rows_base(lapply(parts, `[[`, "effects"),
                             template = empty_effects(fit$design)),
    sets = bind_rows_base(lapply(parts, `[[`, "sets"), template = empty_effect_sets(fit$design)),
    complement = as.character(unlist(lapply(parts, function(p) {
      if (isTRUE(p$complement)) p$response
    })))
  )
}

# Effects of one response on the scale analysed and, for a proportion that the
# stressors alone increase and never decrease, also on its complement, to
# which the guidelines are then applied: there the increases are reductions
# towards the limit 0, and the upper limit of the proportion (1) cannot be
# crossed. `pred_for(flip)` gives the predictions on the scale analysed
# (FALSE) or on its complement (TRUE).
response_effects <- function(f, pred_for, design, baseline, set_df, level, near_limit,
                             sensitive) {
  direct <- effects_from_predictions(f, pred_for(FALSE), design, baseline, set_df, level,
                                     near_limit, sensitive, scale = "response")
  ok <- direct$effects$status == "ok"
  directions <- direct$effects$direction[ok]
  if (!identical(f$kind, "binomial") || !any(directions == "increase") ||
      any(directions == "decrease")) {
    return(c(direct, list(complement = FALSE)))
  }
  flipped <- effects_from_predictions(f, pred_for(TRUE), design, baseline, set_df, level,
                                      near_limit, sensitive, scale = "complement")
  # The type describes the effects on the response as analysed.
  flipped$sets$type <- direct$sets$type
  list(effects = rbind(direct$effects, flipped$effects), sets = flipped$sets,
       complement = TRUE)
}

# Effects of the stressors alone in every identifiable contrast set and time,
# from the predictions `pred` at the cells of f$grid (as in estimate_sets()).
effects_from_predictions <- function(f, pred, design, baseline, set_df, level, near_limit,
                                     sensitive, scale = "response") {
  key <- paste_key(f$grid$.context, f$grid$.cell, f$grid$.time)
  sets <- design$sets[design$sets$net_identifiable, , drop = FALSE]
  effect_rows <- list()
  set_rows <- list()
  for (i in seq_len(nrow(sets))) {
    s <- sets[i, , drop = FALSE]
    sc <- design$set_cells[[s$set_id]]
    sc_context <- set_cell_contexts(sc, s)
    times <- unique(as.character(f$grid$.time[f$grid$.context == s$context]))
    times <- times[!times %in% set_baseline(baseline, design, s)]
    df <- set_df(s)
    use <- c(which(sc$size == 0L), which(sc$size == 1L))
    stressors <- sc$subset[sc$size == 1L]
    stressor_levels <- vapply(seq_along(stressors), function(j) {
      row <- design$cells[match(sc$cell[use[[j + 1L]]], design$cells$cell), , drop = FALSE]
      as.character(row[[stressors[[j]]]][[1L]])
    }, character(1))
    for (tt in times) {
      idx <- match(paste_key(sc_context[use], sc$cell[use], tt), key)
      v <- single_effects(idx, sc$cell[use], pred, df, level, f, sensitive, design,
                          has_time = !is.na(tt))
      cr <- if (identical(v$status, "ok")) {
        effect_criteria(v$estimate, v$lower, v$upper, v$sum, v$sum_lower, v$sum_upper,
                        near_limit, use_points = is.null(level))
      } else {
        list(near_limit = NA, sum = NA_character_, similar = NA, increase = NA,
             opposing = NA, type = NA_character_)
      }
      directions <- if (identical(v$status, "ok")) {
        effect_directions(v$estimate, v$lower, v$upper, use_points = is.null(level))
      } else {
        rep("unknown", length(stressors))
      }
      ids <- data.frame(response = f$response, context = s$context, stringsAsFactors = FALSE)
      for (b in design$spec$by) {
        ids[[b]] <- s[[b]]
      }
      ids <- cbind(ids, data.frame(time = tt, set_id = s$set_id, order = s$order,
                                   stressors = s$stressors, stringsAsFactors = FALSE))
      effect_rows[[length(effect_rows) + 1L]] <- cbind(
        ids[rep(1L, length(stressors)), , drop = FALSE],
        data.frame(
          scale = scale, stressor = stressors, stressor_level = stressor_levels,
          estimate = v$estimate, lower = v$lower, upper = v$upper,
          level = if (is.null(level)) NA_real_ else level, df = v$df,
          direction = directions, status = v$status, reason = v$reason,
          note = ifelse(v$boundary, "boundary", ""), stringsAsFactors = FALSE
        )
      )
      set_rows[[length(set_rows) + 1L]] <- cbind(ids, data.frame(
        levels = s$levels, background = s$background, scale = scale, type = cr$type,
        sum = v$sum, sum_lower = v$sum_lower, sum_upper = v$sum_upper,
        level = if (is.null(level)) NA_real_ else level, near_limit = cr$near_limit,
        sum_vs_limit = cr$sum, increase = cr$increase, opposing = cr$opposing,
        similar_predictions = cr$similar, null_models = NA_character_, status = v$status,
        reason = v$reason, stringsAsFactors = FALSE
      ))
    }
  }
  list(effects = bind_rows_base(effect_rows, template = empty_effects(design)),
       sets = bind_rows_base(set_rows, template = empty_effect_sets(design)))
}

# Effects 1 - y_i / y_0 of the stressors alone and their sum, with intervals.
# `cols` are the grid columns of the reference cell and of each stressor
# alone, `labels` their cell labels.
single_effects <- function(cols, labels, pred, df, level, f, sensitive, design, has_time) {
  m <- length(cols) - 1L
  out <- list(estimate = rep(NA_real_, m), lower = rep(NA_real_, m),
              upper = rep(NA_real_, m), df = rep(NA_real_, m), boundary = rep(FALSE, m),
              sum = NA_real_, sum_lower = NA_real_, sum_upper = NA_real_, status = "ok",
              reason = "")
  fail <- function(status, reason) {
    out$status <- status
    out$reason <- reason
    out
  }
  control <- design$spec$control
  if (anyNA(cols)) {
    return(fail("not_estimable", sprintf(
      "No values of the response%s in %s.", if (has_time) " at this time" else "",
      format_items(quote_items(short_cells(labels[is.na(cols)], control)), max = 3L)
    )))
  }
  problem <- pred$cell_problem[cols]
  if (any(nzchar(problem))) {
    return(fail("not_estimable", cell_problem_reason(problem, labels, design)))
  }
  proportion <- identical(f$kind, "binomial")
  if (sensitive && !proportion) {
    return(fail("not_defined", paste(
      "An increase in this response means a stronger effect, which has no upper",
      "limit: its effects cannot be expressed as shares of an effect limit."
    )))
  }
  if (!proportion && !isTRUE(f$non_negative)) {
    return(fail("not_defined", paste(
      "The response takes negative values, so its effect limit is not zero:",
      "its effects cannot be expressed as shares of an effect limit."
    )))
  }
  mu <- pred$mu[cols]
  mu0 <- mu[[1L]]
  if (!is.finite(mu0) || mu0 <= 0) {
    return(fail("not_defined", paste(
      "The mean of the reference cell is zero: effects relative to it are not",
      "defined."
    )))
  }
  r <- mu[-1L] / mu0
  out$estimate <- 1 - r
  out$sum <- sum(out$estimate)
  # Means fixed at the boundary (only zeros, or only successes) have no
  # variance: the effect of a stressor alone in such a cell gets no interval,
  # and the uncertainty of a reference cell at the boundary is left out.
  alone_boundary <- pred$at_boundary[cols][-1L]
  out$boundary <- alone_boundary | pred$at_boundary[cols][[1L]]
  if (is.null(level)) {
    return(out)
  }
  V <- pred$V_eta[cols, cols, drop = FALSE]
  dmu <- pred$dmu[cols]
  critical <- function(d) {
    p <- 1 - (1 - level) / 2
    if (is.finite(d)) stats::qt(p, d) else stats::qnorm(p)
  }
  contrast_df <- function(g) {
    if (is.na(df)) pred$contrast_df(g, cols) else df
  }
  for (j in seq_len(m)) {
    if (alone_boundary[[j]] || !(r[[j]] > 0)) {
      next
    }
    g <- numeric(m + 1L)
    g[[1L]] <- -dmu[[1L]] / mu0
    g[[j + 1L]] <- dmu[[j + 1L]] / mu[[j + 1L]]
    se <- sqrt(max(drop(crossprod(g, V %*% g)), 0))
    d <- contrast_df(g)
    if (is.na(d) || d < 1) {
      next
    }
    q <- critical(d)
    out$lower[[j]] <- 1 - r[[j]] * exp(q * se)
    out$upper[[j]] <- 1 - r[[j]] * exp(-q * se)
    out$df[[j]] <- d
  }
  if (any(alone_boundary)) {
    return(out)
  }
  # The sum of the effects is m minus the sum of the ratios, whose interval is
  # taken on the log scale, as that of each ratio: the sum stays below m.
  total <- sum(r)
  g <- c(-total * dmu[[1L]] / mu0, dmu[-1L] / mu0) / total
  se <- sqrt(max(drop(crossprod(g, V %*% g)), 0))
  d <- contrast_df(g)
  if (!is.na(d) && d >= 1) {
    q <- critical(d)
    out$sum_lower <- m - total * exp(q * se)
    out$sum_upper <- m - total * exp(-q * se)
  }
  out
}

empty_effects <- function(design) {
  out <- data.frame(response = character(), context = character(), stringsAsFactors = FALSE)
  for (b in design$spec$by) {
    out[[b]] <- character()
  }
  cbind(out, data.frame(
    time = character(), set_id = character(), order = integer(), stressors = character(),
    scale = character(), stressor = character(), stressor_level = character(),
    estimate = numeric(), lower = numeric(), upper = numeric(), level = numeric(),
    df = numeric(), direction = character(), status = character(), reason = character(),
    note = character(), stringsAsFactors = FALSE
  ))
}

empty_effect_sets <- function(design) {
  out <- data.frame(response = character(), context = character(), stringsAsFactors = FALSE)
  for (b in design$spec$by) {
    out[[b]] <- character()
  }
  cbind(out, data.frame(
    time = character(), set_id = character(), order = integer(), stressors = character(),
    levels = character(), background = character(), scale = character(), type = character(),
    sum = numeric(), sum_lower = numeric(), sum_upper = numeric(), level = numeric(),
    near_limit = logical(), sum_vs_limit = character(), increase = logical(),
    opposing = logical(), similar_predictions = logical(), null_models = character(),
    status = character(), reason = character(), stringsAsFactors = FALSE
  ))
}

user_effects_table <- function(effects) {
  if (is.null(effects)) {
    return(NULL)
  }
  data.frame(
    stressor = names(effects), estimate = unname(effects), lower = NA_real_,
    upper = NA_real_, level = NA_real_, df = NA_real_,
    direction = effect_directions(unname(effects), NULL, NULL, use_points = TRUE),
    status = ifelse(is.na(effects), "unknown", "given"), reason = "", note = "",
    stringsAsFactors = FALSE
  )
}

# Tables of the result ----------------------------------------------------------

sp_models_table <- function(reach, criteria, evals, design) {
  codes <- sp_codes()
  reached <- lapply(reach, function(g) unique(g$null_model))
  all_reached <- unique(unlist(reached))
  steps_reached <- sort(unique(unlist(lapply(reach, `[[`, "step"))))
  n_reached <- vapply(codes, function(m) {
    sum(vapply(reached, function(r) m %in% r, logical(1)))
  }, integer(1), USE.NAMES = FALSE)
  status <- ifelse(n_reached == 0L, "excluded",
                   if (length(all_reached) == 1L) "suggested" else "possible")
  # Criteria that depend on the effects, when every evaluation agrees on them.
  effect_criteria_names <- c("near_limit", "sum", "increase", "opposing")
  merged <- lapply(stats::setNames(effect_criteria_names, effect_criteria_names), function(n) {
    v <- unique(evals[[n]])
    if (length(v) == 1L) v else NA
  })
  cr <- c(criteria, merged)
  checked <- design_check(design)
  data.frame(
    null_model = codes,
    name = unname(sp_names()[codes]),
    status = status,
    reached = n_reached,
    evaluations = length(reach),
    reason = vapply(codes, sp_reason, character(1), steps_reached = steps_reached, cr = cr,
                    USE.NAMES = FALSE),
    design = checked$design,
    design_note = checked$design_note,
    note = c(
      "",
      paste(
        "The effect addition of ecotoxicology: \"BL\" (Bliss independence) gives the",
        "same prediction, with the deviation relative to the reference cell."
      ),
      "",
      paste(
        "must assumes a beta distribution of the stress capacities (`sam_shape` in",
        "estimate_interactions())."
      ),
      ""
    ),
    stringsAsFactors = FALSE
  )
}

sp_value_label <- function(x) {
  if (is.logical(x)) {
    ifelse(is.na(x), "unknown", ifelse(x, "yes", "no"))
  } else {
    ifelse(is.na(x), "unknown", x)
  }
}

sp_questions <- function(near_limit) {
  c(
    similar_action = "The stressors act similarly (same mode of action)",
    specific_action = "One stressor acts specifically on physiology",
    mortality = "The effect is mortality",
    stressor_effect = "The stressor-effect relationship is known",
    sensitivity_distribution = "The sensitivity distribution is known",
    correlation = "Correlation of the sensitivities",
    near_limit = sprintf(
      "The effect of one stressor is close to the effect limit (at least %s%%)",
      format(100 * near_limit)
    ),
    sum = "The sum of the effects compared with the effect limit",
    increase = "A stressor alone increases the response",
    opposing = "Stressors alone with opposite effects"
  )
}

sp_criteria_table <- function(criteria, evals, origin, near_limit, complement = character()) {
  user <- vapply(criteria, sp_value_label, character(1))
  effect_value <- function(v) {
    labels <- sp_value_label(v)
    # Estimated, but the interval spans the threshold.
    if (identical(origin, "estimated")) {
      labels[labels == "unknown"] <- "uncertain"
    }
    u <- unique(labels)
    if (length(u) == 1L) {
      return(list(value = u, detail = ""))
    }
    counts <- table(factor(labels, levels = unique(labels)))
    list(value = "varies", detail = sprintf(
      "%s of %d contrast sets and times",
      paste(sprintf("%s in %d", names(counts), as.integer(counts)), collapse = ", "),
      length(labels)
    ))
  }
  # Effects given by the user may leave a criterion unknown (an NA effect).
  effect_source <- function(v) {
    if (identical(origin, "given") && length(v) == 1L && is.na(v)) "unknown" else origin
  }
  names_effects <- c("near_limit", "sum", "increase", "opposing")
  values <- lapply(names_effects, function(n) {
    v <- effect_value(evals[[n]])
    # The guidelines were applied to the complement of these proportions.
    if (length(complement) > 0L) {
      where <- sprintf("for %s, on the complement",
                       format_items(quote_items(complement), max = Inf))
      v$detail <- if (nzchar(v$detail)) paste0(v$detail, "; ", where) else where
    }
    v
  })
  data.frame(
    criterion = names(sp_questions(near_limit)),
    question = unname(sp_questions(near_limit)),
    value = c(unname(user), vapply(values, `[[`, character(1), "value")),
    source = c(ifelse(is.na(unlist(criteria)), "unknown", "given"),
               vapply(names_effects, function(n) effect_source(evals[[n]]), character(1),
                      USE.NAMES = FALSE)),
    detail = c(rep("", length(criteria)), vapply(values, `[[`, character(1), "detail")),
    stringsAsFactors = FALSE
  )
}

# Unknown criteria whose value changes the null models kept, with the models
# kept for each of their values. The criteria that depend on the effects only
# count when there is a single evaluation (effects given or unknown).
sp_decisive <- function(criteria, evals) {
  domains <- sp_domains()
  candidates <- names(criteria)[vapply(criteria, is.na, logical(1))]
  if (nrow(evals) == 1L) {
    effect_names <- c("near_limit", "sum", "increase", "opposing")
    candidates <- c(candidates, effect_names[vapply(effect_names, function(n) {
      is.na(evals[[n]][[1L]])
    }, logical(1))])
  }
  rows <- list()
  for (cn in candidates) {
    kept <- vapply(domains[[cn]], function(v) {
      models <- unlist(lapply(seq_len(nrow(evals)), function(i) {
        cr <- sp_evaluation(criteria, evals, i)
        cr[[cn]] <- v
        sp_scenarios(cr)$null_model
      }))
      # A value that contradicts the known criteria gives no model.
      if (length(models) == 0L) NA_character_ else paste(sp_order(models), collapse = ", ")
    }, character(1), USE.NAMES = FALSE)
    possible <- !is.na(kept)
    if (length(unique(kept[possible])) > 1L) {
      rows[[length(rows) + 1L]] <- data.frame(
        criterion = cn, value = sp_value_label(domains[[cn]])[possible],
        null_models = kept[possible], stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_base(rows, template = data.frame(criterion = character(), value = character(),
                                              null_models = character(),
                                              stringsAsFactors = FALSE))
}

# What the design allows ---------------------------------------------------------

design_check <- function(design) {
  codes <- sp_codes()
  out <- data.frame(design = rep(NA_character_, length(codes)),
                    design_note = rep("", length(codes)), stringsAsFactors = FALSE)
  if (is.null(design)) {
    return(out)
  }
  sets <- design$sets[design$sets$net_identifiable, , drop = FALSE]
  if (nrow(sets) == 0L) {
    out$design <- "no"
    out$design_note <- "The design identifies no contrast set."
    return(out)
  }
  some <- function(ok) if (all(ok)) "yes" else if (any(ok)) "some sets" else "no"
  out$design[codes %in% c("AD", "MU", "SAM")] <- "yes"
  two <- sets$order == 2L
  out$design[codes == "DO"] <- some(two)
  if (!all(two)) {
    out$design_note[codes == "DO"] <- paste(
      "Dominance is defined for two stressors only: sets of three or more",
      "stressors are not estimated."
    )
  }
  ca <- lapply(seq_len(nrow(sets)), function(i) ca_design_check(design, sets[i, , drop = FALSE]))
  ok <- vapply(ca, `[[`, logical(1), "ok")
  out$design[codes == "CA"] <- some(ok)
  notes <- unique(vapply(ca[!ok], `[[`, character(1), "reason"))
  if (any(ok)) {
    alone <- unlist(lapply(ca[ok], `[[`, "levels"))
    fewest <- tapply(alone, names(alone), min)
    few <- fewest[fewest < 5L]
    if (length(few) > 0L) {
      notes <- c(notes, sprintf(
        paste(
          "Sch\u00e4fer and Piggott (2018) suggest at least five levels of a stressor",
          "to establish its stressor-effect relationship; the design has %s alone",
          "(the control included)."
        ),
        and_list(sprintf("%d levels of '%s'", as.integer(few), names(few)))
      ))
    }
  }
  out$design_note[codes == "CA"] <- paste(notes, collapse = " ")
  out
}

# Whether the design allows concentration addition in a set: numeric doses
# above the control and at least two doses of each stressor alone (the other
# stressors at the reference), as in ca_layout(). `levels` counts the levels
# of each stressor alone, the control included.
ca_design_check <- function(design, s) {
  spec <- design$spec
  sc <- design$set_cells[[s$set_id]]
  pool <- pool_of(context_pools(design$cells, spec), s$context)
  cells <- design$cells[design$cells$context %in% pool, , drop = FALSE]
  reference <- cells[match(sc$cell[[1L]], cells$cell), spec$stressors, drop = FALSE]
  full <- cells[match(sc$cell[[nrow(sc)]], cells$cell), spec$stressors, drop = FALSE]
  in_set <- strsplit(s$stressors, " x ", fixed = TRUE)[[1L]]
  fail <- function(reason) list(ok = FALSE, reason = reason, levels = integer())
  levels <- integer()
  for (st in in_set) {
    lv <- spec$levels[[st]]
    num <- suppressWarnings(as.numeric(lv))
    ctrl <- suppressWarnings(as.numeric(spec$control[[st]]))
    if (anyNA(num) || is.na(ctrl)) {
      return(fail(sprintf(
        "Concentration addition needs numeric doses: the levels of '%s' are not numbers.",
        st
      )))
    }
    if (!isTRUE(as.numeric(full[[st]][[1L]]) > ctrl)) {
      return(fail(sprintf(
        "Concentration addition needs doses above the control; '%s' is below it in some sets.",
        st
      )))
    }
    alone <- lv[num > ctrl]
    lab <- vapply(alone, function(l) {
      r <- reference
      r[[st]] <- l
      make_cell_label(r, spec$stressors)
    }, character(1), USE.NAMES = FALSE)
    n_alone <- sum(lab %in% cells$cell)
    if (n_alone < 2L) {
      return(fail(sprintf(
        paste(
          "Concentration addition needs at least two doses of each stressor alone;",
          "the design has %d of '%s' in some sets."
        ),
        n_alone, st
      )))
    }
    levels[[st]] <- n_alone + 1L
  }
  list(ok = TRUE, reason = "", levels = levels)
}

# Use in estimate_interactions() --------------------------------------------------

# What estimate_interactions() keeps of a suggestion in its settings.
null_model_choice <- function(x) {
  list(
    guide = x$settings$guide,
    null_models = x$null_models,
    models = x$models[, c("null_model", "name", "status", "reason"), drop = FALSE],
    criteria = x$criteria,
    decisive = x$decisive,
    near_limit = x$settings$near_limit,
    level = x$settings$level,
    effects = x$settings$effects,
    responses = x$settings$responses,
    sensitivity = x$settings$sensitivity,
    complement = x$settings$complement,
    suggested_sensitivity = x$settings$suggested_sensitivity
  )
}

# The effects behind a suggestion estimated from a fitted model depend on
# `sensitivity` (mortality or survival), and the guidelines may have been
# applied to the complement of a proportion that the stressors increase: the
# interactions must be estimated on the scale of the suggestion.
check_choice_sensitivity <- function(choice, sensitivity, call = sys.call(-1)) {
  if (is.null(choice) || !identical(choice$effects, "estimated")) {
    return(invisible(TRUE))
  }
  wanted <- choice$suggested_sensitivity
  same <- if (is.logical(wanted)) {
    identical(isTRUE(wanted), isTRUE(sensitivity))
  } else {
    setequal(as.character(wanted), as.character(sensitivity))
  }
  if (!same) {
    shown <- function(v) {
      if (is.logical(v)) {
        format(isTRUE(v))
      } else if (length(v) == 0L) {
        "NULL"
      } else if (length(v) == 1L) {
        paste0("\"", v, "\"")
      } else {
        paste0("c(", paste0("\"", v, "\"", collapse = ", "), ")")
      }
    }
    abort_input(sprintf(
      paste(
        "The null models were suggested for the interactions estimated with",
        "`sensitivity` = %s%s, and they would be estimated with %s. Give that",
        "`sensitivity`, or call suggest_null_models() again."
      ),
      shown(wanted),
      if (length(choice$complement) > 0L) {
        sprintf(
          paste(
            " (the guidelines were applied to the complement of %s, which the",
            "stressors alone increase)"
          ),
          format_items(quote_items(choice$complement), max = Inf)
        )
      } else {
        ""
      },
      shown(sensitivity)
    ), call = call)
  }
  invisible(TRUE)
}

#' @export
print.must_null_models <- function(x, ...) {
  s <- x$settings
  cat("<must_null_models> after the guidelines of ", s$guide, "\n", sep = "")
  cr <- x$criteria
  cat("Criteria:\n")
  width <- max(nchar(cr$question))
  for (i in seq_len(nrow(cr))) {
    value <- cr$value[[i]]
    if (identical(cr$source[[i]], "estimated")) {
      value <- paste(value, if (identical(value, "uncertain")) {
        "(estimated; an interval spans the threshold)"
      } else {
        "(estimated)"
      })
    }
    if (nzchar(cr$detail[[i]])) {
      value <- paste0(value, ": ", cr$detail[[i]])
    }
    cat(sprintf("  %-*s  %s\n", width, cr$question[[i]], value))
  }
  if (!is.null(x$sets)) {
    ok <- sum(x$sets$status == "ok")
    n_resp <- length(unique(x$sets$response))
    cat(strwrap(sprintf(
      "Effects of the stressors alone: estimated in %d of %d contrast sets and times%s, %s.",
      ok, nrow(x$sets),
      if (n_resp > 1L) sprintf(" of %d responses", n_resp) else "",
      if (is.null(s$level)) {
        "from the point estimates"
      } else {
        sprintf("with %s%% intervals", format(100 * s$level))
      }
    ), exdent = 2L), sep = "\n")
    similar <- sum(x$sets$similar_predictions %in% TRUE)
    if (similar > 0L) {
      cat(strwrap(sprintf(
        paste(
          "In %d of them the point estimates put one effect close to the limit or",
          "at about ten times the others: the null models predict similar joint",
          "effects there."
        ),
        similar
      ), indent = 2L, exdent = 2L), sep = "\n")
    }
    ok_sets <- x$sets[x$sets$status == "ok", , drop = FALSE]
    if (nrow(ok_sets) > 0L) {
      types <- ifelse(is.na(ok_sets$type), "unknown", ok_sets$type)
      counts <- sort(table(types), decreasing = TRUE)
      cat(strwrap(sprintf(
        "Directions of the effects (types of Piggott et al. 2015): %s.",
        paste(sprintf("%s in %d", names(counts), as.integer(counts)), collapse = ", ")
      ), indent = 2L, exdent = 2L), sep = "\n")
    }
    if (any(ok_sets$increase %in% TRUE)) {
      cat(strwrap(paste(
        "Where a stressor alone increases the response, the increase offsets the",
        "decreases in the sum of the effects, and concentration addition and the",
        "stress addition model are not reached."
      ), indent = 2L, exdent = 2L), sep = "\n")
    }
    if (length(s$complement) > 0L) {
      cat(strwrap(sprintf(
        paste(
          "The stressors alone increase %s, a proportion, and never decrease it:",
          "the guidelines were applied to its complement (1 - p), where the",
          "increases are reductions towards the limit 0."
        ),
        format_items(quote_items(s$complement), max = Inf)
      ), indent = 2L, exdent = 2L), sep = "\n")
    }
  } else if (!is.null(x$effects)) {
    cat(strwrap(sprintf(
      "Effects of the stressors alone (given): %s; type: %s.",
      paste(sprintf("%s = %s (%s)", x$effects$stressor, format(x$effects$estimate),
                    x$effects$direction), collapse = ", "),
      sp_value_label(effect_type(x$effects$direction))
    ), exdent = 2L), sep = "\n")
  } else {
    cat("Effects of the stressors alone: unknown\n")
  }
  cat("Null models:\n")
  m <- x$models
  for (i in seq_len(nrow(m))) {
    head <- sprintf("%s (%s): %s", m$null_model[[i]], m$name[[i]], m$status[[i]])
    if (m$evaluations[[i]] > 1L && m$reached[[i]] > 0L) {
      head <- sprintf("%s in %d of %d", head, m$reached[[i]], m$evaluations[[i]])
    }
    text <- paste0(head, ". ", m$reason[[i]])
    if (!identical(m$status[[i]], "excluded")) {
      if (!is.na(m$design[[i]]) && !identical(m$design[[i]], "yes")) {
        text <- paste(text, sprintf("Allowed by the design: %s.", m$design[[i]]),
                      m$design_note[[i]])
      } else if (nzchar(m$design_note[[i]])) {
        text <- paste(text, m$design_note[[i]])
      }
      if (nzchar(m$note[[i]])) {
        text <- paste(text, m$note[[i]])
      }
    }
    cat(strwrap(text, indent = 2L, exdent = 4L), sep = "\n")
  }
  if (nrow(x$decisive) > 0L) {
    labels <- c(similar_action = "Similar mode of action",
                specific_action = "Specific action on physiology",
                mortality = "Mortality", stressor_effect = "Stressor-effect relationship known",
                sensitivity_distribution = "Sensitivity distribution known",
                correlation = "Correlation of the sensitivities",
                near_limit = "An effect close to the limit",
                sum = "Sum of the effects vs the limit",
                increase = "A stressor alone increases the response",
                opposing = "Opposite effects")
    cat("Unknown criteria that would change the choice:\n")
    for (cn in unique(x$decisive$criterion)) {
      d <- x$decisive[x$decisive$criterion == cn, , drop = FALSE]
      cat(strwrap(sprintf("%s (`%s`): %s", labels[[cn]], cn,
                          paste(sprintf("%s -> %s", d$value, d$null_models), collapse = "; ")),
                  indent = 2L, exdent = 4L), sep = "\n")
    }
  }
  quoted <- function(v) {
    if (length(v) == 1L) {
      paste0("\"", v, "\"")
    } else {
      paste0("c(", paste0("\"", v, "\"", collapse = ", "), ")")
    }
  }
  extra <- ""
  if (length(s$complement) > 0L) {
    extra <- sprintf(", sensitivity = %s", if (is.logical(s$suggested_sensitivity)) {
      format(s$suggested_sensitivity)
    } else if (length(s$suggested_sensitivity) == 0L) {
      "NULL"
    } else {
      quoted(s$suggested_sensitivity)
    })
  }
  cat(sprintf("For estimate_interactions(): null_models = %s%s\n", quoted(x$null_models),
              extra))
  invisible(x)
}
