#' Estimate single, combined and interaction effects under null models
#'
#' From the model selected for each response by [fit_responses()], computes
#' for every identifiable contrast set of the design (see
#' [describe_design()]) the effect of each stressor alone, the effect of
#' their combination and its deviation from what each null model predicts.
#' This is post-estimation inference: one model per response, several null
#' models tested from its predictions (Madge Pimentel et al. 2025).
#'
#' @param fit A `must_fit` object from [fit_responses()].
#' @param null_models Null models to evaluate: any of `"AD"` (simple
#'   addition), `"MU"` (multiplicative) and `"DO"` (dominance). Default: all
#'   three.
#' @param inference `"simulation"` (default; Krinsky–Robb: draws of the fixed
#'   effects from their estimated distribution) or `"delta"` (delta method;
#'   only for `"AD"` and `"MU"`).
#' @param correction `"small_sample"` (default) or `"none"`. With
#'   `"small_sample"`, the covariance of the fixed effects comes from the
#'   REML refit of the selected model and the sampling distribution is a t
#'   with the degrees of freedom of the experimental units (see Details).
#'   `"none"` uses the maximum-likelihood covariance and the normal
#'   distribution, as in Madge Pimentel et al. (2025); with few replicates
#'   its tests reject a true null model too often.
#' @param n_sim Number of draws for `inference = "simulation"`, and for the
#'   tests of combinations with only zeros with either inference. Default
#'   10000.
#' @param level Confidence level of the intervals. Default 0.95.
#' @param p_adjust Method of [stats::p.adjust()] used across contrast sets
#'   and sampling times (see Details). Default `"BH"`.
#' @param overall If `TRUE` (default) and the design has two or more sampling
#'   times after any `baseline`, adds for each set the average of each
#'   effect over those times (`time = "overall"`).
#' @param sensitivity `NULL` (default) or the names of responses for which an
#'   increase means sensitivity (mortality, damage) rather than tolerance.
#'   For binomial responses (`trials` in [fit_responses()]) every null model
#'   is then applied to the complement (1 - p, e.g. survival instead of
#'   mortality); for other responses the MU and DO rows are
#'   `"not_applicable"`, because those null models assume that an increase
#'   means tolerance.
#' @param seed `NULL` (default) or a whole number. The simulation runs with
#'   this seed through [withr::with_seed()], which restores the random
#'   number state of the session afterwards. With `NULL`, a seed is drawn
#'   from the session's random number stream; either way it is stored in
#'   `settings$seed`, so the estimates can be repeated exactly.
#'
#' @return An object of class `must_interactions`, a list with
#'   `estimates`, `design`, `settings` and `diagnostics`. `estimates` is a
#'   data frame with one row per response, contrast set, sampling time, null
#'   model and quantity, and columns:
#'   * `response`, `context` and the `by` columns, `time` (`NA` without
#'     time, `"overall"` for averages over times), `set_id`, `order`,
#'     `stressors`, `levels`, `background` (levels of the stressors held
#'     outside the set): what is estimated.
#'   * `null_model` (`"AD"`, `"MU"`, `"DO"`) and `scale` (`"response"` or
#'     `"log(response)"`): what "no interaction" means and on which scale.
#'   * `component`: `"single"` (one stressor alone vs the reference cell;
#'     its name in `stressor`), `"combined"` (all stressors of the set vs the
#'     reference) or `"interaction"`.
#'   * `estimand` for interactions: `"net"` (deviation of the combination
#'     from the null-model combination of the single effects) or `"pure"`
#'     (highest-order interaction once all lower-order ones are removed;
#'     only for three or more stressors).
#'   * `estimate`, `se` (standard error; with `inference = "simulation"`,
#'     the standard deviation of the simulated estimates before they are
#'     scaled to a t distribution, comparable with the delta method),
#'     `lower`, `upper`, `level`, `df` (degrees of freedom of the t
#'     distribution; `Inf` for the normal), `p_value`, `p_adjusted`, `mc_se`
#'     (Monte Carlo standard error of a simulated p-value).
#'   * `observed` and `expected` (interactions only, response scale): the
#'     fitted mean of the combination and the value the null model predicts.
#'   * `inference`, `status` (`"ok"`, `"boundary"` (tested from the
#'     probability of a combination with only zeros, see Details),
#'     `"not_estimable"`, `"not_identifiable"`, `"not_defined"` or
#'     `"not_applicable"`), `reason`, `note` (`"expected_below_zero"` when
#'     simple addition predicts an impossible negative value for a
#'     non-negative response; `"no_point_expectation"` when dominance gives no
#'     single expected value; `"confounded_with_groups"` for the sets flagged
#'     in `sets$confounded` by [describe_design()]; `"complement"` for
#'     responses analysed as 1 - p; several notes are separated by `"; "`).
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
#' successes) and the model has a log or logit link, the null model still
#' predicts a mean E for it from the other cells. The p-value is twice the
#' probability of observing only zeros if the mean were E, averaged over the
#' draws of E (`status = "boundary"`; Poisson, NB2, Tweedie, binomial and
#' beta-binomial models). Random intercepts are integrated out: observations
#' that share a level of a random effect share one intercept with the
#' variance of all random terms, which is exact with a single grouping (units
#' or rooms, say) and conservative with nested or crossed ones. The same test
#' is used for a single stressor whose cell has only zeros (against the
#' reference cell). The estimate is computed with the observed mean (zero),
#' which makes MU estimates `-Inf`; intervals are not given. If the reference
#' cell, or a single-stressor cell that an interaction needs, has only zeros,
#' that contrast cannot be measured and the row is `"not_estimable"`.
#'
#' **Multiple testing.** p-values are adjusted within each response, null
#' model, component, stressor and estimand, across contrast sets, sampling
#' times and `by` groups; averages over times (`"overall"`) are adjusted
#' separately.
#'
#' @section Assumptions:
#' * Means are conditional predictions (random effects at zero). With a log
#'   link and a random intercept, marginal means differ by a constant
#'   factor: MU estimates and all p-values are unchanged, AD and DO
#'   estimates change scale.
#' * MU and DO only represent their mechanistic assumptions when an
#'   increase in the response means tolerance (survival, abundance) rather
#'   than sensitivity (mortality); declare the latter in `sensitivity`. AD
#'   may predict impossible values below zero, in which case a deviation
#'   from it does not imply a mechanistic interaction (Madge Pimentel et al.
#'   2025).
#' * A non-significant interaction is a failure to detect a deviation, not
#'   evidence of additivity.
#'
#' @section Missing values and zeros:
#' MU needs positive means: if a fitted or simulated mean is not positive
#' (possible with an identity link), the MU rows of that set are
#' `"not_estimable"`. Treatment combinations without values of the response
#' (at a time) make the sets that need them `"not_estimable"` (at that
#' time). Times declared as `baseline` in [describe_design()] are not
#' estimated.
#'
#' @section Reproducibility:
#' With `inference = "simulation"` results depend on the random number
#' stream: the seed used is in `settings$seed`, and `mc_se` gives the Monte
#' Carlo error of each simulated p-value. The same seed may not give
#' identical draws in other versions of R or with other linear algebra
#' libraries.
#'
#' @references
#' Krinsky, I., & Robb, A. L. (1986). On approximating the statistical
#' properties of elasticities. *The Review of Economics and Statistics*,
#' 68(4), 715–719. \doi{10.2307/1924536}
#'
#' Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
#' J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025).
#' Hypothesis-driven research on multiple stressors: An analytical framework
#' for stressor interactions. *Ecology and Evolution*, 15, e71959.
#' \doi{10.1002/ece3.71959}
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin*, 2(6), 110–114.
#' \doi{10.2307/3002019}
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
#' est <- estimate_interactions(fits, n_sim = 2000, seed = 1)
#' est
#' est$estimates[est$estimates$component == "interaction",
#'   c("null_model", "estimate", "lower", "upper", "df", "p_adjusted")]
#' # The procedure of Madge Pimentel et al. (2025), without the correction
#' estimate_interactions(fits, n_sim = 2000, seed = 1, correction = "none")
#' @export
estimate_interactions <- function(fit, null_models = c("AD", "MU", "DO"),
                                  inference = "simulation", correction = "small_sample",
                                  n_sim = 10000L, level = 0.95, p_adjust = "BH",
                                  overall = TRUE, sensitivity = NULL, seed = NULL) {
  if (!inherits(fit, "must_fit")) {
    abort_input(sprintf(
      "`fit` must be the output of fit_responses(), not %s.", describe_value(fit)
    ))
  }
  null_models <- check_null_models(null_models)
  check_choice(inference, c("simulation", "delta"))
  check_choice(correction, c("small_sample", "none"))
  if (identical(inference, "delta") && "DO" %in% null_models) {
    abort_input(paste(
      "The delta method is not valid for the dominance null model, which is",
      "not differentiable where the ranking of the means changes (Madge",
      "Pimentel et al. 2025). Use inference = \"simulation\" or drop \"DO\"",
      "from `null_models`."
    ))
  }
  check_number(n_sim, lower = 100)
  if (n_sim != round(n_sim)) {
    abort_input("`n_sim` must be a whole number.")
  }
  check_number(level, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  check_choice(p_adjust, stats::p.adjust.methods)
  check_flag(overall)
  sensitivity <- check_sensitivity(sensitivity, names(fit$fits))
  check_seed(seed)
  seed <- resolve_seed(seed)
  run <- function() {
    lapply(fit$fits, function(f) {
      estimate_one_response(
        f, design = fit$design, null_models = null_models, inference = inference,
        correction = correction, n_sim = as.integer(n_sim), level = level,
        overall = overall, sensitive = f$response %in% sensitivity
      )
    })
  }
  parts <- withr::with_seed(seed, run())
  estimates <- bind_rows_base(parts, template = empty_estimates(fit$design))
  estimates <- adjust_p_values(estimates, p_adjust)
  failed <- vapply(fit$fits, function(f) !identical(f$status, "ok"), logical(1))
  structure(
    list(
      estimates = estimates,
      design = fit$design,
      settings = list(
        null_models = null_models,
        inference = inference,
        correction = correction,
        n_sim = as.integer(n_sim),
        level = level,
        p_adjust = p_adjust,
        overall = overall,
        sensitivity = sensitivity,
        seed = seed,
        prediction = "conditional (random effects at zero)",
        time_structure = fit$settings$time_structure,
        families = vapply(fit$fits, function(f) {
          if (identical(f$status, "ok")) paste0(f$family, " (", f$link, ")") else NA_character_
        }, character(1))
      ),
      diagnostics = list(
        failed_responses = data.frame(
          response = names(fit$fits)[failed],
          reason = vapply(fit$fits[failed], `[[`, character(1), "reason"),
          stringsAsFactors = FALSE
        ),
        reml = vapply(fit$fits, function(f) !is.null(f$model_reml), logical(1))
      )
    ),
    class = "must_interactions"
  )
}

check_null_models <- function(null_models, call = sys.call(-1)) {
  valid <- c("AD", "MU", "DO")
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
    expected = numeric(), inference = character(), status = character(),
    reason = character(), note = character(), stringsAsFactors = FALSE
  ))
}

estimate_one_response <- function(f, design, null_models, inference, correction, n_sim,
                                  level, overall, sensitive) {
  if (!identical(f$status, "ok")) {
    return(NULL)
  }
  model <- if (identical(correction, "small_sample") && !is.null(f$model_reml)) {
    f$model_reml
  } else {
    f$model
  }
  complement <- sensitive && identical(f$kind, "binomial")
  pred <- prepare_predictions(f, model, n_sim, complement = complement)
  key <- paste_key(f$grid$.context, f$grid$.cell, f$grid$.time)
  sets <- design$sets[design$sets$net_identifiable, , drop = FALSE]
  baseline <- design$spec$baseline
  rows <- list()
  for (i in seq_len(nrow(sets))) {
    s <- sets[i, , drop = FALSE]
    sc <- design$set_cells[[s$set_id]]
    times <- unique(as.character(f$grid$.time[f$grid$.context == s$context]))
    times <- times[!times %in% baseline]
    df <- set_df(f$df_info, s, correction)
    per_time <- lapply(times, function(tt) {
      idx <- match(paste_key(s$context, sc$cell, tt), key)
      set_estimates(s, sc, idx, tt, pred, null_models, inference, level, f, design, df,
                    sensitive)
    })
    rows <- c(rows, lapply(per_time, `[[`, "rows"))
    if (overall && length(times) >= 2L) {
      rows[[length(rows) + 1L]] <- overall_estimates(s, per_time, pred, inference, level, f,
                                                     design)
    }
  }
  bind_rows_base(rows)
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
    expected = values$expected, inference = values$inference, status = values$status,
    reason = values$reason, note = values$note, stringsAsFactors = FALSE
  ))
}

empty_values <- function(status, reason, level, inference) {
  list(estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
       level = level, df = NA_real_, p_value = NA_real_, mc_se = NA_real_,
       observed = NA_real_, expected = NA_real_, inference = inference, status = status,
       reason = reason, note = "", draws = NULL, grad = NULL)
}

set_estimates <- function(s, sc, idx, time_value, pred, null_models, inference, level, f,
                          design, df, sensitive) {
  m <- s$order
  in_set <- strsplit(s$stressors, " x ", fixed = TRUE)[[1L]]
  w <- contrast_weights(m)
  missing_now <- sc$cell[sc$for_net & is.na(idx)]
  has_time <- !is.na(time_value)
  control <- design$spec$control
  rows <- list()
  parts <- list()
  for (nm in null_models) {
    for (sp in component_specs(nm, m, in_set, w)) {
      values <- if (identical(nm, "DO") && m > 2L) {
        empty_values("not_defined", "The dominance null model is defined for two stressors only.",
                     level, inference)
      } else if (identical(sp$estimand, "pure") && !isTRUE(s$pure_identifiable)) {
        lacking <- strsplit(s$missing_for_pure, "; ", fixed = TRUE)[[1L]]
        empty_values("not_identifiable", sprintf(
          paste(
            "The design lacks lower-order combinations (%s): the pure",
            "interaction cannot be separated from them."
          ),
          format_items(quote_items(short_cells(lacking, control)), max = 3L)
        ), level, inference)
      } else if (sensitive && !identical(f$kind, "binomial") && nm %in% c("MU", "DO")) {
        empty_values("not_applicable", paste(
          "An increase in this response means sensitivity: the multiplicative",
          "and dominance null models assume that an increase means tolerance",
          "(Madge Pimentel et al. 2025). Analyse a response where higher means",
          "better, such as survival."
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
      } else if (identical(inference, "simulation") && nzchar(pred$sim_problem())) {
        empty_values("not_estimable", pred$sim_problem(), level, inference)
      } else {
        compute_component(nm, sp, idx, m, pred, inference, level, f, df, sc, design)
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

compute_component <- function(null_model, sp, idx, m, pred, inference, level, f, df, sc,
                              design) {
  use <- if (identical(null_model, "DO")) 1:4 else which(sp$weights != 0)
  cols <- idx[use]
  if (any(pred$at_boundary[cols])) {
    return(boundary_component(null_model, sp, idx, use, m, pred, level, f, df, sc, design))
  }
  mu <- pred$mu[cols]
  weights <- if (identical(null_model, "DO")) dominance_weights(mu) else sp$weights[use]
  df <- contrast_df(df, null_model, weights, cols, pred, f)
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
  if (identical(null_model, "MU")) {
    bad_point <- any(!is.finite(mu) | mu <= 0)
    bad_sims <- !is.null(draws) && any(draws <= 0)
    if (bad_point || bad_sims) {
      return(empty_values(
        "not_estimable",
        "The multiplicative model needs positive means; some fitted or simulated means are not positive.",
        level, inference
      ))
    }
  }
  grad <- NULL
  component_draws <- NULL
  normal_draws <- NULL
  if (identical(null_model, "DO")) {
    estimate <- dominance_interaction(matrix(mu, nrow = 1L))
    component_draws <- dominance_interaction(draws)
    normal_draws <- dominance_interaction(normal)
    stats_list <- summarise_draws(component_draws, level)
    stats_list$se <- stats::sd(normal_draws)
  } else {
    estimate <- sum(weights * transform_means(mu, null_model))
    if (identical(inference, "simulation")) {
      component_draws <- drop(transform_means(draws, null_model) %*% weights)
      normal_draws <- drop(transform_means(normal, null_model) %*% weights)
      stats_list <- summarise_draws(component_draws, level)
      stats_list$se <- stats::sd(normal_draws)
    } else {
      g <- delta_gradient(weights, mu, pred$dmu[cols], null_model)
      se <- sqrt(max(drop(t(g) %*% pred$V_eta[cols, cols, drop = FALSE] %*% g), 0))
      stats_list <- summarise_wald(estimate, se, df, level)
      grad <- numeric(length(pred$mu))
      grad[cols] <- g
    }
  }
  observed <- NA_real_
  expected <- NA_real_
  note <- ""
  if (identical(sp$component, "interaction")) {
    observed <- pred$mu[[idx[[length(idx)]]]]
    if (identical(sp$estimand, "pure")) {
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
  }
  list(
    estimate = estimate, se = stats_list$se, lower = stats_list$lower,
    upper = stats_list$upper, level = level, df = df, p_value = stats_list$p_value,
    mc_se = if (identical(inference, "simulation")) mc_error(stats_list$p_value, pred$n_sim) else NA_real_,
    observed = observed, expected = expected, inference = inference, status = "ok",
    reason = "", note = note, draws = component_draws, normal_draws = normal_draws,
    grad = grad, cols = cols, weights = weights, null_model = null_model
  )
}

# Simulated values of a contrast from draws of the means of its cells.
contrast_draws <- function(draws, null_model, weights) {
  if (identical(null_model, "DO")) {
    dominance_interaction(draws)
  } else {
    drop(transform_means(draws, null_model) %*% weights)
  }
}

# Local weights of the dominance interaction min(yAB, y0) - min(yA, yB) at the
# point estimates: which cells it depends on.
dominance_weights <- function(mu) {
  c(as.numeric(mu[[4L]] >= mu[[1L]]), -as.numeric(mu[[2L]] <= mu[[3L]]),
    -as.numeric(mu[[2L]] > mu[[3L]]), as.numeric(mu[[4L]] < mu[[1L]]))
}

# Degrees of freedom of a contrast: those of its set or, with one dispersion
# per cell (set_df() gives NA), the Welch-Satterthwaite approximation from the
# variance contributions of its cells.
contrast_df <- function(df, null_model, weights, cols, pred, f) {
  if (!is.na(df)) {
    return(df)
  }
  scale <- if (identical(null_model, "MU")) "MU" else "AD"
  g <- delta_gradient(weights, pred$mu[cols], pred$dmu[cols], scale)
  v <- g^2 * diag(pred$V_eta)[cols]
  units <- f$df_info$cell_units[as.character(pred$grid$.cell_time[cols])]
  welch_df(v, as.numeric(units))
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
  full_pos <- length(idx)
  target <- if (identical(sp$component, "single")) which(sp$weights == 1) else full_pos
  at <- use[pred$at_boundary[idx[use]]]
  control <- design$spec$control
  if (!identical(as.integer(at), as.integer(target))) {
    lacking <- sc$cell[setdiff(at, target)]
    return(empty_values("not_estimable", sprintf(
      paste(
        "%s had only zeros: with the response eliminated in the reference cell",
        "or by a stressor alone, this contrast cannot be measured."
      ),
      format_items(quote_items(short_cells(lacking, control)), max = 3L)
    ), level, inference))
  }
  target_col <- idx[[target]]
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
  others <- setdiff(seq_along(idx), target)
  other_weights <- if (identical(null_model, "DO")) rep(1, length(others)) else sp$weights[others]
  df_used <- contrast_df(df, null_model, other_weights, idx[others], pred, f)
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
  draws_other <- pred$mu_draws(idx[others], df_used)
  if (is.null(draws_other)) {
    return(empty_values("not_estimable", pred$sim_problem(), level, inference))
  }
  M <- matrix(NA_real_, nrow = nrow(draws_other), ncol = length(idx))
  M[, others] <- draws_other
  expected_draws <- expected_full(M, sp$component, sp$estimand, null_model, m, sp$weights)
  mu_point <- pred$mu[idx]
  mu_point[[target]] <- NA_real_
  expected <- expected_full(matrix(mu_point, nrow = 1L), sp$component, sp$estimand,
                            null_model, m, sp$weights)
  # The probability is computed on the scale of the data (not the complement).
  mean_data <- if (pred$complement) 1 - expected_draws else expected_draws
  prob <- prob_all_boundary(mean_data, rows, params$family, params$link, params$disp,
                            params$power, params$sd_re, type, factors = params$factors)
  ok <- is.finite(prob)
  p_value <- min(1, 2 * mean(prob[ok]))
  mu_obs <- pred$mu[idx]
  mu_obs[[target]] <- observed
  estimate <- if (identical(null_model, "DO")) {
    dominance_interaction(matrix(mu_obs, nrow = 1L))
  } else {
    sum(sp$weights[use] * transform_means(mu_obs[use], null_model))
  }
  is_interaction <- identical(sp$component, "interaction")
  list(
    estimate = estimate, se = NA_real_, lower = NA_real_, upper = NA_real_, level = level,
    df = df_used, p_value = p_value,
    mc_se = if (p_value < 1) 2 * stats::sd(prob[ok]) / sqrt(sum(ok)) else NA_real_,
    observed = if (is_interaction) observed else NA_real_,
    expected = if (is_interaction) expected[[1L]] else NA_real_,
    inference = inference, status = "boundary",
    reason = sprintf(
      paste(
        "%s had only %s: the p-value is twice the probability of that outcome if",
        "the null model held, averaged over the uncertainty of the expected mean;",
        "no interval is given."
      ),
      quote_items(short_cells(sc$cell[[target]], control)),
      if (identical(type, "zero")) "zeros" else "successes"
    ),
    note = "", draws = NULL, grad = NULL
  )
}

# Averages of each effect over the sampling times of a set.
overall_estimates <- function(s, per_time, pred, inference, level, f, design) {
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
      combine_over_times(parts, pred, inference, level)
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

combine_over_times <- function(parts, pred, inference, level) {
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
      lapply(parts, function(p) {
        contrast_draws(pred$mu_draws(p$cols, df), p$null_model, p$weights)
      })
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
  list(estimate = estimate, se = stats_list$se, lower = stats_list$lower,
       upper = stats_list$upper, level = level, df = df, p_value = stats_list$p_value,
       mc_se = mc, observed = NA_real_, expected = NA_real_, inference = inference,
       status = "ok", reason = "", note = "", draws = NULL, grad = NULL)
}

adjust_p_values <- function(estimates, method) {
  estimates$p_adjusted <- rep(NA_real_, nrow(estimates))
  if (nrow(estimates) > 0L) {
    group <- paste_key(estimates$response, estimates$null_model, estimates$component,
                       estimates$stressor, estimates$estimand,
                       estimates$time %in% "overall")
    ok <- estimates$status %in% c("ok", "boundary") & is.finite(estimates$p_value)
    for (g in unique(group[ok])) {
      sel <- ok & group == g
      estimates$p_adjusted[sel] <- stats::p.adjust(estimates$p_value[sel], method = method)
    }
  }
  first <- setdiff(names(estimates), "p_adjusted")
  at <- match("p_value", first)
  estimates[, c(first[seq_len(at)], "p_adjusted", first[-seq_len(at)]), drop = FALSE]
}

#' @export
print.must_interactions <- function(x, n = 12L, ...) {
  s <- x$settings
  cat(sprintf(
    "<must_interactions> %d response(s); null models: %s\n",
    length(s$families), paste(s$null_models, collapse = ", ")
  ))
  cat(sprintf(
    "Inference: %s, %s; %s%% intervals; p-values adjusted with %s; seed %s\n",
    if (identical(s$inference, "simulation")) sprintf("Krinsky-Robb simulation (%d draws)", s$n_sim) else "delta method",
    if (identical(s$correction, "small_sample")) "small-sample correction (REML, t)" else "no small-sample correction",
    format(100 * s$level), s$p_adjust, format(s$seed)
  ))
  cat("Scales: AD = differences (response scale); MU = log ratios; DO = dominance (response scale)\n")
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
      status = int$status,
      stringsAsFactors = FALSE
    )
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
