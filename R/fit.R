#' Fit a GLMM to each response, choosing its error distribution
#'
#' Fits, for every response, one generalised linear mixed model per candidate
#' error distribution with [glmmTMB::glmmTMB()], and keeps the one that fits
#' the data best. The null models are not involved here: following the
#' post-estimation approach, the model is chosen for the data and the null
#' models are tested afterwards from its predictions
#' ([estimate_interactions()]).
#'
#' @param data The data frame given to [describe_design()] (same treatment
#'   combinations; other columns may be added).
#' @param responses Character vector with the names of one or more numeric
#'   response columns.
#' @param design A `must_design` object from [describe_design()]. Designs
#'   with checks of severity `"problem"` are refused.
#' @param families `NULL` (default) to choose the candidates from each
#'   response (see Details); a named list of family objects used for every
#'   response, e.g. `list(gaussian = gaussian(), gamma_log = Gamma("log"))`;
#'   or a named list with one such list per response.
#' @param response_types `NULL` (default) or a named character vector with
#'   `"count"` or `"continuous"` for some responses, to override the automatic
#'   detection used to choose the default candidates.
#' @param trials `NULL` (default) or a character vector named after
#'   responses that count successes out of a known number of trials
#'   (survivors out of the animals added, emerged out of the eggs), giving
#'   the column with the number of trials of each observation, e.g.
#'   `c(survivors = "n_start")`. These responses are modelled as binomial or
#'   beta-binomial and their means are probabilities.
#' @param offset `NULL` (default) or a character vector named after count
#'   responses giving the column with the sampling effort of each observation
#'   (volume filtered, area or time sampled), e.g. `c(daphnia = "volume_l")`.
#'   The model includes `offset(log(effort))`, only log-link candidates are
#'   used, and means are per unit of effort. Every value of the response
#'   needs an effort.
#' @param time_structure How treatment effects relate to sampling time when
#'   the design has `time`: `"interaction"` (default; a separate mean for
#'   each treatment combination and time, no assumption of constant effects)
#'   or `"additive"` (treatment effects constant over time on the scale of
#'   the link function).
#' @param dispersion Which observations share a dispersion (variance)
#'   parameter: `"common"` (default; one for all), `"by"` (one per `by`
#'   group) or `"cell"` (one per treatment combination and `by` group; the
#'   candidates with a dispersion parameter then need two or more units in
#'   every combination). Ignored for Poisson and binomial candidates, which
#'   have none.
#' @param random `NULL` (default) adds a random intercept for each `groups`
#'   column of the design that can be estimated (see [describe_design()])
#'   and one per experimental unit when some unit has more than one row
#'   (repeated measures or subsamples); `"none"` removes random effects; any
#'   other string is used as the random part of the formula in lme4 syntax,
#'   with column names of `data`, e.g. `"(1 | tank) + (1 | block)"`. Custom
#'   terms can also use the internal columns `.unit` (the units), `.time`
#'   (the sampling times, as a factor) and `.g_<name>` (the groups of each
#'   `groups` column, made unique within `by` groups).
#' @param control `NULL` (default) or the output of
#'   [glmmTMB::glmmTMBControl()], passed to every fit. By default the
#'   iteration and evaluation limits of the optimiser grow with the number of
#'   fixed effects (see Details); the rest of the settings are those of
#'   glmmTMB.
#' @param residual_checks If `TRUE` (default), simulation-based residual
#'   checks of the DHARMa package (uniformity, dispersion and zero inflation)
#'   decide which candidates may be selected. Needs DHARMa installed.
#' @param n_sim_residuals Number of simulations for the residual checks.
#'   Default 500.
#' @param seed `NULL` (default) or a whole number. The residual checks run
#'   with this seed through [withr::with_seed()], which restores the random
#'   number state of the session afterwards. With `NULL`, a seed is drawn
#'   from the session's random number stream; either way it is stored in
#'   `settings$seed`, so the selection can be repeated exactly.
#' @param quiet If `FALSE` (default), a message reports each response as it
#'   is fitted.
#'
#' @return An object of class `must_fit`, a list with:
#'   * `fits`: a list with one element per response, holding `response`,
#'     `status` (`"ok"` or `"failed"`), `reason`, `model` (the selected
#'     `glmmTMB` fit, by maximum likelihood), `model_reml` (the same model
#'     refitted by REML, used for the small-sample correction of
#'     [estimate_interactions()]; `NULL` if the refit failed), `family`
#'     (name of the selected candidate), `link`, `response_type` and
#'     `type_source` (`"detected"`, `"declared"`, `"trials"` or
#'     `"families"`), `selection` (data frame with one row per candidate:
#'     `family`, `link`, `random` (random part used), `converged`, `aic`,
#'     `p_uniformity`, `p_dispersion`, `p_zero_inflation`, `n_tests_passed`,
#'     `residuals_ok`, `selected`, `message`), `checks_passed` (`TRUE`,
#'     `FALSE` if no candidate passed all the residual checks, or `NA` if
#'     they were not run), `formula` (of the selected candidate), `n_obs`,
#'     `n_missing` (rows without a response value) and `grid` (treatment
#'     combinations and times predicted later).
#'   * `design`: the design.
#'   * `settings`: the arguments used, with the seed.
#'
#' @details
#' **Model.** Treatment combinations are coded as one factor (a cell-means
#' model), so any design is handled: full or partial factorials, ray designs
#' or a few combinations of many stressors. With `time_structure =
#' "interaction"` the fixed part is one mean per treatment combination, time
#' and `by` group; with `"additive"` it is one mean per treatment combination
#' and `by` group plus a common effect of time.
#'
#' **Default candidates.** Responses with `trials` are binomial: binomial and
#' beta-binomial, both with a logit link. Otherwise whole non-negative
#' numbers are treated as counts (the `type_source` says the type was
#' detected): Poisson and negative binomial (NB2), both with a log link.
#' Continuous values: Gaussian, Gaussian with log link, Gamma with log link
#' and lognormal if all values are positive; Gaussian, Gaussian with log
#' link and Tweedie if there are zeros; Gaussian only if there are negative
#' values. With an `offset` only the log-link candidates are kept. Use
#' `response_types` when the detection is wrong (a density rounded to whole
#' numbers is not a count) and `families` for other distributions (e.g.
#' [glmmTMB::beta_family()] for proportions without a number of trials, such
#' as cover).
#'
#' **Optimiser.** A cell-means model has one parameter per treatment
#' combination (and time), which can be more than the default limits of the
#' optimiser of glmmTMB (300 iterations, 400 evaluations) allow. By default
#' they are raised to 20 iterations and 30 evaluations per fixed effect
#' (never below glmmTMB's defaults); a fit that reaches them still counts as
#' not converged.
#'
#' **Combinations with only zeros.** On a log link, the mean of a treatment
#' combination (at a time) whose values are all zero is on the boundary of
#' the parameter space (minus infinity on the link scale); on a logit link
#' the same happens when all values are zero or all are successes. For these
#' candidates the mean of such a combination is fixed at that boundary
#' (`map` of glmmTMB), so the rest of the model is estimated normally; their
#' AIC counts the fixed parameters, so that it stays comparable with the
#' identity-link candidates. With `time_structure = "additive"` only
#' combinations with such values at every time are fixed (otherwise their
#' mean is not on the boundary). [estimate_interactions()] tests the
#' contrasts that involve these combinations from the probability of the
#' observed zeros.
#'
#' **Random intercepts at zero.** With few units, the variance of a random
#' intercept is often estimated at zero, which often makes the Hessian
#' singular. If that happens to an automatic random intercept (`random =
#' NULL`), the candidate is refitted without it and the reduced model is kept
#' when it reaches the same maximum likelihood, whether or not the full model
#' converged; the `random` column of `selection` and the `reason` say so.
#' Terms given in `random` are never removed.
#'
#' **Selection.** Candidates that fail to converge or have a non-positive
#' definite Hessian are discarded. With `residual_checks = TRUE`, a candidate
#' passes when the three DHARMa tests give p >= 0.05; if at least one
#' candidate passes, only those compete, and otherwise those passing the
#' most tests (`checks_passed = FALSE`: inspect that model before using
#' it). The candidate with the lowest AIC is selected. Candidates whose AIC
#' differs by less than 1e-4 fit equally well: in a cell-means model without
#' random effects the link does not change the fitted means. The log link is
#' then preferred, because its means stay positive as the multiplicative and
#' dominance null models need, unless some combination has only zeros, in
#' which case the identity link is preferred; the `reason` says which.
#'
#' **REML refit.** The selected model is refitted by restricted maximum
#' likelihood (REML), whose variance estimates account for the number of
#' fixed effects. [estimate_interactions()] uses it for its default
#' small-sample correction; the selection itself uses the maximum-likelihood
#' AIC.
#'
#' @section Assumptions:
#' * A random intercept per unit implies the same correlation between any
#'   two samplings of a unit. Temporal autocorrelation is only modelled if a
#'   structure is given in `random`, e.g. `random = "ar1(.time + 0 | .unit)"`
#'   (first-order autoregressive; it treats the sampling times as equally
#'   spaced).
#' * With `dispersion = "common"`, all treatment combinations and `by`
#'   groups share one dispersion parameter on the scale of the family; with
#'   very different groups (species, populations), use `"by"` or `"cell"`.
#' * The selection among candidates is not accounted for in the uncertainty
#'   of later estimates.
#'
#' @section Missing values and zeros:
#' Rows with a missing response are left out for that response only and
#' counted in `n_missing`. Zeros are kept as observed values and no
#' transformation is applied (see "Combinations with only zeros").
#'
#' @references
#' Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C.
#' W., Nielsen, A., Skaug, H. J., Mächler, M., & Bolker, B. M. (2017).
#' glmmTMB balances speed and flexibility among packages for zero-inflated
#' generalized linear mixed modeling. *The R Journal*, 9(2), 378–400.
#' \doi{10.32614/RJ-2017-066}
#'
#' Madge Pimentel, I., et al. (2025). Hypothesis-driven research on multiple
#' stressors: An analytical framework for stressor interactions. *Ecology
#' and Evolution*, 15, e71959. \doi{10.1002/ece3.71959}
#'
#' @examples
#' jars <- expand.grid(
#'   salinity = c("ambient", "high"),
#'   warming = c("no", "yes"),
#'   replicate = 1:5,
#'   stringsAsFactors = FALSE
#' )
#' jars$jar <- seq_len(nrow(jars))
#' # Algal counts: each stressor halves the expected count
#' mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
#'   ifelse(jars$warming == "yes", 0.5, 1)
#' jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
#' # Survivors out of 20 animals per jar
#' jars$n_start <- 20
#' jars$survivors <- c(18, 15, 16, 9)[match(paste(jars$salinity, jars$warming),
#'   c("ambient no", "high no", "ambient yes", "high yes"))] +
#'   c(-1, 0, 1, 0, 0)[jars$replicate]
#' design <- describe_design(jars, c("salinity", "warming"),
#'   control = c(salinity = "ambient", warming = "no"), unit = "jar")
#' fits <- fit_responses(jars, c("algae", "survivors"), design,
#'   trials = c(survivors = "n_start"), residual_checks = FALSE, seed = 1,
#'   quiet = TRUE)
#' fits
#' fits$fits$algae$selection
#' @export
fit_responses <- function(data, responses, design, families = NULL,
                          response_types = NULL, trials = NULL, offset = NULL,
                          time_structure = "interaction", dispersion = "common",
                          random = NULL, control = NULL, residual_checks = TRUE,
                          n_sim_residuals = 500L, seed = NULL, quiet = FALSE) {
  check_data_frame(data, allow_empty = FALSE)
  data <- as_plain_data_frame(data)
  check_design_object(design)
  check_choice(time_structure, c("interaction", "additive"))
  check_choice(dispersion, c("common", "by", "cell"))
  if (identical(dispersion, "by") && length(design$spec$by) == 0L) {
    abort_input(paste(
      "`dispersion = \"by\"` needs `by` groups in the design; use \"common\" or",
      "\"cell\"."
    ))
  }
  check_flag(residual_checks)
  check_flag(quiet)
  check_number(n_sim_residuals, lower = 50)
  check_seed(seed)
  check_control(control)
  responses <- check_responses(data, responses, design)
  trials <- check_response_map(trials, responses, data, "trials")
  offset <- check_response_map(offset, responses, data, "offset")
  both <- intersect(names(trials), names(offset))
  if (length(both) > 0L) {
    abort_input(sprintf(
      "%s cannot have both `trials` and `offset`: a binomial response has no sampling effort.",
      format_items(quote_items(both), max = Inf)
    ))
  }
  for (r in names(trials)) {
    check_binomial_response(data, r, trials[[r]])
  }
  for (r in names(offset)) {
    check_offset_response(data, r, offset[[r]])
  }
  families <- check_families(families, responses)
  check_family_kinds(families, trials, offset)
  response_types <- check_response_types(response_types, responses)
  problems <- design$checks[design$checks$severity == "problem", , drop = FALSE]
  if (nrow(problems) > 0L) {
    abort_input(
      paste(
        "The design has problems that make the model invalid; fix them first:",
        paste0("\n* ", problems$message, collapse = "")
      )
    )
  }
  frame <- design_model_frame(data, design)
  if (residual_checks) {
    check_installed(
      "DHARMa",
      "run the residual checks (or set `residual_checks = FALSE` to skip them)"
    )
  }
  random_terms <- resolve_random(random, design, frame)
  seed <- resolve_seed(seed)
  run <- function() {
    lapply(responses, function(r) {
      if (!quiet) message("Fitting response: ", r)
      fit_one_response(
        frame = frame,
        response = r,
        families = families[[r]],
        response_type = response_types[[r]],
        trials_col = if (r %in% names(trials)) trials[[r]] else NULL,
        offset_col = if (r %in% names(offset)) offset[[r]] else NULL,
        time_structure = time_structure,
        dispersion = dispersion,
        random = random_terms,
        control = control,
        design = design,
        residual_checks = residual_checks,
        n_sim_residuals = as.integer(n_sim_residuals)
      )
    })
  }
  fits <- withr::with_seed(seed, run())
  names(fits) <- responses
  structure(
    list(
      fits = fits,
      design = design,
      settings = list(
        time_structure = time_structure,
        dispersion = dispersion,
        random = random_text(random_terms$terms),
        random_requested = if (is.null(random)) "automatic" else random,
        trials = trials,
        offset = offset,
        control = if (is.null(control)) "scaled iteration limits" else "user",
        residual_checks = residual_checks,
        n_sim_residuals = as.integer(n_sim_residuals),
        seed = seed
      )
    ),
    class = "must_fit"
  )
}

check_design_object <- function(design, call = sys.call(-1)) {
  if (!inherits(design, "must_design")) {
    abort_input(
      sprintf("`design` must be the output of describe_design(), not %s.", describe_value(design)),
      call = call
    )
  }
  invisible(design)
}

check_seed <- function(seed, call = sys.call(-1)) {
  if (!is.null(seed)) {
    check_number(seed, lower = -.Machine$integer.max, upper = .Machine$integer.max,
                 arg = "seed", call = call)
    if (seed != round(seed)) {
      abort_input("`seed` must be a whole number.", call = call)
    }
  }
  invisible(seed)
}

# The seed actually used: the one given, or one drawn from the session's
# stream (which advances it by one draw, as any random function does).
resolve_seed <- function(seed) {
  if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as.integer(seed)
}

check_control <- function(control, call = sys.call(-1)) {
  if (!is.null(control) && (!is.list(control) || is.null(control$optCtrl))) {
    abort_input(
      sprintf(
        "`control` must be NULL or the output of glmmTMB::glmmTMBControl(), not %s.",
        describe_value(control)
      ),
      call = call
    )
  }
  invisible(control)
}

# Optimiser limits for a model with `n_fixed` fixed effects. glmmTMB's defaults
# (300 iterations, 400 evaluations) are too low for cell-means models with
# many treatment combinations and times.
optimizer_control <- function(n_fixed) {
  glmmTMB::glmmTMBControl(optCtrl = list(
    iter.max = max(300L, 20L * as.integer(n_fixed)),
    eval.max = max(400L, 30L * as.integer(n_fixed))
  ))
}

check_responses <- function(data, responses, design, call = sys.call(-1)) {
  if (!is.character(responses) || length(responses) == 0L || anyNA(responses) ||
      anyDuplicated(responses) > 0L) {
    abort_input("`responses` must be a character vector of distinct column names.", call = call)
  }
  check_columns(data, stats::setNames(responses, rep("responses", length(responses))),
                call = call)
  spec <- design$spec
  clash <- intersect(responses, c(spec$stressors, spec$unit, spec$time, spec$by, spec$groups))
  if (length(clash) > 0L) {
    abort_input(
      sprintf(
        "%s cannot be a response: it is a design column.",
        format_items(quote_items(clash), max = Inf)
      ),
      call = call
    )
  }
  for (r in responses) {
    check_numeric_column(data, r, allow_na = TRUE, call = call)
  }
  responses
}

# `trials` and `offset`: character vectors named after responses, giving a
# numeric column with positive values (whole numbers for trials).
check_response_map <- function(x, responses, data, arg, call = sys.call(-1)) {
  if (is.null(x)) {
    return(stats::setNames(character(), character()))
  }
  example <- if (identical(arg, "trials")) "c(survivors = \"n_start\")" else "c(daphnia = \"volume_l\")"
  if (!is.character(x) || length(x) == 0L || is.null(names(x)) || any(!nzchar(names(x))) ||
      anyNA(x) || anyDuplicated(names(x)) > 0L) {
    abort_input(
      sprintf("`%s` must be a character vector named after the responses, e.g. %s.", arg, example),
      call = call
    )
  }
  unknown <- setdiff(names(x), responses)
  if (length(unknown) > 0L) {
    abort_input(
      sprintf(
        "`%s` names %s, which %s not in `responses`.", arg,
        format_items(quote_items(unknown), max = Inf), if (length(unknown) == 1L) "is" else "are"
      ),
      call = call
    )
  }
  check_columns(data, stats::setNames(unname(x), rep(arg, length(x))), call = call)
  for (col in unique(unname(x))) {
    check_numeric_column(data, col, lower = 0, include_lower = FALSE, allow_na = TRUE,
                         whole_numbers = identical(arg, "trials"), call = call)
  }
  x
}

check_binomial_response <- function(data, response, trials_col, call = sys.call(-1)) {
  y <- data[[response]]
  n <- data[[trials_col]]
  has_y <- !is.na(y)
  bad_n <- which(has_y & is.na(n))
  if (length(bad_n) > 0L) {
    abort_input(
      sprintf("'%s' has no number of trials in %s.", response, format_rows(bad_n)),
      column = trials_col, rows = bad_n, call = call
    )
  }
  bad <- which(has_y & (y < 0 | abs(y - round(y)) > sqrt(.Machine$double.eps) | y > n))
  if (length(bad) > 0L) {
    abort_input(
      sprintf(
        paste(
          "With `trials`, '%s' must hold the number of successes, a whole number",
          "between 0 and '%s' (%s). For proportions without a number of trials",
          "use families = list(beta = glmmTMB::beta_family())."
        ),
        response, trials_col, format_rows(bad)
      ),
      column = response, rows = bad, call = call
    )
  }
  invisible(TRUE)
}

check_offset_response <- function(data, response, effort_col, call = sys.call(-1)) {
  bad <- which(!is.na(data[[response]]) & is.na(data[[effort_col]]))
  if (length(bad) > 0L) {
    abort_input(
      sprintf(
        "'%s' has no sampling effort ('%s') in %s: give it or remove those values.",
        response, effort_col, format_rows(bad)
      ),
      column = effort_col, rows = bad, call = call
    )
  }
  invisible(TRUE)
}

is_family <- function(x) {
  is.list(x) && !is.null(x$family) && is.function(x$linkinv)
}

binomial_family_names <- function() {
  c("binomial", "betabinomial")
}

# Families whose variance has a dispersion parameter estimated from the data.
family_has_dispersion <- function(family_name) {
  !family_name %in% c("poisson", "binomial", "truncated_poisson")
}

check_families <- function(families, responses, call = sys.call(-1)) {
  out <- stats::setNames(vector("list", length(responses)), responses)
  if (is.null(families)) {
    return(out)
  }
  valid_set <- function(fams, label) {
    if (!is.list(fams) || length(fams) == 0L || is.null(names(fams)) ||
        any(!nzchar(names(fams))) || anyDuplicated(names(fams)) > 0L ||
        !all(vapply(fams, is_family, logical(1)))) {
      abort_input(
        sprintf(
          paste(
            "%s must be a named list of family objects, e.g.",
            "list(gaussian = gaussian(), nbinom2 = glmmTMB::nbinom2())."
          ),
          label
        ),
        call = call
      )
    }
    fams
  }
  if (is.list(families) && length(families) > 0L &&
      all(vapply(families, is_family, logical(1)))) {
    fams <- valid_set(families, "`families`")
    for (r in responses) out[[r]] <- fams
    return(out)
  }
  unknown <- setdiff(names(families), responses)
  if (is.null(names(families)) || length(unknown) > 0L) {
    abort_input(
      paste(
        "`families` must be a named list of family objects, or a list named",
        "after the responses with one such list each."
      ),
      call = call
    )
  }
  for (r in names(families)) {
    out[[r]] <- valid_set(families[[r]], sprintf("`families$%s`", r))
  }
  out
}

# Binomial families need trials, trials need binomial families, and an offset
# needs log links.
check_family_kinds <- function(families, trials, offset, call = sys.call(-1)) {
  for (r in names(families)) {
    fams <- families[[r]]
    if (is.null(fams)) {
      next
    }
    is_binomial <- vapply(fams, function(f) f$family %in% binomial_family_names(), logical(1))
    if (r %in% names(trials) && !all(is_binomial)) {
      abort_input(
        sprintf(
          "With `trials`, the candidates for '%s' must be binomial or beta-binomial families.", r
        ),
        call = call
      )
    }
    if (!r %in% names(trials) && any(is_binomial)) {
      abort_input(
        sprintf(
          paste(
            "Binomial families need the number of trials of each observation: give",
            "`trials = c(%s = \"<column>\")` and the number of successes as the response.",
            "For proportions without a number of trials (cover, frequencies) use",
            "glmmTMB::beta_family()."
          ),
          r
        ),
        call = call
      )
    }
    if (r %in% names(offset) &&
        !all(vapply(fams, function(f) identical(f$link, "log"), logical(1)))) {
      abort_input(
        sprintf("An offset acts on the log scale: the candidates for '%s' need a log link.", r),
        call = call
      )
    }
  }
  invisible(TRUE)
}

check_response_types <- function(response_types, responses, call = sys.call(-1)) {
  out <- stats::setNames(rep(NA_character_, length(responses)), responses)
  if (is.null(response_types)) {
    return(out)
  }
  if (!is.character(response_types) || is.null(names(response_types)) ||
      !all(names(response_types) %in% responses) ||
      !all(response_types %in% c("count", "continuous"))) {
    abort_input(
      paste(
        "`response_types` must be a character vector named after responses,",
        "with values \"count\" or \"continuous\"."
      ),
      call = call
    )
  }
  out[names(response_types)] <- response_types
  out
}

# Data with the internal columns used by the models. Rebuilds the treatment
# combinations and stops if they differ from those of the design.
design_model_frame <- function(data, design, call = sys.call(-1)) {
  spec <- design$spec
  roles <- c(spec$stressors, spec$unit, spec$time, spec$by, spec$groups)
  check_columns(data, roles, call = call)
  for (col in roles) {
    if (anyNA(data[[col]])) {
      abort_input(sprintf("Column '%s' has missing values.", col), call = call)
    }
  }
  keys <- design_keys(data, spec)
  observed <- unique(paste_key(keys$context, keys$cell))
  expected <- paste_key(design$cells$context, design$cells$cell)
  if (!setequal(observed, expected)) {
    abort_input(
      paste(
        "`data` does not have the treatment combinations of `design`. Call",
        "describe_design() on the same data."
      ),
      call = call
    )
  }
  group <- if (length(spec$by) == 0L) keys$cell else paste(keys$context, keys$cell, sep = " || ")
  group_levels <- unique(group[order(match(paste_key(keys$context, keys$cell), expected))])
  data$.cell <- keys$cell
  data$.context <- keys$context
  data$.group <- factor(group, levels = group_levels)
  if (!is.null(spec$time)) {
    time_levels <- if (is.numeric(data[[spec$time]])) {
      as.character(sort(unique(data[[spec$time]])))
    } else if (is.factor(data[[spec$time]])) {
      levels(droplevels(data[[spec$time]]))
    } else {
      sort(unique(as.character(data[[spec$time]])))
    }
    data$.time <- factor(as.character(data[[spec$time]]), levels = time_levels)
    cell_time <- paste(as.character(data$.group), as.character(data$.time), sep = " @ ")
    ct_levels <- unique(cell_time[order(as.integer(data$.group), as.integer(data$.time))])
    data$.cell_time <- factor(cell_time, levels = ct_levels)
  } else {
    data$.time <- NA_character_
    data$.cell_time <- data$.group
  }
  data$.unit <- factor(keys$unit_id)
  for (g in spec$groups) {
    value <- as.character(data[[g]])
    data[[group_column(g)]] <- factor(
      if (length(spec$by) == 0L) value else paste(keys$context, value, sep = " || ")
    )
  }
  data
}

# Internal column with the groups of a `groups` column.
group_column <- function(g) {
  if (length(g) == 0L) character() else paste0(".g_", g)
}

# Random part of the models: a named vector of terms (names are the grouping
# columns of the automatic terms) and whether they were chosen automatically.
resolve_random <- function(random, design, frame, call = sys.call(-1)) {
  if (is.null(random)) {
    terms <- character()
    groups <- design$groups
    if (!is.null(groups) && any(groups$random)) {
      use <- group_column(groups$group[groups$random])
      terms <- stats::setNames(sprintf("(1 | %s)", backtick(use)), use)
    }
    if (!is.null(design$spec$unit) && any(table(frame$.unit) > 1L)) {
      terms <- c(terms, .unit = "(1 | .unit)")
    }
    return(list(terms = terms, automatic = TRUE))
  }
  check_string(random, arg = "random", call = call)
  if (identical(random, "none")) {
    return(list(terms = character(), automatic = FALSE))
  }
  list(terms = c(custom = random), automatic = FALSE)
}

random_text <- function(terms) {
  if (length(terms) == 0L) "none" else paste(terms, collapse = " + ")
}

# Formula of a response: `spec` gives the response, its kind ("plain",
# "binomial" or "offset") and the fixed part.
model_formula_text <- function(spec, terms) {
  lhs <- if (identical(spec$kind, "binomial")) {
    sprintf("cbind(%s, .failures)", backtick(spec$response))
  } else {
    backtick(spec$response)
  }
  rhs <- spec$fixed
  if (identical(spec$kind, "offset")) {
    rhs <- paste(rhs, "+ offset(log(.effort))")
  }
  if (length(terms) > 0L) {
    rhs <- paste(rhs, "+", paste(terms, collapse = " + "))
  }
  paste(lhs, "~", rhs)
}

dispersion_formula <- function(dispersion) {
  switch(dispersion, common = "~1", by = "~ .context", cell = "~ .group")
}

# Automatic random terms whose standard deviation is negligible next to the
# fixed effects (variance estimated at the boundary, zero).
zero_variance_terms <- function(model, terms) {
  vc <- tryCatch(glmmTMB::VarCorr(model)$cond, error = function(e) NULL)
  beta <- tryCatch(glmmTMB::fixef(model)$cond, error = function(e) NULL)
  if (is.null(vc) || length(vc) == 0L || is.null(beta)) {
    return(character())
  }
  sds <- vapply(vc, function(v) max(attr(v, "stddev")), numeric(1))
  scale <- max(abs(beta), na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    return(character())
  }
  intersect(names(terms), names(sds)[is.finite(sds) & sds < 1e-3 * scale])
}

is_converged <- function(model) {
  !inherits(model, "error") && isTRUE(model$fit$convergence == 0) &&
    isTRUE(model$sdr$pdHess)
}

# Default candidate families from the values of a response.
default_families <- function(y, response_type = NA_character_, kind = "plain") {
  if (identical(kind, "binomial")) {
    return(list(
      type = "binomial",
      source = "trials",
      families = list(
        binomial = stats::binomial(link = "logit"),
        betabinomial = glmmTMB::betabinomial(link = "logit")
      )
    ))
  }
  y <- y[is.finite(y)]
  whole <- length(y) > 0L && all(y >= 0) && all(abs(y - round(y)) < sqrt(.Machine$double.eps))
  type_source <- if (is.na(response_type)) "detected" else "declared"
  type <- if (is.na(response_type)) {
    if (whole) "count" else "continuous"
  } else {
    response_type
  }
  fams <- if (identical(type, "count")) {
    list(
      poisson = stats::poisson(link = "log"),
      nbinom2 = glmmTMB::nbinom2(link = "log")
    )
  } else if (length(y) > 0L && all(y > 0)) {
    list(
      gaussian = stats::gaussian(),
      gaussian_log = stats::gaussian(link = "log"),
      gamma_log = stats::Gamma(link = "log"),
      lognormal = glmmTMB::lognormal(link = "log")
    )
  } else if (length(y) > 0L && all(y >= 0)) {
    list(
      gaussian = stats::gaussian(),
      gaussian_log = stats::gaussian(link = "log"),
      tweedie = glmmTMB::tweedie(link = "log")
    )
  } else {
    list(gaussian = stats::gaussian())
  }
  if (identical(kind, "offset")) {
    fams <- fams[vapply(fams, function(f) identical(f$link, "log"), logical(1))]
  }
  list(type = type, source = type_source, families = fams)
}

fixed_part <- function(has_time, time_structure) {
  if (!has_time) {
    return("0 + .group")
  }
  if (identical(time_structure, "interaction")) "0 + .cell_time" else "0 + .group + .time"
}

# Treatment combinations (at a time) whose values are all at a boundary of
# the support: all zero, or (binomial) all successes.
boundary_cells <- function(d, response, kind) {
  y <- d[[response]]
  zero <- tapply(y, d$.cell_time, function(v) all(v == 0))
  one <- if (identical(kind, "binomial")) {
    tapply(y >= d$.trials, d$.cell_time, all)
  } else {
    zero & FALSE
  }
  zero <- names(zero)[zero %in% TRUE]
  one <- names(one)[one %in% TRUE]
  c(stats::setNames(rep("zero", length(zero)), zero), stats::setNames(rep("one", length(one)), one))
}

# Combinations whose means are at the boundary of the parameter space of a
# model: only zeros on a log link; only zeros or only successes on a logit
# link. With additive time ("0 + .group + .time") a combination is at the
# boundary only if it has such values at every time.
effective_boundary <- function(boundary, d, fixed, link) {
  if (length(boundary) == 0L || !link %in% c("log", "logit")) {
    return(character())
  }
  if (identical(link, "log")) {
    boundary <- boundary[boundary == "zero"]
  }
  if (identical(fixed, "0 + .group + .time")) {
    ct <- as.character(d$.cell_time)
    g <- as.character(d$.group)
    keep <- vapply(names(boundary), function(lv) {
      times_of_group <- unique(ct[g == g[match(lv, ct)]])
      all(times_of_group %in% names(boundary)[boundary == boundary[[lv]]])
    }, logical(1))
    boundary <- boundary[keep]
  } else if (!fixed %in% c("0 + .cell_time", "0 + .group")) {
    return(character())
  }
  boundary
}

# Start values and map for the fixed effects of a candidate: the means of
# the combinations at a boundary are fixed there (log or logit link).
boundary_map <- function(d, spec, family, boundary, X) {
  link <- family$link
  boundary <- effective_boundary(boundary, d, spec$fixed, link)
  if (length(boundary) == 0L) {
    return(NULL)
  }
  param <- if (identical(spec$fixed, "0 + .cell_time")) ".cell_time" else ".group"
  level <- if (identical(param, ".group")) {
    as.character(d$.group[match(names(boundary), as.character(d$.cell_time))])
  } else {
    names(boundary)
  }
  value <- ifelse(boundary == "zero", -20, 20)
  first <- !duplicated(level)
  # The first columns of the model matrix are the levels of `param`, in order.
  # Columns are found by position: their names can differ from the levels in
  # non-UTF-8 locales (escaped characters).
  lv <- levels(d[[param]])
  if (ncol(X) < length(lv) || !all(startsWith(colnames(X)[seq_along(lv)], param))) {
    abort_internal("unexpected columns in the model matrix of the cell means.")
  }
  hit <- match(level[first], lv)
  if (anyNA(hit)) {
    abort_internal("a combination at the boundary has no column in the model matrix.")
  }
  y <- d[[spec$response]]
  rows_by_level <- split(seq_len(nrow(d)), d[[param]])
  start <- rep(0, ncol(X))
  start[seq_along(lv)] <- vapply(rows_by_level, function(i) {
    if (identical(link, "logit")) {
      p <- sum(y[i]) / sum(d$.trials[i])
      stats::qlogis(min(max(p, 0.01), 0.99))
    } else {
      total <- if (identical(spec$kind, "offset")) sum(y[i]) / sum(d$.effort[i]) else mean(y[i])
      if (is.finite(total) && total > 0) log(total) else -20
    }
  }, numeric(1))
  start[hit] <- value[first]
  codes <- seq_along(start)
  codes[hit] <- NA
  list(start = list(beta = unname(start)), map = list(beta = factor(codes)),
       n_fixed = length(hit))
}

fit_one_response <- function(frame, response, families, response_type, trials_col,
                             offset_col, time_structure, dispersion, random, design,
                             residual_checks, n_sim_residuals, control = NULL) {
  if (is.null(random)) {
    random <- list(terms = character(), automatic = FALSE)
  }
  has_time <- !is.null(design$spec$time)
  y <- frame[[response]]
  keep <- !is.na(y)
  d <- frame[keep, , drop = FALSE]
  d$.group <- droplevels(d$.group)
  d$.cell_time <- droplevels(d$.cell_time)
  if (has_time) {
    d$.time <- droplevels(d$.time)
  }
  if (identical(time_structure, "additive") && has_time && length(design$spec$baseline) > 0L) {
    d <- baseline_as_control(d, design)
  }
  kind <- if (!is.null(trials_col)) "binomial" else if (!is.null(offset_col)) "offset" else "plain"
  if (identical(kind, "binomial")) {
    d$.trials <- d[[trials_col]]
    d$.failures <- d$.trials - d[[response]]
  }
  if (identical(kind, "offset")) {
    d$.effort <- d[[offset_col]]
  }
  base <- list(
    response = response, status = "failed", reason = "", model = NULL, model_reml = NULL,
    family = NA_character_, link = NA_character_, response_type = NA_character_,
    type_source = NA_character_, kind = kind, dispersion = dispersion,
    selection = NULL, checks_passed = NA, formula = NA_character_,
    n_obs = nrow(d), n_missing = sum(!keep), grid = NULL
  )
  if (nrow(d) == 0L) {
    base$reason <- "The response has no non-missing values."
    return(base)
  }
  if (is.null(families)) {
    defaults <- default_families(d[[response]], response_type, kind)
    families <- defaults$families
    base$response_type <- defaults$type
    base$type_source <- defaults$source
    if (length(families) == 0L) {
      base$reason <- "No log-link candidate suits a response with negative values and an offset."
      return(base)
    }
  } else {
    base$response_type <- if (identical(kind, "binomial")) "binomial" else "user"
    base$type_source <- "families"
  }
  few_cells <- character()
  if (identical(dispersion, "cell")) {
    units <- tapply(d$.unit, d$.group, function(u) length(unique(u)))
    few_cells <- names(units)[units < 2L]
    with_dispersion <- vapply(families, function(f) family_has_dispersion(f$family), logical(1))
    if (length(few_cells) > 0L && all(with_dispersion)) {
      base$reason <- cell_dispersion_message(few_cells)
      return(base)
    }
  }
  fixed <- fixed_part(has_time && nlevels(d$.time) > 1L, time_structure)
  spec <- list(response = response, kind = kind, fixed = fixed)
  base$formula <- model_formula_text(spec, random$terms)
  base$fixed <- fixed
  X <- stats::model.matrix(stats::as.formula(paste("~", fixed)), d)
  if (is.null(control)) {
    control <- optimizer_control(ncol(X))
  }
  boundary <- boundary_cells(d, response, kind)
  candidates <- lapply(names(families), function(nm) {
    fit_candidate(spec, random, d, families[[nm]], nm, dispersion, boundary, X,
                  residual_checks, n_sim_residuals, control, few_cells = few_cells)
  })
  selection <- do.call(rbind, lapply(candidates, `[[`, "row"))
  rownames(selection) <- NULL
  chosen <- select_candidate(selection, residual_checks, length(boundary) > 0L)
  if (is.na(chosen$best)) {
    base$selection <- selection
    base$reason <- paste0(
      "No candidate family converged with a positive-definite Hessian.",
      if (length(boundary) > 0L) {
        sprintf(
          " Some treatment combinations have only zeros (%s).",
          format_items(quote_items(describe_cell_times(d, names(boundary), design)), max = 3L)
        )
      } else {
        ""
      }
    )
    return(base)
  }
  best <- chosen$best
  selection$selected <- seq_len(nrow(selection)) == best
  winner <- candidates[[best]]
  grid <- unique(d[, c(".context", ".cell", ".time", ".group", ".cell_time"), drop = FALSE])
  grid <- grid[order(as.integer(grid$.cell_time)), , drop = FALSE]
  rownames(grid) <- NULL
  base$status <- "ok"
  base$model <- winner$model
  base$formula <- winner$formula
  base$family <- selection$family[[best]]
  base$link <- selection$link[[best]]
  family <- families[[best]]
  at_boundary <- effective_boundary(boundary, d, fixed, family$link)
  base$boundary <- at_boundary
  d$.row <- seq_len(nrow(d))
  re_columns <- c(".unit", grep("^\\.g_", names(d), value = TRUE), design$spec$groups,
                  design$spec$unit)
  base$boundary_rows <- d[as.character(d$.cell_time) %in% names(at_boundary),
                          intersect(c(".row", ".cell_time", ".group", ".context", re_columns,
                                      response, ".trials", ".failures", ".effort"), names(d)),
                          drop = FALSE]
  reml <- refit_reml(winner, spec, d, family, dispersion, control)
  base$model_reml <- reml$model
  base$reason <- paste(c(
    chosen$notes,
    if (length(winner$dropped) > 0L) {
      sprintf(
        paste(
          "The random intercept(s) for %s had variance estimated at zero and",
          "were removed from the selected model (same maximum likelihood)."
        ),
        format_items(quote_items(winner$dropped), max = Inf)
      )
    },
    if (length(at_boundary) > 0L) {
      sprintf(
        paste(
          "Some treatment combinations have only zeros%s (%s): their means are",
          "fixed at the boundary of the link function, and contrasts that use",
          "them are tested from the probability of that outcome."
        ),
        if (any(at_boundary == "one")) " or only successes" else "",
        format_items(quote_items(describe_cell_times(d, names(at_boundary), design)), max = 3L)
      )
    },
    reml$note
  ), collapse = " ")
  base$selection <- selection
  base$checks_passed <- chosen$checks_passed
  base$grid <- grid
  base$non_negative <- all(d[[response]] >= 0)
  base$df_info <- df_info(d, X, design, family, winner$terms, dispersion,
                          requested = random$terms)
  base
}

# With additive time the treatment means are shared by every sampling time.
# Before the treatments every unit is untreated, so the baseline rows are
# coded as the control combination of their `by` group: they inform the time
# effects and the differences among units, not the treatment means.
baseline_as_control <- function(d, design) {
  pre <- as.character(d$.time) %in% design$spec$baseline
  if (!any(pre)) {
    return(d)
  }
  cells <- design$cells
  ctrl <- cells[cells$n_stressed == 0L, , drop = FALSE]
  ctrl_group <- if (length(design$spec$by) == 0L) {
    ctrl$cell
  } else {
    paste(ctrl$context, ctrl$cell, sep = " || ")
  }
  hit <- match(d$.context[pre], ctrl$context)
  if (anyNA(hit)) {
    abort_input(paste(
      "With `baseline` and time_structure = \"additive\", every `by` group needs",
      "a control combination, to which the baseline samplings are assigned."
    ))
  }
  group <- as.character(d$.group)
  group[pre] <- ctrl_group[hit]
  d$.cell[pre] <- ctrl$cell[hit]
  d$.group <- droplevels(factor(group, levels = union(levels(d$.group), ctrl_group)))
  cell_time <- paste(as.character(d$.group), as.character(d$.time), sep = " @ ")
  ct_levels <- unique(cell_time[order(as.integer(d$.group), as.integer(d$.time))])
  d$.cell_time <- factor(cell_time, levels = ct_levels)
  d
}

# Readable labels of .cell_time levels: stressors away from control, `by`
# group and time.
describe_cell_times <- function(d, levels, design) {
  first <- match(levels, as.character(d$.cell_time))
  out <- short_cells(d$.cell[first], design$spec$control, d$.context[first])
  if (!is.null(design$spec$time)) {
    out <- paste0(out, " at ", as.character(d$.time[first]))
  }
  out
}

# Fit one model, keeping glmmTMB's warnings as messages.
fit_glmmtmb <- function(spec, terms, d, family, dispersion, control, start = NULL,
                        map = NULL, REML = FALSE) {
  messages <- character()
  disp <- if (family_has_dispersion(family$family)) dispersion_formula(dispersion) else "~1"
  model <- withCallingHandlers(
    tryCatch(
      glmmTMB::glmmTMB(
        stats::as.formula(model_formula_text(spec, terms)), data = d, family = family,
        dispformula = stats::as.formula(disp), start = start, map = map,
        control = control, REML = REML
      ),
      error = function(e) e
    ),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(model = model, messages = messages)
}

cell_dispersion_message <- function(few) {
  sprintf(
    paste(
      "dispersion = \"cell\" estimates one dispersion per treatment combination",
      "and needs two or more units with data in each; %s %s one."
    ),
    format_items(quote_items(few), max = 3L), if (length(few) == 1L) "has" else "have"
  )
}

fit_candidate <- function(spec, random, d, family, name, dispersion, boundary, X,
                          residual_checks, n_sim_residuals,
                          control = glmmTMB::glmmTMBControl(), few_cells = character()) {
  terms <- random$terms
  if (identical(dispersion, "cell") && family_has_dispersion(family$family) &&
      length(few_cells) > 0L) {
    row <- data.frame(
      family = name, link = family$link, random = random_text(terms),
      converged = FALSE, aic = NA_real_,
      p_uniformity = NA_real_, p_dispersion = NA_real_, p_zero_inflation = NA_real_,
      n_tests_passed = NA_integer_, residuals_ok = NA, selected = FALSE,
      message = cell_dispersion_message(few_cells), stringsAsFactors = FALSE
    )
    return(list(row = row, model = NULL, formula = model_formula_text(spec, terms),
                dropped = character(), terms = terms, start = NULL, map = NULL))
  }
  bmap <- boundary_map(d, spec, family, boundary, X)
  messages <- character()
  attempt <- function(terms, start = bmap$start) {
    res <- fit_glmmtmb(spec, terms, d, family, dispersion, control, start = start,
                       map = bmap$map)
    messages <<- c(messages, res$messages)
    res$model
  }
  model <- attempt(terms)
  if (inherits(model, "error") && is.null(bmap) && family$link %in% c("log", "logit")) {
    messages <- c(messages, conditionMessage(model))
    model <- attempt(terms, start = list(beta = unname(start_values(d, spec, family, X))))
  }
  # An automatic random intercept whose variance is estimated at zero adds
  # nothing to the fit and often makes the Hessian singular. Removing it gives
  # the same maximum likelihood; the reduced model is kept only if that is the
  # case. This is done whether or not the full model converged, so that the
  # result does not depend on how a version of glmmTMB handles the boundary.
  dropped <- character()
  if (!inherits(model, "error") && isTRUE(random$automatic)) {
    zero <- zero_variance_terms(model, terms)
    if (length(zero) > 0L) {
      kept <- terms[setdiff(names(terms), zero)]
      reduced <- attempt(kept)
      same_fit <- is_converged(reduced) && is.finite(model$fit$objective) &&
        reduced$fit$objective <= model$fit$objective + 0.05
      if (same_fit) {
        model <- reduced
        dropped <- zero
        terms <- kept
        messages <- c(messages, sprintf(
          "Random intercept(s) for %s removed: variance estimated at zero.",
          format_items(quote_items(zero), max = Inf)
        ))
      }
    }
  }
  row <- data.frame(
    family = name, link = family$link, random = random_text(terms),
    converged = FALSE, aic = NA_real_,
    p_uniformity = NA_real_, p_dispersion = NA_real_, p_zero_inflation = NA_real_,
    n_tests_passed = NA_integer_, residuals_ok = NA, selected = FALSE, message = "",
    stringsAsFactors = FALSE
  )
  formula <- model_formula_text(spec, terms)
  out <- list(row = row, model = NULL, formula = formula, dropped = dropped, terms = terms,
              start = bmap$start, map = bmap$map)
  if (inherits(model, "error")) {
    out$row$message <- conditionMessage(model)
    return(out)
  }
  converged <- is_converged(model)
  out$row$converged <- converged
  # Fixed boundary parameters are counted, so that the AIC stays comparable
  # with candidates that estimate them.
  n_fixed <- if (is.null(bmap)) 0L else bmap$n_fixed
  out$row$aic <- if (converged) stats::AIC(model) + 2 * n_fixed else NA_real_
  if (converged && residual_checks) {
    checks <- residual_tests(model, n_sim_residuals)
    out$row$p_uniformity <- checks$p_uniformity
    out$row$p_dispersion <- checks$p_dispersion
    out$row$p_zero_inflation <- checks$p_zero_inflation
    p <- c(checks$p_uniformity, checks$p_dispersion, checks$p_zero_inflation)
    out$row$n_tests_passed <- sum(p >= 0.05, na.rm = TRUE)
    out$row$residuals_ok <- checks$ok
    messages <- c(messages, checks$message)
  }
  out$row$message <- paste(unique(messages[nzchar(messages)]), collapse = " | ")
  out$model <- model
  out
}

# Start values on the link scale from the mean of the response.
start_values <- function(d, spec, family, X) {
  y <- d[[spec$response]]
  value <- if (identical(family$link, "logit")) {
    stats::qlogis(min(max(sum(y) / sum(d$.trials), 0.01), 0.99))
  } else {
    mu0 <- if (identical(spec$kind, "offset")) sum(y) / sum(d$.effort) else mean(y[y > 0])
    if (is.finite(mu0) && mu0 > 0) log(mu0) else 0
  }
  ifelse(grepl("^\\.time", colnames(X)), 0, value)
}

# Choose among the candidates: converged ones, those that pass the residual
# checks (or pass the most), lowest AIC; ties between links broken on purpose.
select_candidate <- function(selection, residual_checks, has_boundary) {
  usable <- selection$converged & is.finite(selection$aic)
  if (!any(usable)) {
    return(list(best = NA_integer_, checks_passed = NA, notes = character()))
  }
  pool <- usable
  checks_passed <- NA
  notes <- character()
  if (residual_checks) {
    passed <- usable & selection$residuals_ok %in% TRUE
    if (any(passed)) {
      pool <- passed
      checks_passed <- TRUE
    } else {
      n <- selection$n_tests_passed
      n[is.na(n)] <- 0L
      most <- max(n[usable])
      pool <- usable & n == most
      if (!any(pool)) {
        pool <- usable
      }
      checks_passed <- FALSE
      notes <- sprintf(
        paste(
          "No candidate passed all the residual checks; the lowest AIC among those",
          "passing the most (%d of 3) was selected."
        ),
        most
      )
    }
  }
  aic <- ifelse(pool, selection$aic, Inf)
  tied <- which(pool & abs(aic - min(aic)) < 1e-4)
  best <- tied[[1L]]
  if (length(tied) > 1L) {
    prefer <- if (has_boundary) "identity" else "log"
    hit <- tied[selection$link[tied] == prefer]
    if (length(hit) > 0L) {
      best <- hit[[1L]]
    }
    notes <- c(notes, sprintf(
      paste(
        "%s fit equally well (same AIC: in a cell-means model the link does not",
        "change the fitted means); %s was chosen because %s."
      ),
      format_items(quote_items(selection$family[tied]), max = Inf),
      quote_items(selection$family[[best]]),
      if (identical(selection$link[[best]], "log")) {
        "its means stay positive, as the multiplicative and dominance null models need"
      } else if (has_boundary) {
        "some combinations have only zeros, which a log link puts on its boundary"
      } else {
        "it comes first in the list of candidates"
      }
    ))
  }
  list(best = best, checks_passed = checks_passed, notes = notes)
}

# REML refit of the selected candidate (same random terms, same boundary map).
refit_reml <- function(winner, spec, d, family, dispersion, control) {
  has_variance <- family_has_dispersion(family$family) || length(winner$terms) > 0L
  if (!has_variance) {
    return(list(model = winner$model, note = character()))
  }
  res <- fit_glmmtmb(spec, winner$terms, d, family, dispersion, control,
                     start = winner$start, map = winner$map, REML = TRUE)
  if (is_converged(res$model)) {
    return(list(model = res$model, note = character()))
  }
  list(
    model = NULL,
    note = paste(
      "The REML refit did not converge: the small-sample correction uses the",
      "maximum-likelihood covariance."
    )
  )
}

# What the degrees of freedom of the small-sample correction depend on.
# Blocks and groups count when their random intercept was requested
# (`requested`), even if it was then removed because its variance was
# estimated at zero: the degrees of freedom follow the design, as in an
# analysis of variance.
df_info <- function(d, X, design, family, terms, dispersion, requested = terms) {
  repeated <- any(table(droplevels(d$.unit)) > 1L)
  # Units sampled at a single time each (destructive sampling) are nested in
  # treatment combination x time; otherwise in treatment combination.
  times_per_unit <- tapply(as.character(d$.time), droplevels(d$.unit),
                           function(t) length(unique(t)))
  unit_strata <- if (all(times_per_unit == 1L)) d$.cell_time else d$.group
  groups <- design$groups
  in_model <- if (is.null(groups)) {
    logical()
  } else {
    vapply(groups$group, terms_use_group, logical(1), terms = requested)
  }
  block_cols <- if (is.null(groups)) {
    character()
  } else {
    group_column(groups$group[groups$structure == "blocks" & in_model])
  }
  rows_by_context <- split(seq_len(nrow(d)), d$.context)
  context_names <- names(rows_by_context)
  # Indexed by position: the context label is "" without `by`, and [[""]]
  # matches nothing.
  contexts <- do.call(rbind, lapply(seq_along(rows_by_context), function(k) {
    ctx <- context_names[[k]]
    i <- rows_by_context[[k]]
    Xi <- X[i, colSums(abs(X[i, , drop = FALSE])) > 0, drop = FALSE]
    data.frame(
      context = ctx,
      n_obs = length(i),
      rank = qr(Xi)$rank,
      n_units = length(unique(d$.unit[i])),
      n_unit_groups = length(unique(unit_strata[i])),
      block_df = sum(vapply(block_cols, function(b) length(unique(d[[b]][i])) - 1L, numeric(1))),
      stringsAsFactors = FALSE
    )
  }))
  whole <- list()
  if (!is.null(groups)) {
    # Groups to which some stressors are applied whole: whole plots, and
    # partly confounded groups whose random intercept is in the model.
    use <- which(groups$structure == "whole_units" |
                   (groups$structure == "confounded" & in_model))
    for (k in use) {
      g <- group_column(groups$group[[k]])
      between <- strsplit(groups$between[[k]], ", ", fixed = TRUE)[[1L]]
      combo <- make_cell_label(lapply(stats::setNames(between, between), function(s) {
        as.character(d[[s]])
      }), between)
      wc <- do.call(rbind, lapply(seq_along(rows_by_context), function(j) {
        ctx <- context_names[[j]]
        i <- rows_by_context[[j]]
        data.frame(
          context = ctx,
          n_groups = length(unique(d[[g]][i])),
          n_combos = length(unique(combo[i])),
          block_df = sum(vapply(block_cols, function(b) length(unique(d[[b]][i])) - 1L,
                                numeric(1))),
          stringsAsFactors = FALSE
        )
      }))
      whole[[length(whole) + 1L]] <- list(group = groups$group[[k]], between = between,
                                          contexts = wc)
    }
  }
  list(
    has_dispersion = family_has_dispersion(family$family),
    has_random = length(terms) > 0L,
    repeated = repeated,
    dispersion = dispersion,
    rank_total = qr(X)$rank,
    contexts = contexts,
    whole = whole,
    cell_units = tapply(d$.unit, d$.cell_time, function(u) length(unique(u)))
  )
}

# Whether the random terms include an intercept for the groups of `g`:
# automatic terms are named after their column; terms given in `random` are
# searched for "| g)" or "| .g_g)".
terms_use_group <- function(g, terms) {
  if (length(terms) == 0L) {
    return(FALSE)
  }
  if (group_column(g) %in% names(terms)) {
    return(TRUE)
  }
  plain <- function(x) gsub("[[:space:]`]", "", x)
  text <- plain(terms)
  any(grepl(paste0("|", plain(g), ")"), text, fixed = TRUE) |
        grepl(paste0("|", plain(group_column(g)), ")"), text, fixed = TRUE))
}

# DHARMa tests for one fitted model. A test that cannot be computed gives NA
# (and the model is then not counted as passing).
residual_tests <- function(model, n_sim) {
  out <- list(p_uniformity = NA_real_, p_dispersion = NA_real_,
              p_zero_inflation = NA_real_, n_passed = 0L, ok = NA, message = "")
  sim <- tryCatch(
    DHARMa::simulateResiduals(fittedModel = model, n = n_sim, plot = FALSE),
    error = function(e) e
  )
  if (inherits(sim, "error")) {
    out$message <- paste("Residual simulation failed:", conditionMessage(sim))
    return(out)
  }
  get_p <- function(test) {
    p <- tryCatch(as.numeric(test(sim, plot = FALSE)$p.value)[[1L]], error = function(e) NA_real_)
    if (length(p) == 0L) NA_real_ else p
  }
  out$p_uniformity <- get_p(DHARMa::testUniformity)
  out$p_dispersion <- get_p(DHARMa::testDispersion)
  out$p_zero_inflation <- get_p(DHARMa::testZeroInflation)
  p <- c(out$p_uniformity, out$p_dispersion, out$p_zero_inflation)
  out$n_passed <- sum(p >= 0.05, na.rm = TRUE)
  out$ok <- if (anyNA(p)) NA else all(p >= 0.05)
  out
}

#' @export
print.must_fit <- function(x, ...) {
  cat(sprintf(
    "<must_fit> %d response(s); time structure: %s; dispersion: %s; random effects: %s\n",
    length(x$fits), x$settings$time_structure, x$settings$dispersion, x$settings$random
  ))
  cat(sprintf(
    "Residual checks: %s; seed: %s\n",
    if (x$settings$residual_checks) {
      sprintf("DHARMa (%d simulations)", x$settings$n_sim_residuals)
    } else {
      "not run"
    },
    format(x$settings$seed)
  ))
  for (f in x$fits) {
    if (!identical(f$status, "ok")) {
      cat(sprintf("  %s: FAILED. %s\n", f$response, f$reason))
      next
    }
    checks <- if (is.na(f$checks_passed)) {
      "not run"
    } else if (f$checks_passed) {
      "passed"
    } else {
      "NONE PASSED"
    }
    aic <- f$selection$aic[f$selection$selected]
    used <- f$selection$random[f$selection$selected]
    cat(sprintf(
      "  %s: %s (%s link; %s, %s), AIC %.1f, n = %d, missing = %d, residual checks: %s%s\n",
      f$response, f$family, f$link, f$response_type, f$type_source, aic, f$n_obs,
      f$n_missing, checks,
      if (identical(used, x$settings$random)) "" else paste0("; random effects: ", used)
    ))
  }
  invisible(x)
}
