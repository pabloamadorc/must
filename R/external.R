# Models fitted outside must ---------------------------------------------------
#
# The means of the treatment combinations come from the reference grid of
# emmeans (Suggests): their linear predictor, the matrix that maps the
# coefficients of the model to it, the covariance of the coefficients and the
# degrees-of-freedom function. From there the contrasts, draws and tests are
# those of the selected models of fit_responses() (estimate_sets()).

#' @rdname estimate_interactions
#' @export
estimate_interactions.default <- function(fit, design, null_models = "AD",
                                          inference = "simulation",
                                          correction = "small_sample", n_sim = 10000L,
                                          level = 0.95, p_adjust = "BH", overall = TRUE,
                                          sensitivity = FALSE, seed = NULL, sam_shape = 3.2,
                                          margin = NULL, df = NULL, data = NULL, ...) {
  call <- generic_call("estimate_interactions")
  if (missing(design) || !inherits(design, "must_design")) {
    abort_input(sprintf(
      paste(
        "`fit` is %s, not the output of fit_responses(). To analyse a model fitted",
        "outside must, give the output of describe_design() for its data as `design`."
      ),
      describe_value(fit)
    ), call = call)
  }
  opts <- check_estimate_options(null_models, inference, correction, n_sim, level, p_adjust,
                                 overall, seed, sam_shape, margin, call = call)
  check_flag(sensitivity, call = call)
  check_choice_sensitivity(opts$null_model_choice, sensitivity, call = call)
  if (!is.null(df)) {
    check_number(df, lower = 0, include_lower = FALSE, call = call)
    if (identical(correction, "none")) {
      abort_input(paste(
        "`df` is only used with correction = \"small_sample\"; with \"none\" the",
        "normal distribution is used."
      ), call = call)
    }
  }
  if (!is.null(data)) {
    check_data_frame(data, call = call)
  }
  check_installed("emmeans", "estimate interactions from a model fitted outside must",
                  call = call)
  ext <- external_predictions(fit, design, data, opts$n_sim, sensitivity, call, ...)
  seed_given <- !is.null(seed)
  opts$seed <- resolve_seed(seed)
  set_df <- if (identical(correction, "none")) {
    function(s) Inf
  } else if (!is.null(df)) {
    function(s) df
  } else {
    function(s) NA_real_
  }
  parts <- withr::with_seed(opts$seed, list(
    estimate_sets(ext$f, ext$pred, design, opts, sensitivity, baseline = ext$baseline,
                  set_df = set_df)
  ))
  response <- ext$f$response
  new_interactions(
    parts, design, opts, if (sensitivity) response else character(),
    settings = list(
      prediction = paste(
        "emmeans reference grid (random effects at zero; other predictors averaged",
        "over or at their means)"
      ),
      time_structure = NA_character_,
      families = stats::setNames(ext$family, response),
      model = class(fit)[[1L]],
      df = if (!is.null(df)) {
        df
      } else if (!is.null(ext$disp_df) && !identical(correction, "none")) {
        "residual"
      } else {
        "emmeans"
      },
      averaged_over = ext$averaged_over,
      versions = analysis_versions("emmeans"),
      seed_given = seed_given
    ),
    diagnostics = list(
      failed_responses = data.frame(response = character(), reason = character(),
                                    stringsAsFactors = FALSE),
      reml = stats::setNames(NA, response),
      means = ext$means
    )
  )
}

# Predictions of a model fitted outside must at the treatment combinations
# (and times) of `design`, in the form that estimate_sets() uses.
external_predictions <- function(fit, design, data, n_sim, sensitive, call, ...) {
  spec <- design$spec
  vars <- c(spec$stressors, spec$time, spec$by)
  emm_info <- external_grid(fit, design, data, call, ...)
  emm <- emm_info$emm
  g <- emm@grid
  labels <- data.frame(row.names = seq_len(nrow(g)))
  for (s in spec$stressors) {
    labels[[s]] <- grid_labels(g[[s]], spec$levels[[s]])
  }
  for (b in spec$by) {
    labels[[b]] <- grid_labels(g[[b]], unique(design$cells[[b]]))
  }
  keep <- stats::complete.cases(labels)
  time_label <- if (is.null(spec$time)) {
    rep(NA_character_, nrow(g))
  } else {
    as.character(g[[spec$time]])
  }
  # Baseline times of each stressor as they are labelled in the grid.
  baseline <- baseline_times(spec)
  if (!is.null(baseline) && is.numeric(g[[spec$time]])) {
    times <- unique(g[[spec$time]])
    baseline <- lapply(baseline, function(b) {
      as.character(times[match_numeric(suppressWarnings(as.numeric(b)), times)])
    })
  }
  cell <- make_cell_label(labels, spec$stressors)
  context <- make_context_label(labels, spec$by)
  group <- if (length(spec$by) == 0L) cell else paste(context, cell, sep = " || ")
  grid <- data.frame(
    .context = context, .cell = cell, .time = time_label,
    .cell_time = if (is.null(spec$time)) group else paste(group, time_label, sep = " @ "),
    stringsAsFactors = FALSE
  )[keep, , drop = FALSE]
  rownames(grid) <- NULL
  if (anyDuplicated(paste_key(grid$.context, grid$.cell, grid$.time)) > 0L) {
    nesting <- emm@model.info$nesting
    abort_input(sprintf(
      paste(
        "The reference grid of the model has several means for the same treatment",
        "combination%s: check that no stressor is nested in another factor of the",
        "model (for example, applied to whole blocks)."
      ),
      if (length(nesting) > 0L) {
        sprintf(" (emmeans found %s)", paste(
          names(nesting), "%in%", vapply(nesting, paste, character(1), collapse = "*"),
          collapse = ", "
        ))
      } else {
        ""
      }
    ), call = call)
  }
  active <- !is.na(emm@bhat)
  linfct <- emm@linfct[keep, , drop = FALSE]
  X <- linfct[, active, drop = FALSE]
  V <- as.matrix(emm@V)
  if (!identical(dim(V), c(sum(active), sum(active)))) {
    abort_internal("the covariance of emmeans does not match its coefficients.")
  }
  eta <- unname(stats::predict(emm, type = "link"))[keep]
  link <- external_link(emm@misc, call)
  family <- tryCatch(stats::family(fit), error = function(e) NULL)
  proportion <- link$name %in% c("logit", "probit", "cloglog", "cauchit") ||
    (inherits(family, "family") &&
       family$family %in% c("binomial", "quasibinomial", "betabinomial", "beta"))
  pred <- prediction_engine(grid, X, V, eta = eta, link = link, n_sim = n_sim,
                            complement = sensitive && proportion)
  pred$model <- fit
  pred$stable <- stabilizing_scale(
    family,
    theta = if (is.numeric(fit$theta) && length(fit$theta) == 1L) fit$theta else nbinom2_theta(fit),
    transformed = inherits(family, "family") && !identical(link$name, family$link)
  )
  se <- sqrt(pmax(diag(pred$V_eta), 0))
  se[!is.finite(diag(pred$V_eta))] <- NA_real_
  pred$cell_problem <- external_problems(eta, se, link$name)
  cell_df <- summary(emm)$df[keep]
  # emmeans gives a glm whose dispersion is estimated (quasi-Poisson, say)
  # infinite degrees of freedom; summary.glm() uses its residual ones, and so
  # do the contrasts.
  disp_df <- glm_dispersion_df(fit)
  dffun <- emm@dffun
  dfargs <- emm@dfargs
  # Degrees of freedom of the linear combination of coefficients given by the
  # gradient `g` of a contrast, as emmeans computes them for its contrasts;
  # the smallest of its cells if emmeans gives none.
  pred$contrast_df <- function(g, cols) {
    L <- drop(crossprod(g, linfct[cols, , drop = FALSE]))
    d <- tryCatch(suppressMessages(dffun(L, dfargs)), error = function(e) NA_real_)
    d <- suppressWarnings(as.numeric(d)[1L])
    if (is.na(d)) {
      d <- suppressWarnings(min(cell_df[cols], na.rm = TRUE))
    }
    if (!is.finite(d) && !is.null(disp_df)) {
      d <- disp_df
    }
    if (is.finite(d) || identical(d, Inf)) d else Inf
  }
  response <- external_response(fit)
  f <- list(
    response = response,
    status = "ok",
    kind = if (proportion) "binomial" else "model",
    non_negative = external_non_negative(fit, link),
    grid = grid
  )
  means <- data.frame(
    context = grid$.context, cell = grid$.cell, time = grid$.time,
    link_estimate = eta, link_se = se, mean = link$linkinv(eta), df = cell_df,
    problem = pred$cell_problem, stringsAsFactors = FALSE
  )
  family_label <- if (inherits(family, "family")) {
    label <- paste0(family$family, " (", family$link, ")")
    if (!identical(link$name, family$link)) {
      label <- paste0(label, "; response transformed with ", link$name)
    }
    label
  } else {
    paste0("link ", link$name)
  }
  list(f = f, pred = pred, baseline = baseline, means = means, family = family_label,
       disp_df = disp_df,
       averaged_over = setdiff(emm_info$predictors, vars))
}

# Reference grid of the model, with every value of the stressor, time and
# `by` columns, averaged over its other predictors as emmeans does.
external_grid <- function(fit, design, data, call, ...) {
  spec <- design$spec
  vars <- c(spec$stressors, spec$time, spec$by)
  dots <- list(...)
  taken <- intersect(names(dots), c("cov.keep", "specs", "object"))
  if (length(taken) > 0L) {
    abort_input(sprintf(
      "%s cannot be passed on to emmeans: must sets %s from `design`.",
      format_items(paste0("`", taken, "`")), if (length(taken) == 1L) "it" else "them"
    ), call = call)
  }
  rescale <- intersect(names(dots), c("regrid", "transform"))
  if (length(rescale) > 0L) {
    abort_input(sprintf(
      paste(
        "%s cannot be passed on to emmeans: the means are simulated on the scale",
        "of the link function and back-transformed by must."
      ),
      format_items(paste0("`", rescale, "`"))
    ), call = call)
  }
  if (!is.null(dots$at) && any(names(dots$at) %in% vars)) {
    abort_input(sprintf(
      paste(
        "`at` cannot set %s: the values of the stressor, time and `by` columns",
        "come from `design`."
      ),
      format_items(quote_items(intersect(names(dots$at), vars)))
    ), call = call)
  }
  rg <- tryCatch(
    if (is.null(data)) {
      emmeans::ref_grid(fit, cov.keep = vars, ...)
    } else {
      emmeans::ref_grid(fit, cov.keep = vars, data = data, ...)
    },
    error = function(e) {
      abort_input(sprintf(
        paste(
          "emmeans could not build the reference grid of `fit` (%s): %s If the data",
          "cannot be recovered from the model, give them as `data`."
        ),
        describe_value(fit), conditionMessage(e)
      ), call = call)
    }
  )
  roles <- rg@roles
  if (length(roles$multresp) > 0L) {
    abort_input(
      "Models with a multivariate response are not supported: fit one model per response.",
      call = call
    )
  }
  absent <- setdiff(vars, roles$predictors)
  if (length(absent) > 0L) {
    abort_input(sprintf(
      paste(
        "%s of `design` %s not a predictor of the model. The stressor, time and",
        "`by` columns of the design must be predictors of the model fitted to its",
        "data."
      ),
      format_items(quote_items(absent), max = Inf),
      if (length(absent) == 1L) "is" else "are"
    ), call = call)
  }
  for (s in spec$stressors) {
    check_model_levels(spec$levels[[s]], rg@levels[[s]], s, call)
  }
  for (b in spec$by) {
    check_model_levels(unique(design$cells[[b]]), rg@levels[[b]], b, call)
  }
  emm <- emmeans::emmeans(rg, specs = vars)
  list(emm = emm, predictors = roles$predictors)
}

# Every level of `column` in the design must be a value of that predictor in
# the model (numbers are compared as numbers).
check_model_levels <- function(design_levels, model_levels, column, call) {
  found <- if (is.numeric(model_levels)) {
    !is.na(match_numeric(suppressWarnings(as.numeric(design_levels)), model_levels))
  } else {
    design_levels %in% as.character(model_levels)
  }
  if (!all(found)) {
    abort_input(sprintf(
      paste(
        "The model has no level %s of '%s' (its values: %s). Describe the design with",
        "the data the model was fitted to."
      ),
      format_items(quote_items(design_levels[!found])), column,
      format_items(quote_items(as.character(model_levels)), max = 10L)
    ), call = call)
  }
  invisible(TRUE)
}

# Labels of the values `x` of a grid column, as the design writes them
# (`labels`); NA for values that are not in the design.
grid_labels <- function(x, labels) {
  if (is.numeric(x)) {
    labels[match_numeric(x, suppressWarnings(as.numeric(labels)))]
  } else {
    x <- as.character(x)
    ifelse(x %in% labels, x, NA_character_)
  }
}

# Position in `table` of each number in `x`, with a relative tolerance for
# values that went through text.
match_numeric <- function(x, table) {
  vapply(x, function(v) {
    if (is.na(v)) {
      return(NA_integer_)
    }
    hit <- which(abs(table - v) <= 1e-8 * max(1, abs(v)))
    if (length(hit) == 0L) NA_integer_ else hit[[1L]]
  }, integer(1))
}

# Inverse link (linkinv(), mu.eta() and name) from the `misc` slot of an
# emmeans grid: the link function of the model or the transformation of its
# response, with the multiplier and offset that emmeans records for
# transformations such as log(y + 1).
external_link <- function(misc, call) {
  if (!is.null(misc$tran2)) {
    abort_input(paste(
      "Models with both a transformed response and a link function are not",
      "supported: fit the model to the untransformed response."
    ), call = call)
  }
  tran <- misc$tran
  name <- if (is.null(tran)) "identity" else if (is.character(tran)) tran else tran$name
  link <- if (is.null(tran)) {
    stats::make.link("identity")
  } else if (is.list(tran)) {
    tran
  } else {
    named_link(tran)
  }
  if (is.null(link) || !is.function(link$linkinv) || !is.function(link$mu.eta)) {
    abort_input(sprintf(
      paste(
        "The transformation '%s' of the model is not supported: fit the model to",
        "the untransformed response, or with a link function."
      ),
      paste(format(name), collapse = " ")
    ), call = call)
  }
  if (is.null(name)) {
    name <- "custom"
  }
  mult <- if (is.null(misc$tran.mult)) 1 else misc$tran.mult
  offset <- if (is.null(misc$tran.offset)) 0 else misc$tran.offset
  if (mult != 1 || offset != 0) {
    inv <- link$linkinv
    deriv <- link$mu.eta
    link <- list(
      linkinv = function(eta) inv(eta / mult) - offset,
      mu.eta = function(eta) deriv(eta / mult) / mult
    )
    name <- sprintf("%s%s(y%s)", if (mult != 1) paste0(format(mult), "*") else "", name,
                    if (offset != 0) paste0(" + ", format(offset)) else "")
  }
  list(linkinv = link$linkinv, mu.eta = link$mu.eta, name = name)
}

# Links and response transformations that emmeans names, as in its own
# back-transformations.
named_link <- function(name) {
  switch(
    name,
    identity = , log = , logit = , probit = , cauchit = , cloglog = , inverse = ,
    "1/mu^2" = stats::make.link(name),
    sqrt = list(linkinv = function(eta) pmax(0, eta)^2, mu.eta = function(eta) 2 * pmax(0, eta)),
    log10 = list(linkinv = function(eta) 10^eta, mu.eta = function(eta) 10^eta * log(10)),
    log2 = list(linkinv = function(eta) 2^eta, mu.eta = function(eta) 2^eta * log(2)),
    log1p = list(linkinv = expm1, mu.eta = exp),
    NULL
  )
}

# Means that cannot be used: not estimable from the terms of the model
# ("nonestimable"), without a finite standard error ("no_se"), or at the
# boundary of a log, logit, probit, cloglog or cauchit link, where a
# combination with only zeros (or only successes) puts them, with a standard
# error that the model cannot estimate ("boundary").
external_problems <- function(eta, se, link_name) {
  out <- rep("", length(eta))
  low <- link_name %in% c("log", "logit", "probit", "cloglog", "cauchit")
  high <- link_name %in% c("logit", "probit", "cloglog", "cauchit")
  extreme <- (low & eta < -10) | (high & eta > 10)
  out[!is.na(eta) & !is.na(se) & se > 10 & extreme] <- "boundary"
  out[!is.na(eta) & is.na(se)] <- "no_se"
  out[is.na(eta)] <- "nonestimable"
  out
}

# Name of the response: the left-hand side of the formula of the model.
external_response <- function(fit) {
  out <- tryCatch({
    fo <- stats::formula(fit)
    if (length(fo) == 3L) paste(deparse(fo[[2L]], width.cutoff = 500L), collapse = " ") else ""
  }, error = function(e) "")
  if (length(out) != 1L || is.na(out) || !nzchar(out)) "response" else out
}

# Whether the response cannot be negative (for the note on AD expectations
# below zero): a link whose inverse is positive, or a response without
# negative values.
external_non_negative <- function(fit, link) {
  if (link$name %in% c("log", "logit", "probit", "cloglog", "cauchit", "sqrt", "log10",
                       "log2", "inverse", "1/mu^2")) {
    return(TRUE)
  }
  y <- tryCatch(stats::model.response(stats::model.frame(fit)), error = function(e) NULL)
  is.numeric(y) && length(y) > 0L && all(y >= 0, na.rm = TRUE)
}

# Residual degrees of freedom of a glm() whose dispersion is estimated, as
# summary.glm() treats it (every family except Poisson and binomial; not
# MASS::glm.nb(), whose summary uses the normal distribution); NULL otherwise.
glm_dispersion_df <- function(fit) {
  if (!inherits(fit, "glm") || inherits(fit, "negbin") || inherits(fit, "gam")) {
    return(NULL)
  }
  family <- tryCatch(stats::family(fit)$family, error = function(e) NULL)
  if (!is.character(family) || family[[1L]] %in% c("poisson", "binomial")) {
    return(NULL)
  }
  df <- fit$df.residual
  if (is.numeric(df) && length(df) == 1L && is.finite(df) && df > 0) df else NULL
}
