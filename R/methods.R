#' Summarise and plot estimated interactions
#'
#' `summary()` counts, for each response, sampling time, null model and
#' estimand, the interactions that were tested and those that deviate from
#' the null model, tabulates their classes (see [classify_interactions()])
#' and, with two or more null models, how often they agree (see
#' [compare_null_models()]). `plot()` draws the estimated interactions of
#' one response and null model with their intervals.
#'
#' @param object,x A `must_interactions` object from
#'   [estimate_interactions()].
#' @param alpha Significance level applied to the adjusted p-values.
#'   Default 0.05.
#' @param response Name of the response to plot. Default: the first one.
#' @param null_model Null model to plot, among those estimated: `"AD"`
#'   (default), `"MU"`, `"DO"`, `"BL"`, `"SAM"` or `"CA"`.
#' @param estimand `"net"` (default) or `"pure"`.
#' @param time Sampling times to plot (as text, including `"overall"`).
#'   Default: all.
#' @param ... Further arguments passed to [graphics::plot.default()] (for
#'   example `main` or `xlim`); ignored by `summary()`.
#'
#' @return `summary()` returns an object of class
#'   `summary.must_interactions`, a list with `interactions` (one row per
#'   response, time, null model and estimand: `n_sets`, `n_tested` (status
#'   `"ok"` or `"boundary"`), `n_deviating` (adjusted p-value below
#'   `alpha`), `n_positive`, `n_negative`, `n_boundary`, `n_not_estimable`
#'   and `n_other` (not identifiable, not defined or not applicable)),
#'   `classes` (number of sets of each class per response, time and null
#'   model), `agreement` (with two or more null models: number of
#'   interactions per response, time and estimand that are `n_robust`,
#'   `n_model_dependent`, `n_conflicting`, `n_none` or `n_not_compared`;
#'   `NULL` otherwise), `verdicts` (with an equivalence `margin` in
#'   [estimate_interactions()]: number of interactions per response, time,
#'   null model and estimand with each verdict, `n_deviation`,
#'   `n_negligible`, `n_equivalent` and `n_inconclusive`; `NULL` otherwise)
#'   and `settings`. `plot()` returns `x` invisibly.
#'
#' @details
#' In the plot, filled points are interactions whose adjusted p-value is
#' below `alpha`, open points the rest. Interactions tested from a
#' combination with only zeros (`status = "boundary"`) have no interval and
#' are drawn as triangles; infinite estimates (MU with a combination of
#' zeros) are drawn at the edge of the plot. Rows without an estimate are
#' left out.
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
#' est <- estimate_interactions(fits, null_models = c("AD", "MU"), n_sim = 2000, seed = 1)
#' summary(est)
#' plot(est, null_model = "MU")
#' @name must_interactions-methods
NULL

#' @rdname must_interactions-methods
#' @export
summary.must_interactions <- function(object, alpha = 0.05, ...) {
  check_number(alpha, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  est <- object$estimates
  int <- est[est$component == "interaction", , drop = FALSE]
  keys <- c("response", "time", "null_model", "estimand")
  template <- data.frame(
    response = character(), time = character(), null_model = character(),
    estimand = character(), n_sets = integer(), n_tested = integer(),
    n_deviating = integer(), n_positive = integer(), n_negative = integer(),
    n_boundary = integer(), n_not_estimable = integer(), n_other = integer(),
    stringsAsFactors = FALSE
  )
  counts <- template
  if (nrow(int) > 0L) {
    group <- do.call(paste_key, unname(as.list(int[keys])))
    idx <- split(seq_len(nrow(int)), factor(group, levels = unique(group)))
    counts <- bind_rows_base(lapply(idx, function(i) {
      part <- int[i, , drop = FALSE]
      tested <- part$status %in% c("ok", "boundary")
      deviating <- tested & !is.na(part$p_adjusted) & part$p_adjusted < alpha
      data.frame(
        response = part$response[[1L]], time = part$time[[1L]],
        null_model = part$null_model[[1L]], estimand = part$estimand[[1L]],
        n_sets = nrow(part), n_tested = sum(tested), n_deviating = sum(deviating),
        n_positive = sum(deviating & part$estimate > 0),
        n_negative = sum(deviating & part$estimate < 0),
        n_boundary = sum(part$status == "boundary"),
        n_not_estimable = sum(part$status == "not_estimable"),
        n_other = sum(!tested & part$status != "not_estimable"),
        stringsAsFactors = FALSE
      )
    }), template)
  }
  cls <- classify_interactions(object, alpha = alpha)
  cls <- cls[cls$status == "ok", , drop = FALSE]
  class_levels <- c("AD", "+S", "-S", "+A", "-A")
  classes <- data.frame(response = character(), time = character(),
                        null_model = character(), stringsAsFactors = FALSE)
  for (cl in class_levels) {
    classes[[cl]] <- integer()
  }
  if (nrow(cls) > 0L) {
    group <- paste_key(cls$response, cls$time, cls$null_model)
    idx <- split(seq_len(nrow(cls)), factor(group, levels = unique(group)))
    classes <- bind_rows_base(lapply(idx, function(i) {
      row <- data.frame(response = cls$response[[i[[1L]]]], time = cls$time[[i[[1L]]]],
                        null_model = cls$null_model[[i[[1L]]]], stringsAsFactors = FALSE)
      for (cl in class_levels) {
        row[[cl]] <- sum(cls$class[i] == cl, na.rm = TRUE)
      }
      row
    }), classes)
  }
  agreement <- NULL
  if (length(object$settings$null_models) >= 2L) {
    agreement <- agreement_counts(compare_null_models(object, alpha = alpha))
  }
  verdicts <- NULL
  if (!is.null(object$settings$margin)) {
    verdicts <- verdict_counts(int, alpha)
  }
  structure(
    list(interactions = counts, classes = classes, agreement = agreement,
         verdicts = verdicts, settings = object$settings, alpha = alpha),
    class = "summary.must_interactions"
  )
}

# Number of interactions with each verdict (see interaction_verdict()) per
# response, time, null model and estimand.
verdict_counts <- function(int, alpha) {
  template <- data.frame(
    response = character(), time = character(), null_model = character(),
    estimand = character(), n_deviation = integer(), n_negligible = integer(),
    n_equivalent = integer(), n_inconclusive = integer(), stringsAsFactors = FALSE
  )
  if (nrow(int) == 0L) {
    return(template)
  }
  v <- interaction_verdict(int$p_adjusted, int$p_equivalence_adjusted, int$status, alpha,
                           has_margin = TRUE)
  keys <- c("response", "time", "null_model", "estimand")
  group <- do.call(paste_key, unname(as.list(int[keys])))
  idx <- split(seq_len(nrow(int)), factor(group, levels = unique(group)))
  bind_rows_base(lapply(idx, function(i) {
    data.frame(
      response = int$response[[i[[1L]]]], time = int$time[[i[[1L]]]],
      null_model = int$null_model[[i[[1L]]]], estimand = int$estimand[[i[[1L]]]],
      n_deviation = sum(v[i] %in% "deviation"), n_negligible = sum(v[i] %in% "negligible"),
      n_equivalent = sum(v[i] %in% "equivalent"),
      n_inconclusive = sum(v[i] %in% "inconclusive"), stringsAsFactors = FALSE
    )
  }), template)
}

# Number of interactions of each kind of agreement per response, time and
# estimand.
agreement_counts <- function(cmp) {
  template <- data.frame(
    response = character(), time = character(), estimand = character(),
    n_robust = integer(), n_model_dependent = integer(), n_conflicting = integer(),
    n_none = integer(), n_not_compared = integer(), stringsAsFactors = FALSE
  )
  if (nrow(cmp) == 0L) {
    return(template)
  }
  group <- paste_key(cmp$response, cmp$time, cmp$estimand)
  idx <- split(seq_len(nrow(cmp)), factor(group, levels = unique(group)))
  bind_rows_base(lapply(idx, function(i) {
    a <- cmp$agreement[i]
    data.frame(
      response = cmp$response[[i[[1L]]]], time = cmp$time[[i[[1L]]]],
      estimand = cmp$estimand[[i[[1L]]]],
      n_robust = sum(a %in% "robust"), n_model_dependent = sum(a %in% "model_dependent"),
      n_conflicting = sum(a %in% "conflicting"), n_none = sum(a %in% "none"),
      n_not_compared = sum(is.na(a)), stringsAsFactors = FALSE
    )
  }), template)
}

#' @export
print.summary.must_interactions <- function(x, ...) {
  s <- x$settings
  cat(sprintf(
    "Interactions: %s; %s; adjusted p-values (%s) below %s\n",
    if (identical(s$inference, "simulation")) "Krinsky-Robb simulation" else "delta method",
    correction_label(s),
    s$p_adjust, format(x$alpha)
  ))
  tab <- x$interactions
  if (nrow(tab) == 0L) {
    cat("No interaction was estimated.\n")
    return(invisible(x))
  }
  cat("Deviating / tested interactions by sampling time:\n")
  print(summary_wide(tab), row.names = FALSE)
  if (nrow(x$classes) > 0L) {
    cat("\nClasses of the net deviations (Piggott et al. 2015) by sampling time:\n")
    print(classes_wide(x$classes), row.names = FALSE)
  }
  if (!is.null(x$verdicts) && nrow(x$verdicts) > 0L) {
    cat(sprintf(paste0("\nVerdicts with an equivalence margin of +/-%s of the reference cell,",
                       " deviation/negligible/equivalent/inconclusive:\n"),
                format(x$settings$margin)))
    print(verdict_wide(x$verdicts), row.names = FALSE)
  }
  if (!is.null(x$agreement) && nrow(x$agreement) > 0L) {
    cat(sprintf("\nAgreement between null models (%s), robust/model-dependent/conflicting/none:\n",
                paste(x$settings$null_models, collapse = ", ")))
    print(agreement_wide(x$agreement), row.names = FALSE)
  }
  cat("\nCounts of boundary tests and of rows not estimated: `$interactions`;",
      "classes: `$classes`.\n")
  invisible(x)
}

# One row per response, null model and estimand; one column per time with
# "deviation/negligible/equivalent/inconclusive".
verdict_wide <- function(vd) {
  time <- ifelse(is.na(vd$time), "all", vd$time)
  times <- unique(time)
  key <- paste_key(vd$response, vd$null_model, vd$estimand)
  first <- !duplicated(key)
  out <- data.frame(response = vd$response[first], model = vd$null_model[first],
                    estimand = vd$estimand[first], stringsAsFactors = FALSE)
  for (t in times) {
    hit <- match(paste_key(key[first], t), paste_key(key, time))
    out[[t]] <- ifelse(is.na(hit), "", sprintf(
      "%d/%d/%d/%d", vd$n_deviation[hit], vd$n_negligible[hit], vd$n_equivalent[hit],
      vd$n_inconclusive[hit]
    ))
  }
  out
}

# One row per response and estimand; one column per time with
# "robust/model-dependent/conflicting/none".
agreement_wide <- function(ag) {
  time <- ifelse(is.na(ag$time), "all", ag$time)
  times <- unique(time)
  key <- paste_key(ag$response, ag$estimand)
  first <- !duplicated(key)
  out <- data.frame(response = ag$response[first], estimand = ag$estimand[first],
                    stringsAsFactors = FALSE)
  for (t in times) {
    hit <- match(paste_key(key[first], t), paste_key(key, time))
    out[[t]] <- ifelse(is.na(hit), "", sprintf(
      "%d/%d/%d/%d", ag$n_robust[hit], ag$n_model_dependent[hit], ag$n_conflicting[hit],
      ag$n_none[hit]
    ))
  }
  out
}

# One row per response and null model; one column per time with the classes
# found ("3 AD, 1 -S").
classes_wide <- function(cl) {
  levels <- c("AD", "+S", "-S", "+A", "-A")
  time <- ifelse(is.na(cl$time), "all", cl$time)
  times <- unique(time)
  key <- paste_key(cl$response, cl$null_model)
  first <- !duplicated(key)
  out <- data.frame(response = cl$response[first], model = cl$null_model[first],
                    stringsAsFactors = FALSE)
  text <- vapply(seq_len(nrow(cl)), function(i) {
    n <- unlist(cl[i, levels])
    paste(paste(n[n > 0], levels[n > 0]), collapse = ", ")
  }, character(1))
  for (t in times) {
    hit <- match(paste_key(key[first], t), paste_key(key, time))
    out[[t]] <- ifelse(is.na(hit), "", text[hit])
  }
  out
}

# One row per response, null model and estimand; one column per time with
# "deviating/tested".
summary_wide <- function(tab) {
  time <- ifelse(is.na(tab$time), "all", tab$time)
  times <- unique(time)
  key <- paste_key(tab$response, tab$null_model, tab$estimand)
  first <- !duplicated(key)
  out <- data.frame(response = tab$response[first], model = tab$null_model[first],
                    estimand = tab$estimand[first], stringsAsFactors = FALSE)
  for (t in times) {
    hit <- match(paste_key(key[first], t), paste_key(key, time))
    out[[t]] <- ifelse(is.na(hit), "",
                       sprintf("%d/%d", tab$n_deviating[hit], tab$n_tested[hit]))
  }
  out
}

#' @rdname must_interactions-methods
#' @export
plot.must_interactions <- function(x, response = NULL, null_model = "AD", estimand = "net",
                                   time = NULL, alpha = 0.05, ...) {
  est <- x$estimates
  responses <- unique(est$response)
  if (length(responses) == 0L) {
    abort_input("There are no estimates to plot.")
  }
  if (is.null(response)) {
    response <- responses[[1L]]
  }
  check_choice(response, responses)
  check_choice(null_model, x$settings$null_models)
  check_choice(estimand, c("net", "pure"))
  check_number(alpha, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  keep <- est$response == response & est$null_model == null_model &
    est$component == "interaction" & est$estimand %in% estimand & !is.na(est$estimate)
  if (!is.null(time)) {
    keep <- keep & est$time %in% as.character(time)
  }
  d <- est[keep, , drop = FALSE]
  if (nrow(d) == 0L) {
    abort_input(sprintf(
      "No %s interaction (%s) of '%s' has an estimate to plot.", null_model, estimand, response
    ))
  }
  control <- x$design$spec$control
  label <- set_combination(d$levels, d$stressors)
  label <- ifelse(nzchar(d$background), paste0(label, " [", compact_levels(d$background, control), "]"),
                  label)
  if (!all(is.na(d$time))) {
    label <- paste0(label, " @ ", d$time)
  }
  if (length(unique(d$context)) > 1L) {
    label <- paste0(d$context, ": ", label)
  }
  n <- nrow(d)
  y <- rev(seq_len(n))
  finite <- c(d$estimate, d$lower, d$upper)
  finite <- finite[is.finite(finite)]
  xlim <- if (length(finite) > 0L) range(c(finite, 0)) else c(-1, 1)
  if (diff(xlim) == 0) {
    xlim <- xlim + c(-1, 1)
  }
  point <- pmin(pmax(d$estimate, xlim[[1L]]), xlim[[2L]])
  deviating <- !is.na(d$p_adjusted) & d$p_adjusted < alpha
  boundary <- d$status == "boundary"
  pch <- ifelse(boundary, ifelse(deviating, 17L, 2L), ifelse(deviating, 19L, 1L))
  label_cex <- 0.7
  width <- max(graphics::strwidth(label, units = "inches", cex = label_cex))
  left <- min(width + 0.25, 0.55 * graphics::grconvertX(1, from = "ndc", to = "inches"))
  withr::local_par(list(mai = c(0.9, left, 0.8, 0.2)))
  args <- list(
    x = point, y = y, xlim = xlim, ylim = c(0.5, n + 0.5), pch = pch, yaxt = "n",
    ylab = "", xlab = sprintf("%s interaction (%s)", null_model, null_model_scale(null_model)),
    main = sprintf("%s: %s, %s", response, null_model, estimand), bty = "l"
  )
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(graphics::plot.default, args)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  has_ci <- is.finite(d$lower) & is.finite(d$upper)
  if (any(has_ci)) {
    graphics::segments(d$lower[has_ci], y[has_ci], d$upper[has_ci], y[has_ci])
  }
  # mtext() draws every label (axis() drops those that overlap).
  graphics::mtext(label, side = 2, at = y, las = 1, line = 0.3, cex = label_cex)
  # Centred on the figure, just above the plotting region.
  graphics::legend(
    x = graphics::grconvertX(0.5, from = "ndc", to = "user"),
    y = graphics::grconvertY(1, from = "npc", to = "user"), xjust = 0.5, yjust = 0,
    xpd = NA, horiz = TRUE, bty = "n", cex = 0.75,
    legend = c(sprintf("adjusted p < %s", format(alpha)), "not significant",
               "only zeros (no interval)"),
    pch = c(19L, 1L, 17L)
  )
  invisible(x)
}

# Short label of the levels of the stressors of a set: "salinity=high +
# warming=yes".
set_combination <- function(levels, stressors) {
  vapply(seq_along(levels), function(i) {
    parts <- strsplit(levels[[i]], " | ", fixed = TRUE)[[1L]]
    in_set <- strsplit(stressors[[i]], " x ", fixed = TRUE)[[1L]]
    parts <- parts[sub("=.*$", "", parts) %in% in_set]
    paste(parts, collapse = " + ")
  }, character(1))
}
