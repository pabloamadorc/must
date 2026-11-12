#' Plot estimated interactions
#'
#' Draws the interactions estimated by [estimate_interactions()] in one of
#' four ways, chosen with `type`:
#'
#' * `"effects"` (default): for each interaction, the effect of each of its
#'   stressors alone, the combined effect expected under the null model (net
#'   interactions only), the combined effect observed and the interaction,
#'   all on the scale of the null model, with the class of the interaction
#'   at the right (see [classify_interactions()]). It shows how the
#'   interaction arises: it is the distance between the expected and the
#'   observed combined effects.
#' * `"interactions"`: the interactions alone, one row each, on the scale of
#'   the null model or as relative deviations (`scale`), optionally grouped
#'   under headers (`group_by`): an overview of many sets and times.
#' * `"models"`: the net or pure interaction under each null model, as a
#'   deviation relative to the reference cell (the same quantity for every
#'   null model, in percent), with its class: whether the conclusion depends
#'   on the null model (see [compare_null_models()]).
#' * `"time"`: for designs with sampling times, one figure per set with the
#'   effects of the stressors alone and the expected and observed combined
#'   effects over time (top) and the net interaction and its class at each
#'   time (bottom).
#'
#' @param x A `must_interactions` object from [estimate_interactions()].
#' @param type What to draw: `"effects"` (default), `"interactions"`,
#'   `"models"` or `"time"` (see Description).
#' @param response Name of the response to plot. Default: the first one.
#' @param null_model Null model to plot, among those estimated. Default: the
#'   first one estimated (`type = "effects"` and `"time"` need the effects of
#'   the stressors alone, which the dominance model, `"DO"`, does not give,
#'   so they skip it). For `type = "models"`, the null models to compare, in
#'   the order given; default: all those estimated.
#' @param estimand `"net"` (default) or `"pure"` (sets of three or more
#'   stressors only). `type = "time"` draws net interactions.
#' @param time Sampling times to plot (as text, including `"overall"`).
#'   Default: all. `type = "time"` draws the numeric times and leaves out
#'   `"overall"`.
#' @param set_id,context Sets (`set_id`, as in the estimates) and `by`
#'   groups (`context`, as in the estimates, e.g. `"species=T. gravida"`)
#'   to plot. Default: all.
#' @param scale For `type = "interactions"`: `"model"` (default; the scale
#'   of the null model, on which the tests are done) or `"relative"`
#'   (observed minus expected, relative to the reference cell; comparable
#'   between null models).
#' @param group_by For `type = "interactions"`: `"none"` (default; one label
#'   per row), `"set"`, `"time"` or `"context"`, to write what the rows
#'   share once, as a header above them.
#' @param classes How the class of each interaction is written: `"auto"`
#'   (default: names when the plot has a single interaction, codes
#'   otherwise), `"names"` (`type = "effects"` and `"interactions"` only),
#'   `"codes"` or `"none"`. The codes are `S+`, `S-`, `A+` and `A-` for the
#'   classes `+S`, `-S`, `+A` and `-A` of [classify_interactions()]
#'   (synergism or antagonism, positive or negative), `n.d.` for no
#'   deviation detected (class `AD`), `n/a` when must does not classify the
#'   interaction (dominance model and pure interactions) and `NC` when it
#'   could not be classified.
#' @param alpha Significance level applied to the adjusted p-values (filled
#'   points) and passed to [classify_interactions()]. Default 0.05.
#' @param labels Optional named character vector with labels for the sets,
#'   named by `set_id` (for example `c(S01 = "High salinity + warming")`).
#'   Default: the levels of the stressors of each set.
#' @param stressor_labels Optional named character vector with the names to
#'   show for the stressors, named by stressor (for example
#'   `c(temperature = "Warming")`).
#' @param label_style How sets are labelled by default: `"full"` (default;
#'   `"salinity: high + warming: yes"`), `"short"` (stressor names cut to
#'   four letters) or `"levels"` (levels only: shorter, but ambiguous when
#'   stressors share level names).
#' @param label_width Maximum number of characters per line of the labels
#'   in the margin of `type = "effects"` and `"interactions"`; longer labels
#'   are wrapped. Default 42.
#' @param per_page Number of blocks per page: interactions (each with its
#'   effects) for `type = "effects"`, rows or groups of rows under a header
#'   for `type = "interactions"`, and sets at one sampling time for `type =
#'   "models"` (default 4). By default, `type = "effects"` and
#'   `"interactions"` put on each page as many blocks as fit in the height
#'   of the device, up to 4 for `"effects"`. `type = "time"` draws one set
#'   per page. Pages are drawn one after the other; in an interactive
#'   session R asks before each new one.
#' @param xlim,ylim,ylim_interaction Limits of the axis of the estimates:
#'   `xlim` for `type = "effects"`, `"interactions"` and `"models"` (as
#'   proportions for `"models"`: `c(-0.5, 0.5)` is -50% to 50%), and `ylim`
#'   and `ylim_interaction` for the top and bottom panels of `type = "time"`.
#'   By default the limits cover the estimates and their intervals; for
#'   `type = "models"`, intervals that reach beyond the estimates by more
#'   than the range of the estimates are cut. Cut intervals, and estimates
#'   outside given limits, end in an arrow.
#' @param main Title. Default: the response, the null model and the
#'   estimand.
#' @param ... Not used.
#'
#' @return Invisibly, a data frame with one row per estimate drawn and
#'   columns `response`, `context`, `time`, `set_id`, `null_model`,
#'   `estimand`, `component` (`"single"`, `"expected"`, `"combined"` or
#'   `"interaction"`), `stressor`, `estimate`, `lower`, `upper` (the values
#'   drawn: relative deviations for `type = "models"` and `scale =
#'   "relative"`), `p_adjusted`, `status`, `class_code` (for interactions)
#'   and `page`.
#'
#' @details
#' Filled points are estimates whose adjusted p-value is below `alpha`, open
#' points the rest; a diamond (`type = "effects"` and `"interactions"`) or a
#' cross (`type = "models"` and `"time"`) marks an estimate without p-value.
#' Interactions or effects tested from a treatment combination with only
#' zeros (`status = "boundary"`) have no interval and are drawn as
#' triangles in `type = "effects"` and `"interactions"`; infinite estimates
#' (MU with a combination of zeros) are drawn as an arrow at the edge.
#' Estimates that could not be computed, or that the null model does not
#' define, are marked `NE`.
#'
#' The effects of the stressors alone and the combined effects are on the
#' scale of the null model and relative to the reference cell of the set
#' (the control, or the background levels with `background = "all"` in
#' [describe_design()]): differences in response units for AD, differences
#' of logarithms for MU and responses relative to the reference for BL, SAM
#' and CA. The combined effect expected under the null model is the
#' observed combined effect minus the net interaction; its interval is not
#' computed, so it is drawn without one.
#'
#' What the rows of a page share (the `by` group, the sampling time) is
#' written under the title rather than in every label.
#'
#' @seealso [classify_interactions()] for the classes and
#'   [compare_null_models()] for the agreement between null models.
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
#' # How the interaction arises under each null model
#' plot(est)
#' plot(est, null_model = "MU")
#' # Both null models on the same axis
#' plot(est, type = "models")
#' @export
plot.must_interactions <- function(x, type = "effects", response = NULL, null_model = NULL,
                                   estimand = "net", time = NULL, set_id = NULL,
                                   context = NULL, scale = "model", group_by = "none",
                                   classes = "auto", alpha = 0.05, labels = NULL,
                                   stressor_labels = NULL, label_style = "full",
                                   label_width = 42, per_page = NULL, xlim = NULL,
                                   ylim = NULL, ylim_interaction = NULL, main = NULL, ...) {
  call <- generic_call("plot")
  if (...length() > 0L) {
    given <- names(list(...))
    given <- if (is.null(given)) rep("", ...length()) else given
    shown <- ifelse(nzchar(given), quote_items(given), "an unnamed value")
    abort_input(sprintf(
      "Unknown argument(s): %s. The arguments of plot() are listed in ?plot.must_interactions.",
      format_items(shown)
    ), call = call)
  }
  check_choice(type, c("effects", "interactions", "models", "time"), call = call)
  check_choice(estimand, c("net", "pure"), call = call)
  check_choice(scale, c("model", "relative"), call = call)
  check_choice(group_by, c("none", "set", "time", "context"), call = call)
  check_choice(classes, c("auto", "names", "codes", "none"), call = call)
  check_choice(label_style, c("full", "short", "levels"), call = call)
  check_number(alpha, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE,
               call = call)
  check_number(label_width, lower = 8, call = call)
  if (!is.null(per_page)) {
    check_number(per_page, lower = 1, call = call)
    if (per_page != round(per_page)) abort_input("`per_page` must be a whole number.", call = call)
  }
  check_plot_limits(xlim, "xlim", call)
  check_plot_limits(ylim, "ylim", call)
  check_plot_limits(ylim_interaction, "ylim_interaction", call)
  check_string(main, allow_null = TRUE, call = call)
  check_plot_type_args(type, scale, group_by, classes, estimand, per_page, xlim, ylim,
                       ylim_interaction, call)
  est <- x$estimates
  responses <- unique(est$response)
  if (length(responses) == 0L) abort_input("There are no estimates to plot.", call = call)
  if (is.null(response)) response <- responses[[1L]]
  check_choice(response, responses, call = call)
  null_model <- check_plot_null_model(null_model, type, x$settings$null_models, call)
  mine <- est[est$response == response, , drop = FALSE]
  time <- check_plot_filter(time, mine$time, "time", call)
  set_id <- check_plot_filter(set_id, mine$set_id, "set_id", call)
  context <- check_plot_filter(context, mine$context, "context", call)
  labels <- check_plot_labels(labels, "labels", call)
  if (!is.null(labels) && !all(names(labels) %in% est$set_id)) {
    abort_input(paste("The names of `labels` must be set identifiers of the estimates",
                      "(`set_id`, e.g. 'S01')."), call = call)
  }
  stressors <- x$design$spec$stressors
  stressor_labels <- check_plot_labels(stressor_labels, "stressor_labels", call)
  if (!is.null(stressor_labels) && !all(names(stressor_labels) %in% stressors)) {
    abort_input(sprintf("The names of `stressor_labels` must be stressors: %s.",
                        format_items(quote_items(stressors), max = Inf)), call = call)
  }
  opts <- list(
    response = response, null_model = null_model, estimand = estimand, time = time,
    set_id = set_id, context = context, scale = scale, group_by = group_by,
    classes = classes, alpha = alpha, labels = labels, stressor_labels = stressor_labels,
    label_style = label_style, label_width = label_width, per_page = per_page,
    xlim = xlim, ylim = ylim, ylim_interaction = ylim_interaction, main = main,
    stressors = stressors, time_name = x$design$spec$time, level = x$settings$level,
    call = call
  )
  out <- switch(type,
                effects = plot_effects(x, opts),
                interactions = plot_interactions(x, opts),
                models = plot_models(x, opts),
                time = plot_time(x, opts))
  rownames(out) <- NULL
  invisible(out)
}

# Arguments -------------------------------------------------------------------

check_plot_limits <- function(lim, arg, call) {
  if (!is.null(lim) && (!is.numeric(lim) || length(lim) != 2L || any(!is.finite(lim)) ||
                          lim[[1L]] >= lim[[2L]])) {
    abort_input(sprintf("`%s` must be two finite numbers in increasing order.", arg),
                call = call)
  }
  invisible(lim)
}

check_plot_labels <- function(labels, arg, call) {
  if (is.null(labels)) return(NULL)
  if (!is.character(labels) || anyNA(labels) || is.null(names(labels)) ||
      any(!nzchar(names(labels))) || anyDuplicated(names(labels))) {
    abort_input(sprintf("`%s` must be a character vector with unique names.", arg),
                call = call)
  }
  labels
}

# Arguments that only some types use: an error rather than silently ignored.
check_plot_type_args <- function(type, scale, group_by, classes, estimand, per_page, xlim,
                                 ylim, ylim_interaction, call) {
  only <- function(arg, type) {
    abort_input(sprintf("`%s` applies to type = \"%s\" only.", arg, type), call = call)
  }
  if (scale != "model" && type != "interactions") {
    abort_input(paste0(
      "`scale` applies to type = \"interactions\" only: type = \"models\" always draws ",
      "relative deviations, and the other types the scale of the null model."
    ), call = call)
  }
  if (group_by != "none" && type != "interactions") only("group_by", "interactions")
  if (classes == "names" && type %in% c("models", "time")) {
    abort_input(paste0("type = \"models\" and \"time\" write the classes as codes: use ",
                       "classes = \"codes\" or \"none\"."), call = call)
  }
  if (estimand == "pure" && type == "time") {
    abort_input(paste0("type = \"time\" draws net interactions; for pure ones use ",
                       "type = \"effects\" or \"interactions\"."), call = call)
  }
  if (!is.null(per_page) && type == "time") {
    abort_input("`per_page` does not apply to type = \"time\", which draws one set per page.",
                call = call)
  }
  if (!is.null(xlim) && type == "time") {
    abort_input("type = \"time\" has two vertical axes: use `ylim` and `ylim_interaction`.",
                call = call)
  }
  if (!is.null(ylim) && type != "time") only("ylim", "time")
  if (!is.null(ylim_interaction) && type != "time") only("ylim_interaction", "time")
  invisible(TRUE)
}

check_plot_null_model <- function(null_model, type, models, call) {
  if (type == "models") {
    if (is.null(null_model)) return(models)
    if (!is.character(null_model) || length(null_model) == 0L || anyNA(null_model) ||
        !all(null_model %in% models)) {
      abort_input(sprintf("`null_model` must be among the null models estimated: %s.",
                          format_items(quote_items(models), max = Inf)), call = call)
    }
    return(unique(null_model))
  }
  needs_singles <- type %in% c("effects", "time")
  if (is.null(null_model)) {
    usable <- if (needs_singles) setdiff(models, "DO") else models
    if (length(usable) == 0L) {
      abort_input(paste0(
        "Only the dominance model (DO) was estimated, and it gives no effects of the ",
        "stressors alone: use type = \"interactions\" or \"models\"."
      ), call = call)
    }
    return(usable[[1L]])
  }
  check_choice(null_model, models, call = call)
  if (needs_singles && null_model == "DO") {
    abort_input(paste0(
      "The dominance model (DO) gives no effects of the stressors alone: use type = ",
      "\"interactions\" or \"models\", or another null model."
    ), call = call)
  }
  null_model
}

# Values to keep of a column of the estimates (time, set_id, context).
check_plot_filter <- function(value, available, arg, call) {
  if (is.null(value)) return(NULL)
  if (!(is.character(value) || is.numeric(value)) || length(value) == 0L || anyNA(value)) {
    abort_input(sprintf("`%s` must be a vector of values of the column `%s` of the estimates.",
                        arg, arg), call = call)
  }
  value <- as.character(value)
  available <- unique(available[!is.na(available) & nzchar(available)])
  unknown <- setdiff(value, available)
  if (length(unknown) > 0L) {
    abort_input(sprintf(
      "`%s` has values that are not in the estimates of this response: %s. %s",
      arg, format_items(quote_items(unknown)),
      if (length(available) == 0L) {
        sprintf("The estimates have no `%s`.", arg)
      } else {
        sprintf("Available: %s.", format_items(quote_items(available), max = 12L))
      }
    ), call = call)
  }
  value
}

# Selection of rows -----------------------------------------------------------

# Interaction rows of the response, null model(s) and estimand, filtered by
# time, set and context, in plotting order.
plot_interaction_rows <- function(x, opts) {
  est <- x$estimates
  keep <- est$response == opts$response & est$component == "interaction" &
    est$null_model %in% opts$null_model & est$estimand %in% opts$estimand &
    plot_filter(est, opts)
  d <- est[keep, , drop = FALSE]
  if (nrow(d) == 0L) {
    abort_input(sprintf(
      "No %s interaction of '%s' under %s matches the selection%s.",
      opts$estimand, opts$response, paste(opts$null_model, collapse = ", "),
      if (opts$estimand == "pure") {
        " (pure interactions exist for sets of three or more stressors)"
      } else {
        ""
      }
    ), call = opts$call)
  }
  d <- d[plot_order(d, opts$null_model), , drop = FALSE]
  rownames(d) <- NULL
  d
}

plot_filter <- function(est, opts) {
  keep <- rep(TRUE, nrow(est))
  if (!is.null(opts$time)) keep <- keep & est$time %in% opts$time
  if (!is.null(opts$set_id)) keep <- keep & est$set_id %in% opts$set_id
  if (!is.null(opts$context)) keep <- keep & est$context %in% opts$context
  keep
}

# By group, set, time (numeric, then "overall") and null model.
plot_order <- function(d, models) {
  set <- paste_key(d$context, d$set_id)
  order(match(d$context, unique(d$context)), match(set, unique(set)),
        match(d$time, plot_time_order(d$time)), match(d$null_model, models))
}

plot_time_order <- function(time) {
  t <- unique(time)
  num <- suppressWarnings(as.numeric(t))
  c(t[!is.na(num)][order(num[!is.na(num)])], t[is.na(num)])
}

# Rows with the given component of the same estimate as row `i` of `d` (same
# response, group, time, set and null model).
plot_same_set <- function(est, d, i, component) {
  same_time <- if (is.na(d$time[[i]])) {
    is.na(est$time)
  } else {
    !is.na(est$time) & est$time == d$time[[i]]
  }
  est[est$response == d$response[[i]] & est$context == d$context[[i]] & same_time &
        est$set_id == d$set_id[[i]] & est$null_model == d$null_model[[i]] &
        est$component == component, , drop = FALSE]
}

# Labels ----------------------------------------------------------------------

plot_stressor_names <- function(stressor, opts) {
  out <- stressor
  if (!is.null(opts$stressor_labels)) {
    hit <- match(stressor, names(opts$stressor_labels))
    out[!is.na(hit)] <- unname(opts$stressor_labels[hit[!is.na(hit)]])
  }
  out
}

# "salinity: high + warming: yes" from "salinity=high | warming=yes", with the
# stressors in `keep` (all when NULL), in the style of `label_style`.
plot_levels_text <- function(levels, keep, opts) {
  parts <- strsplit(levels, " | ", fixed = TRUE)[[1L]]
  stressor <- sub("=.*$", "", parts)
  level <- sub("^[^=]*=", "", parts)
  if (!is.null(keep)) {
    level <- level[stressor %in% keep]
    stressor <- stressor[stressor %in% keep]
  }
  shown <- plot_stressor_names(stressor, opts)
  if (opts$label_style == "short") shown <- substr(shown, 1L, 4L)
  text <- if (opts$label_style == "levels") level else paste0(shown, ": ", level)
  paste(text, collapse = " + ")
}

# Label of each set (the background, if any, in brackets), or the one given
# in `labels`.
plot_set_labels <- function(d, opts) {
  out <- vapply(seq_len(nrow(d)), function(i) {
    in_set <- strsplit(d$stressors[[i]], " x ", fixed = TRUE)[[1L]]
    text <- plot_levels_text(d$levels[[i]], in_set, opts)
    if (!is.na(d$background[[i]]) && nzchar(d$background[[i]])) {
      text <- paste0(text, " [", plot_levels_text(d$background[[i]], NULL, opts), "]")
    }
    text
  }, character(1))
  if (!is.null(opts$labels)) {
    hit <- match(d$set_id, names(opts$labels))
    out[!is.na(hit)] <- unname(opts$labels[hit[!is.na(hit)]])
  }
  out
}

# "species: T. gravida" from "species=T. gravida" ("" without `by` groups).
plot_context_labels <- function(context) {
  vapply(strsplit(context, " | ", fixed = TRUE), function(parts) {
    paste(sub("=", ": ", parts, fixed = TRUE), collapse = ", ")
  }, character(1))
}

plot_time_name <- function(opts) {
  name <- opts$time_name
  if (is.null(name) || length(name) != 1L || is.na(name) || !nzchar(name)) return("Time")
  paste0(toupper(substr(name, 1L, 1L)), substring(name, 2L))
}

# "Day 7", "Average over sampling times"; NA without times.
plot_time_labels <- function(time, opts) {
  out <- ifelse(time == "overall", "Average over sampling times",
                paste(plot_time_name(opts), time))
  out[is.na(time)] <- NA_character_
  out
}

# Wraps labels at their separators (" + ", " \u00b7 ") rather than inside
# "stressor: level"; `fits(line)` says whether a line fits, and a part too
# long for a line is wrapped at its spaces.
wrap_label <- function(s, fits) {
  if (is.na(s) || !nzchar(s)) return("")
  # A separator that ends a line may overhang it.
  fits_line <- function(l) {
    for (sep in c(" +", " \u00b7")) {
      if (endsWith(l, sep)) return(fits(substr(l, 1L, nchar(l) - 2L)))
    }
    fits(l)
  }
  words <- strsplit(s, " ", fixed = TRUE)[[1L]]
  ends <- words %in% c("+", "\u00b7")
  units <- vapply(split(words, cumsum(c(0, ends[-length(ends)]))), paste, character(1),
                  collapse = " ")
  lines <- character()
  current <- ""
  for (u in units) {
    trial <- if (nzchar(current)) paste(current, u) else u
    if (nzchar(current) && !fits_line(trial)) {
      lines <- c(lines, current)
      current <- u
    } else {
      current <- trial
    }
  }
  lines <- c(lines, current)
  lines <- unlist(lapply(lines, function(l) if (fits_line(l)) l else wrap_words(l, fits_line)))
  paste(lines, collapse = "\n")
}

wrap_words <- function(s, fits) {
  lines <- character()
  current <- ""
  for (w in strsplit(s, " ", fixed = TRUE)[[1L]]) {
    trial <- if (nzchar(current)) paste(current, w) else w
    if (nzchar(current) && !fits(trial)) {
      lines <- c(lines, current)
      current <- w
    } else {
      current <- trial
    }
  }
  c(lines, current)
}

# Wrapped to `width` characters per line.
plot_wrap <- function(text, width) {
  vapply(text, wrap_label, character(1), fits = function(l) nchar(l) <= width,
         USE.NAMES = FALSE)
}

plot_lines <- function(text) {
  lengths(strsplit(ifelse(nzchar(text), text, " "), "\n", fixed = TRUE))
}

plot_interval_note <- function(level) {
  if (is.numeric(level) && length(level) == 1L && is.finite(level)) {
    sprintf("%s%% confidence intervals", format(100 * level, trim = TRUE))
  } else {
    "Intervals of the estimates"
  }
}

# Under the title: what all the rows share (group and, if `time`, sampling
# time) and the intervals.
plot_subtitle <- function(d, opts, time = TRUE, extra = NULL) {
  parts <- character()
  context <- unique(d$context)
  if (length(context) == 1L && nzchar(context)) parts <- plot_context_labels(context)
  times <- unique(d$time)
  if (time && length(times) == 1L && !is.na(times)) {
    parts <- c(parts, plot_time_labels(times, opts))
  }
  paste(c(parts, plot_interval_note(opts$level), extra), collapse = " \u00b7 ")
}

plot_title <- function(opts, default, page, n_pages) {
  if (!is.null(opts$main)) return(opts$main)
  if (n_pages > 1L) default <- sprintf("%s \u00b7 %d/%d", default, page, n_pages)
  default
}

# Classes -----------------------------------------------------------------------

plot_codes <- c("+S" = "S+", "-S" = "S-", "+A" = "A+", "-A" = "A-", AD = "n.d.")

# Code of the class of each interaction row (see `classes` in the help page).
plot_class_codes <- function(x, d, alpha) {
  code <- rep("NC", nrow(d))
  code[d$estimand == "pure" | d$null_model == "DO"] <- "n/a"
  ask <- code == "NC"
  if (any(ask)) {
    cl <- classify_interactions(x, alpha = alpha)
    key <- function(t) {
      paste_key(t$response, t$context, t$time, t$set_id, t$null_model, t$estimand)
    }
    hit <- match(key(d), key(cl))
    found <- cl$class[hit]
    code[ask & !is.na(hit) & cl$status[hit] %in% "not_applicable"] <- "n/a"
    ok <- ask & !is.na(found) & found %in% names(plot_codes)
    code[ok] <- unname(plot_codes[found[ok]])
  }
  code
}

plot_class_names <- function(code) {
  names <- c("S+" = "Positive\nsynergism (S+)", "S-" = "Negative\nsynergism (S-)",
             "A+" = "Positive\nantagonism (A+)", "A-" = "Negative\nantagonism (A-)",
             "n.d." = "No deviation\ndetected (n.d.)", "n/a" = "Not\napplicable (n/a)",
             NC = "Not\nclassified (NC)")
  out <- unname(names[code])
  out[is.na(out)] <- ""
  out
}

plot_class_mode <- function(opts, n) {
  switch(opts$classes, auto = if (n == 1L) "names" else "codes", opts$classes)
}

# Marks and axes ----------------------------------------------------------------

plot_palette <- function() {
  list(
    ink = "#2b2b2a", soft = "#6b6a66", faint = "#9b9a96",
    interaction = "#2a78d6", combined = "#202020",
    models = c(AD = "#2a78d6", MU = "#eb6834", DO = "#1baf7a", BL = "#eda100",
               SAM = "#e87ba4", CA = "#008300"),
    model_pch = c(AD = 21, MU = 22, DO = 24, BL = 23, SAM = 25, CA = 21),
    singles = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"),
    single_pch = c(21, 22, 24, 25, 23)
  )
}

# Limits that cover zero, the finite estimates and their intervals; with a
# finite `extension`, intervals are cut at that many times the range of the
# estimates beyond it.
plot_limits <- function(value, lower, upper, pad, extension = Inf) {
  finite <- function(v) v[is.finite(v)]
  base <- range(c(0, finite(value)))
  if (diff(base) == 0) base <- base + c(-0.5, 0.5)
  full <- range(c(base, finite(lower), finite(upper)))
  if (is.finite(extension)) {
    width <- diff(base)
    full <- c(max(full[[1L]], base[[1L]] - extension * width),
              min(full[[2L]], base[[2L]] + extension * width))
  }
  grDevices::extendrange(full, f = pad)
}

# Interval from `lower` to `upper` at `at` (a y for horizontal intervals, an
# x for vertical ones), cut at the limits with an arrow at each cut end.
draw_interval <- function(lower, upper, at, lim, col, horizontal = TRUE, lwd = 1.5,
                          length = 0.06) {
  if (!is.finite(lower) || !is.finite(upper) || lower > upper) return(invisible())
  tip <- 0.035 * diff(lim)
  a <- max(lower, lim[[1L]])
  b <- min(upper, lim[[2L]])
  line <- function(from, to, arrow) {
    xy <- if (horizontal) c(from, at, to, at) else c(at, from, at, to)
    if (arrow) {
      graphics::arrows(xy[[1L]], xy[[2L]], xy[[3L]], xy[[4L]], length = length, angle = 25,
                       col = col, lwd = lwd)
    } else {
      graphics::segments(xy[[1L]], xy[[2L]], xy[[3L]], xy[[4L]], col = col, lwd = lwd)
    }
  }
  if (a < b) line(a, b, FALSE)
  if (lower < lim[[1L]]) line(lim[[1L]] + tip, lim[[1L]], TRUE)
  if (upper > lim[[2L]]) line(lim[[2L]] - tip, lim[[2L]], TRUE)
  invisible()
}

# An estimate outside the limits (or infinite): an arrow at the edge,
# pointing outwards.
draw_outside <- function(value, at, lim, col, horizontal = TRUE, lwd = 1.5, length = 0.06) {
  low <- value < lim[[1L]]
  edge <- if (low) lim[[1L]] else lim[[2L]]
  from <- edge + (if (low) 1 else -1) * 0.035 * diff(lim)
  if (horizontal) {
    graphics::arrows(from, at, edge, at, length = length, angle = 25, col = col, lwd = lwd)
  } else {
    graphics::arrows(at, from, at, edge, length = length, angle = 25, col = col, lwd = lwd)
  }
  invisible()
}

plot_pages <- function(block, per_page) {
  blocks <- unique(block)
  page <- if (is.null(per_page)) {
    rep(1L, length(blocks))
  } else {
    as.integer(ceiling(seq_along(blocks) / per_page))
  }
  page[match(block, blocks)]
}

# Ask before each new page in an interactive session, as plot.lm() does;
# returns the previous setting for the caller to restore.
ask_new_pages <- function(n_pages) {
  if (n_pages > 1L && grDevices::dev.interactive(orNone = TRUE)) {
    grDevices::devAskNewPage(TRUE)
  } else {
    NULL
  }
}

restore_ask <- function(old) {
  if (!is.null(old)) grDevices::devAskNewPage(old)
  invisible()
}

plot_result <- function(rows) {
  rows[c("response", "context", "time", "set_id", "null_model", "estimand", "component",
         "stressor", "estimate", "lower", "upper", "p_adjusted", "status", "class_code",
         "page")]
}

# type = "effects" ------------------------------------------------------------

plot_effects <- function(x, opts) {
  d <- plot_interaction_rows(x, opts)
  d$class_code <- plot_class_codes(x, d, opts$alpha)
  est <- x$estimates
  header <- plot_set_labels(d, opts)
  if (length(unique(d$context)) > 1L) {
    header <- paste0(plot_context_labels(d$context), " \u00b7 ", header)
  }
  if (length(unique(d$time)) > 1L) {
    header <- paste0(header, " \u00b7 ", plot_time_labels(d$time, opts))
  }
  net <- opts$estimand == "net"
  blocks <- vector("list", nrow(d))
  for (i in seq_len(nrow(d))) {
    singles <- plot_same_set(est, d, i, "single")
    singles <- singles[!duplicated(singles$stressor), , drop = FALSE]
    combined <- plot_same_set(est, d, i, "combined")
    if (nrow(singles) == 0L || nrow(combined) == 0L) {
      abort_input(sprintf(
        "Set %s has no effects of its stressors alone to plot: use type = \"interactions\".",
        d$set_id[[i]]
      ), call = opts$call)
    }
    combined <- combined[1L, , drop = FALSE]
    in_set <- strsplit(d$stressors[[i]], " x ", fixed = TRUE)[[1L]]
    singles <- singles[order(match(singles$stressor, in_set)), , drop = FALSE]
    n <- nrow(singles)
    expected <- combined$estimate - d$estimate[[i]]
    expected_status <- if (is.na(expected)) "not_estimable" else "ok"
    block <- data.frame(
      component = c(rep("single", n), if (net) "expected", "combined", "interaction"),
      stressor = c(singles$stressor, rep(NA_character_, net + 2L)),
      label = c(paste(plot_stressor_names(singles$stressor, opts), "alone"),
                if (net) "Combined, expected (no interval)", "Combined, observed",
                if (net) "Net interaction" else "Pure interaction"),
      estimate = c(singles$estimate, if (net) expected, combined$estimate, d$estimate[[i]]),
      lower = c(singles$lower, if (net) NA_real_, combined$lower, d$lower[[i]]),
      upper = c(singles$upper, if (net) NA_real_, combined$upper, d$upper[[i]]),
      p_adjusted = c(singles$p_adjusted, if (net) NA_real_, combined$p_adjusted,
                     d$p_adjusted[[i]]),
      status = c(singles$status, if (net) expected_status, combined$status, d$status[[i]]),
      class_code = c(rep("", n + net + 1L), d$class_code[[i]]),
      header = header[[i]], block = i, stringsAsFactors = FALSE
    )
    keys <- d[rep(i, nrow(block)), c("response", "context", "time", "set_id", "null_model",
                                     "estimand"), drop = FALSE]
    blocks[[i]] <- cbind(keys, block, row.names = NULL)
  }
  rows <- do.call(rbind, blocks)
  pal <- plot_palette()
  rows$color <- ifelse(rows$component == "interaction", pal$interaction,
                       ifelse(rows$component == "combined", pal$ink, pal$soft))
  xlab <- sprintf("Effects relative to the reference / interaction\n%s (%s)", opts$null_model,
                  paste(unique(d$scale), collapse = ", "))
  rows$page <- draw_margin_forest(rows, d, opts, xlab = xlab,
                                  class_mode = plot_class_mode(opts, nrow(d)), most = 4L)
  plot_result(rows)
}

# type = "interactions" -------------------------------------------------------

plot_interactions <- function(x, opts) {
  d <- plot_interaction_rows(x, opts)
  d$class_code <- plot_class_codes(x, d, opts$alpha)
  several_contexts <- length(unique(d$context)) > 1L
  several_times <- length(unique(d$time)) > 1L
  group <- opts$group_by
  if (group == "time" && !several_times) group <- "none"
  if (group == "context" && !several_contexts) group <- "none"
  if (group == "time") {
    d <- d[order(match(d$time, plot_time_order(d$time)), seq_len(nrow(d))), , drop = FALSE]
  }
  set_label <- plot_set_labels(d, opts)
  context_label <- plot_context_labels(d$context)
  when <- plot_time_labels(d$time, opts)
  shown_time <- if (several_times) when else rep(NA_character_, nrow(d))
  join <- function(a, b) ifelse(is.na(b) | !nzchar(b), a, paste0(a, " \u00b7 ", b))
  in_context <- if (several_contexts) join(context_label, set_label) else set_label
  interaction_name <- paste(if (opts$estimand == "net") "Net" else "Pure", "interaction")
  label <- switch(group,
                  none = join(in_context, shown_time),
                  set = ifelse(is.na(shown_time), interaction_name, shown_time),
                  time = in_context,
                  context = join(set_label, shown_time))
  header <- switch(group,
                   none = rep("", nrow(d)),
                   set = in_context,
                   time = when,
                   context = context_label)
  block <- if (group == "none") seq_len(nrow(d)) else match(header, unique(header))
  relative <- opts$scale == "relative"
  rows <- data.frame(
    d[c("response", "context", "time", "set_id", "null_model", "estimand")],
    component = "interaction", stressor = NA_character_, label = label,
    estimate = if (relative) d$relative else d$estimate,
    lower = if (relative) d$relative_lower else d$lower,
    upper = if (relative) d$relative_upper else d$upper,
    p_adjusted = d$p_adjusted, status = d$status, class_code = d$class_code,
    header = header, block = block, color = plot_palette()$interaction,
    stringsAsFactors = FALSE
  )
  xlab <- if (relative) {
    "Relative deviation: (observed - expected) / reference"
  } else {
    sprintf("%s interaction (%s)", opts$null_model, paste(unique(d$scale), collapse = ", "))
  }
  rows$page <- draw_margin_forest(rows, d, opts, xlab = xlab,
                                  class_mode = plot_class_mode(opts, nrow(d)), most = Inf)
  plot_result(rows)
}

# Forest with the labels in the left margin, block headers in bold and the
# class in a column at the right (types "effects" and "interactions").
# Without `per_page`, each page takes as many blocks as fit in the device,
# up to `most`. Returns the page of each row.
draw_margin_forest <- function(rows, d, opts, xlab, class_mode, most) {
  pal <- plot_palette()
  cex <- 0.8
  text <- plot_wrap(rows$label, opts$label_width)
  head <- plot_wrap(rows$header, opts$label_width)
  class_text <- switch(class_mode, names = plot_class_names(rows$class_code),
                       codes = rows$class_code, none = rep("", nrow(rows)))
  class_text[!nzchar(rows$class_code)] <- ""
  show_class <- class_mode != "none"
  drawable <- rows$status %in% c("ok", "boundary") & !is.na(rows$estimate)
  boundary <- drawable & rows$status == "boundary"
  with_interval <- drawable & !boundary
  xlim <- if (is.null(opts$xlim)) {
    plot_limits(rows$estimate[drawable], rows$lower[with_interval], rows$upper[with_interval],
                pad = 0.06)
  } else {
    opts$xlim
  }
  withr::local_par(list(cex = 1))
  left <- max(c(graphics::strwidth(text, units = "inches", cex = cex),
                graphics::strwidth(head[nzchar(head)], units = "inches", cex = cex, font = 2),
                0)) + 0.25
  left <- max(0.35, min(left, 0.5 * graphics::par("fin")[[1L]]))
  right <- if (show_class) {
    max(c(graphics::strwidth(class_text, units = "inches", cex = cex),
          graphics::strwidth("Class", units = "inches", cex = cex, font = 2))) + 0.25
  } else {
    0.2
  }
  bottom <- if (grepl("\n", xlab, fixed = TRUE)) 0.9 else 0.75
  withr::local_par(list(mai = c(bottom, left, 0.95, right), mgp = c(2.5, 0.7, 0),
                        fg = pal$ink, col.axis = pal$soft, col.lab = pal$ink, cex.axis = 0.85))
  start <- c(TRUE, rows$block[-1L] != rows$block[-nrow(rows)])
  n_text <- plot_lines(text)
  n_head <- ifelse(start & nzchar(head), plot_lines(head), 0L)
  line <- cex * graphics::par("cin")[[2L]]
  height <- graphics::par("pin")[[2L]]
  page <- if (is.null(opts$per_page)) {
    fit_pages(rows$block, n_text, n_head, height, line, most)
  } else {
    plot_pages(rows$block, opts$per_page)
  }
  pages <- unique(page)
  old_ask <- ask_new_pages(length(pages))
  on.exit(restore_ask(old_ask), add = TRUE)
  subtitle <- plot_subtitle(d, opts)
  default_title <- sprintf("%s \u00b7 %s \u00b7 %s", opts$response, opts$null_model,
                           opts$estimand)
  for (pg in pages) {
    k <- which(page == pg)
    r <- rows[k, , drop = FALSE]
    has_head <- n_head[k] > 0L
    pos <- forest_layout(n_text[k], n_head[k], height, line)
    y <- pos$y
    graphics::plot(NA, xlim = xlim, ylim = c(0, pos$top), xaxs = "i", yaxs = "i",
                   yaxt = "n", ylab = "", xlab = xlab, bty = "n")
    graphics::abline(v = 0, lty = 2, col = pal$soft)
    for (i in seq_along(k)) {
      j <- k[[i]]
      if (!drawable[[j]]) {
        graphics::text(mean(xlim), y[[i]], "NE", col = pal$soft, cex = cex)
        next
      }
      col <- r$color[[i]]
      if (with_interval[[j]]) draw_interval(r$lower[[i]], r$upper[[i]], y[[i]], xlim, col)
      v <- r$estimate[[i]]
      if (v < xlim[[1L]] || v > xlim[[2L]]) {
        draw_outside(v, y[[i]], xlim, col)
        if (boundary[[j]]) {
          # The triangle of the tests from only zeros, next to the arrow
          inside <- xlim[[if (v < xlim[[1L]]) 1L else 2L]] +
            (if (v < xlim[[1L]]) 1 else -1) * 0.05 * diff(xlim)
          graphics::points(inside, y[[i]], pch = 24, col = col, bg = col, cex = 1.1)
        }
        next
      }
      p <- r$p_adjusted[[i]]
      expected <- r$component[[i]] == "expected"
      pch <- if (boundary[[j]]) 24 else if (expected || is.na(p)) 23 else 21
      filled <- boundary[[j]] || (!expected && !is.na(p) && p < opts$alpha)
      graphics::points(v, y[[i]], pch = pch, col = col, bg = if (filled) col else "white",
                       cex = 1.1)
    }
    graphics::mtext(text[k], side = 2, at = y, las = 1, line = 0.3, cex = cex, col = pal$ink)
    if (any(has_head)) {
      graphics::mtext(head[k][has_head], side = 2, at = pos$y_head[has_head], las = 1,
                      line = 0.3, cex = cex, font = 2, col = pal$ink)
    }
    if (show_class) {
      x_class <- graphics::grconvertX(graphics::grconvertX(1, from = "npc", to = "inches") + 0.12,
                                      from = "inches", to = "user")
      graphics::text(rep(x_class, length(k)), y, labels = class_text[k], adj = c(0, 0.5),
                     xpd = NA, cex = cex,
                     col = ifelse(r$class_code %in% c("NC", "n/a"), pal$soft, pal$ink))
      graphics::text(x_class, graphics::par("usr")[[4L]], "Class", adj = c(0, 0), xpd = NA,
                     cex = cex, font = 2, col = pal$ink)
    }
    graphics::mtext(subtitle, side = 3, line = 0.3, cex = 0.7, col = pal$soft)
    legend_top <- draw_forest_legend(r, drawable[k], boundary[k], opts$alpha, xlim)
    # The title goes above the legend, whatever its number of rows.
    graphics::title(main = plot_title(opts, default_title, pg, length(pages)), cex.main = 1,
                    line = (legend_top + 0.1) / graphics::par("csi"))
  }
  page
}

# Vertical positions (in rows) of the rows of a page and of the headers above
# the first row of each block, from the number of lines of each label
# (`n_text`) and header (`n_head`, 0 without header), the height of the plot
# and that of a line of text (inches). Rows are one unit apart and blocks
# 0.9 units (more with longer headers), unless the text needs more room: the
# inches per unit are then solved so that no label or header overlaps.
forest_layout <- function(n_text, n_head, height, line) {
  base_step <- pmax(1, 0.55 * n_text)
  base_gap <- ifelse(n_head > 0L, 0.9 + 0.4 * pmax(0, n_head - 1L), 0)
  text_step <- forest_text_step(n_text, line)
  need <- forest_head_need(n_text, n_head, line)
  # Inches taken by the page with `inch` inches per unit
  inches <- function(inch) {
    step <- pmax(base_step * inch, text_step)
    sum(pmax(base_gap * inch + step, need)) + inch
  }
  if (inches(1e-6) < height) {
    inch <- stats::uniroot(function(u) inches(u) - height, c(1e-6, height), tol = 1e-8)$root
    step <- pmax(base_step, text_step / inch)
    gap <- pmax(base_gap + step, need / inch) - step
  } else {
    # Too many rows for the device: the base spacing, which may overlap.
    step <- base_step
    gap <- base_gap
    inch <- height / (sum(gap + step) + 1)
  }
  cursor <- cumsum(gap + step)
  total <- max(cursor) + 0.8
  y <- total - cursor
  list(y = y, y_head = y + ((n_head + n_text) * line / 2 + 0.03) / inch, top = total + 0.2)
}

# Inches that each label needs, and that each row with a header needs from
# the row above it to its own (header and half of both labels).
forest_text_step <- function(n_text, line) n_text * line + 0.06

forest_head_need <- function(n_text, n_head, line) {
  previous <- c(0, n_text[-length(n_text)])
  ifelse(n_head > 0L, n_head * line + (n_text + previous) * line / 2 + 0.15, 0)
}

# Pages with as many consecutive blocks as fit in `height` inches (with 10%
# to spare), up to `most`; a block that does not fit alone gets a page.
fit_pages <- function(block, n_text, n_head, height, line, most) {
  blocks <- unique(block)
  page <- integer(length(blocks))
  current <- 1L
  first <- 1L
  for (b in seq_along(blocks)) {
    rows <- block %in% blocks[first:b]
    needed <- sum(pmax(forest_text_step(n_text[rows], line),
                       forest_head_need(n_text[rows], n_head[rows], line)))
    if (b > first && (b - first + 1L > most || 1.1 * needed > height)) {
      current <- current + 1L
      first <- b
    }
    page[b] <- current
  }
  page[match(block, blocks)]
}

# Only the kinds of point that the page shows.
draw_forest_legend <- function(r, drawable, boundary, alpha, xlim) {
  pal <- plot_palette()
  tested <- drawable & !boundary & r$component != "expected"
  sig <- !is.na(r$p_adjusted) & r$p_adjusted < alpha
  entries <- list(
    list(any(tested & sig), sprintf("adjusted p < %s", format(alpha)), 21, pal$interaction),
    list(any(tested & !sig & !is.na(r$p_adjusted)), sprintf("adjusted p >= %s", format(alpha)),
         21, "white"),
    list(any(tested & is.na(r$p_adjusted)), "p not available", 23, "white"),
    list(any(boundary), "Only zeros: no interval", 24, pal$interaction),
    list(any(!drawable), "NE: no estimate", NA_real_, NA_character_)
  )
  entries <- Filter(function(e) e[[1L]], entries)
  # Height above the plot, in inches, of the top of the legend.
  bottom <- 0.22
  if (length(entries) == 0L) return(bottom)
  n <- length(entries)
  top <- graphics::grconvertY(1, from = "npc", to = "inches")
  y <- graphics::grconvertY(top + bottom, from = "inches", to = "user")
  drawn <- graphics::legend(x = xlim[[2L]], y = y, xjust = 1, yjust = 0, xpd = NA, bty = "n",
                            legend = vapply(entries, `[[`, character(1), 2L),
                            pch = vapply(entries, `[[`, numeric(1), 3L),
                            pt.bg = vapply(entries, `[[`, character(1), 4L),
                            col = pal$interaction, cex = 0.7, text.col = pal$ink,
                            ncol = if (n <= 2L) n else if (n <= 4L) 2L else 3L)
  graphics::grconvertY(drawn$rect$top, from = "user", to = "inches") - top
}

# type = "models" -------------------------------------------------------------

plot_models <- function(x, opts) {
  d <- plot_interaction_rows(x, opts)
  d$class_code <- plot_class_codes(x, d, opts$alpha)
  pal <- plot_palette()
  cex <- 0.8
  show_class <- opts$classes != "none"
  set_label <- plot_set_labels(d, opts)
  several_contexts <- length(unique(d$context)) > 1L
  block <- paste_key(d$context, d$set_id, d$time)
  d$page <- plot_pages(block, if (is.null(opts$per_page)) 4L else opts$per_page)
  xlim <- if (is.null(opts$xlim)) {
    plot_limits(d$relative, d$relative_lower, d$relative_upper, pad = 0.04, extension = 1)
  } else {
    opts$xlim
  }
  pages <- sort(unique(d$page))
  old_ask <- ask_new_pages(length(pages))
  on.exit(restore_ask(old_ask), add = TRUE)
  withr::local_par(list(cex = 1, mai = c(1.25, 1.05, 0.85, 0.2), mgp = c(2.2, 0.6, 0),
                        fg = pal$ink, col.axis = pal$soft, col.lab = pal$ink, cex.axis = 0.85))
  width <- graphics::par("pin")[[1L]]
  # Headers are written inside the plot, wrapped to its width.
  wrap_to_plot <- function(s) {
    wrap_label(s, function(l) {
      graphics::strwidth(l, units = "inches", cex = cex, font = 2) <= 0.96 * width
    })
  }
  # A single sampling time goes under the title; several, above each block.
  several_times <- length(unique(d$time)) > 1L
  subtitle <- plot_subtitle(d, opts, extra = sprintf("filled: adjusted p < %s", format(opts$alpha)))
  default_title <- sprintf("%s \u00b7 %s interaction", opts$response, opts$estimand)
  for (pg in pages) {
    y_row <- rep(NA_real_, nrow(d))
    heads <- list()
    bands <- list()
    cursor <- 0
    last_set <- NULL
    last_context <- NULL
    add_head <- function(text, font, space) {
      heads[[length(heads) + 1L]] <<- list(y = -cursor, text = text, font = font)
      cursor <<- cursor + space
    }
    for (b in unique(block[d$page == pg])) {
      i <- which(block == b)
      first <- i[[1L]]
      if (several_contexts && !identical(last_context, d$context[[first]])) {
        text <- wrap_to_plot(plot_context_labels(d$context[[first]]))
        add_head(text, 2, 0.9 * plot_lines(text) + 0.2)
        last_set <- NULL
      }
      set_key <- paste_key(d$context[[first]], d$set_id[[first]])
      if (!identical(last_set, set_key)) {
        text <- wrap_to_plot(set_label[[first]])
        add_head(text, 2, 0.9 * plot_lines(text) + 0.4)
      }
      if (several_times) add_head(plot_time_labels(d$time[[first]], opts), 1, 1)
      y_row[i] <- -cursor - seq_along(i) + 0.25
      bands[[length(bands) + 1L]] <- range(y_row[i]) + c(-0.45, 0.45)
      cursor <- cursor + length(i) + 0.9
      last_set <- set_key
      last_context <- d$context[[first]]
    }
    graphics::plot(NA, xlim = xlim, ylim = c(-cursor + 0.25, 0.8), xaxs = "i", yaxs = "i",
                   xaxt = "n", yaxt = "n", ylab = "", bty = "n", cex.main = 1,
                   xlab = "Deviation from the null model (% of the reference)",
                   main = plot_title(opts, default_title, pg, length(pages)))
    ticks <- pretty(xlim)
    ticks <- ticks[ticks >= xlim[[1L]] & ticks <= xlim[[2L]]]
    graphics::axis(1, at = ticks, labels = paste0(format(signif(100 * ticks, 5), trim = TRUE), "%"))
    from_left <- function(inches) {
      graphics::grconvertX(graphics::grconvertX(xlim[[1L]], from = "user", to = "inches") + inches,
                           from = "inches", to = "user")
    }
    x_model <- from_left(-0.79)
    x_class <- from_left(-0.25)
    graphics::text(c(x_model, if (show_class) x_class), 0.65, c("Model", if (show_class) "Class"),
                   xpd = NA, cex = cex * 0.9, font = 2, col = pal$ink)
    for (band in bands) graphics::segments(0, band[[1L]], 0, band[[2L]], col = pal$faint, lty = 2)
    for (h in heads) {
      graphics::text(xlim[[1L]], h$y, h$text, adj = c(0, 0.5), cex = cex, font = h$font,
                     col = pal$ink)
    }
    for (i in which(d$page == pg)) {
      nm <- d$null_model[[i]]
      col <- pal$models[[nm]]
      yi <- y_row[[i]]
      graphics::text(x_model, yi, nm, xpd = NA, cex = cex, col = col, font = 2)
      if (show_class) {
        graphics::text(x_class, yi, d$class_code[[i]], xpd = NA, cex = cex,
                       col = if (d$class_code[[i]] %in% c("NC", "n/a")) pal$soft else pal$ink)
      }
      v <- d$relative[[i]]
      if (is.na(v)) {
        graphics::text(mean(xlim), yi, "NE", cex = cex, col = pal$soft)
        next
      }
      draw_interval(d$relative_lower[[i]], d$relative_upper[[i]], yi, xlim, col, lwd = 1.4,
                    length = 0.055)
      if (v < xlim[[1L]] || v > xlim[[2L]]) {
        draw_outside(v, yi, xlim, col, lwd = 1.2, length = 0.055)
        next
      }
      p <- d$p_adjusted[[i]]
      graphics::points(v, yi, pch = if (is.na(p)) 4 else pal$model_pch[[nm]], col = col,
                       bg = if (!is.na(p) && p < opts$alpha) col else "white", cex = 1.05)
    }
    graphics::mtext(subtitle, side = 3, line = 0.25, cex = 0.7, col = pal$soft)
    if (show_class) {
      graphics::mtext(paste("S: synergism; A: antagonism; n.d.: no deviation detected",
                            "(classes of classify_interactions())."),
                      side = 1, line = 3.3, cex = 0.63, col = pal$soft)
    }
    graphics::mtext(paste0(if (show_class) "n/a: class not applicable; NC: not classified; ",
                           "NE: no deviation available; x: p not available."),
                    side = 1, line = if (show_class) 4.2 else 3.3, cex = 0.6, col = pal$soft)
  }
  d$component <- "interaction"
  d$stressor <- NA_character_
  d$estimate <- d$relative
  d$lower <- d$relative_lower
  d$upper <- d$relative_upper
  plot_result(d)
}

# type = "time" ---------------------------------------------------------------

plot_time <- function(x, opts) {
  est <- x$estimates
  numeric_time <- suppressWarnings(as.numeric(est$time))
  keep <- est$response == opts$response & est$null_model == opts$null_model &
    (est$component %in% c("single", "combined") |
       (est$component == "interaction" & est$estimand %in% "net")) &
    is.finite(numeric_time) & plot_filter(est, opts)
  d <- est[keep, , drop = FALSE]
  if (!any(d$component == "interaction")) {
    abort_input(if (all(is.na(est$time))) {
      "type = \"time\" needs a design with sampling times (`time` in describe_design())."
    } else {
      "No sampling time (other than \"overall\") of the selection has net interactions to plot."
    }, call = opts$call)
  }
  d <- d[!duplicated(paste_key(d$context, d$set_id, d$time, d$component, d$stressor)), ,
         drop = FALSE]
  pal <- plot_palette()
  show_class <- opts$classes != "none"
  interactions <- d[d$component == "interaction", , drop = FALSE]
  interactions <- interactions[plot_order(interactions, opts$null_model), , drop = FALSE]
  interactions$class_code <- plot_class_codes(x, interactions, opts$alpha)
  set_key <- paste_key(d$context, d$set_id)
  interaction_set <- paste_key(interactions$context, interactions$set_id)
  sets <- unique(interaction_set)
  old_ask <- ask_new_pages(length(sets))
  on.exit(restore_ask(old_ask), add = TRUE)
  withr::local_par(list(mfrow = c(1L, 1L), cex = 1, las = 1, mar = graphics::par("mar"),
                        oma = c(2.8, 0.2, 3.8, 0.2), mgp = c(2.7, 0.65, 0),
                        col.axis = pal$soft, fg = pal$ink, col.lab = pal$ink))
  footnote <- sprintf("%s; filled: adjusted p < %s; x: p not available.",
                      plot_interval_note(opts$level), format(opts$alpha))
  out <- vector("list", length(sets))
  for (s in seq_along(sets)) {
    z <- d[set_key == sets[[s]], , drop = FALSE]
    int <- interactions[interaction_set == sets[[s]], , drop = FALSE]
    times <- sort(unique(as.numeric(int$time)))
    take <- function(component, stressor = NULL) {
      q <- z[z$component == component, , drop = FALSE]
      if (!is.null(stressor)) q <- q[q$stressor %in% stressor, , drop = FALSE]
      q[match(times, as.numeric(q$time)), , drop = FALSE]
    }
    first <- int[1L, , drop = FALSE]
    in_set <- strsplit(first$stressors, " x ", fixed = TRUE)[[1L]]
    singles <- lapply(in_set, function(st) take("single", st))
    combined <- take("combined")
    interaction <- int[match(times, as.numeric(int$time)), , drop = FALSE]
    if (any(vapply(singles, function(q) all(is.na(q$estimate)), logical(1))) ||
        all(is.na(combined$estimate))) {
      abort_input(sprintf("Set %s has no effects of its stressors alone over time to plot.",
                          first$set_id), call = opts$call)
    }
    expected <- combined$estimate - interaction$estimate
    expected[is.nan(expected)] <- NA_real_
    lim_top <- if (is.null(opts$ylim)) {
      plot_limits(c(unlist(lapply(singles, `[[`, "estimate")), expected, combined$estimate),
                  c(unlist(lapply(singles, `[[`, "lower")), combined$lower),
                  c(unlist(lapply(singles, `[[`, "upper")), combined$upper), pad = 0.09)
    } else {
      opts$ylim
    }
    lim_bottom <- if (is.null(opts$ylim_interaction)) {
      plot_limits(interaction$estimate, interaction$lower, interaction$upper, pad = 0.09)
    } else {
      opts$ylim_interaction
    }
    xlim <- if (length(times) == 1L) {
      times + c(-0.5, 0.5)
    } else {
      grDevices::extendrange(times, f = 0.07)
    }
    title <- plot_set_labels(first, opts)
    if (nzchar(first$context)) {
      title <- paste0(plot_context_labels(first$context), " \u00b7 ", title)
    }
    title <- wrap_label(title, function(l) nchar(l) <= 85)
    # Room above the panels for the title and the set, however many lines.
    extra <- 0.9 * (plot_lines(title) - 1L)
    # Each stressor keeps its colour and symbol in every figure.
    position <- match(in_set, opts$stressors)
    colors <- pal$singles[(position - 1L) %% length(pal$singles) + 1L]
    pch <- pal$single_pch[(position - 1L) %% length(pal$single_pch) + 1L]
    # Both panels share the left margin, wide enough for their tick labels,
    # so that their times line up.
    top_axis <- plot_axis_ticks(lim_top, 5L)
    bottom_axis <- plot_axis_ticks(lim_bottom, 4L)
    tick_width <- max(graphics::strwidth(c(top_axis$labels, bottom_axis$labels),
                                         units = "inches", cex = 0.85))
    ylab_line <- 0.65 + tick_width / graphics::par("csi") + 0.6
    left <- ylab_line + 2.2
    graphics::layout(matrix(1:2, ncol = 1), heights = c(1.15, 1))
    # layout() shrinks the text of two panels; keep it at full size.
    graphics::par(cex = 1, mar = c(0.8, left, 4.5, 0.8), oma = c(2.8, 0.2, 3.8 + extra, 0.2))
    graphics::plot(NA, xlim = xlim, ylim = lim_top, xaxt = "n", yaxt = "n", xlab = "",
                   ylab = "", bty = "n", yaxs = "i")
    graphics::axis(2, at = top_axis$at, labels = top_axis$labels, cex.axis = 0.85)
    graphics::title(ylab = paste0("Effect vs. reference\n", first$scale), line = ylab_line,
                    cex.lab = 0.9)
    graphics::abline(h = 0, lty = 2, col = pal$faint)
    for (j in seq_along(singles)) {
      draw_series(times, singles[[j]]$estimate, singles[[j]]$lower, singles[[j]]$upper,
                  singles[[j]]$p_adjusted, colors[[j]], pch[[j]], 1, lim_top, opts$alpha)
    }
    none <- rep(NA_real_, length(times))
    draw_series(times, expected, none, none, none, pal$soft, 23, 2, lim_top, opts$alpha,
                expected = TRUE)
    draw_series(times, combined$estimate, combined$lower, combined$upper, combined$p_adjusted,
                pal$combined, 21, 1, lim_top, opts$alpha)
    graphics::title("Single and combined effects", line = 3.3, cex.main = 0.85)
    y_legend <- graphics::grconvertY(graphics::grconvertY(1, from = "npc", to = "inches") + 0.1,
                                     from = "inches", to = "user")
    graphics::legend(x = mean(xlim), y = y_legend, xjust = 0.5, yjust = 0, xpd = NA, bty = "n",
                     cex = 0.72, text.col = pal$ink,
                     legend = c(paste(plot_stressor_names(in_set, opts), "alone"),
                                "Combined, expected (no interval)", "Combined, observed"),
                     col = c(colors, pal$soft, pal$combined), pch = c(pch, 23, 21),
                     pt.bg = c(colors, "white", pal$combined),
                     lty = c(rep(1, length(singles)), 2, 1),
                     ncol = min(3L, length(singles) + 2L))
    graphics::par(mar = c(5.5, left, 2.1, 0.8))
    graphics::plot(NA, xlim = xlim, ylim = lim_bottom, xaxt = "n", yaxt = "n", xlab = "",
                   ylab = "", bty = "n", yaxs = "i")
    graphics::axis(2, at = bottom_axis$at, labels = bottom_axis$labels, cex.axis = 0.85)
    graphics::title(ylab = paste0("Net interaction\n", first$scale), line = ylab_line,
                    cex.lab = 0.9)
    graphics::axis(1, at = times, labels = format(times, trim = TRUE), cex.axis = 0.85)
    graphics::abline(h = 0, lty = 2, col = pal$faint)
    draw_series(times, interaction$estimate, interaction$lower, interaction$upper,
                interaction$p_adjusted, pal$interaction, 21, 1, lim_bottom, opts$alpha)
    graphics::title("Deviation of the observed from the expected combined effect",
                    cex.main = 0.85, line = 0.6)
    codes <- interaction$class_code
    if (show_class) {
      graphics::mtext(codes, side = 1, at = times, line = 1.7, cex = 0.85, font = 2,
                      col = ifelse(codes %in% c("NC", "n/a"), pal$soft, pal$ink))
      graphics::mtext("Class:", side = 1, at = graphics::par("usr")[[1L]], line = 1.7, adj = 1,
                      cex = 0.68, col = pal$ink)
    }
    graphics::mtext(plot_time_name(opts), side = 1, line = if (show_class) 3.3 else 2.5,
                    cex = 0.85, col = pal$ink)
    main <- if (is.null(opts$main)) {
      sprintf("%s \u00b7 %s \u00b7 net interaction", opts$response, opts$null_model)
    } else {
      opts$main
    }
    graphics::mtext(main, outer = TRUE, side = 3, line = 2.1 + extra, font = 2, cex = 0.95,
                    col = pal$ink)
    graphics::mtext(title, outer = TRUE, side = 3, line = 0.5, cex = 0.77, col = pal$ink)
    graphics::mtext(footnote, outer = TRUE, side = 1, line = 0.65, cex = 0.63, col = pal$soft)
    if (show_class) {
      graphics::mtext(paste("Classes: S = synergism; A = antagonism; n.d. = no deviation",
                            "detected; NC = not classified."),
                      outer = TRUE, side = 1, line = 0.05, cex = 0.6, col = pal$soft)
    }
    expected_rows <- combined
    expected_rows$component <- "expected"
    expected_rows$estimate <- expected
    expected_rows$lower <- NA_real_
    expected_rows$upper <- NA_real_
    expected_rows$p_adjusted <- NA_real_
    columns <- setdiff(names(interaction), "class_code")
    drawn <- rbind(do.call(rbind, singles)[columns], expected_rows[columns], combined[columns],
                   interaction[columns])
    drawn <- drawn[!is.na(drawn$set_id), , drop = FALSE]
    drawn$estimand <- "net"
    drawn$class_code <- ifelse(drawn$component == "interaction",
                               codes[match(drawn$time, interaction$time)], "")
    drawn$page <- s
    out[[s]] <- plot_result(drawn)
  }
  do.call(rbind, out)
}

# Ticks of a vertical axis within `lim`, with labels as axis() would write
# them but never in scientific notation.
plot_axis_ticks <- function(lim, n) {
  at <- pretty(lim, n = n)
  at <- at[at >= lim[[1L]] & at <= lim[[2L]]]
  list(at = at, labels = format(at, trim = TRUE, scientific = FALSE, drop0trailing = TRUE))
}

# One series over time: line, intervals (cut with arrows) and points.
draw_series <- function(x, v, lower, upper, p, col, pch, lty, lim, alpha, expected = FALSE) {
  graphics::lines(x, v, col = col, lty = lty, lwd = 1.4)
  light <- grDevices::adjustcolor(col, alpha.f = 0.55)
  for (i in seq_along(x)) {
    if (is.na(v[[i]])) next
    draw_interval(lower[[i]], upper[[i]], x[[i]], lim, light, horizontal = FALSE, lwd = 1,
                  length = 0.055)
    if (v[[i]] < lim[[1L]] || v[[i]] > lim[[2L]]) {
      draw_outside(v[[i]], x[[i]], lim, col, horizontal = FALSE, lwd = 1, length = 0.055)
      next
    }
    sig <- !is.na(p[[i]]) && p[[i]] < alpha
    graphics::points(x[[i]], v[[i]], pch = if (is.na(p[[i]]) && !expected) 4 else pch,
                     col = col, bg = if (sig && !expected) col else "white", cex = 1.05)
  }
  invisible()
}
