#' Describe what a multiple-stressor design can identify
#'
#' Inspects the treatment structure of an experiment before any model is
#' fitted: which combinations of stressor levels were observed and which are
#' missing, which interactions each group of cells identifies, which cells are
#' shared by several of those groups, how the experimental units are
#' replicated and whether they were measured repeatedly. [fit_responses()]
#' refuses designs with problems, and [estimate_interactions()] and
#' [classify_interactions()] only work on the contrast sets found here.
#'
#' @param data A data frame with one row per observation. Tibbles are
#'   accepted; grouped data frames are not.
#' @param stressors Character vector with the names of at least two stressor
#'   columns. A column may hold factors, text or numbers; every distinct value
#'   is treated as a level (a category, not a gradient).
#' @param control Named character vector (or list) with, for every stressor,
#'   the level that represents its absence (control or ambient), e.g.
#'   `c(temperature = "20", nitrate = "0")`. There is no default: which level
#'   means "no stressor" is a scientific decision.
#' @param unit Name of the column that identifies the experimental
#'   (replication) unit, or `NULL` (default) if every row is an independent
#'   unit. Required when `time` is given.
#' @param time Name of the column with the sampling time, or `NULL` (default)
#'   for a single sampling.
#' @param by Optional names of columns with other factors (population,
#'   species, experiment). Contrast sets are formed within each combination of
#'   their levels.
#' @param groups Optional names of columns that group experimental units:
#'   blocks, whole plots of a split-plot design, rooms, incubators, header
#'   tanks or batches. For each one the checks report which stressors were
#'   applied to whole groups and whether their effects can be separated from
#'   differences among groups, and [fit_responses()] adds a random intercept
#'   for it when it can be estimated. Group labels must be unique within each
#'   `by` group (number plots 1-32, not 1-8 within each block).
#' @param max_order Largest number of stressors in a contrast set. `NULL`
#'   (default) means all stressors.
#' @param baseline Optional values of `time` taken before the stressors were
#'   applied (pre-treatment samplings), e.g. `baseline = 0`. They are kept in
#'   the models fitted by [fit_responses()], where they inform the
#'   differences among units (with `time_structure = "additive"` they are
#'   coded as untreated, the control combination), but
#'   [estimate_interactions()] does not estimate or test interactions at
#'   them, so they do not enter the multiple-testing adjustment either. If
#'   the stressors were applied at different times, give a named list with
#'   the times before each one, e.g. `baseline = list(drying = c(-36, -1))`
#'   when drying started on day 0 and the other stressors were applied from
#'   the start: interactions that involve drying are then not estimated at
#'   those times, and those among the other stressors are (with
#'   `time_structure = "additive"`, drying is coded at its control there).
#'   Stressors left out of the list were applied before the first sampling.
#'   Each unit keeps the treatment it was assigned in the stressor columns.
#'   Requires `time`.
#' @param background `"control"` (default) or `"all"`. With `"control"`,
#'   interactions are estimated with the stressors outside each contrast set
#'   at their control (interactions "in the absence of the other
#'   stressors"). With `"all"`, sets are also formed with the other stressors
#'   held at each of their observed non-control levels (for example, the
#'   interaction between salinity and warming in the presence of nutrients);
#'   their reference cell is then that background combination instead of the
#'   control.
#' @param shared Optional named character vector that marks rows shared by
#'   every group of a `by` column, e.g. `c(regime = "both")`: the rows whose
#'   `regime` is `"both"` (usually the control, sometimes also the
#'   single-stressor treatments) belong to every other level of `regime`.
#'   Use it when one control serves several `by` groups (two exposure
#'   regimes, two orders of application, upstream and downstream ponds), so
#'   that it enters the analysis once instead of being copied into each
#'   group. Contrast sets are still formed within each `by` group, with the
#'   shared rows supplying the cells that the group lacks; sets of different
#'   groups that use them are correlated. Each name must be a `by` column and
#'   each value one of its levels.
#'
#' @return An object of class `must_design`, a list with:
#'   * `spec`: the arguments used, with the observed levels of each stressor
#'     (control first).
#'   * `cells`: data frame with one row per observed treatment combination
#'     within each `by` group: `context` (label of the `by` group, `""` without
#'     `by`), the `by` columns and the stressor columns (as text), `cell`
#'     (label), `n_stressed` (number of stressors not at their control),
#'     `n_units`, `n_obs` and, with `time`, `n_times`. Ordered by `context`,
#'     `n_stressed` and stressor levels.
#'   * `missing_cells`: combinations of the observed levels that are absent,
#'     within each `by` group (same columns as `cells` without the counts).
#'   * `sets`: one row per candidate contrast set, that is, per observed cell
#'     with two or more stressors away from control (and, with
#'     `background = "all"`, per subset of two or more of those stressors):
#'     `set_id`, `context`, `order` (number of stressors), `stressors`,
#'     `levels` (label of the combined cell), `background` (levels of the
#'     stressors outside the set that are not at control, `""` if none),
#'     `net_identifiable`, `pure_identifiable` (`NA` for two
#'     stressors, where net and pure coincide), `missing_for_net`,
#'     `missing_for_pure`, `min_units` (fewest units among the cells needed for
#'     the net deviation), `shared` (`TRUE` if the set uses a cell that
#'     another identifiable set also uses) and `confounded` (`TRUE` if the set
#'     involves a stressor applied to whole groups of a `groups` column and
#'     some of its combinations were given to a single group).
#'   * `groups`: one row per `groups` column: `group`, `n_groups`, `between`
#'     (stressors applied to whole groups), `structure` (`"blocks"`,
#'     `"whole_units"`, `"confounded"` or `"units"`, see Details),
#'     `groups_per_combination`, `units_per_group`, `cells_per_group` and
#'     `random` (whether [fit_responses()] adds a random intercept for it).
#'   * `set_cells`: for each set, the cells it uses, with the `context` in
#'     which each one was observed (used by [estimate_interactions()]).
#'   * `cell_use`: number of identifiable sets that use each cell.
#'   * `checks`: data frame with `check`, `severity` (`"problem"`,
#'     `"warning"` or `"info"`) and `message`.
#'
#' @details
#' A contrast set is the group of cells needed to estimate how a combination
#' of stressors deviates from a null model. For a combination of the
#' stressors in S at given levels (all other stressors at control):
#'
#' * the **net** deviation needs the control, each stressor of S alone and
#'   their combination. It compares the combined response with the
#'   combination of the single-stressor effects predicted by the null model.
#' * the **pure** interaction of order |S| needs all 2^|S| combinations of
#'   S. It is what remains once all lower-order interactions are removed.
#'
#' With two stressors both are the same. With three or more, the net
#' deviation of ABC includes the pairwise interactions AB, AC and BC, so it is
#' not a pure third-order interaction; designs that lack the pairwise cells
#' cannot separate them.
#'
#' Cells shared by several sets (typically the control in designs with several
#' doses, or the single-stressor cells of a factorial with several levels)
#' make the estimates of those sets correlated.
#'
#' **Rows shared by several `by` groups.** When a control (or another cell)
#' is common to several groups, copying it into each group would count its
#' units several times. With `shared`, those rows form their own group
#' (`context` such as `"regime=both"`), and each of the other groups uses its
#' cells when it lacks them: a set of `"regime=pulse"` takes its control from
#' `"regime=both"`. The shared rows supply cells, not sets: sets are formed
#' from the combinations observed within each group. A design with the
#' order of application as a `by` column (nutrients then salinity, salinity
#' then nutrients) and the control and single stressors shared estimates the
#' interaction for each order.
#'
#' **Groups of units.** A stressor that takes a single level within every
#' group of a `groups` column was applied to whole groups (a whole-plot
#' factor). Each column gets one of four structures:
#' * `"blocks"`: no stressor is applied to whole groups; treatment
#'   combinations are compared within groups (complete or incomplete blocks).
#' * `"whole_units"`: some stressors are applied to whole groups, with at
#'   least two groups per combination of their levels. For effects that
#'   involve them the replication is the number of groups, not of units, and
#'   the random intercept for the groups is needed for valid inference.
#' * `"confounded"`: some combination of the stressors applied to whole
#'   groups was given to a single group, so effects that involve them cannot
#'   be separated from differences among groups (for example, one incubator
#'   per temperature, or one nitrate level per block). The sets affected are
#'   flagged in `sets$confounded`; a random intercept does not remove the
#'   confounding and is not added when it cannot be estimated.
#' * `"units"`: each group holds a single unit, so the column repeats
#'   `unit` and is not used.
#'
#' @section Assumptions:
#' * The `control` level of each stressor represents its absence.
#' * With `background = "control"`, interactions are estimated with the
#'   stressors outside the set at control (interactions "in the absence of
#'   the other stressors"); an interaction can differ in the presence of
#'   another stressor, which `background = "all"` estimates.
#' * Each experimental unit receives a single treatment combination for the
#'   whole experiment. Treatments that change within a unit (sequential
#'   exposures) are reported as a problem, except stressors that start later
#'   than others, declared in `baseline` by stressor.
#' * Units and groups not given in `unit` and `groups` are assumed
#'   independent: the design cannot reveal a shared incubator or water supply
#'   that is not recorded in the data.
#' * Levels are categories. A gradient with many levels and few replicates is
#'   analysed level by level, not as a continuous predictor.
#'
#' @section Missing values and zeros:
#' Missing values in the stressor, `unit`, `time`, `by` or `groups` columns
#' are an error: no treatment combination is guessed. The response is not
#' used here.
#'
#' @references
#' Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
#' J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025).
#' Hypothesis-driven research on multiple stressors: An analytical framework
#' for stressor interactions. *Ecology and Evolution*, 15, e71959.
#' \doi{10.1002/ece3.71959}
#'
#' @examples
#' # Two stressors, two levels each, four replicate jars
#' jars <- expand.grid(
#'   salinity = c("ambient", "high"),
#'   warming = c("no", "yes"),
#'   replicate = 1:4,
#'   stringsAsFactors = FALSE
#' )
#' jars$jar <- seq_len(nrow(jars))
#' design <- describe_design(
#'   jars,
#'   stressors = c("salinity", "warming"),
#'   control = c(salinity = "ambient", warming = "no"),
#'   unit = "jar"
#' )
#' design
#' design$sets
#'
#' # A ray design: the mixture only at matching doses
#' ray <- data.frame(
#'   a = c(0, 2, 20, 0, 0, 2, 20),
#'   b = c(0, 0, 0, 2, 20, 2, 20)
#' )
#' ray <- ray[rep(seq_len(nrow(ray)), each = 3), ]
#' ray_design <- describe_design(ray, c("a", "b"), control = c(a = "0", b = "0"))
#' ray_design$sets
#' ray_design$missing_cells
#'
#' # Temperature applied to whole rooms, one room per temperature: its
#' # effects cannot be separated from differences between rooms
#' rooms <- expand.grid(temperature = c("15", "20"), nutrients = c("no", "yes"),
#'                      jar = 1:3, stringsAsFactors = FALSE)
#' rooms$room <- paste0("R", rooms$temperature)
#' rooms$id <- seq_len(nrow(rooms))
#' room_design <- describe_design(rooms, c("temperature", "nutrients"),
#'   control = c(temperature = "15", nutrients = "no"), unit = "id",
#'   groups = "room")
#' room_design$groups
#' room_design$sets$confounded
#' @export
describe_design <- function(data, stressors, control, unit = NULL, time = NULL,
                            by = NULL, groups = NULL, max_order = NULL, baseline = NULL,
                            background = "control", shared = NULL) {
  check_data_frame(data, allow_empty = FALSE)
  data <- as_plain_data_frame(data)
  spec <- validate_design_spec(data, stressors, control, unit, time, by, groups,
                               max_order, baseline, background, shared)
  keys <- design_keys(data, spec)
  cells <- summarise_cells(data, spec, keys)
  pools <- context_pools(cells, spec)
  missing_cells <- find_missing_cells(spec, cells, pools)
  sets <- find_contrast_sets(spec, cells, pools)
  group_info <- analyse_groups(data, spec, keys)
  sets$sets$confounded <- flag_confounded(sets$sets, sets$set_cells, cells, group_info)
  checks <- design_checks(spec, keys, cells, missing_cells, sets, group_info,
                          nested = nested_group_labels(data, spec, keys),
                          labels = label_checks(data, spec, cells, pools),
                          pools = pools)
  structure(
    list(
      spec = spec,
      cells = cells,
      missing_cells = missing_cells$table,
      sets = sets$sets,
      set_cells = sets$set_cells,
      cell_use = sets$cell_use,
      groups = group_table(group_info),
      checks = checks,
      n_rows = nrow(data)
    ),
    class = "must_design"
  )
}

# Names used for the columns of the output tables.
reserved_names <- function() {
  c("context", "cell", "n_stressed", "n_units", "n_obs", "n_times", "set_id",
    "order", "stressors", "levels", "response")
}

validate_design_spec <- function(data, stressors, control, unit, time, by,
                                 groups, max_order, baseline = NULL,
                                 background = "control", shared = NULL,
                                 call = sys.call(-1)) {
  if (!is.character(stressors) || length(stressors) < 2L || anyNA(stressors) ||
      anyDuplicated(stressors) > 0L) {
    abort_input(
      paste(
        "`stressors` must be a character vector with the names of at least",
        "two different columns."
      ),
      call = call
    )
  }
  check_columns(data, stats::setNames(stressors, rep("stressors", length(stressors))),
                call = call)
  if (!is.null(unit)) {
    check_string(unit, arg = "unit", call = call)
    check_columns(data, c(unit = unit), call = call)
  }
  if (!is.null(time)) {
    check_string(time, arg = "time", call = call)
    check_columns(data, c(time = time), call = call)
    if (is.null(unit)) {
      abort_input(
        paste(
          "`time` is given but `unit` is not. With repeated samplings the",
          "experimental unit is needed to model the dependence between",
          "samplings of the same unit: give its column as `unit`."
        ),
        call = call
      )
    }
  }
  if (!is.null(by)) {
    if (!is.character(by) || length(by) == 0L || anyNA(by) || anyDuplicated(by) > 0L) {
      abort_input("`by` must be NULL or a character vector of column names.", call = call)
    }
    check_columns(data, stats::setNames(by, rep("by", length(by))), call = call)
  }
  if (!is.null(groups)) {
    if (!is.character(groups) || length(groups) == 0L || anyNA(groups) ||
        anyDuplicated(groups) > 0L) {
      abort_input("`groups` must be NULL or a character vector of column names.",
                  call = call)
    }
    check_columns(data, stats::setNames(groups, rep("groups", length(groups))), call = call)
  }
  roles <- c(stressors, unit, time, by, groups)
  if (anyDuplicated(roles) > 0L) {
    abort_input(
      sprintf(
        "A column cannot play two roles (stressor, unit, time, by, groups): %s.",
        format_items(quote_items(unique(roles[duplicated(roles)])), max = Inf)
      ),
      call = call
    )
  }
  clash <- intersect(c(stressors, by), reserved_names())
  if (length(clash) > 0L) {
    abort_input(
      sprintf(
        "Rename the column(s) %s: these names are used in the output tables.",
        format_items(quote_items(clash), max = Inf)
      ),
      call = call
    )
  }
  for (col in roles) {
    na_rows <- which(is.na(data[[col]]))
    if (length(na_rows) > 0L) {
      abort_input(
        sprintf(
          paste(
            "Column '%s' has missing values (%s). Treatment, unit, time, `by`",
            "and `groups` columns cannot be missing: fix or remove those rows first."
          ),
          col, format_rows(na_rows)
        ),
        column = col,
        rows = na_rows,
        call = call
      )
    }
  }
  if (missing(control) || is.null(control)) {
    abort_input(
      paste(
        "`control` is required: give the level that means 'no stressor' for",
        "each stressor, e.g. control = c(temperature = \"20\", nitrate = \"0\")."
      ),
      call = call
    )
  }
  if (is.null(names(control)) || !all(stressors %in% names(control))) {
    abort_input(
      sprintf(
        "`control` must be named and have one element per stressor (%s).",
        format_items(quote_items(stressors), max = Inf)
      ),
      call = call
    )
  }
  control_chr <- character(length(stressors))
  names(control_chr) <- stressors
  levels_list <- list()
  for (s in stressors) {
    value <- control[[s]]
    if (length(value) != 1L || is.na(value)) {
      abort_input(sprintf("`control` must give a single level for '%s'.", s), call = call)
    }
    value <- as.character(value)
    observed <- unique(as.character(data[[s]]))
    if (!value %in% observed) {
      abort_input(
        sprintf(
          "The control level '%s' of '%s' is not in the data. Observed levels: %s.",
          value, s, format_items(quote_items(sort(observed)), max = 20L)
        ),
        call = call
      )
    }
    if (length(observed) < 2L) {
      abort_input(
        sprintf("Stressor '%s' only has its control level: it is not a stressor here.", s),
        call = call
      )
    }
    control_chr[[s]] <- value
    levels_list[[s]] <- stressor_levels(data[[s]], value)
  }
  if (is.null(max_order)) {
    max_order <- length(stressors)
  } else {
    check_number(max_order, lower = 2, upper = length(stressors), arg = "max_order",
                 call = call)
    if (max_order != round(max_order)) {
      abort_input("`max_order` must be a whole number.", call = call)
    }
    max_order <- as.integer(max_order)
  }
  check_choice(background, c("control", "all"), arg = "background", call = call)
  baseline <- check_baseline(baseline, data, time, stressors, call = call)
  shared <- check_shared(shared, data, by, call = call)
  list(
    stressors = stressors,
    control = control_chr,
    unit = unit,
    time = time,
    by = by,
    groups = groups,
    max_order = as.integer(max_order),
    levels = levels_list,
    # Times before every stressor was applied, and times before each one.
    baseline = baseline$common,
    baseline_by = baseline$by_stressor,
    background = background,
    shared = shared
  )
}

# `baseline` is a vector of times before all the stressors were applied, or
# a named list with the times before each stressor was applied. Returns the
# times common to every stressor (`common`, NULL if none) and the times of
# each stressor (`by_stressor`, NULL without baseline).
check_baseline <- function(baseline, data, time, stressors, call = sys.call(-1)) {
  if (is.null(baseline)) {
    return(list(common = NULL, by_stressor = NULL))
  }
  if (is.null(time)) {
    abort_input("`baseline` needs `time`: give the column with the sampling times.",
                call = call)
  }
  observed <- unique(as.character(data[[time]]))
  times_of <- function(x, what) {
    if (!(is.character(x) || is.numeric(x) || is.factor(x)) || length(x) == 0L || anyNA(x)) {
      abort_input(sprintf("%s must be a vector of sampling times.", what), call = call)
    }
    x <- unique(as.character(x))
    absent <- setdiff(x, observed)
    if (length(absent) > 0L) {
      abort_input(
        sprintf(
          "%s has time(s) that are not in column '%s': %s. Observed times: %s.",
          what, time, format_items(quote_items(absent), max = 5L),
          format_items(quote_items(observed), max = 10L)
        ),
        call = call
      )
    }
    x
  }
  if (is.list(baseline)) {
    named <- names(baseline)
    if (is.null(named) || any(!nzchar(named)) || anyDuplicated(named) > 0L) {
      abort_input(
        paste(
          "A list in `baseline` must be named by stressor, with the times before",
          "each one was applied, e.g. baseline = list(drying = c(-36, -1))."
        ),
        call = call
      )
    }
    unknown <- setdiff(named, stressors)
    if (length(unknown) > 0L) {
      abort_input(
        sprintf(
          "The names of `baseline` must be stressors (%s); %s %s not.",
          format_items(quote_items(stressors), max = Inf),
          format_items(quote_items(unknown), max = Inf),
          if (length(unknown) == 1L) "is" else "are"
        ),
        call = call
      )
    }
    by_stressor <- stats::setNames(rep(list(character()), length(stressors)), stressors)
    for (s in named) {
      by_stressor[[s]] <- times_of(baseline[[s]], sprintf("`baseline$%s`", s))
      if (all(observed %in% by_stressor[[s]])) {
        abort_input(
          sprintf("`baseline$%s` includes every sampling time: '%s' was never applied.", s, s),
          call = call
        )
      }
    }
    common <- Reduce(intersect, by_stressor)
  } else {
    common <- times_of(baseline, "`baseline`")
    if (all(observed %in% common)) {
      abort_input("`baseline` cannot include every sampling time.", call = call)
    }
    by_stressor <- stats::setNames(rep(list(common), length(stressors)), stressors)
  }
  list(common = if (length(common) > 0L) common else NULL, by_stressor = by_stressor)
}

# Times before each stressor was applied (NULL without baseline), also for
# designs made before `baseline` could be given by stressor.
baseline_times <- function(spec) {
  if (!is.null(spec$baseline_by)) {
    return(spec$baseline_by)
  }
  if (is.null(spec$baseline)) {
    return(NULL)
  }
  stats::setNames(rep(list(as.character(spec$baseline)), length(spec$stressors)), spec$stressors)
}

# `shared`: named character vector (or list) with, for some `by` columns, the
# level that marks rows shared by every other level of that column.
check_shared <- function(shared, data, by, call = sys.call(-1)) {
  if (is.null(shared)) {
    return(NULL)
  }
  if (is.null(by)) {
    abort_input(
      paste(
        "`shared` needs `by`: it marks the rows of a `by` column that belong to",
        "every group of that column, e.g. by = \"regime\", shared = c(regime = \"both\")."
      ),
      call = call
    )
  }
  if (is.list(shared)) {
    if (any(lengths(shared) != 1L)) {
      abort_input("`shared` must give a single level for each `by` column.", call = call)
    }
    shared <- unlist(shared)
  }
  if (!is.atomic(shared) || length(shared) == 0L || is.null(names(shared)) ||
      any(!nzchar(names(shared))) || anyNA(shared) || anyDuplicated(names(shared)) > 0L) {
    abort_input(
      paste(
        "`shared` must be a named character vector, e.g. c(regime = \"both\"): the",
        "names are `by` columns and the values the level that marks shared rows."
      ),
      call = call
    )
  }
  unknown <- setdiff(names(shared), by)
  if (length(unknown) > 0L) {
    abort_input(
      sprintf(
        "The names of `shared` must be `by` columns (%s); %s %s not.",
        format_items(quote_items(by), max = Inf), format_items(quote_items(unknown), max = Inf),
        if (length(unknown) == 1L) "is" else "are"
      ),
      call = call
    )
  }
  out <- stats::setNames(as.character(shared), names(shared))
  for (b in names(out)) {
    observed <- unique(as.character(data[[b]]))
    if (!out[[b]] %in% observed) {
      abort_input(
        sprintf(
          "The shared level '%s' of '%s' is not in the data. Observed levels: %s.",
          out[[b]], b, format_items(quote_items(sort(observed)), max = 20L)
        ),
        call = call
      )
    }
    if (length(observed) < 2L) {
      abort_input(
        sprintf("'%s' has no group other than the shared one ('%s').", b, out[[b]]),
        call = call
      )
    }
  }
  out
}

# Per-row labels and indicators used by every other step.
design_keys <- function(data, spec) {
  lev <- data.frame(
    lapply(spec$stressors, function(s) as.character(data[[s]])),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(lev) <- spec$stressors
  stressed <- do.call(cbind, lapply(spec$stressors, function(s) lev[[s]] != spec$control[[s]]))
  list(
    levels = lev,
    cell = make_cell_label(lev, spec$stressors),
    context = make_context_label(data, spec$by),
    by_values = if (length(spec$by) > 0L) {
      stats::setNames(lapply(spec$by, function(b) as.character(data[[b]])), spec$by)
    } else {
      NULL
    },
    n_stressed = as.integer(rowSums(stressed)),
    unit_id = if (is.null(spec$unit)) {
      as.character(seq_len(nrow(data)))
    } else {
      as.character(data[[spec$unit]])
    },
    time_value = if (is.null(spec$time)) rep(NA_character_, nrow(data)) else as.character(data[[spec$time]])
  )
}

summarise_cells <- function(data, spec, keys) {
  key <- paste_key(keys$context, keys$cell)
  idx <- split(seq_len(nrow(data)), factor(key, levels = unique(key)))
  first <- vapply(idx, function(i) i[[1L]], integer(1))
  cells <- data.frame(context = keys$context[first], stringsAsFactors = FALSE)
  for (b in spec$by) {
    cells[[b]] <- keys$by_values[[b]][first]
  }
  for (s in spec$stressors) {
    cells[[s]] <- keys$levels[[s]][first]
  }
  cells$cell <- keys$cell[first]
  cells$n_stressed <- keys$n_stressed[first]
  cells$n_units <- vapply(idx, function(i) length(unique(keys$unit_id[i])), integer(1))
  cells$n_obs <- lengths(idx, use.names = FALSE)
  if (!is.null(spec$time)) {
    cells$n_times <- vapply(idx, function(i) length(unique(keys$time_value[i])), integer(1))
  }
  sort_cells(cells, spec)
}

sort_cells <- function(cells, spec) {
  if (nrow(cells) == 0L) {
    return(cells)
  }
  keys <- c(
    list(cells$context, cells$n_stressed),
    lapply(spec$stressors, function(s) match(cells[[s]], spec$levels[[s]]))
  )
  cells <- cells[do.call(order, c(keys, list(method = "radix"))), , drop = FALSE]
  rownames(cells) <- NULL
  cells
}

# The `by` groups whose cells the sets of each group can use: the group
# itself and the groups of rows marked as `shared` that match it in the other
# `by` columns. `pools` is named by the groups that form sets (not shared);
# `shared` lists the groups of shared rows.
context_pools <- function(cells, spec) {
  contexts <- unique(cells$context)
  if (length(spec$shared) == 0L) {
    return(list(pools = stats::setNames(as.list(contexts), contexts), shared = character()))
  }
  values <- cells[!duplicated(cells$context), c("context", spec$by), drop = FALSE]
  is_shared <- Reduce(`|`, lapply(names(spec$shared), function(b) {
    values[[b]] == spec$shared[[b]]
  }))
  pools <- list()
  for (i in which(!is_shared)) {
    covers <- vapply(which(is_shared), function(j) {
      all(vapply(spec$by, function(b) {
        v <- values[[b]][[j]]
        v == values[[b]][[i]] || (b %in% names(spec$shared) && v == spec$shared[[b]])
      }, logical(1)))
    }, logical(1))
    pools[[values$context[[i]]]] <- c(values$context[[i]], values$context[is_shared][covers])
  }
  list(pools = pools, shared = values$context[is_shared])
}

# Groups in the pool of `ctx` (by position: without `by` the context label is
# "", which [[""]] does not match).
pool_of <- function(pools, ctx) {
  i <- match(ctx, names(pools$pools))
  if (is.na(i)) ctx else pools$pools[[i]]
}

# Rows of `cells` with the given cell labels in the pool of `ctx`: the group's
# own cell if it has one, otherwise the first shared group that has it
# (NA if none).
pool_match <- function(labels, ctx, cells, pools) {
  pos <- rep(NA_integer_, length(labels))
  for (p in pool_of(pools, ctx)) {
    hit <- which(cells$context == p)
    m <- match(labels, cells$cell[hit])
    fill <- is.na(pos) & !is.na(m)
    pos[fill] <- hit[m[fill]]
  }
  pos
}

# Combinations of the levels observed in each `by` group (with the rows shared
# with it) that the group lacks. Levels that a group never uses are not counted
# as missing there; label_checks() points them out.
find_missing_cells <- function(spec, cells, pools) {
  template <- cells[0, c("context", spec$by, spec$stressors, "cell", "n_stressed"), drop = FALSE]
  contexts <- names(pools$pools)
  group_levels <- lapply(contexts, function(ctx) {
    in_pool <- cells$context %in% pool_of(pools, ctx)
    lapply(stats::setNames(spec$stressors, spec$stressors), function(s) {
      spec$levels[[s]][spec$levels[[s]] %in% cells[[s]][in_pool]]
    })
  })
  sizes <- vapply(group_levels, function(lv) prod(lengths(lv)), numeric(1))
  n_full <- sum(sizes)
  if (max(sizes) > 1e5) {
    return(list(table = template, n_full = n_full, enumerated = FALSE))
  }
  out <- lapply(seq_along(contexts), function(k) {
    ctx <- contexts[[k]]
    grid <- expand.grid(group_levels[[k]], stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    names(grid) <- spec$stressors
    grid_label <- make_cell_label(grid, spec$stressors)
    # Shared rows are not a group of their own: what is missing is what the
    # group lacks once the shared cells are counted in it.
    absent <- is.na(pool_match(grid_label, ctx, cells, pools))
    if (!any(absent)) {
      return(NULL)
    }
    first <- which(cells$context == ctx)[[1L]]
    rows <- data.frame(context = rep(ctx, sum(absent)), stringsAsFactors = FALSE)
    for (b in spec$by) {
      rows[[b]] <- rep(cells[[b]][[first]], sum(absent))
    }
    for (s in spec$stressors) {
      rows[[s]] <- grid[[s]][absent]
    }
    rows$cell <- grid_label[absent]
    rows$n_stressed <- as.integer(rowSums(do.call(
      cbind, lapply(spec$stressors, function(s) grid[[s]][absent] != spec$control[[s]])
    )))
    rows
  })
  table <- bind_rows_base(out, template)
  list(table = sort_cells(table, spec), n_full = n_full, enumerated = TRUE)
}

# One candidate set per observed cell with two or more stressors away from
# control: the stressors of the set are exactly those that are not at control
# in that cell, and the other stressors stay at control. With
# background = "all", also one set per subset of two or more of those
# stressors, with the others held at their levels in that cell (the
# background), which then define the reference cell.
find_contrast_sets <- function(spec, cells, pools) {
  stressors <- spec$stressors
  candidates <- list()
  # Sets come from the combinations observed in each group; shared rows only
  # supply cells.
  for (r in which(cells$n_stressed >= 2L & !cells$context %in% pools$shared)) {
    level_row <- vapply(stressors, function(s) cells[[s]][[r]], character(1))
    away <- stressors[level_row != spec$control[stressors]]
    choices <- if (identical(spec$background, "all")) {
      subs <- all_subsets(length(away))
      subs <- subs[lengths(subs) >= 2L & lengths(subs) <= spec$max_order]
      lapply(subs, function(t) away[t])
    } else if (length(away) <= spec$max_order) {
      list(away)
    } else {
      list()
    }
    for (in_set in choices) {
      candidates[[length(candidates) + 1L]] <- list(r = r, in_set = in_set, away = away,
                                                     level_row = level_row)
    }
  }
  set_cells <- list()
  rows <- vector("list", length(candidates))
  width <- max(2L, nchar(length(candidates)))
  for (j in seq_along(candidates)) {
    cand <- candidates[[j]]
    r <- cand$r
    in_set <- cand$in_set
    level_row <- cand$level_row
    bg <- setdiff(cand$away, in_set)
    ctx <- cells$context[[r]]
    m <- length(in_set)
    subsets <- all_subsets(m)
    reference <- spec$control
    reference[bg] <- level_row[bg]
    labels <- vapply(subsets, function(t) {
      lv <- reference
      lv[in_set[t]] <- level_row[in_set[t]]
      make_cell_label(as.list(lv), stressors)
    }, character(1))
    position <- pool_match(labels, ctx, cells, pools)
    present <- !is.na(position)
    size <- lengths(subsets)
    for_net <- size == 0L | size == 1L | size == m
    units <- cells$n_units[position]
    set_id <- sprintf("S%0*d", width, j)
    set_cells[[set_id]] <- data.frame(
      subset = vapply(subsets, function(t) {
        if (length(t) > 0L) {
          paste(in_set[t], collapse = "+")
        } else if (length(bg) == 0L) {
          "control"
        } else {
          "reference"
        }
      }, character(1)),
      size = size,
      cell = labels,
      present = present,
      for_net = for_net,
      context = ifelse(present, cells$context[position], ctx),
      stringsAsFactors = FALSE
    )
    rows[[j]] <- data.frame(
      set_id = set_id,
      context = ctx,
      order = m,
      stressors = paste(in_set, collapse = " x "),
      levels = cells$cell[[r]],
      background = if (length(bg) == 0L) "" else make_cell_label(as.list(level_row[bg]), bg),
      net_identifiable = all(present[for_net]),
      pure_identifiable = if (m >= 3L) all(present) else NA,
      missing_for_net = paste(labels[for_net & !present], collapse = "; "),
      missing_for_pure = if (m >= 3L) paste(labels[!present], collapse = "; ") else "",
      min_units = if (all(present[for_net])) min(units[for_net]) else NA_integer_,
      stringsAsFactors = FALSE
    )
    for (b in spec$by) {
      rows[[j]][[b]] <- cells[[b]][[r]]
    }
  }
  template <- data.frame(
    set_id = character(), context = character(), order = integer(),
    stressors = character(), levels = character(), background = character(),
    net_identifiable = logical(), pure_identifiable = logical(),
    missing_for_net = character(), missing_for_pure = character(),
    min_units = integer(), stringsAsFactors = FALSE
  )
  for (b in spec$by) {
    template[[b]] <- character()
  }
  sets <- bind_rows_base(rows, template)
  sets <- sets[, c("set_id", "context", spec$by, setdiff(names(template), c("set_id", "context", spec$by))),
               drop = FALSE]
  cell_use <- count_cell_use(sets, set_cells, cells)
  use_key <- paste_key(cell_use$context, cell_use$cell)
  sets$shared <- vapply(seq_len(nrow(sets)), function(i) {
    if (!sets$net_identifiable[[i]]) {
      return(FALSE)
    }
    used <- set_cells[[sets$set_id[[i]]]]
    used <- used[used$for_net, , drop = FALSE]
    n <- cell_use$n_sets[use_key %in% paste_key(used$context, used$cell)]
    any(n > 1L)
  }, logical(1))
  list(sets = sets, set_cells = set_cells, cell_use = cell_use)
}

count_cell_use <- function(sets, set_cells, cells) {
  use <- cells[, c("context", "cell"), drop = FALSE]
  use$n_sets <- 0L
  use_key <- paste_key(use$context, use$cell)
  ok <- sets$set_id[sets$net_identifiable]
  for (id in ok) {
    sc <- set_cells[[id]]
    sc <- sc[sc$for_net, , drop = FALSE]
    hit <- use_key %in% paste_key(sc$context, sc$cell)
    use$n_sets[hit] <- use$n_sets[hit] + 1L
  }
  use
}

# Structure of each `groups` column: the stressors that are constant within
# its groups (applied to whole groups), the number of groups holding each
# combination of those stressors, and whether a random intercept for the
# groups can be estimated next to the cell means.
analyse_groups <- function(data, spec, keys) {
  cells_in_context <- n_distinct_by(keys$cell, keys$context)
  lapply(spec$groups, function(g) {
    key <- paste_key(keys$context, as.character(data[[g]]))
    constant <- vapply(spec$stressors, function(s) {
      all(n_distinct_by(keys$levels[[s]], key) == 1L)
    }, logical(1))
    between <- spec$stressors[constant]
    combo <- if (length(between) > 0L) {
      paste_key(keys$context, make_cell_label(keys$levels, between))
    } else {
      keys$context
    }
    per_combo <- n_distinct_by(key, combo)
    units <- n_distinct_by(keys$unit_id, key)
    cells <- n_distinct_by(keys$cell, key)
    group_context <- tapply(keys$context, key, function(v) v[[1L]])
    available <- cells_in_context[match(group_context[names(cells)], names(cells_in_context))]
    structure <- if (all(units == 1L)) {
      "units"
    } else if (length(between) == 0L) {
      "blocks"
    } else if (min(per_combo) >= 2L) {
      "whole_units"
    } else {
      "confounded"
    }
    list(
      group = g,
      between = between,
      structure = structure,
      per_combo = per_combo,
      n_groups = length(units),
      units = units,
      cells = cells,
      complete = all(cells == available),
      # A random intercept is estimable only if the groups are not absorbed
      # by the cell means, i.e. some combination is held by several groups.
      random = !identical(structure, "units") && length(units) > length(per_combo)
    )
  })
}

group_table <- function(group_info) {
  template <- data.frame(
    group = character(), n_groups = integer(), between = character(),
    structure = character(), groups_per_combination = character(),
    units_per_group = character(), cells_per_group = character(),
    random = logical(), stringsAsFactors = FALSE
  )
  rows <- lapply(group_info, function(gi) {
    data.frame(
      group = gi$group,
      n_groups = gi$n_groups,
      between = paste(gi$between, collapse = ", "),
      structure = gi$structure,
      groups_per_combination = if (length(gi$between) > 0L) {
        format_range(gi$per_combo)
      } else {
        NA_character_
      },
      units_per_group = format_range(gi$units),
      cells_per_group = format_range(gi$cells),
      random = gi$random,
      stringsAsFactors = FALSE
    )
  })
  bind_rows_base(rows, template)
}

# Sets that involve a stressor applied to whole groups of one `groups` column
# and use a combination of those stressors held by a single group.
confounded_by_group <- function(sets, set_cells, cells, gi) {
  if (nrow(sets) == 0L || length(gi$between) == 0L) {
    return(rep(FALSE, nrow(sets)))
  }
  vapply(seq_len(nrow(sets)), function(i) {
    in_set <- strsplit(sets$stressors[[i]], " x ", fixed = TRUE)[[1L]]
    if (!any(in_set %in% gi$between)) {
      return(FALSE)
    }
    sc <- set_cells[[sets$set_id[[i]]]]
    sc <- sc[sc$present, , drop = FALSE]
    used <- cells[paste_key(cells$context, cells$cell) %in% paste_key(sc$context, sc$cell), ,
                  drop = FALSE]
    combo <- paste_key(used$context, make_cell_label(used, gi$between))
    any(gi$per_combo[match(combo, names(gi$per_combo))] == 1L, na.rm = TRUE)
  }, logical(1))
}

flag_confounded <- function(sets, set_cells, cells, group_info) {
  flags <- rep(FALSE, nrow(sets))
  for (gi in group_info) {
    flags <- flags | confounded_by_group(sets, set_cells, cells, gi)
  }
  flags
}

group_messages <- function(gi, sets, set_cells, cells) {
  g <- gi$group
  one <- length(gi$between) == 1L
  between <- format_items(quote_items(gi$between), max = Inf)
  random_text <- if (gi$random) {
    sprintf("fit_responses() adds a random intercept for '%s'.", g)
  } else {
    sprintf("No random intercept is added for '%s': it cannot be estimated.", g)
  }
  switch(
    gi$structure,
    units = list(severity = "warning", message = sprintf(
      "Each group of '%s' holds a single unit: it repeats `unit` and is not used.", g
    )),
    blocks = list(severity = "info", message = sprintf(
      paste(
        "'%s': %d groups with %s treatment combination(s) each (%s); no",
        "stressor is applied to whole groups. %s"
      ),
      g, gi$n_groups, format_range(gi$cells),
      if (gi$complete) "complete blocks" else "incomplete blocks",
      random_text
    )),
    whole_units = list(severity = "info", message = sprintf(
      paste(
        "%s %s applied to whole groups of '%s' (%s groups per combination of",
        "%s levels): for effects that involve %s, the replication is the",
        "number of groups, not of units. %s"
      ),
      between, if (one) "is" else "are", g, format_range(gi$per_combo),
      if (one) "its" else "their", if (one) "it" else "them", random_text
    )),
    confounded = list(severity = "warning", message = sprintf(
      paste(
        "%s %s applied to whole groups of '%s', and %d combination(s) of %s",
        "levels were given to a single group: effects that involve %s cannot",
        "be separated from differences among groups of '%s' (%d contrast",
        "set(s) flagged in `sets$confounded`). %s"
      ),
      between, if (one) "is" else "are", g, sum(gi$per_combo == 1L),
      if (one) "its" else "their", if (one) "it" else "them", g,
      sum(confounded_by_group(sets, set_cells, cells, gi)), random_text
    ))
  )
}

# Two `groups` columns that only together define groups within which some
# stressor takes a single level, although neither column does alone: the
# groups of one are probably labelled within the other (plots numbered 1-8
# within each block), so that groups with the same label are merged and the
# stressors applied to them go unnoticed.
nested_group_labels <- function(data, spec, keys) {
  out <- list()
  groups <- spec$groups
  if (length(groups) < 2L) {
    return(out)
  }
  constant_within <- function(key) {
    spec$stressors[vapply(spec$stressors, function(s) {
      all(n_distinct_by(keys$levels[[s]], key) == 1L)
    }, logical(1))]
  }
  pairs <- utils::combn(groups, 2L, simplify = FALSE)
  for (pair in pairs) {
    g <- pair[[1L]]
    h <- pair[[2L]]
    key_g <- paste_key(keys$context, as.character(data[[g]]))
    key_h <- paste_key(keys$context, as.character(data[[h]]))
    key_gh <- paste_key(keys$context, as.character(data[[g]]), as.character(data[[h]]))
    if (all(n_distinct_by(keys$unit_id, key_gh) == 1L)) {
      next
    }
    extra <- setdiff(constant_within(key_gh), c(constant_within(key_g), constant_within(key_h)))
    if (length(extra) == 0L) {
      next
    }
    one <- length(extra) == 1L
    out[[length(out) + 1L]] <- data.frame(
      check = "nested_labels",
      severity = "warning",
      message = sprintf(
        paste(
          "Within each combination of '%s' and '%s', %s %s a single level, but",
          "not within the groups of either column alone: the groups of one are",
          "probably labelled within the other (plots numbered 1-8 within each",
          "block), so groups with the same label are merged and %s %s not",
          "recognised as applied to whole groups. Give each nested group a unique",
          "label, e.g. paste(%s, %s), or add that combination as another",
          "`groups` column."
        ),
        g, h, format_items(quote_items(extra), max = Inf), if (one) "takes" else "take",
        if (one) "it" else "they", if (one) "is" else "are", h, g
      ),
      stringsAsFactors = FALSE
    )
  }
  out
}

# Labels that may hide an error in the coding of the treatments: levels that
# differ only in letter case or spaces (a warning), and, with `by`, non-control
# levels of a stressor observed in only some `by` groups, which may be one
# treatment under two labels ("MP" in one population, "PET" in another).
label_checks <- function(data, spec, cells, pools) {
  out <- list()
  for (col in c(spec$stressors, spec$by, spec$groups)) {
    values <- unique(as.character(data[[col]]))
    invalid <- !validUTF8(values)
    if (any(invalid)) {
      out[[length(out) + 1L]] <- data.frame(
        check = "encoding", severity = "warning",
        message = sprintf(
          paste(
            "Levels of '%s' are not valid UTF-8 text (%s): the file was probably",
            "read with the wrong encoding. Read it again with its encoding (for",
            "example fileEncoding = \"latin1\" in read.csv()) so that labels print",
            "correctly."
          ),
          col, format_items(quote_items(iconv(values[invalid], "latin1", "UTF-8")), max = 3L)
        ),
        stringsAsFactors = FALSE
      )
      # Compared as Latin-1, which any byte sequence is.
      values[invalid] <- iconv(values[invalid], "latin1", "UTF-8")
    }
    folded <- tolower(gsub("[[:space:]]+", " ", trimws(values)))
    dup <- unique(folded[duplicated(folded)])
    if (length(dup) == 0L) {
      next
    }
    pieces <- vapply(dup, function(f) {
      format_items(quote_items(sort(values[folded == f], method = "radix")), max = Inf)
    }, character(1))
    out[[length(out) + 1L]] <- data.frame(
      check = "similar_levels", severity = "warning",
      message = sprintf(
        paste(
          "Levels of '%s' that differ only in letter case or spaces: %s. They are",
          "treated as different levels; if they are the same, give them the same label."
        ),
        col, paste(unname(pieces), collapse = "; ")
      ),
      stringsAsFactors = FALSE
    )
  }
  groups <- names(pools$pools)
  if (length(groups) < 2L) {
    return(out)
  }
  for (s in spec$stressors) {
    lev <- setdiff(spec$levels[[s]], spec$control[[s]])
    where <- lapply(lev, function(l) {
      groups[vapply(groups, function(ctx) {
        any(cells[[s]][cells$context %in% pool_of(pools, ctx)] == l)
      }, logical(1))]
    })
    partial <- which(lengths(where) > 0L & lengths(where) < length(groups))
    if (length(partial) == 0L) {
      next
    }
    pieces <- vapply(partial, function(i) {
      sprintf("'%s' (only %s)", lev[[i]], format_items(quote_items(where[[i]]), max = 2L))
    }, character(1))
    out[[length(out) + 1L]] <- data.frame(
      check = "levels_by_group", severity = "info",
      message = sprintf(
        paste(
          "Levels of '%s' observed in only some `by` groups: %s. If one of them is",
          "the same treatment as another level under a different label, give both",
          "the same label; otherwise the combinations they lack are missing."
        ),
        s, format_items(pieces, max = 3L)
      ),
      stringsAsFactors = FALSE
    )
  }
  out
}

design_checks <- function(spec, keys, cells, missing_cells, sets, group_info = list(),
                          nested = list(), labels = list(),
                          pools = context_pools(cells, spec)) {
  out <- list()
  add <- function(check, severity, message) {
    out[[length(out) + 1L]] <<- data.frame(
      check = check, severity = severity, message = message, stringsAsFactors = FALSE
    )
  }
  # Replication unit and treatment constancy
  if (is.null(spec$unit)) {
    add("replication_unit", "info", paste(
      "No `unit` given: every row is treated as an independent experimental",
      "unit. If several rows come from the same tank, plot or jar, give that",
      "column as `unit`."
    ))
  } else {
    unit_cells <- n_distinct_by(paste_key(keys$context, keys$cell), keys$unit_id)
    varying <- names(unit_cells)[unit_cells > 1L]
    if (length(varying) > 0L) {
      add("replication_unit", "problem", sprintf(
        paste(
          "%d unit(s) receive more than one treatment combination or `by`",
          "group (%s). Treatments that change within a unit (sequential",
          "exposures) are not supported: code the treatment as one label per",
          "unit or analyse each phase separately. If unit labels are reused",
          "(tanks numbered 1-4 within each treatment, or the same tank in",
          "several `by` groups such as habitats or species), give each unit a",
          "unique label within its `by` group, e.g. paste(treatment, tank) or",
          "paste(tank, habitat): contrasts are estimated within each `by`",
          "group, so the correlation between `by` groups of the same unit does",
          "not affect them."
        ),
        length(varying), format_items(quote_items(varying))
      ))
    } else {
      add("replication_unit", "info", sprintf(
        "Experimental unit: '%s' (%d units, each in a single treatment combination).",
        spec$unit, length(unique(keys$unit_id))
      ))
    }
  }
  # Replication per cell
  few <- cells$n_units < 2L
  single <- short_cells(cells$cell[few], spec$control, cells$context[few])
  if (length(single) > 0L) {
    add("replication", "warning", sprintf(
      paste(
        "%d treatment combination(s) have a single experimental unit (%s):",
        "their means have no within-combination replication."
      ),
      length(single), format_items(quote_items(single), max = 3L)
    ))
  } else {
    add("replication", "info", sprintf(
      "Units per treatment combination: %s.", format_range(cells$n_units)
    ))
  }
  # Subsamples
  sampling <- if (is.null(spec$time)) keys$unit_id else paste_key(keys$unit_id, keys$time_value)
  rows_per <- table(sampling)
  if (any(rows_per > 1L) && !is.null(spec$unit)) {
    add("subsamples", "info", sprintf(
      paste(
        "Some units have several rows per sampling (up to %d): they are",
        "subsamples, not replicates. Replication is counted in units and",
        "fit_responses() adds a random intercept per unit."
      ),
      max(rows_per)
    ))
  }
  # Time
  if (!is.null(spec$time)) {
    times_per_unit <- n_distinct_by(keys$time_value, keys$unit_id)
    n_times <- length(unique(keys$time_value))
    if (any(times_per_unit > 1L)) {
      add("repeated_measures", "info", sprintf(
        paste(
          "Units are measured repeatedly (%d sampling times; %s per unit).",
          "Samplings of the same unit are not independent: by default",
          "fit_responses() adds a random intercept per unit, which assumes the",
          "same correlation between any two samplings of a unit; temporal",
          "autocorrelation can be modelled through its `random` argument."
        ),
        n_times, format_range(times_per_unit)
      ))
      ctx_times <- tapply(keys$time_value, keys$context, function(v) length(unique(v)))
      unit_ctx <- tapply(keys$context, keys$unit_id, function(v) v[[1L]])
      unit_ctx_times <- ctx_times[match(unit_ctx[names(times_per_unit)], names(ctx_times))]
      lacking <- names(times_per_unit)[times_per_unit < unit_ctx_times]
      if (length(lacking) > 0L) {
        add("unbalanced_time", "warning", sprintf(
          paste(
            "%d unit(s) are not sampled at every time of their group (%s).",
            "Contrasts at a given time use only the units sampled then."
          ),
          length(lacking), format_items(quote_items(lacking))
        ))
      }
    } else {
      add("repeated_measures", "info", sprintf(
        "Each unit is sampled once (%d sampling times): times are independent groups of units.",
        n_times
      ))
    }
  }
  # Pre-treatment samplings
  if (!is.null(spec$time)) {
    times <- unique(keys$time_value)
    per <- spec$baseline_by
    uniform <- !is.null(per) && all(vapply(per, function(t) setequal(t, per[[1L]]), logical(1)))
    if (uniform) {
      add("baseline", "info", sprintf(
        paste(
          "Baseline (pre-treatment) time(s): %s. They are kept in the models,",
          "but interactions are not estimated or tested at them."
        ),
        format_items(quote_items(spec$baseline), max = 5L)
      ))
    } else if (!is.null(per)) {
      with_times <- names(per)[lengths(per) > 0L]
      pieces <- vapply(with_times, function(st) {
        sprintf("'%s' before %s", st, format_items(quote_items(per[[st]]), max = 5L))
      }, character(1))
      add("baseline", "info", sprintf(
        paste(
          "Stressors applied later than others: %s; the other stressors were",
          "applied before the first sampling. Interactions that involve a stressor",
          "are not estimated or tested before it was applied; those samplings are",
          "kept in the models."
        ),
        paste(pieces, collapse = "; ")
      ))
    } else {
      value <- suppressWarnings(as.numeric(times))
      early <- times[!is.na(value) & value <= 0]
      if (!anyNA(value) && length(early) > 0L && length(early) < length(times)) {
        add("baseline", "info", sprintf(
          paste(
            "Sampling time(s) %s may be pre-treatment samplings: if so, declare",
            "them in `baseline` so that interactions are not estimated and",
            "tested there."
          ),
          format_items(quote_items(early), max = 5L)
        ))
      }
    }
  }
  # Missing combinations
  if (!missing_cells$enumerated) {
    add("missing_combinations", "info", sprintf(
      "The full factorial has %s combinations; missing combinations were not listed.",
      format(missing_cells$n_full, big.mark = ",")
    ))
  } else if (nrow(missing_cells$table) > 0L) {
    add("missing_combinations", "info", sprintf(
      paste(
        "%d of the %d combinations of the observed levels are absent (%s).",
        "Interactions that need them are not identifiable (see `sets`)."
      ),
      nrow(missing_cells$table), missing_cells$n_full,
      format_items(quote_items(short_cells(missing_cells$table$cell, spec$control,
                                           missing_cells$table$context)), max = 3L)
    ))
  }
  # Identifiable sets
  ok <- sets$sets[sets$sets$net_identifiable, , drop = FALSE]
  if (nrow(ok) == 0L) {
    add("identifiable_sets", "problem", paste(
      "No group of cells identifies an interaction. Each interaction needs, in",
      "the same `by` group (or among the rows shared with it), the control,",
      "each stressor alone and their combination."
    ))
  } else {
    by_order <- table(ok$order)
    pieces <- sprintf("%d with %s stressors", as.integer(by_order), names(by_order))
    n_bg <- sum(nzchar(ok$background))
    add("identifiable_sets", "info", sprintf(
      "Identifiable contrast sets (net deviation): %s%s.", paste(pieces, collapse = "; "),
      if (n_bg > 0L) {
        sprintf(" (%d of them with other stressors away from control)", n_bg)
      } else {
        ""
      }
    ))
  }
  partial <- ok[ok$order >= 3L & !ok$pure_identifiable %in% TRUE, , drop = FALSE]
  if (nrow(partial) > 0L) {
    add("pure_not_identifiable", "warning", sprintf(
      paste(
        "%d set(s) of three or more stressors lack some lower-order",
        "combinations (%s): their net deviation mixes the highest-order",
        "interaction with lower-order ones and the design cannot split them."
      ),
      nrow(partial),
      format_items(quote_items(short_cells(partial$levels, spec$control, partial$context)),
                   max = 2L)
    ))
  }
  if (any(ok$shared)) {
    n_max <- max(sets$cell_use$n_sets)
    add("shared_cells", "info", sprintf(
      paste(
        "%d identifiable set(s) share cells (one cell is used by up to %d",
        "sets, typically the control): their estimates are correlated, not",
        "independent replicates of an interaction."
      ),
      sum(ok$shared), n_max
    ))
  }
  # Rows shared by several `by` groups
  if (length(pools$shared) > 0L) {
    users <- names(pools$pools)[lengths(pools$pools) > 1L]
    shared_cells <- cells[cells$context %in% pools$shared, , drop = FALSE]
    add("shared", "info", sprintf(
      paste(
        "Rows of %s are shared: their %d treatment combination(s) (%s) complete the",
        "contrast sets of %d other `by` group(s), whose estimates are therefore",
        "correlated."
      ),
      format_items(quote_items(pools$shared), max = 3L), nrow(shared_cells),
      format_items(quote_items(short_cells(shared_cells$cell, spec$control)), max = 3L),
      length(users)
    ))
    unused <- setdiff(pools$shared, unlist(lapply(pools$pools, `[`, -1L)))
    if (length(unused) > 0L) {
      add("shared", "warning", sprintf(
        "No `by` group matches the shared rows of %s in the other `by` columns: no set uses them.",
        format_items(quote_items(unused), max = 3L)
      ))
    }
    both <- character()
    for (ctx in users) {
      own <- cells$cell[cells$context == ctx]
      common <- intersect(own, cells$cell[cells$context %in% pool_of(pools, ctx)[-1L]])
      both <- c(both, short_cells(common, spec$control, rep(ctx, length(common))))
    }
    if (length(both) > 0L) {
      add("shared", "warning", sprintf(
        paste(
          "%d combination(s) are both shared and observed within a `by` group",
          "(%s): the sets of that group use its own and ignore the shared ones."
        ),
        length(both), format_items(quote_items(both), max = 3L)
      ))
    }
  }
  for (gi in group_info) {
    msg <- group_messages(gi, sets$sets, sets$set_cells, cells)
    add("groups", msg$severity, msg$message)
  }
  for (lc in labels) {
    add(lc$check, lc$severity, lc$message)
  }
  for (nl in nested) {
    add(nl$check, nl$severity, nl$message)
  }
  many <- spec$stressors[lengths(spec$levels) > 5L]
  for (s in many) {
    add("many_levels", "info", sprintf(
      paste(
        "Stressor '%s' has %d levels, each treated as a separate category:",
        "a gradient is not modelled as a continuous predictor."
      ),
      s, length(spec$levels[[s]])
    ))
  }
  checks <- bind_rows_base(out)
  severity_rank <- match(checks$severity, c("problem", "warning", "info"))
  checks <- checks[order(severity_rank, method = "radix"), , drop = FALSE]
  rownames(checks) <- NULL
  checks
}

#' @export
print.must_design <- function(x, ...) {
  spec <- x$spec
  cat(sprintf(
    "<must_design> %d rows, %d stressors, %d observed treatment combinations\n",
    x$n_rows, length(spec$stressors), nrow(x$cells)
  ))
  ctrl <- paste0(spec$stressors, " (", spec$control, ")")
  cat("Stressors (control): ", paste(ctrl, collapse = ", "), "\n", sep = "")
  cat(sprintf(
    "Unit: %s | Time: %s | By: %s | Groups: %s\n",
    if (is.null(spec$unit)) "none (rows)" else spec$unit,
    if (is.null(spec$time)) {
      "none"
    } else if (is.null(spec$baseline_by)) {
      spec$time
    } else if (!is.null(spec$baseline) &&
               all(vapply(spec$baseline_by, setequal, logical(1), spec$baseline))) {
      sprintf("%s (baseline: %s)", spec$time, paste(spec$baseline, collapse = ", "))
    } else {
      later <- spec$baseline_by[lengths(spec$baseline_by) > 0L]
      sprintf("%s (baseline: %s)", spec$time, paste(
        vapply(names(later), function(st) {
          sprintf("%s before %s", st, paste(later[[st]], collapse = ", "))
        }, character(1)),
        collapse = "; "
      ))
    },
    if (is.null(spec$by)) {
      "none"
    } else if (is.null(spec$shared)) {
      paste(spec$by, collapse = ", ")
    } else {
      sprintf("%s (shared: %s)", paste(spec$by, collapse = ", "),
              paste(names(spec$shared), "=", spec$shared, collapse = ", "))
    },
    if (is.null(spec$groups)) "none" else paste(spec$groups, collapse = ", ")
  ))
  ok <- x$sets[x$sets$net_identifiable, , drop = FALSE]
  cat(sprintf(
    "Contrast sets: %d candidate(s), %d identifiable%s\n", nrow(x$sets), nrow(ok),
    if (identical(spec$background, "all")) "; other stressors at control or not" else ""
  ))
  cat("Checks:\n")
  print_checks(x$checks)
  invisible(x)
}

print_checks <- function(checks) {
  for (i in seq_len(nrow(checks))) {
    lines <- strwrap(checks$message[[i]], width = 65L)
    cat(sprintf(
      "  %-8s %s\n", checks$severity[[i]],
      paste(lines, collapse = paste0("\n", strrep(" ", 11L)))
    ))
  }
}
