#' Classify interactions (Piggott et al. 2015)
#'
#' Codes the single, combined and interaction effects of each identifiable
#' contrast set as positive, negative or not significant, and assigns the
#' directional class of Piggott et al. (2015, Table 1) under the
#' simple-addition (AD) and multiplicative (MU) null models.
#'
#' @param x A `must_interactions` object from [estimate_interactions()].
#' @param alpha Significance level applied to the adjusted p-values
#'   (`p_adjusted`). Default 0.05.
#'
#' @return A data frame with one row per response, contrast set, sampling
#'   time and null model: `response`, `context` and the `by` columns,
#'   `time`, `set_id`, `order`, `stressors`, `levels`, `background`,
#'   `null_model`, `scale`, `estimand` (`"net"`), `stressor_a`, `code_a`,
#'   `stressor_b`, `code_b` (two-stressor sets only), `single_codes` (the
#'   code of every single effect, `"salinity +; warming 0"`),
#'   `code_combined`, `code_interaction` (`"+"`, `"-"`, `"0"`, or `"++"` /
#'   `"--"` for the strong cases defined by Piggott et al.), `pattern`,
#'   `class` (`"AD"`, `"+S"`, `"-S"`, `"+A"`, `"-A"`), `meaning`,
#'   `borderline` (some adjusted p-value used is within two Monte Carlo
#'   standard errors of `alpha`, so that another seed could change the
#'   class), `shared` (the set shares cells with other sets), `confounded`
#'   (the set is confounded with differences among groups of units, see
#'   [describe_design()]), `min_units`, `status` (`"ok"`,
#'   `"not_applicable"` or `"not_estimable"`) and `reason`.
#'
#' @details
#' Each effect is coded `"0"` if its adjusted p-value is at least `alpha`,
#' and `"+"` or `"-"` after the sign of its estimate otherwise. A
#' significant interaction is coded `"++"` (`"--"`) when the combined effect
#' is larger (smaller) than every single effect, or when the absolute
#' interaction is larger than every absolute single effect (legend of Table
#' 1 in Piggott et al. 2015), comparing point estimates on the scale of the
#' null model. Table 1 then reduces to:
#'
#' | Interaction | Single effects | Class |
#' |---|---|---|
#' | `0` | any | AD |
#' | `++` | any | +S |
#' | `--` | any | -S |
#' | `+` | none negative | +S |
#' | `+` | at least one negative | -A |
#' | `-` | none positive | -S |
#' | `-` | at least one positive | +A |
#'
#' Under MU, +S and -S mean "more positive" or "more negative" than the
#' multiplicative expectation, not than the sum of the single effects.
#'
#' **Three or more stressors.** The net deviation of a set of k stressors is
#' the combined effect minus the sum of the k single effects (AD, or the
#' same on the log scale for MU), the comparison that Table 1 classifies; the
#' same rules are applied with k single effects. This generalization is a
#' design decision of this package: Piggott et al. (2015) define the classes
#' for two stressors. The pure k-way interaction is not classified.
#'
#' **Combinations with only zeros** (`status = "boundary"` in
#' [estimate_interactions()]) are classified from the sign of their
#' estimate (which may be infinite under MU) and their p-value.
#'
#' @section Assumptions:
#' * AD means that no deviation from the null model was detected. It is not
#'   evidence of additivity: with few replicates most interactions will be
#'   classed AD.
#' * The class of a set with `confounded = TRUE` may reflect differences
#'   among groups of units (rooms, blocks) as much as the stressors.
#' * The classification is defined for null models that combine effects by a
#'   sum (AD on the response scale, MU on the log scale). The dominance model
#'   gets `status = "not_applicable"`, as do rows of `"not_applicable"`
#'   estimates (see `sensitivity` in [estimate_interactions()]).
#'
#' @references
#' Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015).
#' Reconceptualizing synergism and antagonism among multiple stressors.
#' *Ecology and Evolution*, 5(7), 1538–1547. \doi{10.1002/ece3.1465}
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
#' classify_interactions(est)
#' @export
classify_interactions <- function(x, alpha = 0.05) {
  if (!inherits(x, "must_interactions")) {
    abort_input(sprintf(
      "`x` must be the output of estimate_interactions(), not %s.", describe_value(x)
    ))
  }
  check_number(alpha, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  est <- x$estimates
  by <- x$design$spec$by
  sets <- x$design$sets
  id_cols <- c("response", "context", by, "time", "set_id", "order", "stressors",
               "levels", "background", "null_model", "scale")
  row_key <- do.call(paste_key, unname(as.list(est[id_cols])))
  is_net <- est$component == "interaction" & est$estimand %in% "net"
  groups <- est[is_net, id_cols, drop = FALSE]
  group_key <- row_key[is_net]
  rows_by_key <- split(seq_len(nrow(est)), row_key)
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    part <- est[rows_by_key[[group_key[[i]]]], , drop = FALSE]
    classify_group(groups[i, , drop = FALSE], part, alpha, sets)
  })
  out <- bind_rows_base(rows, template = classification_template(by))
  rownames(out) <- NULL
  out
}

classification_template <- function(by) {
  out <- data.frame(response = character(), context = character(), stringsAsFactors = FALSE)
  for (b in by) {
    out[[b]] <- character()
  }
  cbind(out, data.frame(
    time = character(), set_id = character(), order = integer(),
    stressors = character(), levels = character(), background = character(),
    null_model = character(), scale = character(), estimand = character(),
    stressor_a = character(), code_a = character(), stressor_b = character(),
    code_b = character(), single_codes = character(), code_combined = character(),
    code_interaction = character(), pattern = character(), class = character(),
    meaning = character(), borderline = logical(), shared = logical(),
    confounded = logical(), min_units = integer(), status = character(),
    reason = character(), stringsAsFactors = FALSE
  ))
}

classify_group <- function(g, part, alpha, sets) {
  set_row <- sets[sets$set_id == g$set_id, , drop = FALSE]
  base <- g
  stressors <- strsplit(g$stressors, " x ", fixed = TRUE)[[1L]]
  two <- length(stressors) == 2L
  base$estimand <- "net"
  base$stressor_a <- if (two) stressors[[1L]] else NA_character_
  base$code_a <- NA_character_
  base$stressor_b <- if (two) stressors[[2L]] else NA_character_
  base$code_b <- NA_character_
  base$single_codes <- NA_character_
  base$code_combined <- NA_character_
  base$code_interaction <- NA_character_
  base$pattern <- NA_character_
  base$class <- NA_character_
  base$meaning <- NA_character_
  base$borderline <- NA
  base$shared <- isTRUE(set_row$shared)
  base$confounded <- isTRUE(set_row$confounded)
  base$min_units <- set_row$min_units
  base$status <- "ok"
  base$reason <- ""
  if (identical(g$null_model, "DO")) {
    base$status <- "not_applicable"
    base$reason <- paste(
      "Table 1 of Piggott et al. (2015) compares the combined effect with a",
      "sum of single effects; the dominance model has no such sum."
    )
    return(base)
  }
  singles <- lapply(stressors, function(st) {
    part[part$component == "single" & part$stressor %in% st, , drop = FALSE]
  })
  combined <- part[part$component == "combined", , drop = FALSE]
  interaction <- part[part$component == "interaction" & part$estimand %in% "net", , drop = FALSE]
  pieces <- c(singles, list(combined, interaction))
  statuses <- unlist(lapply(pieces, `[[`, "status"))
  if (any(vapply(pieces, nrow, integer(1)) != 1L) ||
      !all(statuses %in% c("ok", "boundary"))) {
    base$status <- if (length(statuses) > 0L && all(statuses == "not_applicable")) {
      "not_applicable"
    } else {
      "not_estimable"
    }
    reasons <- unique(unlist(lapply(pieces, function(p) p$reason[nzchar(p$reason)])))
    base$reason <- if (length(reasons) > 0L) paste(reasons, collapse = " ") else "Some effects were not estimated."
    return(base)
  }
  single_estimates <- vapply(singles, `[[`, numeric(1), "estimate")
  single_codes <- vapply(singles, function(p) effect_code(p$estimate, p$p_adjusted, alpha),
                         character(1))
  code_i <- effect_code(interaction$estimate, interaction$p_adjusted, alpha)
  base$code_combined <- effect_code(combined$estimate, combined$p_adjusted, alpha)
  base$code_interaction <- strong_interaction_code(
    code_i, single_estimates, combined$estimate, interaction$estimate
  )
  if (two) {
    base$code_a <- single_codes[[1L]]
    base$code_b <- single_codes[[2L]]
  }
  base$single_codes <- paste(stressors, single_codes, collapse = "; ")
  base$class <- piggott_class(single_codes, base$code_interaction)
  base$meaning <- class_meaning(base$class, g$null_model)
  base$pattern <- paste(c(paste(if (two) c("a", "b") else stressors, single_codes),
                          paste("i", base$code_interaction)), collapse = " | ")
  base$borderline <- any(vapply(pieces, is_borderline, logical(1), alpha = alpha))
  base
}

effect_code <- function(estimate, p_adjusted, alpha) {
  if (is.na(estimate) || !is.finite(p_adjusted)) {
    return(NA_character_)
  }
  if (p_adjusted >= alpha) "0" else if (estimate > 0) "+" else "-"
}

# An adjusted p-value within two Monte Carlo standard errors of alpha (the
# standard error of the raw p-value, scaled by the adjustment).
is_borderline <- function(p, alpha) {
  if (!is.finite(p$mc_se) || !is.finite(p$p_adjusted) || !is.finite(p$p_value) ||
      p$p_value <= 0) {
    return(FALSE)
  }
  abs(p$p_adjusted - alpha) < 2 * p$mc_se * p$p_adjusted / p$p_value
}

# "++" / "--" following the legend of Table 1 in Piggott et al. (2015):
# `singles` are the estimates of the single effects.
strong_interaction_code <- function(code_i, singles, combined, interaction) {
  if (is.na(code_i) || identical(code_i, "0")) {
    return(code_i)
  }
  beyond_singles <- abs(interaction) > max(abs(singles))
  if (identical(code_i, "+") && (combined > max(singles) || beyond_singles)) {
    return("++")
  }
  if (identical(code_i, "-") && (combined < min(singles) || beyond_singles)) {
    return("--")
  }
  code_i
}

# Table 1 of Piggott et al. (2015), reduced to its decision rules;
# `single_codes` has one code per single effect.
piggott_class <- function(single_codes, code_i) {
  if (anyNA(c(single_codes, code_i))) {
    return(NA_character_)
  }
  any_negative <- "-" %in% single_codes
  any_positive <- "+" %in% single_codes
  switch(
    code_i,
    "0" = "AD",
    "++" = "+S",
    "--" = "-S",
    "+" = if (any_negative) "-A" else "+S",
    "-" = if (any_positive) "+A" else "-S"
  )
}

class_meaning <- function(class, null_model) {
  reference <- if (identical(null_model, "MU")) "the multiplicative expectation" else "the sum of the single effects"
  switch(
    class,
    "AD" = paste("No deviation from", reference, "detected (not evidence of additivity)"),
    "+S" = paste("Combined effect more positive than", reference),
    "-S" = paste("Combined effect more negative than", reference),
    "+A" = paste("Combined effect less positive than", reference),
    "-A" = paste("Combined effect less negative than", reference),
    NA_character_
  )
}
