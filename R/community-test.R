# Community tests of the null models MU and AD_taxa on every taxon ------------

#' Test the interactions of a community with multivariate GLMs
#'
#' @description
#' Tests whether the stressors of each contrast set interact on a community
#' as a whole, under a null model applied to every taxon: the multiplicative
#' null model (`"MU"`, Madge Pimentel et al. 2025) or the additive one,
#' truncated at zero where the stressors eliminate the taxon (`"AD_taxa"`,
#' Thompson et al. 2018). Each taxon gets a generalised linear model with and
#' without the interaction: without it, the effects of the stressors add up
#' on the log scale (a log link, MU) or on the response (an identity link,
#' AD). The test of the community is the sum of the statistics of all taxa,
#' with a p-value from resampling whole experimental units under the null
#' model. This is the approach of multivariate generalised linear models
#' (package mvabund) that Szöcs et al. (2015) proposed for mesocosm
#' experiments, applied to the interaction.
#'
#' @param data A data frame with one row per sample: the design columns and
#'   one column per taxon with its counts, densities or biomass.
#' @param taxa Names of the columns of `data` with the taxa.
#' @param design The design, from [describe_design()] on `data`.
#' @param family Distribution of each taxon: `"negative.binomial"`
#'   (default) or `"poisson"` for counts (fitted with [mvabund::manyglm()]
#'   under MU); or `"tweedie"` for amounts that are not counts, with zeros
#'   and positive values, such as densities or biomass (see "Counts,
#'   densities and biomass").
#' @param null_models One or both of `"MU"` (default; multiplicative on each
#'   taxon) and `"AD_taxa"` (additive on each taxon, truncated at zero; see
#'   "Additive null model").
#' @param effort `NULL` (default) or the column of `data` with the sampling
#'   effort of each sample (volume of water, area, fields counted), or one
#'   such column per taxon, in the order of `taxa`. The models then describe
#'   the density, the count divided by the effort, through an offset
#'   `log(effort)`, and the deviations are those of the densities. Efforts
#'   must be positive.
#' @param resampling How the p-values are computed, always resampling whole
#'   experimental units under the null model: `"pit.trap"` (default; the
#'   bootstrap of probability integral transform residuals, Warton et al.
#'   2017), `"perm.resid"` (permutation of residuals) or `"montecarlo"`
#'   (parametric bootstrap), the methods of [mvabund::anova.manyglm()]. With
#'   `family = "tweedie"` or `"AD_taxa"` only `"pit.trap"` is available.
#' @param n_boot Number of resamples (at least 99). Default 999.
#' @param overall If `TRUE` (default) and there are several sampling times,
#'   also tests the interactions of all times together.
#' @param p_adjust Method of [stats::p.adjust()] for the tests of the
#'   community, applied across sets and times within each null model (the
#'   tests of all times apart). Default `"BH"`.
#' @param seed `NULL` or a whole number for the resampling; the seed used is
#'   kept in `settings$seed`.
#' @param quiet If `FALSE` (default), reports each set as it is tested.
#'
#' @details
#' **Models.** For each contrast set identifiable by [describe_design()] and
#' each sampling time (except those in `baseline`), the test uses the
#' samples of the reference cell, of each stressor of the set alone and of
#' their full combination. The null model has one effect for each stressor
#' (`Y ~ s1 + ... + sm`, with `sj` marking the samples where stressor `j`
#' is applied): with the log link the effects add up on the log scale, so
#' the mean of the combination is the product of the mean of the reference
#' and the ratios of the stressors alone, the multiplicative null model of
#' [estimate_interactions()]. The full model has one mean per cell
#' (`Y ~ cell`), and the difference between both, one degree of freedom, is
#' the net interaction under MU. With all times together the models are
#' `Y ~ time * (s1 + ... + sm)` and `Y ~ time * cell`, with one interaction
#' per time.
#'
#' **Additive null model.** Under `"AD_taxa"` the mean of the combination
#' is the mean of the reference plus the differences of the stressors alone,
#' the additive null model on the response scale (an identity link), with
#' the means kept non-negative. Where that sum is zero or less, the stressors
#' are expected to eliminate the taxon, as in the compositional null model of
#' Thompson et al. (2018, eq. 3) that [estimate_community()] uses: when the
#' combination has none of it, the null model fits it (no deviation);
#' otherwise the null model is fitted at the boundary of a zero mean for the
#' combination. For a single response, [estimate_interactions()] with
#' `inference = "resampling"` tests AD without truncation, so that an
#' eliminated response counts as a deviation (flooring; Madge Pimentel et al.
#' 2025). The models are fitted by maximum likelihood with the size of each
#' negative binomial (estimated once per taxon, with one mean per treatment
#' combination and time, from all the samples of the tests) or the power of
#' each Tweedie fixed; the statistic of each taxon is, as with Tweedie
#' models, the drop in deviance divided by the dispersion estimated from
#' the null model, `df0 * (D0 - D1) / D0` (with the size fixed, a likelihood
#' ratio rejected true null models too often), and the p-values come from
#' the PIT-trap of whole units (mvabund only fits log links). The estimate of each taxon is its net
#' interaction under AD_taxa on the response scale, \eqn{\mu_c - \max(0,
#' \mu_0 + \sum_j (\mu_j - \mu_0))}{mu_c - max(0, mu_0 + sum_j (mu_j -
#' mu_0))}.
#'
#' **Test of the community.** The statistic is the sum over taxa of the
#' statistics of the interaction (the taxa are assumed independent in the
#' statistic). Its p-value comes from resampling whole experimental units
#' under the null model, so that the correlation between the taxa of a unit,
#' and between the samples of a unit over time, is kept. Szöcs et al. (2015)
#' tested the community response of mesocosms this way, with negative
#' binomial models and permutations; `resampling = "perm.resid"` permutes the
#' residuals of whole units. With counts, both models are fitted with
#' [mvabund::manyglm()] and tested with [mvabund::anova.manyglm()]
#' (likelihood-ratio statistics, `cor.type = "I"`, `block` the units).
#'
#' **Counts, densities and biomass.** Densities computed from counts
#' (individuals per litre, cells per millilitre) are best tested from the
#' counts, with the amount that converts them (volume of the sample, of the
#' subsample or of the chamber counted) as `effort`: the counts carry the
#' information on their variance that the densities lose. When the counts
#' are not available, or the response is not a count (biomass, biovolume),
#' `family = "tweedie"` models each taxon with a Tweedie distribution of
#' power between 1 and 2, a Poisson sum of gamma amounts that has zeros and
#' positive values, as [fit_responses()] can. The power of each taxon is
#' estimated once by maximum likelihood with glmmTMB, from a model with one
#' mean per treatment combination and time fitted to all the samples of the
#' tests (without the combinations and times where the taxon is absent,
#' whose mean is zero), and kept fixed. Each test fits the null model with
#' [stats::glm.fit()] and the full model by its closed form. The statistic
#' of a taxon is the drop in deviance divided by the dispersion of the null
#' model estimated from its deviance, `df0 * (D0 - D1) / D0`, which does not
#' depend on the unit of the taxon and is at most `df0`, the residual
#' degrees of freedom of the null model. The p-values come from the PIT-trap
#' of whole units, as in mvabund: the probability integral transform
#' residuals of the null models of all taxa of a unit are resampled
#' together and transformed back with the Tweedie quantile function,
#' computed with the series of the compound Poisson-gamma distribution
#' (Dunn & Smyth 2005).
#'
#' **Calibration.** In simulated experiments with four units per combination
#' and eight taxa that followed MU, with or without a factor of each unit
#' shared by all its taxa, the tests rejected the null model, at the 5%
#' level, in 2–4% of the experiments with counts and `"pit.trap"` (also with
#' an effort that varied between samples), in 4–5% with `"perm.resid"` and
#' in 1–4% with Tweedie densities. When the combination left two of the
#' eight taxa at 30% of what MU expects, the tests detected it in 30% of the
#' experiments with negative binomial counts and in 61% with Tweedie
#' densities (less variable, in that simulation, than the counts). Under
#' `"AD_taxa"`, with six taxa that followed the additive model truncated at
#' zero (two of them eliminated by the sum of the effects), the tests
#' rejected it in 3–6% of the experiments with negative binomial counts and
#' in 5–6% with Tweedie densities; when the combination left two taxa at 30%
#' of what AD expects, they detected it in 20% and 28% of the experiments.
#'
#' **Taxa.** Each test uses the taxa found in its samples (`n_taxa`). For
#' each, `taxa` gives its statistic, its p-value adjusted for testing every
#' taxon with mvabund's step-down resampling procedure (which controls the
#' family-wise error rate within the test) and, under MU, its net
#' interaction on the log scale, as in [estimate_interactions()]:
#' \eqn{\log\mu_c - \sum_j \log\mu_j + (m - 1)\log\mu_0}{log mu_c - sum_j log
#' mu_j + (m - 1) log mu_0}, from the means of the full model (densities with
#' `effort`): positive when the combination has more than MU expects (`NA`
#' when some cell has none; with all times, the average of the times).
#'
#' **Limits.** The models are generalised linear models without random
#' effects, with the same distribution for every taxon (not the one that
#' [fit_responses()] would choose for each). Sets with a stressor applied to
#' whole groups of units (`groups` of [describe_design()]) are not tested,
#' because the groups, not the units, are then the replicates. Only MU and AD_taxa are
#' tested: they are the null models that a GLM can express (a log and an
#' identity link); dominance and the others are not. Tests under AD_taxa
#' take longer than under MU, because each taxon of each resample needs a
#' constrained fit: in the Tweedie example below, about 20 s with seven taxa
#' and 199 resamples, against 3 s under MU. Only the net
#' interaction of the stressors of a set is tested. Unlike the properties of
#' [estimate_community()], the test says whether the community deviates,
#' through any of its taxa, not how a property changes.
#'
#' @return An object of class `must_community_test`, a list with
#'   * `tests`: one row per contrast set, time (`time = "overall"` for all
#'     times together) and null model, with the columns of the set (as in
#'     [estimate_interactions()]), `null_model`, `n_taxa`,
#'     `statistic` (sum of the statistics of the taxa), `df`, `p_value`,
#'     `p_adjusted`, `status` (`"ok"` or `"not_estimable"`), `reason` and
#'     `note` (warnings of the fits, if any).
#'   * `taxa`: one row per test and taxon used: `set_id`, `time`, `taxon`,
#'     `null_model`, `scale` (`"log(response)"` under MU, `"response"` under
#'     AD_taxa), `estimate`, `statistic` and `p_adjusted`.
#'   * `design`, `settings` (the options, with the statistic of each null
#'     model, the seed and the versions) and `diagnostics`: `excluded` (the
#'     taxa left out of each test, with the null model and the reason), with
#'     `family = "tweedie"`, `power` (the power of each taxon, `NA` when it
#'     could not be estimated, and the warnings of its fit) and, for counts
#'     under AD_taxa, `size` (the size of the negative binomial of each
#'     taxon).
#'
#' @references
#' Dunn, P. K. & Smyth, G. K. (2005). Series evaluation of Tweedie
#' exponential dispersion model densities. *Statistics and Computing*, 15,
#' 267–280. \doi{10.1007/s11222-005-4070-y}
#'
#' Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
#' J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025).
#' Hypothesis-driven research on multiple stressors: An analytical framework
#' for stressor interactions. *Ecology and Evolution*, 15, e71959.
#' \doi{10.1002/ece3.71959}
#'
#' Szöcs, E., Van den Brink, P. J., Lagadic, L., Caquet, T., Roucaute, M.,
#' Auber, A., Bayona, Y., Liess, M., Ebke, P., Ippolito, A., ter Braak, C. J.
#' F., Brock, T. C. M. & Schäfer, R. B. (2015). Analysing chemical-induced
#' changes in macroinvertebrate communities in aquatic mesocosm experiments:
#' a comparison of methods. *Ecotoxicology*, 24, 760–769.
#' \doi{10.1007/s10646-015-1421-0}
#'
#' Thompson, P. L., MacLennan, M. M. & Vinebrooke, R. D. (2018). An improved
#' null model for assessing the net effects of multiple stressors on
#' communities. *Global Change Biology*, 24, 517–525.
#' \doi{10.1111/gcb.13852}
#'
#' Wang, Y., Naumann, U., Eddelbuettel, D., Wilshire, J. & Warton, D.
#' mvabund: Statistical Methods for Analysing Multivariate Abundance Data. R
#' package. <https://CRAN.R-project.org/package=mvabund>
#'
#' Warton, D. I., Thibaut, L. & Wang, Y. A. (2017). The PIT-trap—A
#' "model-free" bootstrap procedure for inference about regression models
#' with discrete, multivariate responses. *PLoS ONE*, 12(7), e0181790.
#' \doi{10.1371/journal.pone.0181790}
#'
#' @seealso [estimate_community()] for the effects on the total abundance and
#'   the richness, [estimate_interactions()] for each taxon on its own.
#'
#' @examplesIf requireNamespace("mvabund", quietly = TRUE)
#' \donttest{
#' macro <- polazzo_macroinvertebrates[polazzo_macroinvertebrates$day == 15 &
#'   polazzo_macroinvertebrates$herbicide == "absent", ]
#' design <- describe_design(macro, c("insecticide", "nutrients"),
#'   control = c(insecticide = "absent", nutrients = "ambient"), unit = "mesocosm")
#' taxa <- c("chironomini", "orthocladiinae", "physella_acuta", "cloeon", "hydrobiidae",
#'   "lumbriculidae")
#' community <- test_community(macro, taxa, design, n_boot = 199, seed = 1)
#' community
#' community$taxa
#' }
#' @examples
#' \donttest{
#' # Densities without their counts: Tweedie models.
#' zoo <- polazzo_zooplankton[polazzo_zooplankton$day == 14 &
#'   polazzo_zooplankton$herbicide == "absent", ]
#' design <- describe_design(zoo, c("insecticide", "nutrients"),
#'   control = c(insecticide = "absent", nutrients = "ambient"), unit = "mesocosm")
#' taxa <- c("cyclopoida", "nauplii", "simocephalus", "lecanidae", "ceriodafnia",
#'   "polyarthra", "daphnia")
#' zoo_test <- test_community(zoo, taxa, design, family = "tweedie",
#'   null_models = c("MU", "AD_taxa"), n_boot = 199, seed = 1)
#' zoo_test
#' zoo_test$diagnostics$power
#' }
#' @export
test_community <- function(data, taxa, design, family = "negative.binomial",
                           null_models = "MU", effort = NULL, resampling = "pit.trap",
                           n_boot = 999L, overall = TRUE, p_adjust = "BH", seed = NULL,
                           quiet = FALSE) {
  check_data_frame(data, allow_empty = FALSE)
  data <- as_plain_data_frame(data)
  check_design_object(design)
  if (!is.character(taxa) || length(taxa) < 2L || anyNA(taxa) || anyDuplicated(taxa) > 0L) {
    abort_input("`taxa` must give the columns of two or more taxa, without repetitions.")
  }
  check_columns(data, stats::setNames(taxa, rep("taxa", length(taxa))))
  spec <- design$spec
  roles <- c(spec$stressors, spec$unit, spec$time, spec$by, spec$groups)
  clash <- intersect(taxa, roles)
  if (length(clash) > 0L) {
    abort_input(sprintf("%s cannot be a taxon: it is a design column.",
                        format_items(quote_items(clash), max = Inf)))
  }
  check_choice(family, c("negative.binomial", "poisson", "tweedie"))
  null_models <- check_test_null_models(null_models)
  effort <- check_effort(data, effort, taxa, roles)
  check_choice(resampling, c("pit.trap", "perm.resid", "montecarlo"))
  own <- identical(family, "tweedie") || "AD_taxa" %in% null_models
  if (own && !identical(resampling, "pit.trap")) {
    abort_input(paste(
      'With `family = "tweedie"` or the null model "AD_taxa" the resampling is',
      '`"pit.trap"`: the other methods are those of mvabund, which only fits',
      "the multiplicative model of counts."
    ))
  }
  check_number(n_boot, lower = 99)
  if (n_boot != round(n_boot)) {
    abort_input("`n_boot` must be a whole number.")
  }
  check_flag(overall)
  check_choice(p_adjust, stats::p.adjust.methods)
  check_seed(seed)
  check_flag(quiet)
  check_taxa_values(data, taxa, family)
  problems <- design$checks[design$checks$severity == "problem", , drop = FALSE]
  if (nrow(problems) > 0L) {
    abort_input(paste(
      "The design has problems that make the models invalid; fix them first:",
      paste0("\n* ", problems$message, collapse = "")
    ))
  }
  counts <- !identical(family, "tweedie")
  uses_mvabund <- counts && "MU" %in% null_models
  if (uses_mvabund) {
    check_installed("mvabund", "test the interactions of a community")
  }
  frame <- design_model_frame(data, design)
  seed_given <- !is.null(seed)
  seed <- resolve_seed(seed)
  opts <- list(family = family, null_models = null_models, effort = effort,
               resampling = resampling, n_boot = as.integer(n_boot), overall = overall)
  run <- withr::with_seed(seed, community_tests(frame, taxa, design, opts, quiet))
  parts <- run$parts
  tests <- bind_rows_base(lapply(parts, `[[`, "test"), template = community_test_template(design))
  tests$p_adjusted <- rep(NA_real_, nrow(tests))
  for (nm in null_models) {
    for (all_times in c(FALSE, TRUE)) {
      sel <- tests$status == "ok" & is.finite(tests$p_value) & tests$null_model == nm &
        (tests$time %in% "overall") == all_times
      tests$p_adjusted[sel] <- stats::p.adjust(tests$p_value[sel], method = p_adjust)
    }
  }
  rest <- setdiff(names(tests), "p_adjusted")
  at <- match("p_value", rest)
  tests <- tests[, c(rest[seq_len(at)], "p_adjusted", rest[-seq_len(at)]), drop = FALSE]
  taxa_rows <- bind_rows_base(lapply(parts, `[[`, "taxa"), template = data.frame(
    set_id = character(), time = character(), taxon = character(), null_model = character(),
    scale = character(), estimate = numeric(), statistic = numeric(), p_adjusted = numeric(),
    stringsAsFactors = FALSE
  ))
  excluded <- bind_rows_base(lapply(parts, `[[`, "excluded"), template = data.frame(
    set_id = character(), time = character(), null_model = character(), taxon = character(),
    reason = character(), stringsAsFactors = FALSE
  ))
  rownames(tests) <- NULL
  rownames(taxa_rows) <- NULL
  rownames(excluded) <- NULL
  statistic <- vapply(null_models, function(nm) {
    if (counts && identical(nm, "MU")) {
      "sum of likelihood-ratio statistics"
    } else {
      "sum of drops in deviance scaled by the deviance of the null model"
    }
  }, character(1))
  versions <- c(if (uses_mvabund) "mvabund", if (!counts) "glmmTMB")
  structure(
    list(
      tests = tests,
      taxa = taxa_rows,
      design = design,
      settings = c(opts, list(
        taxa = taxa, p_adjust = p_adjust, seed = seed, seed_given = seed_given,
        statistic = statistic, versions = analysis_versions(versions)
      )),
      diagnostics = list(excluded = excluded, power = run$power, size = run$size)
    ),
    class = "must_community_test"
  )
}

check_test_null_models <- function(null_models, call = sys.call(-1)) {
  valid <- c("MU", "AD_taxa")
  if (is.character(null_models) && "AD" %in% null_models) {
    abort_input(paste(
      'The additive null model of test_community() is "AD_taxa": the sum of the',
      "effects of the stressors on each taxon, truncated at zero where it",
      "eliminates the taxon (Thompson et al. 2018), as in estimate_community()."
    ), call = call)
  }
  check_community_choices(null_models, valid, "null_models", call = call)
}

# Sampling effort: NULL, or one column (for all taxa) or one per taxon, with
# positive values. Returns NULL or one column name per taxon.
check_effort <- function(data, effort, taxa, roles, call = sys.call(-1)) {
  if (is.null(effort)) {
    return(NULL)
  }
  if (!is.character(effort) || !(length(effort) %in% c(1L, length(taxa))) || anyNA(effort)) {
    abort_input(
      paste(
        "`effort` must be NULL, the name of one column of `data`, or one column",
        "per taxon in the order of `taxa`."
      ),
      call = call
    )
  }
  check_columns(data, stats::setNames(effort, rep("effort", length(effort))), call = call)
  clash <- intersect(effort, c(taxa, roles))
  if (length(clash) > 0L) {
    abort_input(sprintf("%s cannot be the effort: it is a taxon or a design column.",
                        format_items(quote_items(clash), max = Inf)), call = call)
  }
  for (col in unique(effort)) {
    check_numeric_column(data, col, lower = 0, include_lower = FALSE, call = call)
  }
  rep_len(effort, length(taxa))
}

# Values of the taxa: numbers, not negative and not missing; whole numbers
# for the count families.
check_taxa_values <- function(data, taxa, family, call = sys.call(-1)) {
  for (taxon in taxa) {
    check_numeric_column(data, taxon, lower = 0, call = call)
  }
  if (identical(family, "tweedie")) {
    return(invisible(data))
  }
  for (taxon in taxa) {
    x <- data[[taxon]]
    bad <- which(abs(x - round(x)) > 1e-8)
    if (length(bad) > 0L) {
      abort_input(sprintf(
        paste(
          "Column '%s' must hold counts (whole numbers) with family = \"%s\";",
          "%s %s not. Densities can be tested from the counts they come from, with",
          "the sampling effort as `effort`, or with family = \"tweedie\"."
        ),
        taxon, family, format_rows(bad), if (length(bad) == 1L) "is" else "are"
      ), call = call)
    }
  }
  invisible(data)
}

community_test_template <- function(design) {
  out <- data.frame(context = character(), stringsAsFactors = FALSE)
  for (b in design$spec$by) {
    out[[b]] <- character()
  }
  cbind(out, data.frame(
    time = character(), set_id = character(), order = integer(), stressors = character(),
    levels = character(), background = character(), null_model = character(),
    n_taxa = integer(), statistic = numeric(), df = numeric(), p_value = numeric(),
    status = character(), reason = character(), note = character(), stringsAsFactors = FALSE
  ))
}

# The tests to run: the samples of each set at each time, and at all times.
community_plan <- function(frame, design, opts) {
  sets <- design$sets[design$sets$net_identifiable, , drop = FALSE]
  baseline <- baseline_times(design$spec)
  has_time <- !is.null(design$spec$time)
  plan <- list()
  for (i in seq_len(nrow(sets))) {
    s <- sets[i, , drop = FALSE]
    sc <- design$set_cells[[s$set_id]]
    contexts <- set_cell_contexts(sc, s)
    use <- c(which(sc$size == 0L), which(sc$size == 1L), which(sc$size == s$order))
    in_set <- strsplit(s$stressors, " x ", fixed = TRUE)[[1L]]
    in_cells <- paste_key(frame$.context, frame$.cell) %in%
      paste_key(contexts[use], sc$cell[use])
    times <- if (has_time) {
      levels(droplevels(frame$.time[in_cells]))
    } else {
      NA_character_
    }
    times <- times[!times %in% set_baseline(baseline, design, s)]
    for (tt in times) {
      rows <- which(in_cells & (if (has_time) as.character(frame$.time) == tt else TRUE))
      plan[[length(plan) + 1L]] <- list(s = s, in_set = in_set, time = tt, rows = rows,
                                        all_times = FALSE)
    }
    if (opts$overall && has_time && length(times) >= 2L) {
      rows <- which(in_cells & as.character(frame$.time) %in% times)
      plan[[length(plan) + 1L]] <- list(s = s, in_set = in_set, time = "overall", rows = rows,
                                        all_times = TRUE)
    }
  }
  plan
}

community_tests <- function(frame, taxa, design, opts, quiet) {
  plan <- community_plan(frame, design, opts)
  effort <- if (is.null(opts$effort)) {
    NULL
  } else {
    matrix(vapply(opts$effort, function(col) as.numeric(frame[[col]]), numeric(nrow(frame))),
           nrow = nrow(frame), dimnames = list(NULL, taxa))
  }
  rows <- sort(unique(unlist(lapply(plan, `[[`, "rows"), use.names = FALSE)))
  power <- NULL
  size <- NULL
  if (identical(opts$family, "tweedie")) {
    power <- tweedie_powers(frame, taxa, rows, effort)
  } else if (identical(opts$family, "negative.binomial") && "AD_taxa" %in% opts$null_models) {
    size <- nb2_sizes(frame, taxa, rows, effort)
  }
  parts <- vector("list", length(plan))
  last <- ""
  for (i in seq_along(plan)) {
    p <- plan[[i]]
    if (!quiet && !identical(p$s$set_id, last)) {
      message("Testing set ", p$s$set_id, ": ", p$s$stressors)
    }
    last <- p$s$set_id
    parts[[i]] <- community_test_one(
      frame[p$rows, , drop = FALSE], taxa, p$s, p$in_set, p$time, design, opts, p$all_times,
      if (is.null(effort)) NULL else effort[p$rows, , drop = FALSE], power, size
    )
  }
  parts <- list(
    test = lapply(parts, function(x) bind_rows_base(lapply(x, `[[`, "test"))),
    taxa = lapply(parts, function(x) bind_rows_base(lapply(x, `[[`, "taxa"))),
    excluded = lapply(parts, function(x) bind_rows_base(lapply(x, `[[`, "excluded")))
  )
  list(parts = lapply(seq_along(parts$test), function(i) {
    list(test = parts$test[[i]], taxa = parts$taxa[[i]], excluded = parts$excluded[[i]])
  }), power = power, size = size)
}

# One test of the community for each null model: the samples `d` of a set at
# one time (or at all times), the null model and the cell-means model.
community_test_one <- function(d, taxa, s, in_set, time_value, design, opts, all_times,
                               effort, power, size) {
  control <- design$spec$control
  row <- data.frame(context = s$context, stringsAsFactors = FALSE)
  for (b in design$spec$by) {
    row[[b]] <- s[[b]]
  }
  row <- cbind(row, data.frame(
    time = as.character(time_value), set_id = s$set_id, order = s$order,
    stressors = s$stressors, levels = s$levels, background = s$background,
    null_model = NA_character_, n_taxa = 0L, statistic = NA_real_, df = NA_real_,
    p_value = NA_real_, status = "not_estimable", reason = "", note = "",
    stringsAsFactors = FALSE
  ))
  not_estimable <- function(reason) {
    lapply(opts$null_models, function(nm) {
      r <- row
      r$null_model <- nm
      r$reason <- reason
      list(test = r, taxa = NULL, excluded = NULL)
    })
  }
  whole <- whole_group_problem(design, in_set)
  if (nzchar(whole)) {
    return(not_estimable(whole))
  }
  samples <- rs_samples(d, in_set, control, all_times)
  if (nzchar(samples$problem)) {
    return(not_estimable(samples$problem))
  }
  env <- samples$env
  away <- samples$away
  layout <- samples$layout
  unit <- samples$unit
  Y_all <- as.matrix(d[taxa])
  storage.mode(Y_all) <- "double"
  lapply(opts$null_models, function(nm) {
    r <- row
    r$null_model <- nm
    reason <- unname(ifelse(colSums(Y_all) > 0, NA_character_,
                            "Not found in the samples of this test."))
    if (!is.null(power)) {
      no_power <- is.na(reason) & is.na(power$power[match(taxa, power$taxon)])
      reason[no_power] <- "The power of its Tweedie distribution could not be estimated."
    }
    if (!is.null(size) && identical(nm, "AD_taxa")) {
      no_size <- is.na(reason) & is.na(size$size[match(taxa, size$taxon)])
      reason[no_size] <- "The size of its negative binomial distribution could not be estimated."
    }
    use <- is.na(reason)
    excluded <- if (any(!use)) {
      data.frame(set_id = s$set_id, time = as.character(time_value), null_model = nm,
                 taxon = taxa[!use], reason = reason[!use], stringsAsFactors = FALSE)
    } else {
      NULL
    }
    if (!any(use)) {
      r$reason <- "No taxon can be tested with the samples of this test."
      return(list(test = r, taxa = NULL, excluded = excluded))
    }
    Y <- Y_all[, use, drop = FALSE]
    E <- if (is.null(effort)) NULL else effort[, use, drop = FALSE]
    fitted <- if (identical(nm, "MU") && !identical(opts$family, "tweedie")) {
      mvabund_community_test(Y, E, env, away, unit, all_times, opts)
    } else {
      fams <- lapply(colnames(Y), function(taxon) {
        switch(opts$family,
               tweedie = rs_family("tweedie", power = power$power[power$taxon == taxon]),
               poisson = rs_family("poisson"),
               rs_family("nbinom2", theta = size$size[size$taxon == taxon]))
      })
      tryCatch(
        rs_test(Y, E, layout, unit, fams, nm, opts$n_boot),
        error = function(e) {
          list(error = paste("The models could not be fitted or tested:", conditionMessage(e)))
        }
      )
    }
    if (!is.null(fitted$error)) {
      r$reason <- fitted$error
      return(list(test = r, taxa = NULL, excluded = excluded))
    }
    r$n_taxa <- ncol(Y)
    r$statistic <- fitted$statistic
    r$df <- fitted$df
    r$p_value <- fitted$p_value
    r$status <- "ok"
    r$note <- paste(unique(fitted$warnings), collapse = "; ")
    estimate <- if (identical(nm, "MU")) {
      mu_net_deviation(fitted$rate, Y, env, length(in_set), all_times)
    } else {
      ad_net_deviation(fitted$rate, layout, truncate = TRUE)
    }
    taxa_rows <- data.frame(
      set_id = s$set_id, time = as.character(time_value), taxon = colnames(Y), null_model = nm,
      scale = if (identical(nm, "MU")) null_model_scale("MU") else "response",
      estimate = unname(estimate), statistic = unname(fitted$taxa_statistic),
      p_adjusted = unname(fitted$taxa_p), stringsAsFactors = FALSE
    )
    list(test = r, taxa = taxa_rows, excluded = excluded)
  })
}

# Counts under MU: both models with mvabund::manyglm() and the test with
# mvabund::anova.manyglm(), resampling whole units.
mvabund_community_test <- function(Y, effort, env, away, unit, all_times, opts) {
  main <- paste(away, collapse = " + ")
  if (all_times) {
    null_formula <- stats::as.formula(paste("Y ~ .time * (", main, ")"))
    full_formula <- Y ~ .time * cell
  } else {
    null_formula <- stats::as.formula(paste("Y ~", main))
    full_formula <- Y ~ cell
  }
  offset <- if (is.null(effort)) NULL else log(effort)
  warnings <- character()
  keep_warning <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  out <- tryCatch(withCallingHandlers({
    Y <- mvabund::mvabund(Y)
    if (is.null(offset)) {
      null_fit <- mvabund::manyglm(null_formula, data = env, family = opts$family)
      full_fit <- mvabund::manyglm(full_formula, data = env, family = opts$family)
    } else {
      null_fit <- mvabund::manyglm(null_formula, data = env, family = opts$family,
                                   offset = offset)
      full_fit <- mvabund::manyglm(full_formula, data = env, family = opts$family,
                                   offset = offset)
    }
    anova_fun <- stats::anova
    test <- anova_fun(null_fit, full_fit, resamp = opts$resampling, test = "LR",
                      p.uni = "adjusted", nBoot = opts$n_boot, block = unit,
                      show.time = "none")
    list(full = full_fit, test = test)
  }, warning = keep_warning), error = function(e) e)
  if (inherits(out, "error")) {
    return(list(error = paste("mvabund could not fit or test the models:",
                              conditionMessage(out))))
  }
  tab <- out$test$table
  rate <- as.matrix(out$full$fitted.values)
  if (!is.null(effort)) {
    rate <- rate / effort
  }
  list(statistic = tab[2L, 3L], df = tab[2L, 2L], p_value = tab[2L, 4L],
       taxa_statistic = out$test$uni.test[2L, ], taxa_p = out$test$uni.p[2L, ], rate = rate,
       warnings = warnings)
}

# Net interaction of each taxon under MU on the log scale, from the means of
# the full model (`rate`, one row per sample: densities with an effort):
# log mu_c - sum_j log mu_j + (m - 1) log mu_0, the average over times with
# all times. NA where a cell has nothing.
mu_net_deviation <- function(rate, Y, env, m, all_times) {
  away <- vapply(seq_len(m), function(j) env[[paste0(".s", j)]] == "TRUE",
                 logical(nrow(env)))
  away <- matrix(away, nrow = nrow(env))
  n_away <- rowSums(away)
  weight <- ifelse(n_away == m, 1, ifelse(n_away == 1L, -1, ifelse(n_away == 0L, m - 1, 0)))
  groups <- if (all_times) split(seq_len(nrow(env)), env$.time) else list(seq_len(nrow(env)))
  per_time <- vapply(groups, function(rows) {
    cells <- split(rows, droplevels(env$cell[rows]))
    first <- vapply(cells, `[[`, integer(1), 1L)
    totals <- vapply(cells, function(r) colSums(Y[r, , drop = FALSE]), numeric(ncol(Y)))
    totals <- matrix(totals, nrow = ncol(Y))
    logs <- log(rate[first, , drop = FALSE])
    est <- colSums(weight[first] * logs)
    est[rowSums(totals == 0) > 0] <- NA_real_
    est
  }, numeric(ncol(Y)))
  per_time <- matrix(per_time, nrow = ncol(Y))
  rowMeans(per_time)
}

# Shape parameters of the taxa ------------------------------------------------

# Power of each taxon by maximum likelihood (glmmTMB), from one mean per
# treatment combination and time in the samples `rows`, leaving out the
# combinations and times without the taxon.
tweedie_powers <- function(frame, taxa, rows, effort) {
  group <- droplevels(factor(paste(paste_key(frame$.context, frame$.cell),
                                   as.character(frame$.time))[rows]))
  out <- data.frame(taxon = taxa, power = NA_real_, note = "", stringsAsFactors = FALSE)
  for (j in seq_along(taxa)) {
    y <- as.numeric(frame[[taxa[[j]]]][rows])
    total <- tapply(y, group, sum)
    keep <- as.character(group) %in% names(total)[total > 0]
    if (sum(y[keep] > 0) < 2L) {
      out$note[[j]] <- "Found in fewer than two samples."
      next
    }
    dat <- data.frame(y = y[keep], g = droplevels(group[keep]),
                      off = if (is.null(effort)) 0 else log(effort[rows[keep], j]))
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(
      glmmTMB::glmmTMB(y ~ 0 + g + offset(off), family = glmmTMB::tweedie(), data = dat,
                       se = FALSE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ), error = function(e) e)
    if (inherits(fit, "error")) {
      out$note[[j]] <- conditionMessage(fit)
      next
    }
    p <- tryCatch(unname(glmmTMB::family_params(fit)[[1L]]), error = function(e) NA_real_)
    if (is.finite(p) && p > 1 && p < 2) {
      out$power[[j]] <- p
    }
    out$note[[j]] <- paste(unique(warnings), collapse = "; ")
  }
  out
}

# Size (theta) of the negative binomial of each taxon by maximum likelihood,
# with the means of every treatment combination and time at their
# estimates; limited to 1e8 (Poisson counts).
nb2_sizes <- function(frame, taxa, rows, effort) {
  group <- droplevels(factor(paste(paste_key(frame$.context, frame$.cell),
                                   as.character(frame$.time))[rows]))
  out <- data.frame(taxon = taxa, size = NA_real_, stringsAsFactors = FALSE)
  for (j in seq_along(taxa)) {
    y <- as.numeric(frame[[taxa[[j]]]][rows])
    if (sum(y > 0) < 2L) {
      next
    }
    e <- if (is.null(effort)) rep(1, length(y)) else effort[rows, j]
    out$size[[j]] <- nb2_size(y, e, group)
  }
  out
}

nb2_size <- function(y, e, group) {
  loglik <- function(log_theta) {
    theta <- exp(log_theta)
    rate <- rs_rates(y, e, group, rs_family("nbinom2", theta = theta))[as.integer(group)]
    sum(stats::dnbinom(y, size = theta, mu = e * rate, log = TRUE))
  }
  opt <- stats::optimize(loglik, c(log(1e-4), log(1e8)), maximum = TRUE)
  exp(opt$maximum)
}

#' @export
print.must_community_test <- function(x, n = 12L, ...) {
  s <- x$settings
  models <- c(MU = "multiplicative (MU)", AD_taxa = "additive truncated at zero (AD_taxa)")
  cat(sprintf(
    "<must_community_test> %d taxa; null model%s tested on every taxon: %s\n",
    length(s$taxa), if (length(s$null_models) > 1L) "s" else "",
    paste(models[s$null_models], collapse = ", ")
  ))
  cat(sprintf(
    "Models: %s GLM for each taxon (%s)%s; %s resampling of whole units (%d); p-values adjusted with %s; seed %s\n",
    s$family,
    paste(c(MU = "log link for MU", AD_taxa = "identity link for AD_taxa")[s$null_models],
          collapse = ", "),
    if (is.null(s$effort)) "" else ", offset log(effort)",
    s$resampling, s$n_boot, s$p_adjust, format(s$seed)
  ))
  cat("Statistic: ", paste(sprintf("%s, %s", names(s$statistic), s$statistic), collapse = "; "),
      "\n", sep = "")
  tests <- x$tests
  if (nrow(tests) == 0L) {
    cat("No test could be done.\n")
    return(invisible(x))
  }
  deviating <- vapply(seq_len(nrow(tests)), function(i) {
    rows <- x$taxa$set_id == tests$set_id[[i]] & x$taxa$time %in% tests$time[[i]] &
      x$taxa$null_model == tests$null_model[[i]]
    sum(x$taxa$p_adjusted[rows] < 0.05, na.rm = TRUE)
  }, integer(1))
  show <- data.frame(
    time = tests$time,
    combination = set_combination(tests$levels, tests$stressors),
    background = compact_levels(tests$background, x$design$spec$control),
    model = tests$null_model,
    n_taxa = tests$n_taxa,
    statistic = signif(tests$statistic, 3),
    df = tests$df,
    p_adj = signif(tests$p_adjusted, 2),
    taxa_deviating = deviating,
    status = tests$status,
    stringsAsFactors = FALSE
  )
  if (all(is.na(show$time))) {
    show$time <- NULL
  }
  if (!any(nzchar(show$background))) {
    show$background <- NULL
  }
  if (length(s$null_models) == 1L) {
    show$model <- NULL
  }
  cat(sprintf("Tests of the community (%d rows%s):\n", nrow(show),
              if (nrow(show) > n) sprintf(", first %d shown", n) else ""))
  print(utils::head(show, n), row.names = FALSE)
  cat("taxa_deviating: taxa whose adjusted p-value (within the test) is below 0.05; see `$taxa`.\n")
  invisible(x)
}
