# Community properties (Thompson et al. 2018) -----------------------------------

#' Estimate non-additive effects on community properties
#'
#' @description
#' Estimates the effects of the stressors, alone and combined, on two
#' properties of a community, the total abundance (or biomass) of its taxa
#' and their richness, and whether the combined effect deviates from four
#' null models: the compositional null model of Thompson et al. (2018),
#' which adds the effects of the stressors on each taxon and computes the
#' property of the community it predicts, its multiplicative version, and
#' the additive and multiplicative null models applied to the property
#' itself. The effects on each taxon come from the models that
#' [fit_responses()] fits to each one, generalised linear (mixed) models of
#' the counts, densities or biomasses of each population, as Szöcs and
#' Schäfer (2015) recommend for counts in ecotoxicology and as the
#' multivariate generalised linear models compared by Szöcs et al. (2015)
#' do for the communities of mesocosm experiments.
#'
#' @param fit A `must_fit` object from [fit_responses()] with one response
#'   per taxon (counts, densities or biomasses; not proportions).
#' @param taxa `NULL` (default) to use every response of `fit`, or the names
#'   of the responses that make up the community (two or more). Every taxon
#'   needs a fitted model.
#' @param properties Properties of the community: one or both of
#'   `"abundance"` (the sum of the mean abundances of the taxa) and
#'   `"richness"` (the expected number of taxa found in a sample). Default
#'   both.
#' @param null_models One or more of `"AD_taxa"` (compositional, additive
#'   on each taxon), `"MU_taxa"` (compositional, multiplicative on each
#'   taxon), `"AD"` (additive on the property) and `"MU"` (multiplicative on
#'   the property). Default all four.
#' @param inference `"simulation"` (default; Krinsky-Robb draws of the
#'   fitted models, with the taxa linked as `taxa_correlation` says) or
#'   `"resampling"` (resampling whole experimental units with all their
#'   taxa; less well calibrated with few units per combination, see
#'   "Resampling").
#' @param taxa_correlation How the estimates of different taxa are linked
#'   in the simulation: `"residual"` (default), through the correlation of
#'   their residuals shrunk towards zero (see Details), or `"none"`
#'   (independent taxa). Not used with `inference = "resampling"`.
#' @param correction `"small_sample"` (default) or `"none"`, as in
#'   [estimate_interactions()].
#' @param n_sim Number of simulation draws (at least 100). Default 10000.
#' @param n_boot Number of resamples with `inference = "resampling"` (at
#'   least 99). Default 999.
#' @param level Confidence level of the intervals. Default 0.95.
#' @param p_adjust Method of [stats::p.adjust()] for the p-values, applied
#'   within each property, null model and component across sets and times
#'   (the averages over times apart). Default `"BH"`.
#' @param overall If `TRUE` (default), also estimates the average of each
#'   effect over the sampling times of each set.
#' @param seed `NULL` or a whole number for the simulation; the seed used is
#'   kept in `settings$seed`.
#'
#' @details
#' **Properties.** In a treatment combination (at a sampling time), let
#' \eqn{x_i} be the mean of taxon \eqn{i} given by its model (with the
#' random effects at zero; zero in a combination where the taxon had only
#' zeros, whose mean is fixed at the boundary by [fit_responses()]). The
#' total abundance is \eqn{R = \sum_i x_i}{R = sum(x_i)}. The richness is
#' the expected number of taxa found in a sample, \eqn{R = \sum_i P(Y_i > 0
#' \mid x_i)}{R = sum(P(Y_i > 0 | x_i))}, with the probability of absence of
#' the distribution fitted to each taxon and the dispersion fitted for that
#' combination: Poisson, negative binomial (NB1 or NB2) or Tweedie. Richness
#' is not defined when some taxon was fitted with a distribution without
#' zeros (Gaussian, Gamma, lognormal).
#'
#' **Effects and null models.** The effects of each stressor alone
#' (`component = "single"`) and of their combination (`"combined"`) are the
#' relative changes of the property from the reference cell, \eqn{\Delta R
#' = R / R_0 - 1} (Thompson et al. 2018, eq. 2). The null models expect for
#' the combination (`"expected"`), from the cells of the stressors alone:
#' * `"AD_taxa"`: the community predicted by adding the effects of the
#'   stressors on each taxon, \eqn{\hat{x}_i = \max(0, x_{i0} + \sum_j
#'   (x_{ij} - x_{i0}))}{x_i = max(0, x_i0 + sum_j (x_ij - x_i0))} (eq. 3),
#'   with its property relative to the reference.
#' * `"MU_taxa"`: the same with multiplicative effects on each taxon,
#'   \eqn{\hat{x}_i = x_{i0} \prod_j x_{ij} / x_{i0}}{x_i = x_i0 prod_j
#'   (x_ij / x_i0)}, a model that Thompson et al. (2018) suggest for species
#'   whose responses multiply. It is not defined when a taxon with only
#'   zeros in the reference cell appears with some stressor alone; a taxon
#'   absent from the reference cell and from every stressor alone is
#'   expected to be absent.
#' * `"AD"`: \eqn{\Delta\hat{R} = \max(-1, \sum_j \Delta R_j)}{max(-1, sum_j
#'   dR_j)} (eq. 4, limited to the loss of the whole property, as Thompson
#'   et al. did).
#' * `"MU"`: \eqn{\Delta\hat{R} = \prod_j (1 + \Delta R_j) - 1}{prod_j (1 +
#'   dR_j) - 1} (eq. 5).
#'
#' The interaction is the observed minus the expected combined effect,
#' \eqn{\Delta R_c - \Delta\hat{R}_c}{dR_c - dR_expected}, on the same scale
#' (`scale = "response / control"`): positive when the combination keeps
#' more of the property than the null model expects. Thompson et al. (2018)
#' measured the error of a prediction as \eqn{|\Delta R_c| -
#' |\Delta\hat{R}_c|}, which can be computed from the rows `combined` and
#' `expected`. The interactions can be classified with
#' [classify_interactions()] and the null models compared with
#' [compare_null_models()].
#'
#' **Uncertainty.** With `inference = "simulation"`, the intervals and
#' p-values come from Krinsky-Robb simulation: draws of the fixed effects of
#' each taxon's model (a
#' multivariate t with the degrees of freedom of that model with
#' `correction = "small_sample"`, as in [estimate_interactions()]), from
#' which the means of the taxa, the properties and their effects are
#' computed draw by draw. Taxa that vary together, such as taxa that are all
#' more abundant in some units, have correlated estimates, and that
#' correlation changes the uncertainty of their total. With
#' `taxa_correlation = "residual"` the draws of the taxa are linked by a
#' Gaussian copula with the correlation of their residuals: the deviations
#' of their values from the mean of their treatment combination at each
#' time, standardised within that combination. The correlation is shrunk
#' towards zero with the intensity of Schäfer and Strimmer (2005), whose
#' variances of the correlations are computed from the sums of the products
#' of the residuals within each experimental unit, so that the repeated
#' observations of a unit are not taken as independent. The dispersion
#' parameters (and the power of the Tweedie) are fixed at their estimates.
#' The intervals are percentile intervals of the draws, and the p-values
#' twice the smaller proportion of draws on either side of zero. The
#' expected combined effect has an interval but is not tested.
#'
#' **Resampling.** With `inference = "resampling"` the taxa are not linked
#' through a model of their correlation: whole experimental units are
#' resampled with all their taxa (and all their times), so that whatever
#' makes the taxa of a unit vary together is kept. The mean of each taxon in
#' each treatment combination and time is the mean of its samples (the
#' maximum-likelihood mean of a GLM of its distribution without random
#' effects; with an offset, per unit of effort), which can differ from the
#' conditional means of the fitted models when units vary. The resamples
#' come from the PIT-trap (Warton et al. 2017) under these means: the
#' probability integral transform residuals of each unit are resampled
#' together, with replacement, and turned back into values of each taxon
#' with the quantile function of its distribution (the family of its fitted
#' model, with the dispersion of the cell means); the means, the properties
#' and their effects are computed again in each resample. The richness uses
#' the dispersion parameters of the fitted models. The intervals are
#' percentile intervals of the resamples and the p-values twice the smaller
#' proportion of them on either side of zero, as with the simulation. It
#' needs every unit sampled the same number of times at every sampling time,
#' distributions with a quantile function (Poisson, negative binomial,
#' Gaussian, gamma, lognormal and Tweedie), and no stressor of the set
#' applied to whole groups of units (whose groups are its replicates). With
#' few units per combination the resamples vary less than the data (the
#' samples deviate less from the means fitted to them than from the true
#' means): in simulated experiments with four units per combination and six
#' taxa, the resampling tests of the total abundance rejected true null
#' models in 11–17% of the experiments, and its intervals covered the true
#' combined effect in 83–93% of them, against 4–8% and 87–99% with the
#' simulation, which therefore stays the default.
#'
#' **Richness under `"AD_taxa"`.** A taxon whose expected mean is close to
#' zero, compared with its uncertainty, is expected to persist in some draws
#' and to be eliminated in others, and its expected presence jumps between
#' zero and nearly one. The expected richness is then biased, and in
#' simulated experiments with four units per combination its tests rejected
#' true null models in about one experiment in ten (up to one in five) at
#' the 5% level. When the sum of the effects of the stressors on some taxon
#' is below zero in 20% to 80% of the draws, the expected richness and its
#' interaction under `"AD_taxa"` are therefore not estimated (`status =
#' "not_estimable"`, with the taxa in `reason`); in the other experiments the
#' tests kept their nominal error (the same rule applies to the resamples).
#' With few units per combination and rare
#' taxa this is frequent. The total abundance, the multiplicative
#' compositional model and the null models of the property are not
#' affected.
#'
#' **Limits.** The compositional additive model fails when a stressor that
#' favours a taxon is combined with one that eliminates it alone: the
#' expected effect of the latter cannot exceed the abundance of the taxon in
#' the reference cell, while that of the former has no limit, so the loss of
#' the taxon is not expected (Thompson et al. 2018, Fig. 3). The null models
#' of the property ignore which taxa each stressor affects: taxa lost to
#' either stressor are counted twice, taxa that survive each stressor alone
#' but not both are missed, and opposite effects on the same taxa are mixed,
#' so that additive effects on every taxon can appear as interactions
#' (Thompson et al. 2018, Fig. 1). Properties that aggregate the community
#' detect fewer effects than analyses of all its taxa: in eleven mesocosm
#' studies, tests of the total abundance, the richness and a diversity index
#' missed about half of the effects that principal response curves detected
#' (Szöcs et al. 2015); [estimate_interactions()] gives the effects on each
#' taxon. The copula only links the taxa through the correlation of their
#' residuals, shrunk towards zero: in simulated experiments whose taxa varied
#' strongly together, it halved the false deviations of the total abundance
#' (8% instead of 16% at the 5% level) but its intervals still covered the
#' true combined effect less often than they should (82% instead of 95%). The
#' uncertainty of the dispersion parameters is not propagated.
#'
#' @return An object of class `must_community`, which inherits from
#'   `must_interactions`: a list with
#'   * `estimates`: a data frame with the columns of the estimates of
#'     [estimate_interactions()], with one row per property (column
#'     `response`), contrast set, time, null model and component (`single`,
#'     `expected`, `combined` and `interaction`, estimand `"net"`). For the
#'     interactions, `observed` and `expected` are the property in the
#'     combination and the property expected there by the null model; `df`
#'     is the smallest number of degrees of freedom among the taxa (`NA`
#'     with resampling).
#'   * `design`: the design of `fit`.
#'   * `settings`: the options used, with `taxa`, `families` (family and
#'     link of each taxon), `inference` and `n_boot`, `shrinkage` (the
#'     intensity applied to the correlation of the taxa, or `NA`), `seed` and
#'     `versions`.
#'   * `diagnostics`: `correlation` (the correlation used between the
#'     taxa), `correlation_raw` (before shrinkage), `n_obs` and `n_units`
#'     (observations and experimental units used to estimate it), `richness`
#'     (`""` or why richness is not defined) and `failed_responses` (empty:
#'     every taxon has a fitted model).
#'
#' @references
#' Thompson, P. L., MacLennan, M. M. & Vinebrooke, R. D. (2018). An improved
#' null model for assessing the net effects of multiple stressors on
#' communities. *Global Change Biology*, 24, 517–525.
#' \doi{10.1111/gcb.13852}
#'
#' Szöcs, E. & Schäfer, R. B. (2015). Ecotoxicology is not normal: a
#' comparison of statistical approaches for analysis of count and proportion
#' data in ecotoxicology. *Environmental Science and Pollution Research*,
#' 22, 13990–13999. \doi{10.1007/s11356-015-4579-3}
#'
#' Szöcs, E., Van den Brink, P. J., Lagadic, L., Caquet, T., Roucaute, M.,
#' Auber, A., Bayona, Y., Liess, M., Ebke, P., Ippolito, A., ter Braak, C. J.
#' F., Brock, T. C. M. & Schäfer, R. B. (2015). Analysing chemical-induced
#' changes in macroinvertebrate communities in aquatic mesocosm experiments:
#' a comparison of methods. *Ecotoxicology*, 24, 760–769.
#' \doi{10.1007/s10646-015-1421-0}
#'
#' Schäfer, J. & Strimmer, K. (2005). A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.
#' *Statistical Applications in Genetics and Molecular Biology*, 4, 32.
#' \doi{10.2202/1544-6115.1175}
#'
#' Warton, D. I., Thibaut, L. & Wang, Y. A. (2017). The PIT-trap—A
#' "model-free" bootstrap procedure for inference about regression models
#' with discrete, multivariate responses. *PLoS ONE*, 12(7), e0181790.
#' \doi{10.1371/journal.pone.0181790}
#'
#' @seealso [estimate_interactions()] for the effects on each taxon,
#'   [classify_interactions()], [compare_null_models()].
#'
#' @examples
#' \donttest{
#' # Mesocosms without the herbicide, sampled one and two weeks after the
#' # insecticide was applied
#' zoo <- polazzo_zooplankton[polazzo_zooplankton$herbicide == "absent" &
#'   polazzo_zooplankton$day %in% c(7, 14), ]
#' design <- describe_design(zoo, c("insecticide", "nutrients"),
#'   control = c(insecticide = "absent", nutrients = "ambient"),
#'   unit = "mesocosm", time = "day")
#' # One model per taxon: Tweedie models of the densities give the
#' # probability that a taxon is absent from a sample
#' fits <- fit_responses(zoo, c("daphnia", "simocephalus", "cyclopoida", "nauplii"),
#'   design, families = list(tweedie = glmmTMB::tweedie()),
#'   residual_checks = FALSE, seed = 1, quiet = TRUE)
#' community <- estimate_community(fits, n_sim = 2000, seed = 1)
#' community
#' classify_interactions(community)[, c("response", "time", "null_model", "class")]
#' # Resampling whole mesocosms instead of modelling the correlation of the taxa
#' estimate_community(fits, properties = "abundance", inference = "resampling",
#'   n_boot = 499, seed = 1)
#' }
#' @export
estimate_community <- function(fit, taxa = NULL, properties = c("abundance", "richness"),
                               null_models = c("AD_taxa", "MU_taxa", "AD", "MU"),
                               inference = "simulation", taxa_correlation = "residual",
                               correction = "small_sample", n_sim = 10000L, n_boot = 999L,
                               level = 0.95, p_adjust = "BH", overall = TRUE, seed = NULL) {
  if (!inherits(fit, "must_fit")) {
    abort_input(sprintf(
      "`fit` must be the output of fit_responses(), with one response per taxon, not %s.",
      describe_value(fit)
    ))
  }
  taxa <- check_taxa(taxa, fit)
  properties <- check_community_choices(properties, community_properties(), "properties")
  null_models <- check_community_models(null_models)
  check_choice(inference, c("simulation", "resampling"))
  resampling <- identical(inference, "resampling")
  check_choice(taxa_correlation, c("residual", "none"))
  check_choice(correction, c("small_sample", "none"))
  check_number(n_sim, lower = 100)
  if (n_sim != round(n_sim)) {
    abort_input("`n_sim` must be a whole number.")
  }
  check_number(n_boot, lower = 99)
  if (n_boot != round(n_boot)) {
    abort_input("`n_boot` must be a whole number.")
  }
  check_number(level, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  check_choice(p_adjust, stats::p.adjust.methods)
  check_flag(overall)
  check_seed(seed)
  fits <- fit$fits[taxa]
  check_community_fits(fits)
  opts <- list(null_models = null_models, properties = properties, correction = correction,
               n_sim = as.integer(n_sim), level = level, p_adjust = p_adjust,
               overall = overall, inference = inference,
               n_draws = as.integer(if (resampling) n_boot else n_sim))
  seed_given <- !is.null(seed)
  opts$seed <- resolve_seed(seed)
  info <- lapply(fits, community_taxon, correction = correction)
  corr <- if (identical(taxa_correlation, "residual") && !resampling) {
    residual_correlation(fits)
  } else {
    list(matrix = diag(length(taxa)), raw = NULL, lambda = NA_real_, n_obs = NA_integer_,
         n_units = NA_integer_)
  }
  dimnames(corr$matrix) <- list(taxa, taxa)
  if (!is.null(corr$raw)) {
    dimnames(corr$raw) <- list(taxa, taxa)
  }
  richness <- richness_problem(info)
  parts <- withr::with_seed(opts$seed, {
    boot <- if (resampling) community_bootstrap(fits, info[[1L]]$grid, n_boot) else NULL
    community_sets(info, corr$matrix, fit$design, opts, richness, boot)
  })
  estimates <- bind_rows_base(parts, template = empty_estimates(fit$design))
  estimates <- adjust_p_values(estimates, p_adjust)
  structure(
    list(
      estimates = estimates,
      design = fit$design,
      settings = c(
        opts[c("null_models", "properties", "correction", "n_sim", "level", "p_adjust",
               "overall", "seed")],
        list(
          taxa = taxa, n_boot = if (resampling) as.integer(n_boot) else NULL,
          taxa_correlation = if (resampling) "resampled" else taxa_correlation,
          shrinkage = corr$lambda, inference = inference, sensitivity = character(),
          margin = NULL,
          prediction = if (resampling) {
            "means of the samples of each treatment combination and time"
          } else {
            "conditional (random effects at zero)"
          },
          time_structure = fit$settings$time_structure,
          families = vapply(info, function(x) paste0(x$family, " (", x$link_name, ")"),
                            character(1)),
          versions = analysis_versions(), seed_given = seed_given
        )
      ),
      diagnostics = list(
        correlation = corr$matrix, correlation_raw = corr$raw, n_obs = corr$n_obs,
        n_units = corr$n_units, richness = richness,
        failed_responses = data.frame(response = character(), reason = character(),
                                      stringsAsFactors = FALSE)
      )
    ),
    class = c("must_community", "must_interactions")
  )
}

community_properties <- function() {
  c("abundance", "richness")
}

community_null_models <- function() {
  c("AD_taxa", "MU_taxa", "AD", "MU")
}

# Families whose probability of a zero gives the richness.
zero_families <- function() {
  c("poisson", "nbinom1", "nbinom2", "tweedie")
}

check_community_models <- function(null_models, call = sys.call(-1)) {
  check_community_choices(null_models, community_null_models(), "null_models", call = call)
}

# One or more options from a closed list, without repetitions, returned in
# the order of the list.
check_community_choices <- function(x, valid, arg, call = sys.call(-1)) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) || !all(x %in% valid) ||
      anyDuplicated(x) > 0L) {
    abort_input(sprintf(
      "`%s` must contain one or more of %s, without repetitions.",
      arg, format_items(quote_items(valid), max = Inf)
    ), call = call)
  }
  valid[valid %in% x]
}

check_taxa <- function(taxa, fit, call = sys.call(-1)) {
  available <- names(fit$fits)
  if (is.null(taxa)) {
    taxa <- available
  }
  if (!is.character(taxa) || anyNA(taxa) || anyDuplicated(taxa) > 0L ||
      !all(taxa %in% available)) {
    abort_input(sprintf(
      "`taxa` must be NULL or names of responses of `fit`, without repetitions (%s).",
      format_items(quote_items(available), max = 8L)
    ), call = call)
  }
  if (length(taxa) < 2L) {
    abort_input(paste(
      "A community needs two or more taxa: fit one response per taxon with",
      "fit_responses() and give them in `taxa`."
    ), call = call)
  }
  taxa
}

# Every taxon needs a fitted model of non-negative values on the same cells.
check_community_fits <- function(fits, call = sys.call(-1)) {
  failed <- names(fits)[!vapply(fits, function(f) identical(f$status, "ok"), logical(1))]
  if (length(failed) > 0L) {
    abort_input(sprintf(
      paste(
        "%s %s no fitted model, and the community needs every taxon. Fit %s again",
        "(e.g. with other `families`), or leave %s out of `taxa`."
      ),
      format_items(quote_items(failed), max = 8L), if (length(failed) == 1L) "has" else "have",
      if (length(failed) == 1L) "it" else "them", if (length(failed) == 1L) "it" else "them"
    ), call = call)
  }
  binomial <- names(fits)[vapply(fits, function(f) identical(f$kind, "binomial"), logical(1))]
  if (length(binomial) > 0L) {
    abort_input(sprintf(
      paste(
        "%s %s proportions (`trials`): the properties of the community need counts,",
        "densities or biomasses of each taxon."
      ),
      format_items(quote_items(binomial), max = 8L), if (length(binomial) == 1L) "is" else "are"
    ), call = call)
  }
  negative <- names(fits)[!vapply(fits, function(f) isTRUE(f$non_negative), logical(1))]
  if (length(negative) > 0L) {
    abort_input(sprintf(
      "%s %s negative values: abundances cannot be negative.",
      format_items(quote_items(negative), max = 8L), if (length(negative) == 1L) "has" else "have"
    ), call = call)
  }
  first <- fits[[1L]]
  same <- vapply(fits, function(f) {
    identical(as.character(f$grid$.cell_time), as.character(first$grid$.cell_time)) &&
      identical(f$fixed, first$fixed)
  }, logical(1))
  if (!all(same)) {
    abort_input(sprintf(
      paste(
        "The taxa must be observed in the same treatment combinations and times;",
        "%s %s not observed in the same as %s."
      ),
      format_items(quote_items(names(fits)[!same]), max = 8L),
      if (sum(!same) == 1L) "was" else "were", quote_items(names(fits)[[1L]])
    ), call = call)
  }
  invisible(fits)
}

# What the simulation needs of each taxon's model.
community_taxon <- function(f, correction) {
  model <- if (identical(correction, "small_sample") && !is.null(f$model_reml)) {
    f$model_reml
  } else {
    f$model
  }
  beta <- glmmTMB::fixef(model)$cond
  V <- conditional_covariance(model, beta)
  X <- stats::model.matrix(stats::as.formula(paste("~", f$fixed)), f$grid)
  if (!identical(colnames(X), names(beta))) {
    abort_internal("the prediction matrix does not match the fixed effects of the model.")
  }
  free <- !is.na(diag(V))
  V[!free, ] <- 0
  V[, !free] <- 0
  root <- tryCatch(chol(V[free, free, drop = FALSE]), error = function(e) NULL)
  fam <- stats::family(model)
  family <- fam$family
  cell_time <- as.character(f$grid$.cell_time)
  list(
    response = f$response, family = family, link = fam, link_name = fam$link,
    X = X, eta = drop(X %*% beta), free = free, root = root,
    var_eta = rowSums((X %*% V) * X),
    at_boundary = cell_time %in% names(f$boundary),
    disp = cell_dispersion(f, model, family, cell_time),
    power = if (identical(family, "tweedie")) {
      unname(glmmTMB::family_params(model)[[1L]])
    } else {
      NA_real_
    },
    df_info = f$df_info,
    units = as.numeric(f$df_info$cell_units[cell_time]),
    grid = f$grid
  )
}

# Dispersion parameter of each cell of the grid (NA for families without
# one), in the parameterisation of sigma().
cell_dispersion <- function(f, model, family, cell_time) {
  if (!family %in% c("nbinom1", "nbinom2", "tweedie")) {
    return(rep(NA_real_, length(cell_time)))
  }
  obs <- tryCatch(stats::predict(model, type = "disp"), error = function(e) NULL)
  rows <- observation_table(f)
  if (is.null(obs) || length(obs) != nrow(rows)) {
    return(rep(stats::sigma(model), length(cell_time)))
  }
  by_cell <- tapply(obs, rows$cell_time, `[`, 1L)
  out <- unname(by_cell[cell_time])
  out[is.na(out)] <- stats::sigma(model)
  out
}

# Data row, experimental unit and treatment combination (and time) of each
# observation of a taxon's model. Fits made before `observations` was kept
# take them from the model frame, in its order.
observation_table <- function(f) {
  if (all(c("row", "unit", "cell_time") %in% names(f$observations))) {
    return(f$observations)
  }
  frame <- f$model$frame
  n <- nrow(frame)
  cell <- if (".cell_time" %in% names(frame)) frame$.cell_time else frame$.group
  unit <- if (".unit" %in% names(frame)) frame$.unit else seq_len(n)
  data.frame(row = seq_len(n), unit = as.character(unit), cell_time = as.character(cell),
             stringsAsFactors = FALSE)
}

# Why richness cannot be computed ("" if it can).
richness_problem <- function(info) {
  bad <- names(info)[!vapply(info, function(x) x$family %in% zero_families(), logical(1))]
  if (length(bad) == 0L) {
    return("")
  }
  sprintf(
    paste(
      "Richness needs the probability that each taxon is absent from a sample, which",
      "the distribution fitted to %s (%s) does not give. Fit %s with a distribution",
      "that has zeros (e.g. `families = list(tweedie = glmmTMB::tweedie())` for",
      "densities with zeros), or leave %s out of `taxa`."
    ),
    format_items(quote_items(bad), max = 5L),
    paste(unique(vapply(info[bad], `[[`, character(1), "family")), collapse = ", "),
    if (length(bad) == 1L) "it" else "them", if (length(bad) == 1L) "it" else "them"
  )
}

# Probability that a taxon is present in a sample, from its mean (a vector).
presence_probability <- function(mu, family, disp, power) {
  1 - prob_boundary_obs(family, mu, n = NA, disp = disp, power = power, type = "zero")
}

# Correlation of the taxa ------------------------------------------------------
#
# Residuals of each taxon: deviations of its values from the mean of their
# treatment combination at each time, divided by the standard deviation of
# that combination (zero where it has no variation). They are paired between
# taxa by data row; their correlation (with the residuals scaled to a mean
# square of one) is shrunk towards zero with the intensity of Schäfer &
# Strimmer (2005, eq. 10 and Appendix A), with the variance of each
# correlation from the sums of products within each experimental unit
# (clusters).
residual_correlation <- function(fits) {
  obs <- lapply(fits, taxon_residuals)
  rows <- Reduce(intersect, lapply(obs, `[[`, "row"))
  n_taxa <- length(fits)
  if (length(rows) < 3L) {
    return(list(matrix = diag(n_taxa), raw = NULL, lambda = 1, n_obs = length(rows),
                n_units = NA_integer_))
  }
  D <- vapply(obs, function(o) o$residual[match(rows, o$row)], numeric(length(rows)))
  D <- matrix(D, nrow = length(rows))
  unit <- obs[[1L]]$unit[match(rows, obs[[1L]]$row)]
  shrink_correlation(D, unit)
}

taxon_residuals <- function(f) {
  frame <- f$model$frame
  y <- as.numeric(frame[[f$response]])
  offset <- grep("^offset\\(", names(frame), value = TRUE)
  if (length(offset) > 0L) {
    y <- y / exp(frame[[offset[[1L]]]])
  }
  rows <- observation_table(f)
  g <- rows$cell_time
  centred <- y - stats::ave(y, g)
  n <- stats::ave(y, g, FUN = length)
  s <- sqrt(stats::ave(centred^2, g, FUN = sum) / pmax(n - 1, 1))
  residual <- ifelse(n > 1 & s > 0, centred / s, 0)
  list(row = rows$row, unit = rows$unit, residual = residual)
}

# Shrinkage of the correlation of the columns of D (residuals, centred within
# cells) towards zero; `unit` gives the cluster of each row.
shrink_correlation <- function(D, unit) {
  n <- nrow(D)
  scale <- sqrt(colMeans(D^2))
  X <- sweep(D, 2L, ifelse(scale > 0, scale, 1), "/")
  R <- crossprod(X) / n
  diag(R) <- 1
  units <- split(seq_len(n), unit)
  n_units <- length(units)
  off <- upper.tri(R)
  lambda <- 1
  if (n_units >= 2L && sum(R[off]^2) > 0) {
    var_r <- matrix(0, ncol(D), ncol(D))
    for (u in units) {
      W <- crossprod(X[u, , drop = FALSE])
      var_r <- var_r + (W - length(u) * R)^2
    }
    var_r <- var_r * n_units / (n_units - 1) / n^2
    lambda <- min(1, max(0, sum(var_r[off]) / sum(R[off]^2)))
  }
  shrunk <- (1 - lambda) * R
  diag(shrunk) <- 1
  list(matrix = shrunk, raw = R, lambda = lambda, n_obs = n, n_units = n_units)
}

# Upper factor C of a correlation matrix (t(C) %*% C = R), by Cholesky or,
# for a singular matrix, from its eigen decomposition.
correlation_root <- function(R) {
  root <- tryCatch(chol(R), error = function(e) NULL)
  if (!is.null(root)) {
    return(root)
  }
  e <- eigen(R, symmetric = TRUE)
  t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow = length(e$values)))
}

# Contrast sets -----------------------------------------------------------------

community_sets <- function(info, R, design, opts, richness, boot = NULL) {
  grid <- info[[1L]]$grid
  key <- paste_key(grid$.context, grid$.cell, grid$.time)
  sets <- design$sets[design$sets$net_identifiable, , drop = FALSE]
  baseline <- baseline_times(design$spec)
  root <- if (all(R[upper.tri(R)] == 0)) NULL else correlation_root(R)
  rows <- list()
  for (i in seq_len(nrow(sets))) {
    s <- sets[i, , drop = FALSE]
    sc <- design$set_cells[[s$set_id]]
    times <- unique(as.character(grid$.time[grid$.context == s$context]))
    times <- times[!times %in% set_baseline(baseline, design, s)]
    if (length(times) == 0L) {
      next
    }
    rows[[length(rows) + 1L]] <- community_set(s, sc, times, info, root, design, opts, key,
                                               richness, boot)
  }
  rows
}

# Estimates of one contrast set at each of its times (and their average).
community_set <- function(s, sc, times, info, root, design, opts, key, richness, boot = NULL) {
  m <- s$order
  in_set <- strsplit(s$stressors, " x ", fixed = TRUE)[[1L]]
  contexts <- set_cell_contexts(sc, s)
  # Cells: the reference, each stressor alone and the full combination.
  use <- c(which(sc$size == 0L), which(sc$size == 1L), which(sc$size == m))
  labels <- sc$cell[use]
  idx <- vapply(times, function(tt) match(paste_key(contexts[use], sc$cell[use], tt), key),
                integer(length(use)))
  idx <- matrix(idx, nrow = length(use))
  complete <- colSums(is.na(idx)) == 0L
  cols <- sort(unique(idx[!is.na(idx)]))
  # Degrees of freedom of each taxon (none with resampling).
  dfs <- if (!is.null(boot)) rep(NA_real_, length(info)) else vapply(info, function(x) {
    df <- set_df(x$df_info, s, opts$correction, contexts = unique(contexts))
    if (is.na(df)) {
      keep <- intersect(cols, which(!x$at_boundary))
      df <- welch_df(x$var_eta[keep], x$units[keep])
    }
    df
  }, numeric(1))
  problem <- ""
  no_root <- names(info)[vapply(info, function(x) is.null(x$root), logical(1))]
  if (!is.null(boot)) {
    problem <- boot$problem
    if (!nzchar(problem)) {
      problem <- whole_group_problem(design, in_set)
    }
  } else if (length(no_root) > 0L) {
    problem <- sprintf(
      "The covariance matrix of the fixed effects of %s is not positive definite.",
      format_items(quote_items(no_root), max = 5L)
    )
  } else if (any(is.na(dfs) | dfs < 1)) {
    problem <- sprintf(
      paste(
        "No degrees of freedom are left to estimate the variance of %s (too few",
        "replicated units)."
      ),
      format_items(quote_items(names(info)[is.na(dfs) | dfs < 1]), max = 5L)
    )
  }
  layout <- list(m = m, idx = idx, cols = cols, complete = complete,
                 pos = matrix(match(idx, cols), nrow = nrow(idx)))
  # Richness is only computed when every taxon has a probability of absence.
  computed <- if (nzchar(richness)) setdiff(opts$properties, "richness") else opts$properties
  sim <- NULL
  point <- NULL
  if (length(computed) > 0L && any(complete) && !(is.null(boot$point) && !is.null(boot))) {
    point_mu <- array(0, dim = c(1L, length(cols), length(info)))
    for (i in seq_along(info)) {
      point_mu[, , i] <- if (is.null(boot)) {
        taxon_means(matrix(info[[i]]$eta[cols], nrow = 1L), info[[i]], cols)
      } else {
        boot$point[cols, i]
      }
    }
    point <- community_statistics(point_mu, layout, info, computed, opts$null_models,
                                  point = TRUE)
    if (!nzchar(problem)) {
      sim <- if (is.null(boot)) {
        simulate_community(info, root, layout, dfs, opts, computed)
      } else {
        resampled_community(boot, layout, info, opts, computed)
      }
    }
  }
  out <- list()
  for (prop in opts$properties) {
    out[[length(out) + 1L]] <- property_rows(prop, s, in_set, times, layout, sim, point[[prop]],
                                             info, design, opts, dfs, problem, richness, labels)
  }
  bind_rows_base(out)
}

# Statistics of a draw (or of the point estimates) -------------------------------
#
# `mu` is an array [draws, cells of the set, taxa] of means. For each time
# with every cell, returns per property a matrix [draws, statistics] with the
# effects of the stressors alone, the combined effect and the expected
# combined effect of each null model, plus (for the point estimates) the
# property in the reference cell, in the combination and expected there.
community_statistics <- function(mu, layout, info, properties, null_models, point = FALSE) {
  m <- layout$m
  n_times <- ncol(layout$pos)
  families <- vapply(info, `[[`, character(1), "family")
  power <- vapply(info, `[[`, numeric(1), "power")
  disp_at <- function(col) vapply(info, function(x) x$disp[[layout$cols[[col]]]], numeric(1))
  presence <- function(x, col) {
    disp <- disp_at(col)
    out <- x
    for (i in seq_len(ncol(x))) {
      out[, i] <- presence_probability(x[, i], families[[i]], disp[[i]], power[[i]])
    }
    out
  }
  stats <- lapply(properties, function(p) vector("list", n_times))
  names(stats) <- properties
  nonpositive <- vector("list", n_times)
  for (t in seq_len(n_times)) {
    if (!layout$complete[[t]]) {
      next
    }
    pos <- layout$pos[, t]
    x0 <- mu[, pos[[1L]], , drop = TRUE]
    xs <- lapply(seq_len(m), function(j) mu[, pos[[1L + j]], , drop = TRUE])
    xc <- mu[, pos[[m + 2L]], , drop = TRUE]
    if (!is.matrix(x0)) {
      # One draw: keep [draws, taxa].
      x0 <- matrix(x0, nrow = 1L)
      xs <- lapply(xs, matrix, nrow = 1L)
      xc <- matrix(xc, nrow = 1L)
    }
    expected_taxa <- list(
      AD_taxa = if ("AD_taxa" %in% null_models) expected_taxa_means(x0, xs, "AD_taxa"),
      MU_taxa = if ("MU_taxa" %in% null_models) expected_taxa_means(x0, xs, "MU_taxa")
    )
    if ("AD_taxa" %in% null_models) {
      # Draws in which the sum of the effects on each taxon eliminates it.
      sums <- x0 + Reduce(`+`, lapply(xs, function(x) x - x0))
      nonpositive[[t]] <- colSums(sums <= 0)
    }
    for (p in properties) {
      value <- if (identical(p, "abundance")) {
        function(x, col) rowSums(x)
      } else {
        function(x, col) rowSums(presence(x, col))
      }
      r0 <- value(x0, pos[[1L]])
      rs <- vapply(seq_len(m), function(j) value(xs[[j]], pos[[1L + j]]), numeric(nrow(x0)))
      rs <- matrix(rs, nrow = nrow(x0))
      rc <- value(xc, pos[[m + 2L]])
      single <- rs / r0 - 1
      combined <- rc / r0 - 1
      expected <- vapply(null_models, function(nm) {
        switch(
          nm,
          AD_taxa = ,
          MU_taxa = value(expected_taxa[[nm]], pos[[m + 2L]]) / r0 - 1,
          AD = pmax(-1, rowSums(single)),
          MU = apply(1 + single, 1L, prod) - 1
        )
      }, numeric(nrow(x0)))
      expected <- matrix(expected, nrow = nrow(x0), dimnames = list(NULL, null_models))
      res <- cbind(single, combined = combined, expected)
      if (point) {
        attr(res, "reference") <- r0
        attr(res, "observed") <- rc
        if (!is.null(expected_taxa$MU_taxa)) {
          attr(res, "undefined") <- names(info)[is.nan(expected_taxa$MU_taxa[1L, ])]
        }
      }
      stats[[p]][[t]] <- res
    }
  }
  attr(stats, "nonpositive") <- nonpositive
  stats
}

# Means of the taxa expected in the combination from the reference (`x0`,
# draws x taxa) and the stressors alone (`xs`, a list): Thompson et al.
# (2018, eq. 3) for "AD_taxa"; their product for "MU_taxa", NaN for a taxon
# absent from the reference but present with a stressor alone.
expected_taxa_means <- function(x0, xs, null_model) {
  if (identical(null_model, "AD_taxa")) {
    return(pmax(x0 + Reduce(`+`, lapply(xs, function(x) x - x0)), 0))
  }
  m <- length(xs)
  out <- Reduce(`*`, xs) / x0^(m - 1L)
  absent <- x0 == 0
  all_absent <- Reduce(`&`, lapply(xs, function(x) x == 0))
  out[absent & all_absent] <- 0
  out[absent & !all_absent] <- NaN
  out
}

# Simulation draws ---------------------------------------------------------------
#
# Draws of the fixed effects of every taxon, in chunks: for each parameter, the
# standard normal deviates of the taxa are correlated by `root` (the upper
# factor of their correlation) and then turned into each taxon's deviations
# by its own covariance (Cholesky factor of its free parameters), and scaled
# to a multivariate t with each taxon's degrees of freedom. Returns, per
# property and time, the statistics of the t draws and of the normal draws
# (for the standard errors).
simulate_community <- function(info, root, layout, dfs, opts, properties) {
  n_sim <- opts$n_sim
  n_taxa <- length(info)
  n_par <- ncol(info[[1L]]$X)
  cols <- layout$cols
  chunk <- max(100L, min(n_sim, floor(2e6 / (length(cols) * n_taxa))))
  starts <- seq(1L, n_sim, by = chunk)
  X_cols <- lapply(info, function(x) x$X[cols, , drop = FALSE])
  keep <- list(t = NULL, normal = NULL)
  for (start in starts) {
    n <- min(chunk, n_sim - start + 1L)
    Z <- matrix(stats::rnorm(n * n_par * n_taxa), nrow = n * n_par)
    if (!is.null(root)) {
      Z <- Z %*% root
    }
    mu_t <- array(0, dim = c(n, length(cols), n_taxa))
    mu_n <- mu_t
    for (i in seq_len(n_taxa)) {
      x <- info[[i]]
      U <- matrix(Z[, i], nrow = n)
      dev_beta <- matrix(0, nrow = n, ncol = n_par)
      dev_beta[, x$free] <- U[, x$free, drop = FALSE] %*% x$root
      dev_eta <- dev_beta %*% t(X_cols[[i]])
      scale <- if (is.finite(dfs[[i]])) sqrt(dfs[[i]] / stats::rchisq(n, df = dfs[[i]])) else 1
      eta <- matrix(x$eta[cols], nrow = n, ncol = length(cols), byrow = TRUE)
      mu_t[, , i] <- taxon_means(eta + dev_eta * scale, x, cols)
      mu_n[, , i] <- taxon_means(eta + dev_eta, x, cols)
    }
    part_t <- community_statistics(mu_t, layout, info, properties, opts$null_models)
    part_n <- community_statistics(mu_n, layout, info, properties, opts$null_models)
    keep$t <- append_draws(keep$t, part_t)
    keep$normal <- append_draws(keep$normal, part_n)
    counts <- attr(part_t, "nonpositive")
    keep$nonpositive <- if (is.null(keep$nonpositive)) {
      counts
    } else {
      Map(function(a, b) if (is.null(a)) NULL else a + b, keep$nonpositive, counts)
    }
  }
  keep
}

# Statistics of the resampled means of the cells of a set, in the shape of
# simulate_community() (the same draws serve for the standard errors).
resampled_community <- function(boot, layout, info, opts, properties) {
  draws <- community_statistics(boot$draws[, layout$cols, , drop = FALSE], layout, info,
                                properties, opts$null_models)
  list(t = draws, normal = draws, nonpositive = attr(draws, "nonpositive"))
}

# Means of a taxon from values of its linear predictor (draws x cells):
# zero at the boundary and never negative.
taxon_means <- function(eta, x, cols) {
  mu <- x$link$linkinv(eta)
  mu[, x$at_boundary[cols]] <- 0
  pmax(mu, 0)
}

append_draws <- function(acc, part) {
  if (is.null(acc)) {
    return(part)
  }
  for (p in names(part)) {
    for (t in seq_along(part[[p]])) {
      if (!is.null(part[[p]][[t]])) {
        acc[[p]][[t]] <- rbind(acc[[p]][[t]], part[[p]][[t]])
      }
    }
  }
  acc
}

# Rows of the estimates ----------------------------------------------------------

property_rows <- function(prop, s, in_set, times, layout, sim, point, info, design, opts,
                          dfs, problem, richness, labels) {
  m <- layout$m
  level <- opts$level
  control <- design$spec$control
  f <- list(response = prop)
  df_min <- min(dfs)
  rows <- list()
  per_time <- list()
  for (t in seq_along(times)) {
    status <- "ok"
    reason <- ""
    if (identical(prop, "richness") && nzchar(richness)) {
      status <- "not_defined"
      reason <- richness
    } else if (!layout$complete[[t]]) {
      missing <- labels[is.na(layout$idx[, t])]
      status <- "not_estimable"
      reason <- sprintf(
        if (is.na(times[[t]])) "No values of the taxa in %s." else "No values of the taxa at this time in %s.",
        format_items(quote_items(short_cells(missing, control)), max = 3L)
      )
    } else if (nzchar(problem)) {
      status <- "not_estimable"
      reason <- problem
    } else if (!(attr(point[[t]], "reference") > 0)) {
      status <- "not_estimable"
      reason <- sprintf(
        "The %s is zero in the reference cell: its relative changes are not defined.",
        if (identical(prop, "abundance")) "total abundance" else "richness"
      )
    }
    uncertain <- if (identical(status, "ok") && identical(prop, "richness")) {
      uncertain_persistence(sim$nonpositive[[t]], opts$n_draws, names(info))
    } else {
      character()
    }
    entries <- list()
    for (nm in opts$null_models) {
      for (sp in community_specs(m, in_set)) {
        values <- if (identical(status, "ok")) {
          community_values(sp, nm, point[[t]], sim$t[[prop]][[t]], sim$normal[[prop]][[t]],
                           level, df_min, opts$n_draws, opts$inference)
        } else {
          empty_values(status, reason, level, opts$inference)
        }
        if (identical(nm, "AD_taxa") && sp$component %in% c("expected", "interaction") &&
            identical(values$status, "ok") && length(uncertain) > 0L) {
          values <- persistence_values(uncertain, level, opts$inference)
        }
        row <- estimate_row(s, f, design, times[[t]], nm, sp, values)
        row$scale <- "response / control"
        rows[[length(rows) + 1L]] <- row
        entries[[length(entries) + 1L]] <- list(null_model = nm, sp = sp, values = values)
      }
    }
    per_time[[t]] <- entries
  }
  if (opts$overall && length(times) >= 2L) {
    rows <- c(rows, overall_community(s, f, design, per_time, opts, df_min))
  }
  bind_rows_base(rows)
}

community_specs <- function(m, in_set) {
  singles <- lapply(seq_len(m), function(j) {
    list(component = "single", stressor = in_set[[j]], estimand = NA_character_, column = j)
  })
  c(singles, list(
    list(component = "expected", stressor = NA_character_, estimand = NA_character_),
    list(component = "combined", stressor = NA_character_, estimand = NA_character_),
    list(component = "interaction", stressor = NA_character_, estimand = "net")
  ))
}

# Values of one row from the point estimates and the draws of its time.
community_values <- function(sp, null_model, point, draws_t, draws_n, level, df, n_sim,
                             inference = "simulation") {
  pick <- function(stats) {
    switch(
      sp$component,
      single = stats[, sp$column],
      combined = stats[, "combined"],
      expected = stats[, null_model],
      interaction = stats[, "combined"] - stats[, null_model]
    )
  }
  estimate <- unname(pick(point))
  if (!is.finite(estimate)) {
    undefined <- attr(point, "undefined")
    if (identical(null_model, "MU_taxa") && length(undefined) > 0L) {
      return(empty_values("not_defined", sprintf(
        paste(
          "The multiplicative compositional model is not defined: %s had only zeros",
          "in the reference cell but not with some stressor alone."
        ),
        format_items(quote_items(undefined), max = 5L)
      ), level, inference))
    }
    return(empty_values("not_estimable", "This quantity is not finite.", level, inference))
  }
  interaction <- identical(sp$component, "interaction")
  community_summary(
    estimate, pick(draws_t), pick(draws_n), sp, level, df, n_sim, inference = inference,
    observed = if (interaction) attr(point, "observed") else NA_real_,
    expected = if (interaction) {
      unname(attr(point, "reference") * (1 + point[, null_model]))
    } else {
      NA_real_
    }
  )
}

# Summary of the draws of a row; draws that are not finite (a reference
# without any taxon in some draw) are left out and noted.
community_summary <- function(estimate, draws, normal, sp, level, df, n_sim,
                              inference = "simulation", observed = NA_real_,
                              expected = NA_real_) {
  ok <- is.finite(draws)
  note <- ""
  if (!all(ok)) {
    note <- "nonfinite_draws_dropped"
  }
  if (sum(ok) < 2L) {
    return(empty_values("not_estimable", sprintf("The %s values of this quantity are not finite.",
                                                 if (identical(inference, "resampling")) {
                                                   "resampled"
                                                 } else {
                                                   "simulated"
                                                 }),
                        level, inference))
  }
  s <- summarise_draws(draws[ok], level)
  p <- if (identical(sp$component, "expected")) NA_real_ else s$p_value
  list(
    estimate = estimate, se = stats::sd(normal[is.finite(normal)]), lower = s$lower,
    upper = s$upper, level = level, df = df, p_value = p,
    mc_se = mc_error(p, n_sim), observed = observed, expected = expected,
    relative = NA_real_, relative_lower = NA_real_, relative_upper = NA_real_,
    p_equivalence = NA_real_, rescaled = NA_real_, rescaled_class = NA_character_,
    inference = inference, status = "ok", reason = "", note = note,
    draws = draws, normal_draws = normal
  )
}

# Averages over the times of a set, draw by draw.
overall_community <- function(s, f, design, per_time, opts, df_min) {
  out <- list()
  for (k in seq_along(per_time[[1L]])) {
    first <- per_time[[1L]][[k]]
    parts <- lapply(per_time, function(p) p[[k]]$values)
    statuses <- vapply(parts, `[[`, character(1), "status")
    values <- if (all(statuses == "ok")) {
      estimate <- mean(vapply(parts, `[[`, numeric(1), "estimate"))
      draws <- rowMeans(do.call(cbind, lapply(parts, `[[`, "draws")))
      normal <- rowMeans(do.call(cbind, lapply(parts, `[[`, "normal_draws")))
      community_summary(estimate, draws, normal, first$sp, opts$level, df_min, opts$n_draws,
                        inference = opts$inference)
    } else if (any(statuses == "not_defined")) {
      # Not defined at some time: neither is the average.
      bad <- which(statuses == "not_defined")[[1L]]
      empty_values("not_defined", parts[[bad]]$reason, opts$level, opts$inference)
    } else {
      empty_values("not_estimable", "Not estimable (with an interval) at every sampling time.",
                   opts$level, opts$inference)
    }
    row <- estimate_row(s, f, design, "overall", first$null_model, first$sp, values)
    row$scale <- "response / control"
    out[[length(out) + 1L]] <- row
  }
  out
}

# Taxa whose elimination by the sum of the effects of the stressors (the
# additive compositional model) is uncertain: it happens in more than 20%
# and less than 80% of the draws. Their expected presence then jumps between
# zero and nearly one from draw to draw, and the expected richness is biased:
# in simulations its tests rejected true null models in about 11% of these
# cases at the 5% level, and in at most 4% of the others.
uncertain_persistence <- function(counts, n_sim, taxa) {
  if (is.null(counts)) {
    return(character())
  }
  share <- counts / n_sim
  taxa[share > 0.2 & share < 0.8]
}

persistence_values <- function(uncertain, level, inference = "simulation") {
  empty_values("not_estimable", sprintf(
    paste(
      "Whether %s would persist under the additive compositional model is too",
      "uncertain (the sum of the effects of the stressors on %s is below zero in",
      "20%% to 80%% of the draws): the expected richness is then biased",
      "and its test unreliable. The total abundance is not affected."
    ),
    format_items(quote_items(uncertain), max = 5L),
    if (length(uncertain) == 1L) "it" else "each"
  ), level, inference)
}

#' @export
print.must_community <- function(x, n = 12L, ...) {
  s <- x$settings
  cat(sprintf(
    "<must_community> %d taxa; properties: %s; null models: %s\n",
    length(s$taxa), paste(s$properties, collapse = ", "), paste(s$null_models, collapse = ", ")
  ))
  if (identical(s$inference, "resampling")) {
    cat(sprintf(
      paste0(
        "Inference: resampling of whole units with all their taxa (PIT-trap, %d resamples) ",
        "from the means of the samples; %s%% intervals; p-values adjusted with %s; seed %s\n"
      ),
      s$n_boot, format(100 * s$level), s$p_adjust, format(s$seed)
    ))
  } else {
    cat(sprintf(
      "Inference: Krinsky-Robb simulation (%d draws), %s; %s%% intervals; p-values adjusted with %s; seed %s\n",
      s$n_sim, correction_label(s, detail = TRUE), format(100 * s$level), s$p_adjust,
      format(s$seed)
    ))
  }
  cat(sprintf(
    "Taxa: %s\n", if (identical(s$taxa_correlation, "resampled")) {
      "resampled together within each unit (their correlation is kept, not modelled)"
    } else if (identical(s$taxa_correlation, "none")) {
      "independent (`taxa_correlation = \"none\"`)"
    } else if (isTRUE(s$shrinkage >= 1)) {
      "independent (their residual correlation was shrunk to zero)"
    } else {
      sprintf("correlated (residual correlation shrunk with intensity %s)",
              format(signif(s$shrinkage, 2)))
    }
  ))
  models <- c(AD_taxa = "AD_taxa = additive on each taxon (Thompson et al. 2018)",
              MU_taxa = "MU_taxa = multiplicative on each taxon",
              AD = "AD = additive on the property", MU = "MU = multiplicative on the property")
  cat("Null models: ", paste(models[s$null_models], collapse = "; "), "\n", sep = "")
  cat("Scale: relative changes of the property from the reference cell\n")
  est <- x$estimates
  int <- est[est$component == "interaction", , drop = FALSE]
  if (nrow(int) == 0L) {
    cat("No interaction could be estimated.\n")
  } else {
    show <- data.frame(
      property = int$response,
      time = int$time,
      combination = set_combination(int$levels, int$stressors),
      background = compact_levels(int$background, x$design$spec$control),
      model = int$null_model,
      observed = signif(int$observed, 3),
      expected = signif(int$expected, 3),
      estimate = signif(int$estimate, 3),
      lower = signif(int$lower, 3),
      upper = signif(int$upper, 3),
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
  if (nzchar(x$diagnostics$richness) && "richness" %in% s$properties) {
    cat("Richness not defined: ", x$diagnostics$richness, "\n", sep = "")
  }
  invisible(x)
}
