# Uncertainty of the contrasts (internal) --------------------------------------
#
# Two corrections for small samples, applied together with
# correction = "small_sample":
#   * the covariance of the fixed effects comes from the REML refit, whose
#     variance estimates account for the fixed effects estimated;
#   * the sampling distribution is a t with the degrees of freedom of the
#     experimental units, instead of a normal: Krinsky-Robb draws come from a
#     multivariate t (normal draws scaled by sqrt(df / chi-square(df))) and
#     the delta method uses t quantiles.
# correction = "none" reproduces Madge Pimentel et al. (2025): maximum
# likelihood covariance and normal draws.

# Degrees of freedom of the error that applies to a contrast set, after the
# strata of the design ("between-within" rule):
# * without repeated measures: observations minus estimated means;
# * with repeated measures or subsamples: experimental units minus the
#   treatment combinations they are nested in (combination x time when each
#   unit is sampled at a single time);
# * minus one per extra block of each blocking column in the model;
# * for groups of units to which some stressors are applied whole (whole
#   plots), the groups minus the combinations of those stressors (minus
#   blocks) form their own stratum: sets with such a stressor take at most
#   those degrees of freedom, and the other sets lose them.
# Families without a dispersion parameter and without random effects
# (Poisson, binomial) keep the normal distribution, as in glm().
set_df <- function(info, s, correction) {
  if (identical(correction, "none") || is.null(info) ||
      (!info$has_dispersion && !info$has_random)) {
    return(Inf)
  }
  if (identical(info$dispersion, "cell") && info$has_dispersion) {
    return(NA_real_)
  }
  # One dispersion per `by` group only matters for families that have one.
  per_context <- identical(info$dispersion, "by") && info$has_dispersion
  rows <- if (per_context) info$contexts[info$contexts$context == s$context, , drop = FALSE] else info$contexts
  base <- if (info$repeated) {
    sum(rows$n_units) - sum(rows$n_unit_groups)
  } else {
    (if (per_context) sum(rows$n_obs) - sum(rows$rank) else sum(rows$n_obs) - info$rank_total)
  }
  base <- base - sum(rows$block_df)
  in_set <- strsplit(s$stressors, " x ", fixed = TRUE)[[1L]]
  limits <- numeric()
  for (w in info$whole) {
    wr <- if (per_context) w$contexts[w$contexts$context == s$context, , drop = FALSE] else w$contexts
    whole_df <- sum(wr$n_groups) - sum(wr$n_combos) - sum(wr$block_df)
    if (any(in_set %in% w$between)) {
      limits <- c(limits, whole_df)
    } else {
      base <- base - max(whole_df, 0)
    }
  }
  min(c(base, limits))
}

# Welch-Satterthwaite degrees of freedom of a contrast whose cells have their
# own dispersion: `v` are the variance contributions of the cells and `m` their
# numbers of units.
welch_df <- function(v, m) {
  keep <- is.finite(v) & v > 0
  if (!any(keep)) {
    return(Inf)
  }
  if (any(m[keep] < 2)) {
    return(0)
  }
  sum(v[keep])^2 / sum(v[keep]^2 / (m[keep] - 1))
}

# Predictions of the selected model at the grid of treatment combinations and
# times, with what the estimators need. Simulation draws are made lazily, on
# the first request, so that the delta method uses no random numbers unless a
# boundary test needs them.
prepare_predictions <- function(f, model, n_sim, complement = FALSE) {
  beta <- glmmTMB::fixef(model)$cond
  V <- as.matrix(stats::vcov(model)$cond)
  grid <- f$grid
  X <- stats::model.matrix(stats::as.formula(paste("~", f$fixed)), grid)
  if (!identical(colnames(X), names(beta))) {
    abort_internal("the prediction matrix does not match the fixed effects of the model.")
  }
  # Parameters fixed at a boundary (map) have no variance.
  fixed_par <- is.na(diag(V))
  V[fixed_par, ] <- 0
  V[, fixed_par] <- 0
  fam <- stats::family(model)
  eta <- drop(X %*% beta)
  env <- new.env(parent = emptyenv())
  env$dev_eta <- NULL
  env$scales <- list()
  env$problem <- ""
  draws_dev <- function() {
    if (is.null(env$dev_eta) && !nzchar(env$problem)) {
      free <- which(!fixed_par)
      L <- tryCatch(chol(V[free, free, drop = FALSE]), error = function(e) NULL)
      if (is.null(L)) {
        env$problem <- "The covariance matrix of the fixed effects is not positive definite."
      } else {
        Z <- matrix(stats::rnorm(n_sim * length(free)), nrow = n_sim)
        dev_beta <- matrix(0, nrow = n_sim, ncol = length(beta))
        dev_beta[, free] <- Z %*% L
        env$dev_eta <- dev_beta %*% t(X)
      }
    }
    env$dev_eta
  }
  scale_for <- function(df) {
    if (!is.finite(df)) {
      return(1)
    }
    key <- format(df, digits = 10)
    if (is.null(env$scales[[key]])) {
      env$scales[[key]] <- sqrt(df / stats::rchisq(n_sim, df = df))
    }
    env$scales[[key]]
  }
  to_mean <- function(eta_values) {
    mu <- fam$linkinv(eta_values)
    if (complement) 1 - mu else mu
  }
  list(
    grid = grid,
    eta = eta,
    mu = to_mean(eta),
    dmu = if (complement) -fam$mu.eta(eta) else fam$mu.eta(eta),
    V_eta = X %*% V %*% t(X),
    n_sim = n_sim,
    family = fam,
    model = model,
    complement = complement,
    at_boundary = as.character(grid$.cell_time) %in% names(f$boundary),
    boundary_type = f$boundary[match(as.character(grid$.cell_time), names(f$boundary))],
    sim_problem = function() {
      draws_dev()
      env$problem
    },
    # Draws of the means of the grid columns `cols`, as a matrix with one row
    # per draw, from a multivariate t with `df` degrees of freedom (normal if
    # df is infinite).
    mu_draws = function(cols, df) {
      dev <- draws_dev()
      if (is.null(dev)) {
        return(NULL)
      }
      s <- scale_for(df)
      eta_d <- dev[, cols, drop = FALSE] * s
      eta_d <- eta_d + matrix(eta[cols], nrow = nrow(eta_d), ncol = length(cols), byrow = TRUE)
      to_mean(eta_d)
    }
  )
}

# Gradient of sum(w * f(mu)) with respect to the linear predictor.
delta_gradient <- function(w, mu, dmu, null_model) {
  if (identical(null_model, "MU")) w * dmu / mu else w * dmu
}

# Wald summaries with a t distribution (normal if df is infinite).
summarise_wald <- function(estimate, se, df, level) {
  crit <- if (is.finite(df)) stats::qt(1 - (1 - level) / 2, df) else stats::qnorm(1 - (1 - level) / 2)
  p <- if (se > 0) {
    if (is.finite(df)) 2 * stats::pt(-abs(estimate / se), df) else 2 * stats::pnorm(-abs(estimate / se))
  } else {
    NA_real_
  }
  list(se = se, lower = estimate - crit * se, upper = estimate + crit * se, p_value = p)
}
