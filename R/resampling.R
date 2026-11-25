# Resampling tests of the null models AD and MU (internal) ---------------------
#
# A test compares, for one response or for each taxon of a community, a null
# model with the full model (one mean per treatment combination, and time)
# and resamples whole experimental units under the null model with the
# PIT-trap (Warton et al. 2017): the probability integral transform
# residuals of all the responses of a unit, at all its times, are resampled
# together and turned back into values with the quantile function of the
# null model at the samples that receive them.
#
# The null models constrain the means of the cells of a contrast set at
# each time: the reference (role 0), each stressor alone (roles 1..m) and
# their full combination (role m + 1).
# * "MU": log mu = b0 + sum(bj sj), a GLM with a log link (stats::glm.fit()).
# * "AD": mu_c = mu_0 + sum(mu_j - mu_0), an identity link. For responses
#   that cannot be negative the means are kept non-negative: when AD expects
#   a negative mean for the combination, the null model is fitted at the
#   boundary, and an absent combination deviates from it (Madge Pimentel et
#   al. 2025, Fig. 1: flooring).
# * "AD_taxa": AD truncated at zero (Thompson et al. 2018, eq. 3): when the
#   means of the reference and the stressors alone add up to zero or less and
#   the combination has nothing, the null model fits it.
# The means of each time are fitted by maximum quasi-likelihood with the
# variance function of the family, with its shape parameters fixed (the
# power of a Tweedie, the size of a negative binomial); the dispersion of
# the null model (Pearson) gives the distribution of the PIT-trap.

# Families ---------------------------------------------------------------------

rs_family_names <- function() {
  c("poisson", "nbinom2", "nbinom1", "gaussian", "gamma", "lognormal", "tweedie")
}

# Variance mu^p (p = 0 Gaussian, 1 Poisson and NB1, 2 gamma and lognormal,
# the Tweedie power) or mu + mu^2 / theta (NB2). `continuous`: a dispersion
# parameter of a continuous distribution (its quantiles are the means when
# it is zero). The statistic of every family is the drop in deviance scaled
# by the deviance of the null model (see rs_statistic()).
rs_family <- function(name, power = NA_real_, theta = NA_real_) {
  list(
    name = name,
    var_power = switch(name, gaussian = 0, poisson = , nbinom1 = 1, gamma = , lognormal = 2,
                       tweedie = power, nbinom2 = NA_real_),
    theta = theta,
    support = switch(name, gaussian = "real", gamma = , lognormal = "positive", "nonnegative"),
    discrete = name %in% c("poisson", "nbinom1", "nbinom2"),
    continuous = name %in% c("gaussian", "gamma", "lognormal", "tweedie")
  )
}

rs_variance <- function(mu, fam) {
  if (identical(fam$name, "nbinom2")) mu + mu^2 / fam$theta else mu^fam$var_power
}

rs_unit_deviance <- function(y, mu, fam) {
  mu <- rep_len(mu, length(y))
  if (identical(fam$name, "nbinom2")) {
    th <- fam$theta
    d <- 2 * (y + th) * log((mu + th) / (y + th))
    pos <- y > 0
    d[pos] <- d[pos] + 2 * y[pos] * log(y[pos] / mu[pos])
    return(pmax(d, 0))
  }
  p <- fam$var_power
  if (p == 0) {
    return((y - mu)^2)
  }
  if (p == 1) {
    d <- 2 * mu
    pos <- y > 0
    d[pos] <- 2 * (y[pos] * log(y[pos] / mu[pos]) - (y[pos] - mu[pos]))
    return(pmax(d, 0))
  }
  if (p == 2) {
    return(pmax(2 * (log(mu / y) + (y - mu) / mu), 0))
  }
  tweedie_unit_deviance(y, mu, p)
}

# Family object of stats::glm.fit() with a log link and the variance and
# deviance of `fam`.
rs_glm_family <- function(fam) {
  force(fam)
  structure(list(
    family = paste0("must_", fam$name),
    link = "log",
    linkfun = function(mu) log(mu),
    linkinv = function(eta) pmax(exp(eta), .Machine$double.eps),
    mu.eta = function(eta) pmax(exp(eta), .Machine$double.eps),
    valideta = function(eta) TRUE,
    variance = function(mu) rs_variance(mu, fam),
    validmu = function(mu) all(is.finite(mu)) && all(mu > 0),
    dev.resids = function(y, mu, wt) wt * rs_unit_deviance(y, mu, fam),
    aic = function(y, n, mu, wt, dev) NA_real_,
    initialize = expression({
      n <- rep(1, nobs)
      mustart <- pmax(y, 0) + 0.1 * (y <= 0)
    }),
    power = fam$var_power
  ), class = "family")
}

# Dispersion of the null model for the PIT-trap: the variance (Gaussian),
# Pearson's estimate (gamma, lognormal as its coefficient of variation
# squared, Tweedie), the excess over Poisson (NB1); NB2 keeps its size.
rs_dispersion <- function(y, mu, fam, df) {
  if (identical(fam$name, "poisson")) {
    return(NA_real_)
  }
  if (identical(fam$name, "nbinom2")) {
    return(fam$theta)
  }
  if (identical(fam$name, "gaussian")) {
    return(sum((y - mu)^2) / df)
  }
  ok <- mu > 0
  x2 <- sum((y[ok] - mu[ok])^2 / rs_variance(mu[ok], fam)) / df
  if (identical(fam$name, "nbinom1")) max(x2 - 1, 0) else x2
}

rs_cdf_discrete <- function(q, mu, disp, fam) {
  out <- switch(
    fam$name,
    poisson = stats::ppois(q, mu),
    nbinom2 = stats::pnbinom(q, size = fam$theta, mu = mu),
    nbinom1 = if (disp > 0) {
      stats::pnbinom(q, size = pmax(mu, 1e-300) / disp, mu = mu)
    } else {
      stats::ppois(q, mu)
    }
  )
  out[mu <= 0] <- as.numeric(q[mu <= 0] >= 0)
  out
}

# Probability integral transform residuals (randomised for counts and for
# the zeros of a Tweedie), kept inside (0, 1).
rs_pit <- function(y, mu, disp, fam) {
  n <- length(y)
  u <- switch(
    fam$name,
    gaussian = stats::pnorm(y, mu, sqrt(disp)),
    gamma = stats::pgamma(y, shape = 1 / disp, scale = mu * disp),
    lognormal = {
      s2 <- log1p(disp)
      stats::plnorm(y, log(mu) - s2 / 2, sqrt(s2))
    },
    tweedie = tweedie_pit(y, mu, disp, fam$var_power),
    {
      hi <- rs_cdf_discrete(y, mu, disp, fam)
      lo <- numeric(n)
      pos <- y > 0
      lo[pos] <- rs_cdf_discrete(y[pos] - 1, mu[pos], disp, fam)
      lo + stats::runif(n) * (hi - lo)
    }
  )
  pmin(pmax(u, 1e-10), 1 - 1e-10)
}

rs_quantile <- function(u, mu, disp, fam) {
  if (fam$continuous && !(disp > 0)) {
    return(mu)
  }
  switch(
    fam$name,
    gaussian = stats::qnorm(u, mu, sqrt(disp)),
    gamma = stats::qgamma(u, shape = 1 / disp, scale = mu * disp),
    lognormal = {
      s2 <- log1p(disp)
      stats::qlnorm(u, log(mu) - s2 / 2, sqrt(s2))
    },
    tweedie = qtweedie_series(u, mu, disp, fam$var_power),
    poisson = stats::qpois(u, mu),
    nbinom2 = stats::qnbinom(u, size = fam$theta, mu = mu),
    nbinom1 = {
      out <- if (disp > 0) {
        stats::qnbinom(u, size = pmax(mu, 1e-300) / disp, mu = mu)
      } else {
        stats::qpois(u, mu)
      }
      out[mu <= 0] <- 0
      out
    }
  )
}

# Fits ---------------------------------------------------------------------------

# Layout of the samples of a test: the role of each sample, its time (1
# without times), the cells of the full model (time and role) and the
# design of the multiplicative model. `env` holds the indicators `away`
# (factors with levels FALSE and TRUE) and, with all times, `.time`.
rs_layout <- function(env, away, all_times) {
  S <- vapply(away, function(nm) as.character(env[[nm]]) == "TRUE", logical(nrow(env)))
  S <- matrix(S, nrow = nrow(env))
  m <- ncol(S)
  n_away <- rowSums(S)
  role <- ifelse(n_away == 0L, 0L, ifelse(n_away == m, m + 1L, max.col(S, ties.method = "first")))
  time <- if (all_times) as.integer(droplevels(factor(env$.time))) else rep(1L, nrow(env))
  main <- paste(away, collapse = " + ")
  x0 <- stats::model.matrix(stats::as.formula(
    if (all_times) paste("~ .time * (", main, ")") else paste("~", main)
  ), env)
  key <- (time - 1L) * (m + 2L) + role
  list(m = m, role = as.integer(role), time = time, n_times = max(time), x0 = x0,
       group = factor(key, levels = sort(unique(key))), s = S)
}

# Samples of a test: `d` holds the samples of the cells of a set at one time
# (or at all its times), with `.context`, `.cell`, `.time`, `.unit` and the
# stressor columns. Returns the indicators of the stressors of the set
# (`env`), the layout, the units, and why the units cannot be resampled
# whole ("" when they can).
# Why the units of a set are not independent replicates to resample: some
# stressor of the set is applied to whole groups of units (whole plots, or
# groups confounded with it), so that the groups are the replicates.
whole_group_problem <- function(design, in_set) {
  g <- design$groups
  if (is.null(g) || nrow(g) == 0L) {
    return("")
  }
  for (k in which(g$structure %in% c("whole_units", "confounded"))) {
    hit <- intersect(strsplit(g$between[[k]], ", ", fixed = TRUE)[[1L]], in_set)
    if (length(hit) > 0L) {
      return(sprintf(
        paste("%s %s applied to whole groups of %s, so the groups, not the units, are",
              "the replicates."),
        format_items(quote_items(hit), max = Inf), if (length(hit) == 1L) "is" else "are",
        g$group[[k]]
      ))
    }
  }
  ""
}

rs_samples <- function(d, in_set, control, all_times) {
  cell <- paste_key(d$.context, d$.cell)
  time <- if (all_times) droplevels(factor(d$.time)) else NULL
  unit <- droplevels(factor(as.character(d$.unit)))
  out <- list(problem = "", unit = unit)
  units_per_cell <- tapply(as.character(unit), if (all_times) paste(cell, time) else cell,
                           function(u) length(unique(u)))
  if (any(units_per_cell < 2L)) {
    out$problem <- "Some combination of the set has fewer than two units at this time."
    return(out)
  }
  samples <- if (all_times) table(unit, time) else table(unit)
  if (length(unique(as.integer(samples))) > 1L) {
    out$problem <- if (all_times) {
      paste(
        "Units sampled a different number of times: resampling whole units over",
        "all times needs every unit sampled at every time."
      )
    } else {
      "Units with a different number of samples at this time cannot be resampled whole."
    }
    return(out)
  }
  env <- data.frame(cell = factor(cell, levels = unique(cell)), stringsAsFactors = FALSE)
  away <- character()
  for (j in seq_along(in_set)) {
    nm <- paste0(".s", j)
    env[[nm]] <- factor(as.character(d[[in_set[[j]]]]) != control[[in_set[[j]]]],
                        levels = c(FALSE, TRUE))
    away <- c(away, nm)
  }
  if (all_times) {
    env$.time <- time
  }
  out$env <- env
  out$away <- away
  out$layout <- rs_layout(env, away, all_times)
  out
}

# Maximum likelihood rate (mean per unit of effort) of each group.
rs_rates <- function(y, e, group, fam) {
  g <- as.integer(group)
  k <- nlevels(group)
  if (identical(fam$name, "nbinom2")) {
    n <- tabulate(g, k)
    e_mean <- as.vector(rowsum(e, g, reorder = TRUE)) / n
    out <- as.vector(rowsum(y, g, reorder = TRUE)) / (n * e_mean)
    varying <- which(as.vector(rowsum(abs(e - e_mean[g]), g, reorder = TRUE)) >
                       1e-12 * e_mean)
    for (j in varying) {
      out[[j]] <- nb2_rate(y[g == j], e[g == j], fam$theta)
    }
    return(out)
  }
  p <- fam$var_power
  as.vector(rowsum(y * e^(1 - p), g, reorder = TRUE)) /
    as.vector(rowsum(e^(2 - p), g, reorder = TRUE))
}

nb2_rate <- function(y, e, theta) {
  if (all(y == 0)) {
    return(0)
  }
  score <- function(r) sum((y - e * r) / (e * r + theta))
  stats::uniroot(score, c(0, max(y / e) + 1), tol = 1e-12)$root
}

# Full model: one rate per group; rates per sample and the deviance.
rs_full <- function(y, e, group, fam) {
  rate <- rs_rates(y, e, group, fam)[as.integer(group)]
  list(rate = rate, deviance = sum(rs_unit_deviance(y, e * rate, fam)))
}

# Null model. Returns the means, the deviance and where its fit starts in a
# resample (`start`).
rs_null <- function(y, e, layout, fam, null_model, start = NULL) {
  if (identical(null_model, "MU")) {
    fit <- stats::glm.fit(layout$x0, y, offset = log(e), family = rs_glm_family(fam),
                          etastart = start, control = stats::glm.control(maxit = 100L))
    mu <- fit$fitted.values
    mu[y == 0 & mu <= 1e-10 * max(mu)] <- 0
    return(list(mu = mu, deviance = sum(rs_unit_deviance(y, mu, fam)),
                start = fit$linear.predictors, converged = fit$converged))
  }
  mu <- numeric(length(y))
  pars <- vector("list", layout$n_times)
  for (t in seq_len(layout$n_times)) {
    i <- layout$time == t
    fit <- rs_ad_time(y[i], e[i], layout$role[i], layout$m, fam,
                      truncate = identical(null_model, "AD_taxa"), start = start[[t]])
    mu[i] <- fit$mu
    pars[[t]] <- fit$par
  }
  list(mu = mu, deviance = sum(rs_unit_deviance(y, mu, fam)), start = pars, converged = TRUE)
}

# Additive null model at one time. Gaussian means are unconstrained (least
# squares). Otherwise the rates of the stressors alone are exp(a_j) and the
# reference s * sum(exp(a)) / (m - 1) with s = plogis(b), so that the
# combination, (1 - s) * sum(exp(a)), is never negative. When the effort is
# the same within each cell, the deviance is minimised from the cell means
# (for these families, sum(d(y, mu)) = n d(mean(y), mu) + constant).
rs_ad_time <- function(y, e, role, m, fam, truncate, start = NULL) {
  idx <- role + 1L
  if (identical(fam$support, "real")) {
    applied <- outer(role, seq_len(m), function(r, j) r == j | r == m + 1L)
    x <- e * cbind(1, applied * 1)
    coef <- stats::lm.fit(x, y)$coefficients
    return(list(mu = drop(x %*% coef), par = NULL))
  }
  free <- rs_rates(y, e, factor(role, levels = 0:(m + 1L)), fam)
  n <- tabulate(idx, m + 2L)
  empty <- as.vector(rowsum(as.numeric(y != 0), idx, reorder = TRUE)) == 0
  if (truncate && empty[[m + 2L]] && sum(free[2:(m + 1L)]) - (m - 1) * free[[1L]] <= 0) {
    return(list(mu = e * free[idx], par = NULL))
  }
  e_cell <- as.vector(rowsum(e, idx, reorder = TRUE)) / n
  if (all(abs(e - e_cell[idx]) <= 1e-12 * e_cell[idx])) {
    yy <- as.vector(rowsum(y, idx, reorder = TRUE)) / n
    ee <- e_cell
    ww <- n
    ii <- seq_len(m + 2L)
  } else {
    yy <- y
    ee <- e
    ww <- rep(1, length(y))
    ii <- idx
  }
  if (is.null(start)) {
    tiny <- 1e-8 * max(free, 1e-300)
    singles <- pmax(free[2:(m + 1L)], tiny)
    start <- c(log(singles), stats::qlogis(min(max((m - 1) * free[[1L]] / sum(singles), 1e-6),
                                               1 - 1e-6)))
  }
  opt <- stats::optim(start, ad_deviance, ad_gradient, method = "BFGS",
                      control = list(maxit = 1000L, reltol = 1e-10),
                      m = m, fam = fam, yy = yy, ee = ee, ww = ww, ii = ii)
  best <- list(rates = ad_rates(opt$par, m), value = opt$value)
  # Means of cells without individuals that the iterations drive towards zero
  # only approach it: the boundary (those means at zero) is fitted as well.
  near <- empty & best$rates <= 1e-3 * max(best$rates)
  if (any(near)) {
    boundary <- ad_boundary(opt$par, near, m, fam, yy, ee, ww, ii)
    if (!is.null(boundary) && boundary$value <= best$value + 1e-10 * (1 + abs(best$value))) {
      best <- boundary
    }
  }
  rates <- best$rates
  rates[rates <= 1e-10 * max(rates) & empty] <- 0
  list(mu = e * rates[idx], par = opt$par)
}

# Additive null model with the means of the cells `near` fixed at zero: a
# stressor alone (its rate), the reference (s = 0) or the combination
# (s = 1). Returns the rates and the deviance, or NULL.
ad_boundary <- function(prm, near, m, fam, yy, ee, ww, ii) {
  free_a <- which(!near[2:(m + 1L)])
  if (length(free_a) == 0L) {
    return(NULL)
  }
  modes <- c(if (near[[m + 2L]]) "one", if (near[[1L]]) "zero")
  if (length(modes) == 0L) {
    modes <- "free"
  }
  fits <- lapply(modes, function(mode) {
    rates_of <- function(p) {
      singles <- numeric(m)
      singles[free_a] <- exp(p[seq_along(free_a)])
      s <- switch(mode, one = 1, zero = 0, free = stats::plogis(p[[length(free_a) + 1L]]))
      total <- sum(singles)
      c(s * total / (m - 1), singles, (1 - s) * total)
    }
    deviance_of <- function(p) {
      mu <- ee * rates_of(p)[ii]
      pos <- yy > 0
      mu[pos] <- pmax(mu[pos], 1e-300)
      sum(ww * rs_unit_deviance(yy, mu, fam))
    }
    opt <- stats::optim(prm[c(free_a, if (identical(mode, "free")) m + 1L)], deviance_of,
                        method = "BFGS", control = list(maxit = 1000L, reltol = 1e-12))
    list(rates = rates_of(opt$par), value = opt$value)
  })
  fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
}

ad_rates <- function(prm, m) {
  singles <- exp(prm[seq_len(m)])
  s <- stats::plogis(prm[[m + 1L]])
  total <- sum(singles)
  c(s * total / (m - 1), singles, (1 - s) * total)
}

ad_deviance <- function(prm, m, fam, yy, ee, ww, ii) {
  mu <- ee * ad_rates(prm, m)[ii]
  pos <- yy > 0
  mu[pos] <- pmax(mu[pos], 1e-300)
  sum(ww * rs_unit_deviance(yy, mu, fam))
}

ad_gradient <- function(prm, m, fam, yy, ee, ww, ii) {
  singles <- exp(prm[seq_len(m)])
  s <- stats::plogis(prm[[m + 1L]])
  total <- sum(singles)
  mu <- pmax(ee * c(s * total / (m - 1), singles, (1 - s) * total)[ii], 1e-300)
  g <- -2 * ww * ee * (yy - mu) / rs_variance(mu, fam)
  if (length(g) != m + 2L) {
    g <- as.vector(rowsum(g, ii, reorder = TRUE))
  }
  c(singles * (g[[1L]] * s / (m - 1) + g[2:(m + 1L)] + g[[m + 2L]] * (1 - s)),
    s * (1 - s) * total * (g[[1L]] / (m - 1) - g[[m + 2L]]))
}

# The drop in deviance scaled by the deviance of the null model, df0 (D0 -
# D1) / D0, for every family: a likelihood ratio divided by the dispersion
# estimated from the null model, which does not depend on the dispersion
# assumed (the size of a negative binomial, fixed in the resamples, or the
# absence of overdispersion of a Poisson) nor on the unit of the response.
# The resamples of residuals of fitted models are less variable than the
# data, and in simulations a likelihood ratio with a fixed dispersion
# rejected true null models too often.
rs_statistic <- function(d0, d1, df0, fam) {
  scaled_drop(d0, d1, df0)
}

scaled_drop <- function(d0, d1, df0) {
  if (!is.finite(d0) || d0 <= 0) {
    return(0)
  }
  max(0, df0 * (d0 - d1) / d0)
}

# Test -----------------------------------------------------------------------------

# Resampling test of one null model on the columns of Y (one response, or the
# taxa of a community) with families `fams`, efforts `effort` (NULL or a
# matrix like Y), the layout of the samples and their units (a factor).
# Returns the statistic (sum over columns), its degrees of freedom and
# p-value, the statistics and step-down p-values of the columns, the rates
# of the full model (per sample, in the units of Y) and warnings.
rs_test <- function(Y, effort, layout, unit, fams, null_model, n_boot) {
  n <- nrow(Y)
  k <- ncol(Y)
  E <- if (is.null(effort)) matrix(1, n, k) else effort
  # Continuous non-negative responses are divided by their mean positive
  # value: these families are closed under scaling, which leaves the
  # statistics and residuals unchanged and keeps the means away from the
  # lower bound of glm.fit() whatever the unit of the response.
  scale <- vapply(seq_len(k), function(j) {
    if (fams[[j]]$discrete || identical(fams[[j]]$support, "real")) 1 else mean(Y[Y[, j] > 0, j])
  }, numeric(1))
  Y <- sweep(Y, 2L, scale, `/`)
  # Positions within units: the samples of a unit in time order.
  ord <- order(as.integer(unit), layout$time)
  position <- integer(n)
  position[ord] <- stats::ave(seq_len(n), unit[ord], FUN = seq_along)
  n_units <- nlevels(unit)
  sample_at <- matrix(NA_integer_, n_units, max(position))
  sample_at[cbind(as.integer(unit), position)] <- seq_len(n)
  from <- as.vector(t(sample_at[, position, drop = FALSE]))
  target <- rep(seq_len(n), times = n_units)
  n_par0 <- layout$n_times * (layout$m + 1L)
  df0 <- n - n_par0
  warnings <- character()
  keep_warning <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  statistic <- numeric(k)
  rate <- matrix(NA_real_, n, k, dimnames = list(NULL, colnames(Y)))
  null <- vector("list", k)
  boot_y <- vector("list", k)
  for (j in seq_len(k)) {
    fam <- fams[[j]]
    fit0 <- withCallingHandlers(rs_null(Y[, j], E[, j], layout, fam, null_model),
                                warning = keep_warning)
    full <- rs_full(Y[, j], E[, j], layout$group, fam)
    statistic[[j]] <- rs_statistic(fit0$deviance, full$deviance, df0, fam)
    rate[, j] <- full$rate * scale[[j]]
    disp <- rs_dispersion(Y[, j], fit0$mu, fam, df0)
    u <- rs_pit(Y[, j], fit0$mu, disp, fam)
    # Element [r, k]: the residual of unit k at the position of sample r,
    # given the null mean of sample r (computed once for each mean).
    key <- (match(fit0$mu, unique(fit0$mu))[target] - 1) * n + from
    first <- which(!duplicated(key))
    q <- rs_quantile(u[from[first]], fit0$mu[target[first]], disp, fam)
    boot_y[[j]] <- matrix(q[match(key, key[first])], nrow = n, ncol = n_units)
    null[[j]] <- fit0
  }
  boot <- matrix(0, n_boot, k)
  failed <- 0L
  for (b in seq_len(n_boot)) {
    drawn <- sample.int(n_units, n_units, replace = TRUE)[as.integer(unit)]
    for (j in seq_len(k)) {
      y <- boot_y[[j]][cbind(seq_len(n), drawn)]
      if (!any(y != 0)) {
        next
      }
      fam <- fams[[j]]
      fit0 <- tryCatch(
        suppressWarnings(rs_null(y, E[, j], layout, fam, null_model, start = null[[j]]$start)),
        error = function(e) NULL
      )
      if (is.null(fit0)) {
        failed <- failed + 1L
        next
      }
      full <- rs_full(y, E[, j], layout$group, fam)
      boot[b, j] <- rs_statistic(fit0$deviance, full$deviance, df0, fam)
    }
  }
  if (failed > 0L) {
    warnings <- c(warnings, sprintf("%d null models of resamples could not be fitted.", failed))
  }
  total <- sum(statistic)
  list(statistic = total, df = layout$n_times,
       p_value = (1 + sum(rowSums(boot) >= total - 1e-8)) / (n_boot + 1),
       taxa_statistic = stats::setNames(statistic, colnames(Y)),
       taxa_p = step_down_p(stats::setNames(statistic, colnames(Y)), boot), rate = rate,
       warnings = unique(warnings))
}

# Step-down adjusted p-values of the columns (free step-down resampling, as
# mvabund's p.uni = "adjusted"): successive maxima of the resampled
# statistics, from the smallest observed one up, and monotone p-values.
step_down_p <- function(observed, boot) {
  k <- length(observed)
  ord <- order(observed, decreasing = TRUE)
  q <- boot[, ord, drop = FALSE]
  if (k > 1L) {
    for (j in (k - 1L):1L) {
      q[, j] <- pmax(q[, j], q[, j + 1L])
    }
  }
  exceed <- colSums(q >= matrix(observed[ord] - 1e-8, nrow(q), k, byrow = TRUE))
  p <- cummax((exceed + 1) / (nrow(boot) + 1))
  out <- numeric(k)
  out[ord] <- p
  stats::setNames(out, names(observed))
}

# Net interaction of each column under AD or AD_taxa on the response scale,
# from the rates of the full model: mu_c - (mu_0 + sum(mu_j - mu_0)), with the
# expectation truncated at zero for AD_taxa; the average over times.
ad_net_deviation <- function(rate, layout, truncate) {
  m <- layout$m
  per_time <- vapply(seq_len(layout$n_times), function(t) {
    first <- vapply(0:(m + 1L), function(r) which(layout$time == t & layout$role == r)[[1L]],
                    integer(1))
    x <- rate[first, , drop = FALSE]
    expected <- x[1L, ] + colSums(x[2:(m + 1L), , drop = FALSE] -
                                    matrix(x[1L, ], m, ncol(x), byrow = TRUE))
    if (truncate) {
      expected <- pmax(expected, 0)
    }
    x[m + 2L, ] - expected
  }, numeric(ncol(rate)))
  rowMeans(matrix(per_time, nrow = ncol(rate)))
}

# Bootstrap of a community -----------------------------------------------------
#
# For estimate_community(inference = "resampling"): the PIT-trap under the
# full model of every taxon (one mean per treatment combination and time, the
# means of its samples), resampling whole units with all their taxa and
# times. Returns the means of each taxon at every cell of the grid (`point`,
# cells x taxa) and in each resample (`draws`, resamples x cells x taxa), or
# why it cannot be done (`problem`).
community_bootstrap <- function(fits, grid, n_boot) {
  cells <- as.character(grid$.cell_time)
  taxa <- names(fits)
  obs <- lapply(fits, observation_table)
  rows <- sort(unique(unlist(lapply(obs, `[[`, "row"), use.names = FALSE)))
  unit_of <- character(length(rows))
  cell_of <- character(length(rows))
  for (o in obs) {
    i <- match(o$row, rows)
    unit_of[i] <- o$unit
    cell_of[i] <- o$cell_time
  }
  # Any order of the times aligns the samples of the units, which are all
  # sampled at the same times.
  time_rank <- as.integer(factor(as.character(grid$.time)))[match(cell_of, cells)]
  time_rank[is.na(time_rank)] <- 1L
  unit <- factor(unit_of)
  out <- list(problem = "", point = NULL, draws = NULL)
  per_unit <- table(unit, time_rank)
  if (length(unique(as.integer(per_unit))) > 1L) {
    out$problem <- paste(
      "Resampling whole units needs every unit sampled the same number of times at",
      "every sampling time."
    )
    return(out)
  }
  fams <- lapply(fits, resampling_family_of)
  unsupported <- taxa[!vapply(fams, function(x) is.null(x$problem), logical(1))]
  if (length(unsupported) > 0L) {
    out$problem <- sprintf("%s cannot be resampled: %s",
                           format_items(quote_items(unsupported), max = 5L),
                           fams[[unsupported[[1L]]]]$problem)
    return(out)
  }
  n <- length(rows)
  ord <- order(as.integer(unit), time_rank)
  position <- integer(n)
  position[ord] <- stats::ave(seq_len(n), unit[ord], FUN = seq_along)
  n_units <- nlevels(unit)
  sample_at <- matrix(NA_integer_, n_units, max(position))
  sample_at[cbind(as.integer(unit), position)] <- seq_len(n)
  n_cells <- length(cells)
  k <- length(fits)
  point <- matrix(NA_real_, n_cells, k, dimnames = list(cells, taxa))
  draws <- array(NA_real_, dim = c(n_boot, n_cells, k))
  taxon <- vector("list", k)
  for (j in seq_len(k)) {
    f <- fits[[j]]
    frame <- f$model$frame
    i <- match(obs[[j]]$row, rows)
    y <- as.numeric(frame[[f$response]])
    offset <- grep("^offset\\(", names(frame), value = TRUE)
    e <- if (length(offset) > 0L) exp(frame[[offset[[1L]]]]) else rep(1, length(y))
    group <- factor(obs[[j]]$cell_time, levels = cells)
    fam <- fams[[j]]$family
    if (identical(fam$name, "nbinom2") && is.na(fam$theta)) {
      fam$theta <- nb2_size(y, e, droplevels(group))
    }
    rate <- rates_by_level(y, e, group, fam)
    point[, j] <- rate
    mu <- e * rate[as.integer(group)]
    disp <- rs_dispersion(y, mu, fam, length(y) - sum(!is.na(rate)))
    u <- rep(NA_real_, n)
    u[i] <- rs_pit(y, mu, disp, fam)
    # Element [r, k]: the residual of unit k at the position of sample r of
    # this taxon (NA where that unit has no value of the taxon).
    from <- as.vector(t(sample_at[, position[i], drop = FALSE]))
    target <- rep(seq_along(i), times = n_units)
    q <- rep(NA_real_, length(from))
    ok <- !is.na(u[from])
    key <- (match(mu, unique(mu))[target] - 1) * n + from
    first <- which(ok & !duplicated(key))
    q[ok] <- rs_quantile(u[from[first]], mu[target[first]], disp, fam)[
      match(key[ok], key[first])
    ]
    taxon[[j]] <- list(q = matrix(q, nrow = length(i), ncol = n_units), e = e, group = group,
                       fam = fam, unit = as.integer(unit)[i])
  }
  for (b in seq_len(n_boot)) {
    drawn <- sample.int(n_units, n_units, replace = TRUE)
    for (j in seq_len(k)) {
      tx <- taxon[[j]]
      y <- tx$q[cbind(seq_len(nrow(tx$q)), drawn[tx$unit])]
      ok <- !is.na(y)
      draws[b, , j] <- rates_by_level(y[ok], tx$e[ok], tx$group[ok], tx$fam)
    }
  }
  out$point <- point
  out$draws <- draws
  out
}

# Rates of the levels of `group` (NA for levels without samples).
rates_by_level <- function(y, e, group, fam) {
  present <- droplevels(group)
  out <- rep(NA_real_, nlevels(group))
  out[match(levels(present), levels(group))] <- rs_rates(y, e, present, fam)
  out
}
