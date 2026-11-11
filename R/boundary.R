# Contrasts that involve a combination whose values are all at a boundary
# (internal) -----------------------------------------------------------------
#
# On a log link, a treatment combination with only zeros has its mean at the
# boundary of the parameter space, and Wald-type or simulation intervals of
# contrasts that use it are not valid. When that combination is the full
# combination of a contrast set (and the control and the single stressors
# are not at the boundary), the null model still predicts a positive mean for
# it, E, from the other cells. The probability of observing only zeros if the
# true mean were E, averaged over the uncertainty of E (the simulation draws
# of the other cells), is a p-value for the null model: small when the zeros
# are far fewer than the null model allows to occur by chance. It is doubled,
# like the two-sided p-values of the other contrasts. The same holds for a
# binomial combination with only zeros or only successes on a logit link.

boundary_families <- function() {
  c("poisson", "nbinom2", "tweedie", "binomial", "betabinomial")
}

# Gauss-Hermite nodes and weights (Golub-Welsch) for integrals against
# exp(-x^2).
gauss_hermite <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- sqrt(i / 2)
  J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = sqrt(pi) * e$vectors[1L, ]^2)
}

# Probability that one observation is at the boundary given its mean `m`
# (count families: mean count; binomial families: probability).
prob_boundary_obs <- function(family, m, n, disp, power, type) {
  switch(
    family,
    poisson = exp(-m),
    nbinom2 = (disp / (disp + m))^disp,
    tweedie = exp(-m^(2 - power) / (disp * (2 - power))),
    binomial = if (identical(type, "zero")) (1 - m)^n else m^n,
    betabinomial = {
      a <- m * disp
      b <- (1 - m) * disp
      if (identical(type, "zero")) {
        exp(lbeta(a, b + n) - lbeta(a, b))
      } else {
        exp(lbeta(a + n, b) - lbeta(a, b))
      }
    },
    abort_internal(sprintf("no boundary probability for family '%s'.", family))
  )
}

# Groups of observations linked by a random-effect level: observations that
# share the level of any grouping factor in `factors` are in the same group
# (connected components). A grouping factor missing from `rows` links them
# all.
random_components <- function(rows, factors) {
  n <- nrow(rows)
  if (length(factors) == 0L) {
    return(seq_len(n))
  }
  if (!all(factors %in% names(rows))) {
    return(rep(1L, n))
  }
  comp <- seq_len(n)
  repeat {
    previous <- comp
    for (g in factors) {
      comp <- stats::ave(comp, as.character(rows[[g]]), FUN = min)
    }
    if (identical(comp, previous)) {
      break
    }
  }
  as.integer(factor(comp))
}

# Probability that all observations of a combination are at the boundary
# when its mean is `mean` (one value per draw, on the response scale: rate
# per unit of effort, or probability). Observations linked by a random-effect
# level (`factors`) share one intercept with the standard deviation of all
# random terms, `sd_re`, integrated by Gauss-Hermite quadrature; unlinked
# groups are independent. This is exact with one grouping factor (units, or
# rooms) and conservative with nested or crossed ones, where it overstates
# the dependence and hence the probability.
prob_all_boundary <- function(mean, rows, family, link, disp, power, sd_re, type,
                              factors = ".unit") {
  out <- rep(NA_real_, length(mean))
  certain <- if (identical(link, "log")) {
    mean <= 0
  } else if (identical(type, "zero")) {
    mean <= 0
  } else {
    mean >= 1
  }
  impossible <- if (!identical(link, "logit")) {
    rep(FALSE, length(mean))
  } else if (identical(type, "zero")) {
    mean >= 1
  } else {
    mean <= 0
  }
  out[certain] <- 1
  out[impossible & !certain] <- 0
  todo <- which(is.na(out) & is.finite(mean))
  if (length(todo) == 0L) {
    return(out)
  }
  gh <- if (sd_re > 0) gauss_hermite(20L) else list(nodes = 0, weights = sqrt(pi))
  b <- sqrt(2) * sd_re * gh$nodes
  w <- gh$weights / sqrt(pi)
  m0 <- mean[todo]
  log_p <- numeric(length(todo))
  components <- if (sd_re > 0) random_components(rows, factors) else seq_len(nrow(rows))
  for (u in split(seq_len(nrow(rows)), components)) {
    acc <- matrix(1, nrow = length(todo), ncol = length(b))
    for (j in u) {
      m <- if (identical(link, "log")) {
        effort <- if (is.null(rows$.effort)) 1 else rows$.effort[[j]]
        outer(m0 * effort, exp(b))
      } else {
        stats::plogis(outer(stats::qlogis(m0), b, "+"))
      }
      n <- if (is.null(rows$.trials)) NA_real_ else rows$.trials[[j]]
      acc <- acc * prob_boundary_obs(family, m, n, disp[[j]], power, type)
    }
    log_p <- log_p + log(pmax(drop(acc %*% w), .Machine$double.xmin))
  }
  out[todo] <- exp(log_p)
  out
}

# Dispersion, Tweedie power and random-effect standard deviation needed by
# prob_all_boundary() for the rows of a combination.
boundary_parameters <- function(model, rows) {
  family <- stats::family(model)$family
  disp <- rep(NA_real_, nrow(rows))
  if (family %in% c("nbinom2", "tweedie", "betabinomial")) {
    all_disp <- tryCatch(stats::predict(model, type = "disp"), error = function(e) NULL)
    disp <- if (is.null(all_disp)) rep(stats::sigma(model), nrow(rows)) else all_disp[rows$.row]
  }
  power <- if (identical(family, "tweedie")) {
    unname(glmmTMB::family_params(model)[[1L]])
  } else {
    NA_real_
  }
  vc <- tryCatch(glmmTMB::VarCorr(model)$cond, error = function(e) NULL)
  sd_re <- if (is.null(vc) || length(vc) == 0L) {
    0
  } else {
    sqrt(sum(vapply(vc, function(v) max(attr(v, "stddev"))^2, numeric(1))))
  }
  # Grouping factors of the random terms; a factor used twice gets a suffix
  # (".unit.1").
  factors <- if (is.null(vc)) character() else unique(sub("\\.[0-9]+$", "", names(vc)))
  list(family = family, link = stats::family(model)$link, disp = disp, power = power,
       sd_re = sd_re, factors = factors)
}

# Expected mean of the cell tested in a contrast under each null model, from
# a matrix of means (one row per draw, one column per subset in
# all_subsets() order; the column of the tested cell is not used). The tested
# cell is the single-stressor cell for single effects and the full
# combination otherwise; for single and combined effects the null
# hypothesis is equality with the reference cell.
expected_full <- function(M, component, estimand, null_model, m, weights, shape = NULL) {
  y0 <- M[, 1L]
  singles <- M[, 1L + seq_len(m), drop = FALSE]
  full <- ncol(M)
  if (component %in% c("combined", "single")) {
    return(y0)
  }
  if (identical(estimand, "pure")) {
    if (identical(null_model, "BL")) {
      # No emergent interaction: w_XYZ = w_X w_YZ + w_Y w_XZ + w_Z w_XY - 2 w_X w_Y w_Z.
      W <- M / y0
      return(y0 * (W[, 2L] * W[, 7L] + W[, 3L] * W[, 6L] + W[, 4L] * W[, 5L] -
                     2 * W[, 2L] * W[, 3L] * W[, 4L]))
    }
    rest <- M[, -full, drop = FALSE]
    w <- weights[-full]
    return(if (identical(null_model, "MU")) exp(-drop(log(rest) %*% w)) else -drop(rest %*% w))
  }
  switch(
    null_model,
    AD = y0 + rowSums(singles - y0),
    MU = ,
    BL = y0 * apply(singles / y0, 1L, prod),
    DO = pmin(apply(singles, 1L, min), y0),
    SAM = y0 * sam_expected(singles / y0, shape),
    abort_internal(sprintf("no expected boundary mean for null model '%s'.", null_model))
  )
}
